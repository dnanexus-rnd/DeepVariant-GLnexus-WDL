# retrieve a BAM slice using htsget
task htsget {
    # htsget endpoint
    String? endpoint = "https://htsnexus.rnd.dnanex.us/v1/reads/BroadHiSeqX_hg38"

    # accession/sample ID e.g. NA12878
    String accession

    # genomic range to retrieve, typically the full length of a chromosome.
    String range

    # default BAM; set to CRAM if the endpoint can only serve this accession as CRAM
    String? format
    
    # optional tar of ~/.cache/hts-ref for CRAM decoding; e.g.
    #   tar czvf my_hts_ref.tar.gz -C ~/.cache hts-ref/
    # if not provided, then samtools will download any necessary reference
    # sequences from the CRAM Reference Registry.
    File? hts_ref_tar

    command <<<
        set -ex -o pipefail
        export SHELL=/bin/bash
        apt-get update -qq && apt-get install -y -qq parallel pigz

        # install samtools from bioconda to get an up-to-date version with
        # htsget support (1.7+)
        conda config --add channels bioconda
        conda install -y samtools==1.7=1

        # formulate the htsget URL & query
        range_tsv=$(echo "${range}" | tr :- '\t' | tr -d ,)
        referenceName=$(echo "$range_tsv" | cut -f1)
        start=$(echo "$range_tsv" | cut -f2)
        end=$(echo "$range_tsv" | cut -f3)
        htsget_url=$(printf "%s/%s?referenceName=%s&start=%s&end=%s" "${endpoint}" "${accession}" "$referenceName" "$start" "$end")
        if [ -n "${format}" ]; then
            htsget_url=$(printf "%s&format=%s" "$htsget_url" "${format}")
        fi

        # figure out output name
        mkdir -p output/
        output_name="${accession}.$(echo "${range}" | tr :- _ | tr -d ,)"

        # extract any given CRAM reference sequences
        if [ -n "${hts_ref_tar}" ]; then
            mkdir -p ~/.cache
            pigz -dc "${hts_ref_tar}" | tar xv -C ~/.cache
        fi

        # get the BAM slice using samtools, using parallel just for retry logic
        parallel --retries 3 "samtools view -b -1 '$htsget_url' > 'output/{}.bam'" ::: "$output_name"
    >>>

    runtime {
        disks: "local-disk 64 HDD"
        docker: "continuumio/miniconda"
    }

    output {
        File bam = glob("output/*.bam")[0]
    }
}

# DeepVariant make_examples
task make_examples {
    File ref_fasta_gz

    # bam & bai (bai can be omitted in which case it will be generated temporarily)
    File bam
    File? bai

    # Genomic range to run on. Either range or ranges_bed is required.
    String? range
    File? ranges_bed

    Int? gvcf_gq_binsize
    
    String deepvariant_docker

    parameter_meta {
        ref_fasta_gz: "stream"
    }

    command <<<
        range_arg="${range}"
        if [ -n "${ranges_bed}" ]; then
            range_arg="${ranges_bed}"
        fi
        if [ -z "$range_arg" ]; then
            1>&2 echo "Either range or ranges_bed is required"
            exit 1
        fi
        binsize_arg=""
        if [ -n "${gvcf_gq_binsize}" ]; then
            binsize_arg="--gvcf_gq-binsize ${gvcf_gq_binsize}"
        fi

        set -ex -o pipefail
        export SHELL=/bin/bash
        apt-get update -qq && apt-get install -y -qq samtools
        (zcat "${ref_fasta_gz}" > ref.fasta && samtools faidx ref.fasta) &
        output_name=$(basename "${bam}" .bam)
        mv "${bam}" input.bam
        if [ -n "${bai}" ]; then
            mv "${bai}" input.bam.bai
        else
            samtools index "input.bam"
        fi
        wait -n
        ls -lh

        mkdir examples/ gvcf/ logs/
        output_fn="examples/$output_name.tfrecord@$(nproc).gz"
        gvcf_fn="gvcf/$output_name.gvcf.tfrecord@$(nproc).gz"

        seq 0 $(( `nproc` - 1 )) | NO_GCE_CHECK=True parallel --halt 2 -t --results logs/ \
            "/opt/deepvariant/bin/make_examples --mode calling --ref ref.fasta --reads input.bam --examples '$output_fn' --gvcf '$gvcf_fn' --task {} --regions '$range_arg' $binsize_arg 2>&1" > /dev/null

        mkdir tar/
        tar -cvzf "tar/$output_name.logs.tar.gz" -C logs/ .
        tar -cvf "tar/$output_name.examples.tar" -C examples/ .
        tar -cvf "tar/$output_name.gvcf_tfrecords.tar" -C gvcf/ .
    >>>

    runtime {
        docker: "${deepvariant_docker}"
        cpu: "32"
    }

    output {
        File examples_tar = glob("tar/*.examples.tar")[0]
        File gvcf_tfrecords_tar = glob("tar/*.gvcf_tfrecords.tar")[0]
        File logs_tar_gz = glob("tar/*.logs.tar.gz")[0]
    }
}

# DeepVariant call_variants (CPU)
task call_variants {
    # DeepVariant model files (with no folder component)
    File model_tar
    File examples_tar
    String deepvariant_docker

    parameter_meta {
        examples_tar: "stream"
        model_tar: "stream"
    }

    command <<<
        set -ex -o pipefail

        tar xvf "${model_tar}" &
        mkdir examples output
        tar xvf "${examples_tar}" -C examples/
        wait -n

        n_examples=$(find examples/ -type f | wc -l)
        output_name=$(basename "${examples_tar}" .examples.tar)

        NO_GCE_CHECK=True /opt/deepvariant/bin/call_variants \
            --outfile "output/$output_name.call_variants.tfrecord.gz" \
            --examples "examples/$output_name.tfrecord@$n_examples.gz" \
            --checkpoint model.ckpt
    >>>

    runtime {
        docker: "${deepvariant_docker}"
        cpu: "32"
    }

    output {
        File call_variants_output = glob("output/*.gz")[0]
    }
}

# DeepVariant postprocess_variants
task postprocess_variants {
    File ref_fasta_gz
    File gvcf_tfrecords_tar
    File call_variants_output
    String deepvariant_docker

    parameter_meta {
        ref_fasta_gz: "stream"
        gvcf_tfrecords_tar: "stream"
    }

    command <<<
        set -ex -o pipefail
        apt-get update -qq && apt-get install -y -qq samtools wget
        # download a multithreaded version of the tabix bgzip utility
        wget --quiet -O bgzip "https://github.com/dnanexus-rnd/GLnexus/blob/master/cli/dxapplet/resources/usr/local/bin/bgzip?raw=true"
        chmod +x bgzip

        zcat "${ref_fasta_gz}" > ref.fasta
        samtools faidx ref.fasta

        mkdir gvcf output
        tar xvf "${gvcf_tfrecords_tar}" -C gvcf/
        n_gvcf_tfrecords=$(find gvcf/ -type f | wc -l)
        output_name=$(basename "${call_variants_output}" .call_variants.tfrecord.gz)
        NO_GCE_CHECK=True /opt/deepvariant/bin/postprocess_variants \
            --ref ref.fasta --infile "${call_variants_output}" \
            --nonvariant_site_tfrecord_path "gvcf/$output_name.gvcf.tfrecord@$n_gvcf_tfrecords.gz" \
            --outfile "output/$output_name.vcf" \
            --gvcf_outfile "output/$output_name.gvcf"

        ./bgzip -@ $(nproc) output/*.vcf &
        ./bgzip -@ $(nproc) output/*.gvcf
        wait -n
    >>>

    runtime {
        docker: "${deepvariant_docker}"
        cpu: "8"
    }

    output {
        File vcf_gz = glob("output/*.vcf.gz")[0]
        File gvcf_gz = glob("output/*.gvcf.gz")[0]
    }
}

task bcftools_concat {
    Array[File]+ vcfs_gz
    String output_filename

    command <<<
        set -ex -o pipefail
        apt-get update -qq && apt-get install -y -qq bcftools
        mkdir out
        cat "${write_lines(vcfs_gz)}" | xargs -n 9999 bcftools concat -o "${output_filename}" -O z
    >>>

    runtime {
        docker: "ubuntu:xenial"
        disks: "local-disk 64 HDD"
    }

    output {
        File vcf_gz = "${output_filename}"
    }
}

task GLnexus {
    Array[File]+ gvcf
    Array[String]+ ranges
    String config
    String output_name

    command <<<
        set -ex -o pipefail
        export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.1
        cat "${write_lines(ranges)}" | tr :- '\t' | tr -d , > ranges.bed
        numactl --interleave=all glnexus_cli --config "${config}" --list --bed ranges.bed "${write_lines(gvcf)}" | bcftools view - | bgzip -@ 4 -c > "${output_name}.vcf.gz"
    >>>

    runtime {
        docker: "quay.io/mlin/glnexus:v0.3.3-24-g3db0285"
        disks: "local-disk 64 HDD"
        cpu: "16"
    }

    output {
        File pvcf_gz = "${output_name}.vcf.gz"
    }
}

