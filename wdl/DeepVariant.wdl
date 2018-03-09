# Use DeepVariant to generate VCF & gVCF for one sample

workflow DeepVariant {
    # reference genome
    File ref_fasta_gz

    # Genomic range(s) to call. Either range or ranges_bed is required
    String? range
    File? ranges_bed

    # Read alignments - bam & bai (bai will be generated if not provided)
    # The output vcf/gvcf filename is derived from the bam's.
    File bam
    File? bai

    # DeepVariant model files (with no folder component)
    File model_tar

    # Advanced setting for DeepVariant make_examples
    Int? gvcf_gq_binsize

    # DeepVariant docker image tag
    # c.f. https://github.com/google/deepvariant/blob/master/docs/deepvariant-docker.md
    String? deepvariant_docker = "gcr.io/deepvariant-docker/deepvariant:0.5.2"

    call make_examples { input:
        ref_fasta_gz = ref_fasta_gz,
        range = range,
        ranges_bed = ranges_bed,
        bam = bam,
        bai = bai,
        gvcf_gq_binsize = gvcf_gq_binsize,
        deepvariant_docker = deepvariant_docker
    }

    call call_variants { input:
        model_tar = model_tar,
        examples_tar = make_examples.examples_tar,
        deepvariant_docker = deepvariant_docker
    }

    call postprocess_variants { input:
        ref_fasta_gz = ref_fasta_gz,
        call_variants_output = call_variants.call_variants_output,
        gvcf_tfrecords_tar = make_examples.gvcf_tfrecords_tar,
        deepvariant_docker = deepvariant_docker
    }

    output {
        File vcf_gz = postprocess_variants.vcf_gz
        File gvcf_gz = postprocess_variants.gvcf_gz
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
