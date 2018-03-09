# Use DeepVariant to generate VCF & gVCF for one sample. Scatters the variant
# calling across a given list of genomic ranges (typically chromosomes),
# fetching the necessary BAM slices using htsget.

import "DeepVariant.wdl" as dv

workflow htsget_DeepVariant {
    # reference genome
    File ref_fasta_gz

    # htsget endpoint e.g. https://htsnexus.rnd.dnanex.us/v1/reads/BroadHiSeqX_b37
    String htsget_endpoint

    # sample to be fetched using htsget e.g. NA12878
    String accession

    # htsget advanced settings
    String? htsget_format
    File? htsget_ref_tar

    # "chr1:1-100000" etc.
    Array[String]+ ranges

    # DeepVariant model files (tar with no folder component)
    File model_tar

    # DeepVariant advanced settings
    Int? gvcf_gq_binsize
    String? deepvariant_docker

    scatter (range in ranges) {
        call htsget { input:
            range = range,
            accession = accession,
            endpoint = htsget_endpoint,
            format = htsget_format,
            hts_ref_tar = htsget_ref_tar
        }
        call dv.DeepVariant { input:
            ref_fasta_gz = ref_fasta_gz,
            range = range,
            bam = htsget.bam,
            model_tar = model_tar,
            gvcf_gq_binsize = gvcf_gq_binsize,
            deepvariant_docker = deepvariant_docker
        }
    }

    # concatenate pieces to get final VCF & gVCF
    call bcftools_concat as concat_vcf { input:
        vcfs_gz = DeepVariant.vcf_gz,
        output_filename = "${accession}.vcf.gz"
    }
    call bcftools_concat as concat_gvcf { input:
        vcfs_gz = DeepVariant.gvcf_gz,
        output_filename = "${accession}.gvcf.gz"
    }

    output {
        File vcf_gz = concat_vcf.vcf_gz
        File gvcf_gz = concat_gvcf.vcf_gz
    }
}

# retrieve a BAM slice using htsget
task htsget {
    # htsget endpoint e.g. https://htsnexus.rnd.dnanex.us/v1/reads/BroadHiSeqX_b37
    String endpoint

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

