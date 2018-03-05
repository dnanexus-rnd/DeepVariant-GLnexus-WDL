# Use DeepVariant to generate VCF & gVCF for one sample. Scatters the variant
# calling across a given list of genomic ranges (typically chromosomes),
# fetching the necessary BAM slices using htsget.

import "tasks.wdl" as t

workflow htsget_DeepVariant {
    # reference genome
    File ref_fasta_gz

    # sample to be fetched using htsget e.g. NA12878
    String accession

    # htsget advanced settings
    String? htsget_endpoint
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
        call t.htsget { input:
            range = range,
            accession = accession,
            endpoint = htsget_endpoint,
            format = htsget_format,
            hts_ref_tar = htsget_ref_tar
        }
        call t.make_examples { input:
            ref_fasta_gz = ref_fasta_gz,
            bam = htsget.bam,
            range = range,
            deepvariant_docker = deepvariant_docker,
            gvcf_gq_binsize = gvcf_gq_binsize
        }
        call t.call_variants { input:
            model_tar = model_tar,
            examples_tar = make_examples.examples_tar,
            deepvariant_docker = deepvariant_docker
        }
        call t.postprocess_variants { input:
            ref_fasta_gz = ref_fasta_gz,
            call_variants_output = call_variants.call_variants_output,
            gvcf_tfrecords_tar = make_examples.gvcf_tfrecords_tar,
            deepvariant_docker = deepvariant_docker
        }
    }

    # concatenate pieces to get final VCF & gVCF
    call t.bcftools_concat as concat_vcf { input:
        vcfs_gz = postprocess_variants.vcf_gz,
        output_filename = "${accession}.vcf.gz"
    }
    call t.bcftools_concat as concat_gvcf { input:
        vcfs_gz = postprocess_variants.gvcf_gz,
        output_filename = "${accession}.gvcf.gz"
    }

    output {
        File vcf_gz = concat_vcf.vcf_gz
        File gvcf_gz = concat_gvcf.vcf_gz
    }
}

