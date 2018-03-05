# Use DeepVariant to generate VCF & gVCF for one sample. Scatters the variant
# calling across a given list of genomic ranges (typically chromosomes),
# fetching the necessary BAM slices using htsget.

import "DeepVariant.wdl" as sub
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
        call sub.DeepVariant { input:
            ref_fasta_gz = ref_fasta_gz,
            range = range,
            bam = htsget.bam,
            model_tar = model_tar,
            gvcf_gq_binsize = gvcf_gq_binsize,
            deepvariant_docker = deepvariant_docker
        }
    }

    # concatenate pieces to get final VCF & gVCF
    call t.bcftools_concat as concat_vcf { input:
        vcfs_gz = DeepVariant.vcf_gz,
        output_filename = "${accession}.vcf.gz"
    }
    call t.bcftools_concat as concat_gvcf { input:
        vcfs_gz = DeepVariant.gvcf_gz,
        output_filename = "${accession}.gvcf.gz"
    }

    output {
        File vcf_gz = concat_vcf.vcf_gz
        File gvcf_gz = concat_gvcf.vcf_gz
    }
}

