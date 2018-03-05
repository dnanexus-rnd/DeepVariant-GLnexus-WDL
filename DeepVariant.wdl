# Use DeepVariant to generate VCF & gVCF for one sample
import "tasks.wdl" as t

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
    String? deepvariant_docker = "gcr.io/deepvariant-docker/deepvariant:0.5.1"

    call t.make_examples { input:
        ref_fasta_gz = ref_fasta_gz,
        range = range,
        ranges_bed = ranges_bed,
        bam = bam,
        bai = bai,
        gvcf_gq_binsize = gvcf_gq_binsize,
        deepvariant_docker = deepvariant_docker
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

    output {
        File vcf_gz = postprocess_variants.vcf_gz
        File gvcf_gz = postprocess_variants.gvcf_gz
    }
}

