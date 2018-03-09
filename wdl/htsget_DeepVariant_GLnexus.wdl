# Scatter the htsget_DeepVariant workflow on several sample accessions, then
# use GLnexus to merge the resulting gVCFs.
#
#               +-----------------------------------------------------------+
#               |                                                           |
#               |  htsget_DeepVariant_GLnexus.wdl                           |
#               |                                                           |
#               |       +--------------------------+                        |
#               |       |                          |   sample gVCF          |
#               |   +--->  htsget_DeepVariant.wdl  |----+                   |
#               |   |   |                          |    |                   |
#               |   |   +--------------------------+    |    +-----------+  |
#               |   |                                   +---->           |  |
# sample IDs -------+---> ...                      ...  ----->  GLnexus  +----> project VCF
#               |   |                                   +---->           |  |
#               |   |   +--------------------------+    |    +-----------+  |
#               |   |   |                          |    |                   |
#               |   +--->  htsget_DeepVariant.wdl  |----+                   |
#               |       |                          |   sample gVCF          |
#               |       +--------------------------+                        |
#               |                                                           |
#               +-----------------------------------------------------------+

import "htsget_DeepVariant.wdl" as swf

workflow htsget_DeepVariant_GLnexus {
    # sample accessions
    Array[String]+ accessions

    # genomic ranges to scatter on in htsget_DeepVariant.wdl
    Array[String]+ ranges

    # htsget settings
    String htsget_endpoint
    String? htsget_format
    File? htsget_ref_tar

    # reference genome
    File ref_fasta_gz

    # DeepVariant model files (tar with no folder component)
    File model_tar

    # DeepVariant advanced settings
    Int? gvcf_gq_binsize
    String? deepvariant_docker

    # pVCF output name
    String output_name

    scatter (accession in accessions) {
        call swf.htsget_DeepVariant as hgdv { input:
            accession = accession,
            ranges = ranges,
            ref_fasta_gz = ref_fasta_gz,
            htsget_endpoint = htsget_endpoint,
            htsget_format = htsget_format,
            htsget_ref_tar = htsget_ref_tar,
            model_tar = model_tar,
            gvcf_gq_binsize = gvcf_gq_binsize,
            deepvariant_docker = deepvariant_docker
        }
    }

    call GLnexus { input:
        gvcf = hgdv.gvcf_gz,
        ranges = ranges,
        config = "DeepVariant",
        output_name = output_name
    }

    output {
        Array[File] vcf_gz = hgdv.vcf_gz
        Array[File] gvcf_gz = hgdv.gvcf_gz
        File pvcf_gz = GLnexus.pvcf_gz
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

