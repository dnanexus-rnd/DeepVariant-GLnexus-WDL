# Scatter the htsget_DeepVariant workflow on several sample accessions, then
# use GLnexus to merge the resulting gVCFs.
import "htsget_DeepVariant.wdl" as sub

workflow htsget_DeepVariant_GLnexus {
    Array[String]+ accessions
    Array[String]+ ranges

    scatter (accession in accessions) {
        call sub.htsget_DeepVariant as hgdv { input:
            accession = accession,
            ranges = ranges
        }
    }

    call GLnexus { input:
        gvcf = hgdv.gvcf_gz,
        ranges = ranges,
        config = "DeepVariant"
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

