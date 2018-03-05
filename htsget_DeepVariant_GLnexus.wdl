# Scatter the htsget_DeepVariant workflow on several sample accessions, then
# use GLnexus to merge the resulting gVCFs.
import "htsget_DeepVariant.wdl" as sub
import "tasks.wdl" as t

workflow htsget_DeepVariant_GLnexus {
    Array[String]+ accessions
    Array[String]+ ranges

    scatter (accession in accessions) {
        call sub.htsget_DeepVariant as hgdv { input:
            accession = accession,
            ranges = ranges
        }
    }

    call t.GLnexus { input:
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

