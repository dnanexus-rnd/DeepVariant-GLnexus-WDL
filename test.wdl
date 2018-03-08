import "htsget_DeepVariant_GLnexus.wdl" as main

workflow DVGLx_test {  
  call main.htsget_DeepVariant_GLnexus
  call check_outputs {
    input:
      pvcf = htsget_DeepVariant_GLnexus.pvcf
  }
}

task check_outputs {
  File pvcf

  command {
    set -ex -o pipefail
    [ -s "${pvcf}" ] || exit 1
  }
}

