import "htsget_DeepVariant_GLnexus.wdl" as main

workflow DVGLx_test {  
  call main.htsget_DeepVariant_GLnexus
  call check_outputs {
    input:
      pvcf_gz = htsget_DeepVariant_GLnexus.pvcf_gz
  }
}

task check_outputs {
  File pvcf_gz

  command {
    set -ex -o pipefail
    actual=$(zcat "${pvcf_gz}" | grep -v \# | cut -f1 | uniq -c | tr -d ' ' | tr '\n' ,)
    expected="19chr12,185chr17,"
    if [ "$actual" != "$expected" ]; then
      exit 1
    fi
  }
}

