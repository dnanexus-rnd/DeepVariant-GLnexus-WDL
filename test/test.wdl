import "htsget_DeepVariant_GLnexus.wdl" as main

workflow DVGLx_test {  
  Array[String]+ accessions
  Array[String]+ ranges
  String htsget_endpoint
  String? htsget_format
  File ref_fasta_gz
  File model_tar
  String deepvariant_docker
  String output_name
  String expected

  call main.htsget_DeepVariant_GLnexus {
    input:
      accessions = accessions,
      ranges = ranges,
      htsget_endpoint = htsget_endpoint,
      htsget_format = htsget_format,
      ref_fasta_gz = ref_fasta_gz,
      model_tar = model_tar,
      deepvariant_docker = deepvariant_docker,
      output_name = output_name
  }
  call check_outputs {
    input:
      pvcf_gz = htsget_DeepVariant_GLnexus.pvcf_gz,
      expected = expected
  }
}

task check_outputs {
  File pvcf_gz
  String expected

  command {
    set -ex -o pipefail
    actual=$(zcat "${pvcf_gz}" | grep -v \# | cut -f1 | uniq -c | tr -d ' ' | tr '\n' ,)
    if [ "$actual" != "${expected}" ]; then
      exit 1
    fi
  }
}

