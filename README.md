# DeepVariant+GLnexus workflows

These [portable WDL](http://openwdl.org/) workflows use [DeepVariant](https://github.com/google/deepvariant) to call variants from WGS read alignments, followed by [GLnexus](https://github.com/dnanexus-rnd/GLnexus) to merge the resulting Genome VCF (gVCF) files for several samples into a Project VCF (pVCF). The [`wdl/`](https://github.com/dnanexus-rnd/DeepVariant-GLnexus-workflow/tree/master/wdl) directory has three nested workflows:

### [DeepVariant.wdl](https://github.com/dnanexus-rnd/DeepVariant-GLnexus-workflow/blob/master/wdl/DeepVariant.wdl)

Based on the [DeepVariant docs](https://github.com/google/deepvariant/blob/r0.5/docs/deepvariant-gvcf-support.md), the sequential workflow to generate gVCF from a given BAM file and genomic range.

                 +----------------------------------------------------------------------------+
                 |                                                                            |
                 |  DeepVariant.wdl                                                           |
                 |                                                                            |
                 |  +-----------------+    +-----------------+    +------------------------+  |
    sample.bam   |  |                 |    |                 |    |                        |  |
     genome.fa ----->  make_examples  |---->  call_variants  |---->  postprocess_variants  |-----> gVCF
         range   |  |                 |    |                 |    |                        |  |
                 |  +-----------------+    +--------^--------+    +------------------------+  |
                 |                                  |                                         |
                 |                                  |                                         |
                 +----------------------------------|-----------------------------------------+
                                                    |
                                           DeepVariant Model
    

`make_examples` and `call_variants` internally parallelize across CPUs on the machine they run on. The tasks use the docker image [published](https://github.com/google/deepvariant/blob/r0.5/docs/deepvariant-docker.md) by the DeepVariant team.

### [htsget_DeepVariant.wdl](https://github.com/dnanexus-rnd/DeepVariant-GLnexus-workflow/blob/master/wdl/htsget_DeepVariant.wdl)

To further parallelize WGS calling accross several machines, scatters DeepVariant.wdl across several genomic ranges (typically full-length chromosomes). For each range, fetches a BAM slice using the [GA4GH htsget](https://blog.dnanexus.com/2017-10-17-introducing-htsget-a-new-ga4gh-protocol-for-genomic-data-delivery/) client in samtools 1.7+, given an [htsget server endpoint](http://samtools.github.io/hts-specs/htsget_interop.html) and sample ID. Finally, concatenates the per-range gVCFs to the  complete product.

                 +--------------------------------------------------------------------------------+
                 |                                                                                |
                 |  htsget_DeepVariant.wdl                                                        |
                 |                                                                                |
                 |       +-----------------+    +-------------------+                             |
                 |       |                 |    |                   |  range gVCF                 |
                 |   +--->  htsget client  |---->  DeepVariant.wdl  |---+                         |
                 |   |   |  (samtools)     |    |                   |   |                         |
                 |   |   |                 |    +-------------------+   |                         |
    sample ID    |   |   +-----------------+                            |  +-------------------+  |
                 |   |                                                  +-->                   |  |
       ranges -------+---> ...                  ...                 ... --->  bcftools concat  +-----> sample gVCF
        (e.g.    |   |                                                  +-->                   |  |
         chr1    |   |   +-----------------+                            |  +-------------------+  |
         chr2    |   |   |                 |    +-------------------+   |                         |
         ...)    |   +--->  htsget client  |    |                   |   |                         |
                 |       |  (samtools)     |---->  DeepVariant.wdl  |---+                         |
                 |       |                 |    |                   |  range gVCF                 |
                 |       +------------^----+    +-------------------+                             |
                 |            |       |                                                           |
                 |            |       |                                                           |
                 +------------|-------|-----------------------------------------------------------+
                              |       |
                   sample ID  |       |
                       range  |       |  range BAM
                              |       |
                         +----v------------+
                         |                 |
                         |  htsget server  |
                         |                 |
                         +-----------------+

By using htsget, the workflow scatters across the ranges without first having to download and slice up a monolithic BAM file.

### [htsget_DeepVariant_GLnexus.wdl](https://github.com/dnanexus-rnd/DeepVariant-GLnexus-workflow/blob/master/wdl/htsget_DeepVariant_GLnexus.wdl)

Scatters htsget_DeepVariant.wdl across several samples to generate an array of gVCF files, then feeds these to GLnexus to merge them into a pVCF.

                  +-----------------------------------------------------------+
                  |                                                           |
                  |  htsget_DeepVariant_GLnexus.wdl                           |
                  |                                                           |
                  |       +--------------------------+                        |
                  |       |                          |   sample gVCF          |
                  |   +--->  htsget_DeepVariant.wdl  |----+                   |
                  |   |   |                          |    |                   |
                  |   |   +--------------------------+    |    +-----------+  |
                  |   |                                   +---->           |  |
    sample IDs -------+---> ...                      ...  ----->  GLnexus  +----> project VCF
                  |   |                                   +---->           |  |
                  |   |   +--------------------------+    |    +-----------+  |
                  |   |   |                          |    |                   |
                  |   +--->  htsget_DeepVariant.wdl  |----+                   |
                  |       |                          |   sample gVCF          |
                  |       +--------------------------+                        |
                  |                                                           |
                  +-----------------------------------------------------------+
                  
Here's an example inputs JSON providing everything required to launch this top-level workflow with dxWDL or Cromwell:

```
{
    "htsget_DeepVariant_GLnexus.accessions": ["NA12878","NA12891","NA12892"],
    "htsget_DeepVariant_GLnexus.htsget_endpoint": "https://htsnexus.rnd.dnanex.us/v1/reads/BroadHiSeqX_b37",
    "htsget_DeepVariant_GLnexus.ranges": ["12:112204691-112247789","17:41196312-41277500"],
    "htsget_DeepVariant_GLnexus.ref_fasta_gz": (REFERENCE GENOME FILE),
    "htsget_DeepVariant_GLnexus.model_tar": (DEEPVARIANT MODEL FILES),
    "htsget_DeepVariant_GLnexus.output_name": "b37_CEUtrio_ALDH2_BRCA1",
}
```
