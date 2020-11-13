# SEmRNA-seq-nf

A Nextflow pipeline used for quantification and quanlity control for single-end RNA-seq data of C. elegans.

## Usage

Clone this pipeline use command below to [quest](http://www.it.northwestern.edu/research/user-services/quest/index.html):
```
git clone https://github.com/gaotian52/SEmRNA-seq-nf.git

```
and `cd SEmRNA-seq-nf`


```

nextflow kallisto_SE.nf --fqs=test.tsv

```



## Pipeliine parameters

* --fqs

We use a sample sheet as the input of sequences here, see the example in `test.tsv`.

Each column represent `strain_name` `sample_name` `raw_FASTQ`
Note that `sample_name` should be unique for each sequence.

* --ref

Reference genome. Default = [c_elegans.PRJNA13758.WS276.genomic.fa.gz](https://wormbase.org/)

* --vcf

Variant Call Format (VCF) file. Default = [WI.20200815.hard-filter.vcf.gz](https://www.elegansvariation.org/data/release/latest)


* --gff3

GFF3 file. Default = [c_elegans.PRJNA13758.WS276.annotations.gff3.gz](https://wormbase.org/)


* -out

Used to specify the output directory. Default = "RNAseq_SE-${date}" 

* --fragment_len

Estimated average fragment length. see the [kallisto document](https://pachterlab.github.io/kallisto/manual) for details. Default = "70" 

* --fragment_sd

Estimated standard deviation of fragment length. see the [kallisto document](https://pachterlab.github.io/kallisto/manual) for details. Default = "50" 

* --bootstrap

Number of bootstrap samples used by kallisto. Default = "100" 



## Output

This pipeline will generate two folders, `kallisto` and `multiqc_report` in your working directory.


`kallisto/`: RNA-Seq mapping results 


`multiqc_report/`:
```
 
├── multiqc_pre_trim_fastqc.html     # Summary of FastQC results on raw FASTQ files
├── multiqc_post_trim_fastqc.html    # Summary of FastQC results on FASTQ files trimmed by fastp
└── multiqc_kallisto.html       				 # Summary of kallisto log
```

## Dependencies

* [Nextflow](https://github.com/nextflow-io/nextflow/)

* [bcftools](https://samtools.github.io/bcftools/bcftools.html/)

* [gffread](https://github.com/gpertea/gffread/)

* [fastp](https://github.com/OpenGene/fastp/)

* [Kallisto](https://pachterlab.github.io/kallisto/)

* [FastQC](https://pachterlab.github.io/kallisto/)

* [Multiqc](https://github.com/s-andrews/FastQC/)




