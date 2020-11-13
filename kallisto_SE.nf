#!/usr/bin/env nextflow



/*
~ ~ ~ > * USER INPUT PARAMETERS 
*/


params.fqs                    = null
params.ref                    = "/projects/b1059/data/genomes/c_elegans/WS276/c_elegans.PRJNA13758.WS276.genomic.fa.gz"
params.vcf                    = "/projects/b1059/analysis/WI-20200815/variation/WI.20200815.hard-filter.vcf.gz"
params.gff3                   = "/projects/b1059/projects/Gaotian/resource/WS276/c_elegans.PRJNA13758.WS276.annotations.gff3.gz"
params.fragment_len           = '70'
params.fragment_sd            = '50'
params.bootstrap              = '100'



/*
~ ~ ~ > * OUTPUT DIRECTORY 
*/


date = new Date().format( 'yyyyMMdd' )

params.out = "RNAseq_SE-${date}"



/*
~ ~ ~ > * INITIATE INPUT Fastq channel 
*/


File fq_file = new File(params.fqs)

fileinfo = Channel
    .from(fq_file.collect { it.tokenize("\t") })
    .map { strain, SM, reads -> [ strain, SM, file("${reads}") ] }
    .into { fq_file_1; fq_file_2 }

Channel
    .from(fq_file.collect { it.tokenize("\t") })
    .map { ST, SM, fq -> ST  }
    .unique()
    .set {Strains}




/*
~ ~ ~ > * INITIATE reference genome channel 
*/



File reference = new File("${params.ref}")
reference_handle = reference.getAbsolutePath()



/*
~ ~ ~ > * INITIATE Variant VCF file channel 
*/



File vcf = new File("${params.vcf}")
vcf_handle = vcf.getAbsolutePath()


/*
~ ~ ~ > * INITIATE GFF3 file channel 
*/


File gff3 = new File("${params.gff3}")
gff3_handle = gff3.getAbsolutePath()







/*
   ==================================
    Trim raw fastq files
   ==================================
*/ 


process fastp_Trim {

   

    tag "reads: ${strain}_${SM}"

    

    cpus 8

    input:
       set val(strain), val(SM), file(reads) from fq_file_1

    
    output:

       set val(strain), val(SM), file ("${SM}_fastp.fq.gz") into fq_file_trim1, fq_file_trim2

       



    script:
    //
    // fastp
    //
  
        """
        fastp --thread 8 --length_required 20 -i ${reads} -o ${SM}_fastp.fq.gz
    
        """
    }








/*
   ==================================
    Generate Strain specific Transcriptome file
   ==================================
*/ 







process generate_specific_trancriptome {

   
 
    tag "${ST}"
    
    input:

      val(ST) from Strains
      
    output:
      
      set val(ST), file("${ST}.transcriptome2.fa") into sample_transcriptome 
 


    """

        bcftools consensus -f ${reference_handle} -i "TYPE='snp'" --sample ${ST} ${vcf_handle} > ${ST}.fa

        gunzip -c ${gff3_handle} | grep "WormBase\\|##" > ce_gff.gff

        
        gffread -w ${ST}.transcriptome.fa -g ${ST}.fa ce_gff.gff       


        cat ${ST}.transcriptome.fa | sed 's/ CDS=.*//' | sed 's/CDS://' | sed 's/Transcript://' | sed 's/Pseudogene://' > ${ST}.transcriptome2.fa
       
    """

}


/* 
   ==================================
   Kallisto index
   ==================================
*/ 


process kal_index {
    
    input:
        set val(ST), file(transcriptome_file) from sample_transcriptome


    output:
        set val(ST), file ("${ST}.transcriptome.index") into transcriptome_index

    script:
        //
        // Kallisto mapper index
        //
        """
        kallisto index -i ${ST}.transcriptome.index ${transcriptome_file}
        
        """
}


fq_file_trim1
	.combine(transcriptome_index, by: 0)
	.into{mapping_data_set;
		  print_mapping}




/* 
   ==================================
   Kallisto mappping
   ==================================
*/ 



process kal_mapping {

    publishDir "${params.out}/kallisto", mode: 'copy' , pattern: "kallisto_*"


    tag "reads: ${SM}"

    cpus 8

    input:
       set val(strain), val(SM), file(reads), file(kalIndex) from mapping_data_set
       


    output:
    
       set val(strain), val(SM), file ("kallisto_${SM}") into kallisto_out_dirs

       file("*_log.txt") into kallisto_log


    script:
    //
    // Kallisto tools mapper
    //
  
        """
        mkdir kallisto_${SM}

        kallisto quant --single \\
        -l ${params.fragment_len} \\
        -s ${params.fragment_sd} \\
        --bootstrap ${params.bootstrap} \\
        -i ${kalIndex} \\
        -o kallisto_${SM} \\
        --threads=8 \\
        ${reads}  &> ${SM}_kallisto_log.txt
        
        
    
        """
    }










/* 
   ==================================
                FastQC
   ==================================
*/ 

    


process pre_trim_fastqc {

 
    tag "${SM}"

    cpus 4

    input:
        set ST, SM, file(reads) from fq_file_2

    output:


       file("${SM}_prelog") into prefastqc_ch

    script:

        """

        
        mkdir -p ${SM}_prelog

        fastqc -t 4 -o ${SM}_prelog -f fastq -q ${reads}

       
        """
}


process post_trim_fastqc {

 
    tag "${SM}"

    cpus 4

    input:
        set ST, SM, file(reads) from fq_file_trim2

    output:


        file("${SM}_postlog") into postfastqc_ch

    script:

        """

        
        mkdir -p ${SM}_postlog

        fastqc -t 4 -o ${SM}_postlog -f fastq -q ${reads}

       
        """
}





/* 
   ==================================
                MultiQC
   ==================================
*/ 

    



prefastqc_ch
    .mix(postfastqc_ch)
    .mix(kallisto_log)
    .collect()
    .into{qc_data;
    qc_data_2;
    qc_data_3}






process summary_multi_qc {


    publishDir "${params.out}/multiqc_report", mode: 'copy'

    cpus 4
    memory '32 GB'

    input:

        file('*') from qc_data_3


    output:

         set file("multiqc_post_trim_fastqc.html"), file("multiqc_pre_trim_fastqc.html"), file("multiqc_kallisto.html"),file("multiqc_*/*")

    script:
    
        """



        multiqc *_postlog/*_fastqc.zip --filename multiqc_post_trim_fastqc.html --interactive

        multiqc *_prelog/*_fastqc.zip --filename multiqc_pre_trim_fastqc.html --interactive

        multiqc *_kallisto_log.txt --filename multiqc_kallisto.html --interactive


        """
}
