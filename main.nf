#!/usr/bin/env nextflow

params.fqs = "$baseDir/test.tsv"
params.transcriptome = "/home/ywq9361/RNA-seq/c.elegans.cdna.ncrna.fa"
params.output = "results"
params.multiqc = "$baseDir/multiqc"
params.fragment_len = '250'
params.fragment_sd = '50'
params.bootstrap = '100'
params.experiment = "$baseDir/experiment_info.txt"
File fq_file = new File(params.fqs)

log.info """\
         R N A S E Q - N F   P I P E L I N E  (Kallisto plus QC)
         ===================================
         transcriptome: ${params.transcriptome}
         fqs          : ${params.fqs}
         output       : ${params.outdir}
         fragment_len : ${params.fragment_len}
         fragment_sd  : ${params.fragment_sd}
         bootstrap    : ${params.bootstrap}
         experiment   : ${params.experiment}

         """
         .stripIndent()


transcriptome_file = file(params.transcriptome)
multiqc_file = file(params.multiqc)
exp_file = file(params.experiment)
/*
 * Make sure files exist
 */

if( !transcriptome_file.exists() ) exit 1, "Missing transcriptome file: ${transcriptome_file}"

if( !exp_file.exists() ) exit 1, "Missing Experiment parameters file: ${exp_file}"

Channel
    .from(fq_file.collect { it.tokenize("\t")})
    .map { strain, SM, ID, NB, fq, folder, sub_folder, uniq_label -> [ strain, SM, ID, NB, file("${fq}"), folder, sub_folder, uniq_label ] }
    .into { read_1_ch; read_2_ch; read_3_ch }

process qc_index {

    tag "$transcriptome_file.simpleName"

    input:
        file transcriptome from transcriptome_file

    output:
        file 'index' into index_ch

        """
        salmon index -t $transcriptome -i index
        """
}

process kal_index {

    input:
        file transcriptome_file

    output:
        file "transcriptome.index" into transcriptome_index

    script:
        //
        // Kallisto mapper index
        //
        """
        kallisto index -i transcriptome.index ${transcriptome_file}
        """
}

process kal_mapping {

    tag "reads: $SM"

    input:
        file index from transcriptome_index
        set strain, SM, ID, NB, fq, folder, sub_folder, uniq_label from read_1_ch

    output:
        file "kallisto_${SM}" into kallisto_out_dirs

    script:
    //
    // Kallisto tools mapper
    //
    def single = fq instanceof Path
    if( !single ){
        """
        mkdir kallisto_${SM}
        kallisto quant --bootstrap ${params.bootstrap} -i ${index} -t ${task.cpus} -o kallisto_${SM} ${fq}
        """
    }
    else {
        """
        mkdir kallisto_${SM}
        kallisto quant --single -l ${params.fragment_len} -s ${params.fragment_sd} -bootstrap ${params.bootstrap} -i ${index} -t ${task.cpus} -o kallisto_${SM} ${fq}
        """
    }
}

process quant {

    tag "${ SM }"

    input:
        file index from index_ch
        set strain, SM, ID, NB, fq, folder, sub_folder, uniq_label from read_2_ch

    output:
        file(ID) into quant_ch

        """
           salmon quant -p 10 -i index -l U -r ${fq} -o ${ID}
        """
}

process fastqc {

    tag "${ SM }"

    input:
        set strain, SM, ID, NB, fq, folder, sub_folder, uniq_label from read_3_ch

    output:
        file("${SM}_log") into fastqc_ch

    script:
        """
        mkdir -p ${SM}_log
        fastqc -o ${SM}_log -f fastq -q ${fq}
        """
}

process multiqc {

    input:
        file('*') from quant_ch.mix(fastqc_ch).collect()
        file(config) from multiqc_file

    output:
        file('multiqc_report.html')

    script:
        """
        cp $config/* .
        echo "custom_logo: \$PWD/logo.png" >> multiqc_config.yaml
        multiqc .
        """
}

process sleuth {

    input:
        file 'kallisto/*' from kallisto_out_dirs.collect()   
        file exp_file

    output: 
        file 'sleuth_object.so'
        file 'gene_table_results.txt'

    script:
        //
        // Setup sleuth R dependancies and environment
        //
     
        """
        sleuth.R kallisto ${exp_file}
        """
}

workflow.onComplete {
    println ( workflow.success ? "\nDone! Open all the results --> $params.output\n" : "Oops .. something went wrong" )
}
