#!/usr/bin/env nextflow

params.fqs = "$baseDir/test.tsv"
params.transcriptome = "/home/ywq9361/RNA-seq/c.elegans.cdna.ncrna.fa"
params.outdir = "results"
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
         outdir       : ${params.outdir}
         fragment_len : ${params.fragment_len}
         fragment_sd  : ${params.fragment_sd}
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

Channel.from(fq_file.collect { it.tokenize("\t")})
             .map { strain, SM, ID, NB, fq, folder, sub_folder, uniq_label -> [ strain, SM, ID, NB, file("${fq}"), folder, sub_folder, uniq_label ] }
             .into { read_1_ch; read_2_ch; read_3_ch }

process index {
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
process mapping {
    tag "reads: $name"

    input:
    file index from transcriptome_index
    set strain, SM, ID, NB, fq, folder, sub_folder, uniq_label from read_1_ch

    output:
    file "kallisto_${name}" into kallisto_out_dirs

    script:
    //
    // Kallisto tools mapper
    //
    def single = fq instanceof Path
    if( !single ){
        """
        mkdir kallisto_${name}
        kallisto quant -b ${params.bootstrap} -i ${index} -t ${task.cpus} -o kallisto_${name} ${fq}
        """
    }
    else {
        """
        mkdir kallisto_${name}
        kallisto quant --single -l ${params.fragment_len} -s ${params.fragment_sd} -b ${params.bootstrap} -i ${index} -t ${task.cpus} -o 
        kallisto_${name} ${fq}
        """
    }
}
process quant {

        tag "${ uniq_label }"

        input:
           file index from index_ch
           set strain, SM, ID, NB, fq, folder, sub_folder, uniq_label from read_2_ch

        output:
           file(uniq_label) into quant_ch

        """
           salmon quant -p 10 -i index -l U -r ${fq} -o ${uniq_label}
        """
    }

process fastqc {

        tag "${ uniq_label }"

        input:
            set strain, SM, ID, NB, fq, folder, sub_folder, uniq_label from read_3_ch

        output:
            file("${uniq_label}_log") into fastqc_ch

        script:
            """
            mkdir -p ${uniq_label}_log
            fastqc -o ${uniq_label}_log -f fastq -q ${fq}
            """
        }

process multiqc {
        publishDir params.outdir, mode:'copy'

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
        println ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
        }
