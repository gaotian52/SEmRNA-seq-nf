#!/usr/bin/env nextflow

params.fqs = "$baseDir/test_data/**.gz"
params.transcriptome = "$baseDir/test_data/c.elegans.cdna.ncrna.fa"
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

Channel
             .fromFilePairs( params.fqs, size: -1 )
             .ifEmpty { error "Cannot find any reads matching: ${params.fqs}" }
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
process kal_mapping {
    tag "reads: $name"

    input:
    file index from transcriptome_index
    set val(name), file(fq) from read_1_ch

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
        kallisto quant --single -l ${params.fragment_len} -s ${params.fragment_sd} -b ${params.bootstrap} -i ${index} -t ${task.cpus} -o kallisto_${name} ${fq}
        """
    }
}
process quant {

        tag "${ name }"

        input:
           file index from index_ch
           set val(name), file( fq ) from read_2_ch

        output:
           file(name) into quant_ch

        """
           salmon quant -p 10 -i index -l U -r ${fq} -o ${name}
        """
    }

process fastqc {

        tag "${ name }"

        input:
            set val(name), file(fq) from read_3_ch

        output:
            file("${name}_log") into fastqc_ch

        script:
            """
            mkdir -p ${name}_log
            fastqc -o ${name}_log -f fastq -q ${fq}
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
        println ( workflow.success ? "\nDone! Open all the results --> $params.output" : "Oops .. something went wrong" )
        }
