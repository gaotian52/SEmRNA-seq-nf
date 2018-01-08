#!/usr/bin/env nextflow

params.fqs = "/home/ywq9361/RNA-seq/sequences_check_for_each_folder/test_HS280_and_HS282.tsv"
params.transcriptome = "/home/ywq9361/RNA-seq/c.elegans.cdna.ncrna.fa"
params.outdir = "results"
params.multiqc = "$baseDir/multiqc"
File fq_file = new File(params.fqs)

log.info """\
         R N A S E Q - N F   P I P E L I N E
         ===================================
         transcriptome: ${params.transcriptome}
         fqs          : ${params.fqs}
         outdir       : ${params.outdir}

         """
         .stripIndent()


transcriptome_file = file(params.transcriptome)
multiqc_file = file(params.multiqc)

Channel.from(fq_file.collect { it.tokenize("\t")})
             .map { strain, SM, ID, NB, fq, folder, sub_folder, uniq_label -> [ strain, SM, ID, NB, file("${fq}"), folder, sub_folder, uniq_label] }
             .into { read_1_ch; read_2_ch }

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

process quant {

        tag "${ uniq_label }"

        input:
           file index from index_ch
           set strain, SM, ID, NB, fq, folder, sub_folder, uniq_label from read_1_ch

        output:
           file(uniq_label) into quant_ch

        """
           salmon quant -p 10 -i index -l U -r ${fq} -o ${uniq_label}
        """
    }

process fastqc {

        tag "${ uniq_label }"

        input:
            set strain, SM, ID, NB, fq, folder, sub_folder, uniq_label from read_2_ch

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

workflow.onComplete {
        println ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
        }
