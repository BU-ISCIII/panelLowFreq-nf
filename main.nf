#!/usr/bin/env nextflow

/*
========================================================================================
                          LOW  FREQUENCY  PANEL
========================================================================================
 #### Homepage / Documentation
 https://github.com/BU-ISCIII/panelLowFreq-nf
 @#### Authors
 Sara Monzon <smonzon@isciii.es>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Pipeline overview:
 - 1. : Preprocessing
 	- 1.1: FastQC for raw sequencing reads quality control
 	- 1.2: Trimmomatic
 - 2. : Mapping
 	- 2.1 : BWA alignment against reference genome
 	- 2.2 : Picard metrics for the analysis of target-capture sequencing experiments
	- 2.3 : Bedtools for calculating exons with less than 20x of depth coverage
 - 3. : Variant Calling
 	- 3.1 : Samtools for generate a pileup for the BAM files
	- 3.2 : VarScan for variant calling
 - 4. : KGGSeq 
 	- 4.1 : Post-Analysis variant annotation and filtering
 - 5. : MultiQC
 - 6. : Output Description HTML
 ----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     BU-ISCIII/panelLowFreq-nf : Low frequency panel v${version}
    =========================================
    Usage:
	
    The typical command for running the pipeline is as follows:
	
    nextflow run BU-ISCIII/panelLowFreq-nf --reads '*_R{1,2}.fastq.gz' --fasta hg38.fullAnalysisSet.fa --step preprocessing
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes).
      --fasta                       Path to Fasta reference
    References
      --bwa_index                   Path to BWA index
      --gtf							Path to GTF reference file. (Mandatory if step = assembly)
      --saveReference				Save reference file and indexes.
	Steps available:
	  --step [str]					Select which step to perform (preprocessing|mapping|variantCalling|annotFilter)
    Options:
      --singleEnd                   Specifies that the input is single end reads
    Trimming options
      --notrim                      Specifying --notrim will skip the adapter trimming step.
      --saveTrimmed                 Save the trimmed Fastq files in the the Results directory.
      --trimmomatic_adapters_file   Adapters index for adapter removal
      --trimmomatic_adapters_parameters Trimming parameters for adapters. <seed mismatches>:<palindrome clip threshold>:<simple clip threshold>. Default 2:30:10
      --trimmomatic_window_length   Window size. Default 4
      --trimmomatic_window_value    Window average quality requiered. Default 20
      --trimmomatic_mininum_length  Minimum length of reads
    Mapping options
	  --keepduplicates				Keep duplicate reads. Picard MarkDuplicates step skipped.
	  --saveAlignedIntermediates	Save intermediate bam files.
    PlasmidID options
      --plasmidid_database          Plasmids database
      --plasmidid_config            PlasmidID annotation config file
    Strain Characterization options
      --srst2_resistance            Fasta file/s for gene resistance databases
      --srst2_virulence             Fasta file/s for gene virulence databases
      --srst2_db_mlst               Fasta file of MLST alleles
      --srst2_def_mlst              ST definitions for MLST scheme
      --srst2_db_sero               Fasta file of serogroup
      --srst2_def_sero              ST definitions for serogroup scheme
    OutbreakSNP options
      --outbreaker_config			Config needed by wgs-outbreaker.
	OutbreakMLST options
    Other options:
      --outdir                      The output directory where the results will be saved
    """.stripIndent()
}