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
      --fasta                       Path to human Fasta reference
    References
      --bwa_index                   Path to BWA index
      --gtf                         Path to GTF reference file. (Mandatory if step = assembly)
      --saveReference               Save reference file and indexes.
    Steps available:
      --step [str]                  Select which step to perform (preprocessing|mapping|variantCalling|annotFilter)
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
      --saveAlignedIntermediates    Save intermediate bam files.
    Annotation options
      --resourceDatasets            Path to resource datasets.
    Other options:
      --outdir                      The output directory where the results will be saved
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = '1.0'

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

/*
 * Default and custom value for configurable variables
 */

params.fasta = false
if( params.fasta ){
    fasta_file = file(params.fasta)
    if( !fasta_file.exists() ) exit 1, "Fasta file not found: ${params.fasta}."
}

// bwa index
params.bwa_index = false

if( params.bwa_index ){
    bwa_file = file(params.bwa_index)
    if( !fasta_file.exists() ) exit 1, "BWAIndex file not found: ${params.bwa_index}."
}

// gtf file
params.gtf = false

if( params.gtf ){
    gtf_file = file(params.gtf)
    if( !gtf_file.exists() ) exit 1, "GTF file not found: ${params.gtf}."
}

// Steps
params.step = "preprocessing"
if ( ! (params.step =~ /(preprocessing|mapping|assembly|plasmidID|outbreakSNP|outbreakMLST|strainCharacterization|mapAnnotation)/) ) {
	exit 1, 'Please provide a valid --step option [preprocessing,mapping,assembly,plasmidID,outbreakSNP,outbreakMLST,strainCharacterization,mapAnnotation]'
}

// SingleEnd option
params.singleEnd = false


// Trimming default
params.notrim = false

// Default trimming options
params.trimmomatic_adapters_file = "\$TRIMMOMATIC_PATH/adapters/NexteraPE-PE.fa"
params.trimmomatic_adapters_parameters = "2:30:10"
params.trimmomatic_window_length = "4"
params.trimmomatic_window_value = "20"
params.trimmomatic_mininum_length = "50"

// Annotation options
params.resourceDatasets = false
if( params.resourceDatasets ){
    resourceDatasets_file = file(params.resourceDatasets)
    if( !resourceDatasets_file.exists() ) exit 1, "Resource dataset file not found: ${params.resourceDatasets}."
}

// Output files options
params.saveReference = false
params.saveTrimmed = false
params.saveAlignedIntermediates = false

// Validate  mandatory inputs

if (! params.reads ) exit 1, "Missing reads: $params.reads. Specify path with --reads"

if( ! params.fasta ) exit 1, "Missing Reference human genome: '$params.fasta'. Specify path with --fasta"

if ( ! params.gtf ){
    exit 1, "GTF file not provided for assembly step, please declare it with --gtf /path/to/gtf_file"
}

if ( ! params.resourceDatasets ){
    exit 1, "Resource dataset not provided for annotation step, please declare it with --resourceDatasets /path/to/resource_datasets_file"
}


/*
 * Create channel for input files
 */

// Create channel for input reads.
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { raw_reads_fastqc; raw_reads_trimming }

// Create channel for bwa_index if supplied
if( params.bwa_index ){
    bwa_index = Channel
        .fromPath(params.bwa_index)
        .ifEmpty { exit 1, "BWA index not found: ${params.bwa_index}" }
}


// Header log info
log.info "========================================="
log.info " BU-ISCIII/panelLowFreq-nf : Low frequency panel v${version}"
log.info "========================================="
def summary = [:]
summary['Reads']               = params.reads
summary['Data Type']           = params.singleEnd ? 'Single-End' : 'Paired-End'
if(params.bwa_index)  summary['BWA Index'] = params.bwa_index
else if(params.fasta) summary['Fasta Ref'] = params.fasta
if(params.gtf)  summary['GTF File'] = params.gtf
summary['Step']                = params.step
summary['Container']           = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']        = "$HOME"
summary['Current user']        = "$USER"
summary['Current path']        = "$PWD"
summary['Working dir']         = workflow.workDir
summary['Output dir']          = params.outdir
summary['Script dir']          = workflow.projectDir
summary['Save Reference']      = params.saveReference
summary['Save Trimmed']        = params.saveTrimmed
summary['Save Intermeds']      = params.saveAlignedIntermediates
if( params.notrim ){
    summary['Trimming Step'] = 'Skipped'
} else {
    summary['Trimmomatic adapters file'] = params.trimmomatic_adapters_file
    summary['Trimmomatic adapters parameters'] = params.trimmomatic_adapters_parameters
    summary["Trimmomatic window length"] = params.trimmomatic_window_length
    summary["Trimmomatic window value"] = params.trimmomatic_window_value
    summary["Trimmomatic minimum length"] = params.trimmomatic_mininum_length
}
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "===================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}
