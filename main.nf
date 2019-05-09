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
     - 2.2 : Samtools
     - 2.3 : Picard metrics for the analysis of target-capture sequencing experiments
     - 2.4 : Samtools for generate a pileup for the depud BAM files
 - 3. : Variant Calling
    - 3.1 : VarScan for variant calling
 - 4. : KGGSeq 
     - 4.1 : Post-Analysis variant annotation and filtering
     - 4.2 : R for table merging
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
    
    nextflow run BU-ISCIII/panelLowFreq-nf --reads '*_R{1,2}.fastq.gz' --fasta hg38.fullAnalysisSet.fa
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes).
      --fasta                       Path to human Fasta reference
    References
      --bwa_index                   Path to BWA index
      --saveReference               Save reference file and indexes.
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
      --keepduplicates                Keep duplicate reads. Picard MarkDuplicates step skipped.
    Veriant calling
      --maxDepth                    Maximum number of reads per input file to read at a position. Default 20000
      --minBaseQ                    Minimum base quality for a base to be considered. Default 0.
      --minVarFreq                  Minimum variant allele frequency threshold. Default 0.05
      --pValue                      Default p-value threshold for calling variants. Default 0.99
    Annotation options
      --resourceDatasets            Path to resource datasets.
      --saveDuplicates              Save duplicated sequences.
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


// SingleEnd option
params.singleEnd = false

//Mapping-duplicates defaults
params.keepduplicates = false

// Trimming default
params.notrim = false

// Default trimming options
params.trimmomatic_adapters_file = "\$TRIMMOMATIC_PATH/adapters/NexteraPE-PE.fa"
params.trimmomatic_adapters_parameters = "2:30:10"
params.trimmomatic_window_length = "4"
params.trimmomatic_window_value = "20"
params.trimmomatic_mininum_length = "50"

//Default variant calling options
params.maxDepth = "20000"
params.minBaseQ = "0"
params.minVarFreq = "0.05"
params.pValue = "0.99"

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
params.saveDuplicates = false

// Validate  mandatory inputs

if (! params.reads ) exit 1, "Missing reads: $params.reads. Specify path with --reads"

if( ! params.fasta ) exit 1, "Missing Reference human genome: '$params.fasta'. Specify path with --fasta"

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
summary['Keep Duplicates'] = params.keepduplicates
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
summary['Save Duplicates']     = params.saveDuplicatesbed
summary['VarScan p-value threshold'] = params.pValue
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

/*
 * Build BWA index
 */
if(!params.bwa_index && fasta_file){
    process makeBWAindex {
        tag "${fasta.baseName}"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta_file

        output:
        file "${fasta}*" into bwa_index

        script:
        """
        mkdir BWAIndex
        bwa index -a bwtsw $fasta
        """
    }
}


/*
 * STEP 1.1 - FastQC
 */


process fastqc {
    tag "$prefix"
    publishDir "${params.outdir}/01-fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results
    file '.command.out' into fastqc_stdout

    script:

    prefix = name - ~/(_S[0-9]{2})?(_L00[1-9])?(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(_00*)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    fastqc --nogroup -t 1 -k 8 $reads
    """
}


/*
 * STEP 1.2 - Trimming
 */

process trimming {
    tag "$prefix"
    publishDir "${params.outdir}/02-trimming", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf(".log") > 0) "logs/$filename"
   else if (filename.indexOf(".fastq.gz") > 0) "trimmed/$filename"
            else params.saveTrimmed ? filename : null
    }

    input:
    set val(name), file(reads) from raw_reads_trimming

    output:
    file '*_paired_*.fastq.gz' into trimmed_paired_reads,trimmed_paired_reads_bwa
    file '*_unpaired_*.fastq.gz' into trimmed_unpaired_reads
    file '*_fastqc.{zip,html}' into trimmomatic_fastqc_reports
    file '*.log' into trimmomatic_results

    script:
    prefix = name - ~/(_S[0-9]{2})?(_L00[1-9])?(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(_00*)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    trimmomatic PE -phred33 $reads -threads 1 $prefix"_paired_R1.fastq" $prefix"_unpaired_R1.fastq" $prefix"_paired_R2.fastq" $prefix"_unpaired_R2.fastq" ILLUMINACLIP:${params.trimmomatic_adapters_file}:${params.trimmomatic_adapters_parameters} SLIDINGWINDOW:${params.trimmomatic_window_length}:${params.trimmomatic_window_value} MINLEN:${params.trimmomatic_mininum_length} 2> ${name}.log
    gzip *.fastq
    fastqc -q *_paired_*.fastq.gz
    """
}

/*
 * STEP 2.1 - BWA alignment
 */

process bwa {
    tag "$prefix"
    publishDir path: { params.saveAlignedIntermediates ? "${params.outdir}/03-bwa" : params.outdir }, mode: 'copy',
            saveAs: {filename -> params.saveAlignedIntermediates ? filename : null }

    input:
    file reads from trimmed_paired_reads_bwa
    file index from bwa_index

    output:
    file '*.bam' into bwa_bam

    script:
    prefix = reads[0].toString() - ~/(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    bwa mem -M $index $reads > ${prefix}.sam
    samtools view -bS ${prefix}.sam -o ${prefix}.bam
    """
}


/*
 * STEP 2.2 - Samtools post-alignment processing
 */

process samtools {
    tag "${bam.baseName}"
    publishDir path: { params.saveAlignedIntermediates ? "${params.outdir}/03-bwa" : params.outdir }, mode: 'copy',
            saveAs: { filename ->
                    if (filename.indexOf("_stats.txt") > 0) "stats/$filename"
                    else params.saveAlignedIntermediates ? filename : null
            }

    input:
    file bam from bwa_bam

    output:
        file '*_sorted.bam' into bam_for_mapped, bam_picard, bam_samtolls
        file '*_sorted.bam.bai' into bwa_bai, bai_picard,bai_for_mapped
        file '*_stats.txt' into samtools_stats

    script:
    """
    samtools sort $bam -T ${bam.baseName}_sorted -o ${bam.baseName}_sorted.bam
    samtools index ${bam.baseName}_sorted.bam
    samtools stats ${bam.baseName}_sorted.bam > ${bam.baseName}_stats.txt
    """
}


/*
 * STEP 2.3 - Picard
 */
if (!params.keepduplicates){
    process picard {
        tag "$prefix"
        publishDir "${params.outdir}/04-picard", mode: 'copy'
            
            
        input:
        file bam from bam_picard

        output:
        file '*_dedup_sorted.bam' into bam_dedup_spp, bam_dedup_ngsplot, bam_dedup_deepTools, bam_dedup_macs, bam_dedup_saturation, bam_dedup_epic, bam_dedup_mpileup
        file '*_dedup_sorted.bam.bai' into bai_dedup_deepTools, bai_dedup_spp, bai_dedup_ngsplot, bai_dedup_macs, bai_dedup_saturation, bai_dedup_epic
        file '*_picardDupMetrics.txt' into picard_reports

        script:
        prefix = bam[0].toString() - ~/(_sorted)?(\.bam)?$/
        """
        picard MarkDuplicates \\
            INPUT=$bam \\
            OUTPUT=${prefix}_dedup.bam \\
            ASSUME_SORTED=true \\
            REMOVE_DUPLICATES=true \\
            METRICS_FILE=${prefix}_picardDupMetrics.txt \\
            VALIDATION_STRINGENCY=LENIENT \\
            PROGRAM_RECORD_ID='null'
        samtools sort ${prefix}_dedup.bam -o ${prefix}_dedup_sorted.bam -T ${prefix}
        samtools index ${prefix}_dedup_sorted.bam
        """
    }
/*
 * STEP 2.4 - Samtools pileup picard
 */

    process mpileup {
        tag "$prefix"
        publishDir "${params.outdir}/05-samtools_NC", mode: 'copy'

        input:
        file dedup_bam from bam_dedup_mpileup
        file fasta from fasta_file

        output:
        file '*.pileup' into pileup_results

        script:
        prefix = dedup_bam[0].toString() - ~/(_dedup)?(\.sorted)?(\.bam)?$/
        """
        samtools mpileup -A -d ${params.maxDepth} -Q ${params.minBaseQ} -f $fasta $dedup_bam > $prefix".pileup"
        """
    }
    
/*
 * STEP 2.4 - Samtools pileup keepduplicates
 */

} else {
    process mpileup {
        tag "$prefix"
        publishDir "${params.outdir}/05-samtools_NC", mode: 'copy'

        input:
        file dup_bam from bam_samtolls
        file fasta from fasta_file

        output:
        file '*.pileup' into pileup_results

        script:
        prefix = dup_bam[0].toString() - ~/(_dedup)?(\.sorted)?(\.bam)?$/
        """
        samtools mpileup -A -d ${params.maxDepth} -Q ${params.minBaseQ} -f $fasta $dup_bam > $prefix".pileup"
        """
    }

}



/*
 * STEP 3.1 - VarScan
 */

process varscan {
    tag "${pileup.baseName}"
    publishDir "${params.outdir}/06-VarScan", mode: 'copy'

    input:
    file pileup from pileup_results

    output:
    file '*.vcf' into vcf_file

    script:
    """
    varscan mpileup2cns $pileup --min-var-freq ${params.minVarFreq} --p-value ${params.pValue} --variants --output-vcf 1 > ${pileup.baseName}.vcf
    """
}

/*
 * STEP 4.1 - KGGSeq
 */

process varscan {
    tag "${vcf.baseName}"
    publishDir "${params.outdir}/07-Annotation", mode: 'copy',
            saveAs: { filename ->
                    if (filename.indexOf(".table") > 0) "bcftools/$filename"
                    else if (filename.indexOf("_all_annotated.tab") > 0) "R_merge/$filename"
                    else if (filename.indexOf("_annot.txt*") > 0) "kggseq/$filename"
                    else if (filename.indexOf("_header.table") > 0) "header_table/$filename"
            }

    input:
    file vcf from vcf_file

    output:
    file '*.table' into bcftools_tables
    file '*_all_annotated.tab' into r_merged_tables
    file '*_annot.txt.flt.txt' into kggseq_flt_file
    file '*_annot.txt.log' into kggseq_annot_log
    file '*_header.table' into header_table

    script:
    """
    bcftools query -H $vcf -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT\t%DP\t%RD\t%AD\t%FREQ\t%PVAL\t%RBQ\t%ABQ\t%RDF\t%RDR\t%ADF\t%ADR]\n' > ${vcf.baseName}.table
    kggseq --no-resource-check --no-lib-check --buildver hg38 --vcf-file $vcf --db-gene refgene --db-score dbnsfp --genome-annot --db-filter ESP5400,dbsnp141,1kg201305,exac --rare-allele-freq 1 --mendel-causing-predict best --omim-annot --out ${vcf.baseName}_annot.txt
    gunzip *_annot.txt.flt.txt.gz
    cp header ${vcf.baseName}_header.table && tail -n +2 ${vcf.baseName}.table >> ${vcf.baseName}_header.table
    """
}