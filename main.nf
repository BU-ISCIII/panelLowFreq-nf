#!/usr/bin/env nextflow

/*
========================================================================================
                          LOW  FREQUENCY  PANEL
========================================================================================
 #### Homepage / Documentation
 https://github.com/BU-ISCIII/panelLowFreq-nf
 @#### Authors
 Sarai Varona <s.varona@isciii.es>
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
 - 5. : Stats
     - 5.2 : Bamstats
     - 5.3 : Picard CalculateHsMetrics
     - 5.1 : MultiQC
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
    Statistics:
      --bamstatsTargets             Path to capture targets bed file
      --picardstatsTargets          Path to capture header targets file
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
    if( !fasta_file.exists() ) exit 1, "Fasta reference file not found: ${params.fasta}."
}


// SingleEnd option
params.singleEnd = false

//Mapping-duplicates defaults
params.keepduplicates = true
params.indexFiles = false

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

params.multiqc_config = "${baseDir}/conf/multiqc_config.yaml"

if (params.multiqc_config){
    multiqc_config = file(params.multiqc_config)
}

// Output files options
params.saveReference = false
params.saveTrimmed = false
params.saveAlignedIntermediates = false

// Stats options
params.bamstatsTargets = false
if( params.bamstatsTargets ){
    bamstatsTargets_file = file(params.bamstatsTargets)
    if( !bamstatsTargets_file.exists() ) exit 1, "bamstats Targets file not found: ${params.bamstatsTargets}."
}

params.picardstatsTargets = false
if( params.picardstatsTargets ){
    picardstatsTargets_file = file(params.picardstatsTargets)
    if( !picardstatsTargets_file.exists() ) exit 1, "bamstats Targets file not found: ${params.picardstatsTargets}."
}



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

// Create channel for fasta reference if supplied
if( params.fasta ){
    Channel
        .fromPath(params.fasta)
        .ifEmpty { exit 1, "Fasta reference not found: ${params.fasta}" }
        .into { fasta_file; fasta_bwamem; fasta_file_pileup_picard; fasta_file_pileup; fasta_bwa_index }
    //Channel
         //.fromPath("${params.fasta}.*")
         //.into { bwa_index_path }
}

// Create channel for picard stat targets
if( params.picardstatsTargets ){
    Channel
        .fromPath(params.picardstatsTargets)
        .ifEmpty { exit 1, "Picard stats file not found: ${params.picardstatsTargets}" }
        .into { picardstatsTargets_file; picardstatsTargets_file_picard }
}

// Create channel for bamstats stat targets
if( params.bamstatsTargets ){
    Channel
        .fromPath(params.bamstatsTargets)
        .ifEmpty { exit 1, "Picard stats file not found: ${params.bamstatsTargets}" }
        .into { bamstatsTargets_file; bamstatsTargets_file_picard }
}

// Create channel for resource Datasets
if( params.resourceDatasets ){
    Channel
        .fromPath(params.resourceDatasets)
        .ifEmpty { exit 1, "resource Datasets file not found: ${params.resourceDatasets}" }
        .into { resourceDatasets_file_picard; resourceDatasets_file }
}

//Create multiQC config chanel
if (params.multiqc_config) {
    Channel
        .fromPath(params.multiqc_config, checkIfExists: true)
        .set { ch_config_for_multiqc }
}

// Create channel for reference index files
if( params.indexFiles ){
    Channel
        .fromPath("${params.fasta}.*")
		.into { bwa_index; bwa_index_pileup_picard; bwa_index_pileup }
        //if( !bwa_index.exists() ) exit 1, "Index files not found: ${params.indexFiles}."
}

// Header log info
log.info "========================================="
log.info " BU-ISCIII/panelLowFreq-nf : Low frequency panel v${version}"
log.info "========================================="
def summary = [:]
summary['Reads']               = params.reads
summary['Data Type']           = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Fasta Ref']            = params.fasta
summary['Keep Duplicates'] = params.keepduplicates
summary['Container']           = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']        = "$HOME"
summary['Current user']        = "$USER"
summary['Current path']        = "$PWD"
summary['Index files']         = params.indexFiles
summary['Working dir']         = workflow.workDir
summary['Output dir']          = params.outdir
summary['Script dir']          = workflow.projectDir
summary['Save Reference']      = params.saveReference
summary['Save Trimmed']        = params.saveTrimmed
summary['Save Intermeds']      = params.saveAlignedIntermediates
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

if ( !params.indexFiles ){
/*
 * Build BWA index
 */
    process makeBWAindex {
        tag "${fasta.baseName}"
        publishDir path: { params.saveReference ? "${params.outdir}/../REFERENCES" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta_file

        output:
        file "${fasta}*" into bwa_index, bwa_index_pileup_picard, bwa_index_pileup

        script:
        """
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
    file '*_fastqc.{zip,html}' into fastqc_results, fastqc_results_picard
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
    file '*_unpaired_*.fastq.gz' into trimmed_unpaired_reads, trimmed_unpaired_reads_picard
    file '*_fastqc.{zip,html}' into trimmomatic_fastqc_reports, trimmomatic_fastqc_reports_picard
    file '*.log' into trimmomatic_results, trimmomatic_results_picard

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
    file fasta from fasta_bwamem
    file bwa_index from bwa_index.collect()

    output:
    file '*.bam' into bwa_bam

    script:
    prefix = reads[0].toString() - ~/(.R1)?(_1)?(_R1)?(_trimmed)?(_paired)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    bwa mem -M $fasta $reads | samtools view -bT $fasta - > ${prefix}.bam
    """
}

/*
 * STEP 2.2 - Samtools post-alignment processing
 */

process samtools {
    tag "${bam.baseName}"
    publishDir path: "${params.outdir}/03-bwa", mode: 'copy',
            saveAs: { filename ->
                    if (filename.indexOf("_stats.txt") > 0) "stats/$filename"
                    else params.saveAlignedIntermediates ? filename : null
            }

    input:
    file bam from bwa_bam

    output:
    file '*_sorted.bam' into bam_for_mapped, bam_picard, bam_samtolls, bam_stats, picard_stats
    file '*_sorted.bam.bai' into bwa_bai, bai_picard,bai_samtools, bai_bamstats, bai_picard_stats
    file '*_stats.txt' into samtools_stats

    script:
    """
    samtools sort $bam -T ${bam.baseName}_sorted -o ${bam.baseName}_sorted.bam
    samtools index ${bam.baseName}_sorted.bam
    samtools stats ${bam.baseName}_sorted.bam > ${bam.baseName}_stats.txt
    """
}


if (!params.keepduplicates){
/*
 * STEP 2.3 - Picard
 */
    process picard {
        tag "$prefix"
        publishDir "${params.outdir}/04-picard", mode: 'copy'
            
        input:
        file bam from bam_picard
        file bai from bai_picard

        output:
        file '*_dedup_sorted.bam' into bam_dedup_mpileup, dedup_bam_stats, dedup_picard_stats
        file '*_dedup_sorted.bam.bai' into bai_dedup_stats, bai_dedup_picard_stats, bai_dedup_mpileup
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

    process mpileup_picard {
        tag "$prefix"
        publishDir "${params.outdir}/05-samtools_NC", mode: 'copy'

        input:
        file dedup_bam from bam_dedup_mpileup
        file bai_file from bai_dedup_mpileup
        file fasta from fasta_file_pileup_picard
		file bwa_index from bwa_index_pileup_picard.collect()

        output:
        file '*.pileup' into pileup_results_picard

        script:
        prefix = dedup_bam[0].toString() - ~/(_dedup)?(\.sorted)?(\.bam)?$/
        """
        samtools mpileup -A -d ${params.maxDepth} -Q ${params.minBaseQ} -f $fasta $dedup_bam > ${prefix}.pileup
        """
    }
    
/*
 * STEP 3.1 - VarScan picard
 */

    process varscan_picard {
        tag "${pileup.baseName}"
        publishDir "${params.outdir}/06-VarScan", mode: 'copy'

        input:
        file pileup from pileup_results_picard

        output:
        file '*.vcf' into vcf_file_picard

        script:
        """
        varscan mpileup2cns $pileup --min-var-freq ${params.minVarFreq} --p-value ${params.pValue} --variants --output-vcf 1 > ${pileup.baseName}.vcf
        """
    }

/*
 * STEP 4.1 - KGGSeq picard
 */

    process kggseq_picard {
        tag "${vcf.baseName}"
        publishDir "${params.outdir}/07-Annotation", mode: 'copy'

        input:
        file vcf from vcf_file_picard
        file resource from resourceDatasets_file_picard

        output:
        file '*.table' into bcftools_tables_picard
        file '*_annot.txt.flt.txt' into kggseq_flt_file_picard
        file '*_annot.txt.log' into kggseq_annot_log_picard
        file '*_header.table' into header_table_picard

        script:
        """
        bcftools query -H $vcf -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT\t%DP\t%RD\t%AD\t%FREQ\t%PVAL\t%RBQ\t%ABQ\t%RDF\t%RDR\t%ADF\t%ADR]\n' > ${vcf.baseName}.table
        kggseq --no-lib-check --buildver hg38 --vcf-file $vcf --resource $resource --db-gene refgene --db-score dbnsfp --genome-annot --db-filter ESP5400,dbsnp141,1kg201305,exac --rare-allele-freq 1 --mendel-causing-predict best --omim-annot --out ${vcf.baseName}_annot.txt
        gunzip *_annot.txt.flt.txt.gz
        cp ${baseDir}/assets/header ${vcf.baseName}_header.table && tail -n +2 ${vcf.baseName}.table >> ${vcf.baseName}_header.table
        """
    }

/*
 * STEP 4.2 - R merge picard
 */

    process rmerge_picard {
        tag "$prefix"
        publishDir "${params.outdir}/07-Annotation/R_merge", mode: 'copy'

        input:
        file header_table from header_table_picard
        file flt from kggseq_flt_file_picard

        output:
        file '*_all_annotated.tab' into r_merged_tables_picard

        script:
        prefix = header_table.baseName - ~/(_header)?(\.table)?$/
        """
        Rscript $baseDir/bin/merge_parse.R $prefix
        """
    }

/*
 * STEP 5.2 - Bamstats picard
 */    
     process bamstats_picard {
        tag "$prefix"
        publishDir "${params.outdir}/08-stats/bamstats", mode: 'copy'

        input:
        file region_list from bamstatsTargets_file_picard
        file sorted_bam from dedup_bam_stats
        file bai_file from bai_dedup_stats


        output:
        file '*_bamstat.txt' into bamstats_result_picard

        script:
        prefix = sorted_bam[0].toString() - ~/(_paired)?(_sorted)?(\.bam)?$/

        """
        bam stats --regionList $region_list --in $sorted_bam --baseSum --basic 2> ${prefix}_bamstat.txt
        """
    }

/*
 * STEP 5.2 - Picard CalculateHsMetrics picard
 */
 
     process picardmetrics_picard {
        tag "$prefix"
        publishDir "${params.outdir}/08-stats/picardmetrics", mode: 'copy'

        input:
        file picard_targer from picardstatsTargets_file_picard
        file sorted_bam from dedup_picard_stats
        file bai_file from bai_dedup_picard_stats


        output:
        file '*_hsMetrics.out' into picardstats_result_picard
        file 'hsMetrics_all.out' into picardstats_all_result_picard

        script:
        prefix = sorted_bam[0].toString() - ~/(_sorted)?(\.bam)(_paired)(_dedup)?$/
        """
        picard CalculateHsMetrics BI=$picard_targer TI=$picard_targer I=$sorted_bam O=${prefix}_hsMetrics.out VALIDATION_STRINGENCY='LENIENT'
        echo "SAMPLE","MEAN TARGET COVERAGE", "PCT USABLE BASES ON TARGET","FOLD ENRICHMENT","PCT TARGET BASES 10X","PCT TARGET BASES 20X","PCT TARGET BASES 30X","PCT TARGET BASES 40X","PCT TARGET BASES 50X" > hsMetrics_all.out
        grep '^RB' ${prefix}_hsMetrics.out | awk 'BEGIN{FS="\t";OFS=","}{print "${prefix}",\$22,\$24,\$25,\$29,\$30,\$31,\$32,\$33}' >> hsMetrics_all.out
        """
    }


/*
 * STEP 5.1 - MultiQC picard
 */

    process multiqc_picard {
        tag "$prefix"
        publishDir "${params.outdir}/08-stats/multiQC", mode: 'copy'

        input:
        file multiqc_config
        file (fastqc:'fastqc/*') from fastqc_results_picard.collect()
        file ('trimommatic/*') from trimmomatic_results_picard.collect()
        file ('trimommatic/*') from trimmomatic_fastqc_reports_picard.collect()
        file ('bamstats/*') from bamstats_result_picard.collect()
        file ('picardstats/*') from picardstats_result_picard.collect()

        output:
        file '*multiqc_report.html' into multiqc_report_picard
        file '*_data' into multiqc_data_picard
        file '.command.err' into multiqc_stderr_picard
        val prefix into multiqc_prefix_picard

        script:
        prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'

        """
        multiqc -d . --config $multiqc_config
        """
    }

}


if (params.keepduplicates){

/*
 * STEP 2.4 - Samtools pileup keepduplicates
 */

    process mpileup {
        tag "$prefix"
        publishDir "${params.outdir}/05-samtools_NC", mode: 'copy'

        input:
        file dup_bam from bam_samtolls
        file bai_file from bai_samtools
        file fasta from fasta_file_pileup
		file bwa_index from bwa_index_pileup.collect()

        output:
        file '*.pileup' into pileup_results

        script:
        prefix = dup_bam[0].toString() - ~/(_dedup)?(\_sorted)?(\.bam)?$/
        """
        samtools mpileup -A -d ${params.maxDepth} -Q ${params.minBaseQ} -f $fasta $dup_bam > ${prefix}.pileup
        """
    }

/*
 * STEP 3.1 - VarScan keepduplicates
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

    process kggseq {
        tag "${vcf.baseName}"
        publishDir "${params.outdir}/07-Annotation", mode: 'copy'

        input:
        file vcf from vcf_file
        file resource from resourceDatasets_file

        output:
        file '*.table' into bcftools_tables
        file '*_annot.txt.flt.txt' into kggseq_flt_file
        file '*_annot.txt.log' into kggseq_annot_log
        file '*_header.table' into header_table

        script:
        """
        bcftools query -H $vcf -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT\t%DP\t%RD\t%AD\t%FREQ\t%PVAL\t%RBQ\t%ABQ\t%RDF\t%RDR\t%ADF\t%ADR]\n' > ${vcf.baseName}.table
        kggseq --no-lib-check --buildver hg38 --vcf-file $vcf --resource $resource --db-gene refgene --db-score dbnsfp --genome-annot --db-filter ESP5400,dbsnp141,1kg201305,exac --rare-allele-freq 1 --mendel-causing-predict best --omim-annot --out ${vcf.baseName}_annot.txt
        gunzip *_annot.txt.flt.txt.gz
        cp ${baseDir}/assets/header ${vcf.baseName}_header.table && tail -n +2 ${vcf.baseName}.table >> ${vcf.baseName}_header.table
        """
    }

/*
 * STEP 4.2 - R merge
 */

    process rmerge {
        tag "$prefix"
        publishDir "${params.outdir}/07-Annotation/R_merge", mode: 'copy'

        input:
        file header_table from header_table
        file flt from kggseq_flt_file

        output:
        file '*_all_annotated.tab' into r_merged_tables

        script:
        prefix = header_table.baseName - ~/(_header)?(\.table)?$/
        """
        Rscript ${baseDir}/bin/merge_parse.R $prefix
        """
    }

/*
 * STEP 5.2 - Bamstats
 */
     process samstats {
        tag "$prefix"
        publishDir "${params.outdir}/08-stats/bamstats", mode: 'copy'

        input:
        file region_list from bamstatsTargets_file
        file sorted_bam from bam_stats
        file bai_file from bai_bamstats


        output:
        file '*_bamstat.txt' into bamstats_result

        script:
        prefix = sorted_bam[0].toString() - ~/(_paired)?(_sorted)?(\.bam)?$/

        """
        bam stats --regionList $region_list --in $sorted_bam --baseSum --basic > ${prefix}_bamstat.txt
        """
    }

/*
 * STEP 5.2 - Picard CalculateHsMetrics
 */
 
     process picardmetrics {
        tag "$prefix"
        publishDir "${params.outdir}/08-stats/picardmetrics", mode: 'copy'

        input:
        file picard_targer from picardstatsTargets_file
        file sorted_bam from picard_stats
        file bai_file from bai_picard_stats


        output:
        file '*_hsMetrics.out' into picardstats_result
        file 'hsMetrics_all.out' into picardstats_all_result

        script:
        prefix = sorted_bam[0].toString() - ~/(_sorted)?(\.bam)(_paired)(_dedup)?$/

        """
        picard CalculateHsMetrics BI=$picard_targer TI=$picard_targer I=$sorted_bam O=${prefix}_hsMetrics.out VALIDATION_STRINGENCY='LENIENT'
        echo "SAMPLE","MEAN TARGET COVERAGE", "PCT USABLE BASES ON TARGET","FOLD ENRICHMENT","PCT TARGET BASES 10X","PCT TARGET BASES 20X","PCT TARGET BASES 30X","PCT TARGET BASES 40X","PCT TARGET BASES 50X" > hsMetrics_all.out
        grep '^RB' ${prefix}_hsMetrics.out | awk 'BEGIN{FS="\t";OFS=","}{print "${prefix}",\$22,\$24,\$25,\$29,\$30,\$31,\$32,\$33}' >> hsMetrics_all.out
        """
    }
 
/*
 * STEP 5.1 - MultiQC
 */

    process multiqc {
        tag "$prefix"
        publishDir "${params.outdir}/08-stats/multiQC", mode: 'copy'

        input:
        file multiqc_config
        file (fastqc:'fastqc/*') from fastqc_results.collect()
        file ('trimommatic/*') from trimmomatic_results.collect()
        file ('trimommatic/*') from trimmomatic_fastqc_reports.collect()
        file ('bamstats/*') from bamstats_result.collect()
        file ('picardstats/*') from picardstats_result.collect()

        output:
        file '*multiqc_report.html' into multiqc_report
        file '*_data' into multiqc_data
        file '.command.err' into multiqc_stderr
        val prefix into multiqc_prefix

        script:
        prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'

        """
        multiqc -d . --config $multiqc_config
        """
    }
}

workflow.onComplete {
    log.info "BU-ISCIII - Pipeline complete"
}