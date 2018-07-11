# Output description for panelLowFreq pipeline

**panelLowFreq** is a bioinformatics best-practice variant calling analysis pipeline used for WES-Seq (whole exome sequencing) or target sequencing. The pipeline focused in variant calling and annotation of candidate low frequency variants.

This document describes the output produced by the pipeline and location of output files.

## Pipeline overview:
The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [Trimmomatic](#trimming) - adapter and low quality trimming
* [BWA](#bwa) - mapping against reference genome
* [Picard](#picard) - enrichment and alignment metrics
* [SAMtools](#samtools) - alignment result processing and variant calling.
* [VarScan](#varscan) - variant calling.
* [KGGSeq](#kggseq) - variant annotation.
* [MultiQC](#multiqc) - quality statistics summary

## Preprocessing
### FastQC
Quality control is performed using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). FastQC gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.
For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Results directory**: ANALYSIS/{ANALYSIS_ID}/01-fastqc
- There is one folder per sample.
- Files:
   - `{sample_id}/{sample_id}_R[12]_fastqc.html`: html report. This file can be opened in your favourite web browser (Firefox/chrome preferable) and it contains the different graphs that fastqc calculates for QC.
   - `{sample_id}/{sample_id}_R[12]_fastqc` : folder with fastqc output in plain text.
   - `{sample_id}/{sample_id}_R[12]_fastqc.zip`: zip compression of above folder.

### Trimming
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is used for removal of adapter contamination and trimming of low quality regions. 
Parameters included for trimming are:
-  Nucleotides with phred quality < 10 in 3'end.
-  Mean phred quality < 15 in a 4 nucleotide window.
-  Read lenght < 70

MultiQC reports the percentage of bases removed by trimming in bar plot showing percentage or reads trimmed in forward and reverse.

**Note**:The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the ANALYSIS/{ANALYSIS_ID}/03-preprocQC directory.

**Results directory**: ANALYSIS/{ANALYSIS_ID}/02-preprocessing
- There is one folder per sample.
- Files:
   - `{sample_id}/{sample_id}_R[12]_filtered.fastq.gz`: contains high quality reads with both forward and reverse tags surviving.
   - `{sample_id}/{sample_id}_R[12]_unpaired.fastq.gz`: contains high quality reads with only forward or reverse tags surviving.

**NOTE:** This results are not delivered to the researcher by default due to disk space issues. If you are interesested in using them, please contact us and we will add them to your delivery.

## Mapping
### BWA
[BWA](http://bio-bwa.sourceforge.net/), or Burrows-Wheeler Aligner, is designed for mapping low-divergent sequence reads against reference genomes. The result alignment files are further processed with [SAMtools](http://samtools.sourceforge.net/), sam format is converted to bam, sorted and an index *.bai* is generated.

**Results directory**: ANALYSIS/{ANALYSIS_ID}/04-mapping.
- There is one folder per sample.
- This files can be used in [IGV](https://software.broadinstitute.org/software/igv/) for alignment visualization.
- Files:
   - `{sample_id}/{sample_id}_sorted.bam` : sorted aligned bam file.
   - `{sample_id}/{sample_id}_sorted.bam.bai`: index file for soreted aligned bam.
   
## Variant Calling
### Samtools
Samtools mpileup command is used for generate a pileup for one the BAM files. In the pileup format each line represents a genomic position, consisting of chromosome name, 1-based coordinate, reference base, the number of reads covering the site, read bases, base qualities and alignment mapping qualities. Information on match, mismatch, indel, strand, mapping quality and start and end of a read are all encoded at the read base column. This information is used by [VarScan](#varscan) for doing the proper variant calling step.

**Results directory**: ANALYSIS/{ANALYSIS_ID}/06-samtools
- There is a folder per sample.
- Files:
   - `{sample_id}/{sample_id}.pileup`: pileup format file.
   
**NOTE:** This results are not delivered to the researcher by default due to disk space issues. If you are interesested in using them, please contact us and we will add them to your delivery.

### VarScan
[VarScan](http://varscan.sourceforge.net/using-varscan.html) is used for variant calling using the command mpileup2cns with the following parameters:
- --min-var-freq 0.05: output variants with minimum 0.05 alternate allele frequency (this paramter allow the detection of low frequency variants)
- --p-value 0.99 : p-value filter is removed for posterior manual filtering.

**Results directory:**: ANALYSIS/{ANALYSIS_ID}/07-VarScan
- There is one folder per sample.
- File:
   - `{sample_id}/{sample_id}.vcf` : file with variants detected by VarScan in vcf format.

## Post-Analysis: annotation and filtering
### KGGSeq

Para el post análisis de las variantes obtenidas se utiliza el software KGGSeq (Li, Gui, Kwan, Bao, & Sham, 2012), una herramienta diseñada para priorizar variantes en el estudio de enfermedades mendelianas. En resumen, se trata de una herramienta de anotación de variantes. Permite incluir información de efecto, gen, tránscrito, enfermedades conocidas relacionadas, artículos de pubmed, etc. 

Además de las anotaciones funcionales y predictivas se realizan una serie de filtros: 
- Depth < 4
- GQ < 10.0
- PL < 20
- Sequencing quality < 50.0
- Population frequency in ANY database (ESP5400,dbsnp141,1kg201305,exac) > 0.005

Los ficheros finales anotados con variantes se pueden encontrar en ..\ANALYSIS\09-annotation/{sample_id}_all_annotated.tab
La explicación del significado de cada uno de los campos del excell se puede encontrar en: fields_description.xlsx


## Bibliography

1. *Li, M.-X., Gui, H.-S., Kwan, J. S. H., Bao, S.-Y., & Sham, P. C. (2012). A comprehensive framework for prioritizing variants in exome sequencing studies of Mendelian diseases. Nucleic acids research, 40(7), e53. doi:10.1093/nar/gkr1257*

2. *McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., … DePristo, M. a. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297–1303. doi:10.1101/gr.107524.110.20*
