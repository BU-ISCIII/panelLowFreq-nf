# Output description for panelLowFreq pipeline

**panelLowFreq** is a bioinformatics best-practice variant calling analysis pipeline used for WES-Seq (whole exome sequencing) or target sequencing. The pipeline focused in variant calling and annotation of candidate low frequency variants.

This document describes the output produced by the pipeline and location of output files.

## Pipeline overview:
The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [Trimmomatic](#trimmomatic) - adapter and low quality trimming
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

**Results directory**: ANALYSIS/{ANALYSIS_ID}/02-preprocessing
- There is one folder per sample.
- Files:
   - `{sample_id}/{sample_id}_R[12]_filtered.fastq.gz`: contains high quality reads with both forward and reverse tags surviving.
   - `{sample_id}/{sample_id}_R[12]_unpaired.fastq.gz`: contains high quality reads with only forward or reverse tags surviving.

**NOTE:** This results are not delivered to the researcher by default due to disk space issues. If you are interesested in using them, please contact us and we will add them to your delivery.

## Alineamiento
Se realiza el mapping con bwa. Los ficheros bam resultantes se pueden consultar en ../ANALYSIS/04-mapping. Las estadísticas del mapeo y cobertura obtenidos se pueden consultar en el mismo html report que para el control de calidad del fastq. 
Se realiza un análisis de cobertura, para ver el porcentaje de exoma cubierto a distintas profundidades con picard HsMetrics.

## Variant Calling
La llamada a variantes se realiza utilizando VarScan con el parámetro mpileup2cns, emitiendo sólo las posiciones variantes y con los siguientes parámetros --min-var-freq 0.05 --p-value 0.99, permitimos una frecuencia del alelo alternativo mínima de 0.05 y eliminamos el filtro del p-value para poder realizar los filtros luego manualmente y buscar la mayor cantidad de variantes posibles.

## Post-Análisis: Anotación y filtrado
Para el post análisis de las variantes obtenidas se utiliza el software KGGSeq (Li, Gui, Kwan, Bao, & Sham, 2012), una herramienta diseñada para priorizar variantes en el estudio de enfermedades mendelianas. En resumen, se trata de una herramienta de anotación de variantes. Permite incluir información de efecto, gen, tránscrito, enfermedades conocidas relacionadas, artículos de pubmed, etc. 
Además de las anotaciones funcionales y predictivas se realizan una serie de filtros: 
- Depth < 4
- GQ < 10.0
- PL < 20
- Sequencing quality < 50.0
- Population frequency in ANY database (ESP5400,dbsnp141,1kg201305,exac) > 0.005

Los ficheros finales anotados con variantes se pueden encontrar en ..\ANALYSIS\09-annotation/{sample_id}_all_annotated.tab
La explicación del significado de cada uno de los campos del excell se puede encontrar en: fields_description.xlsx


## Bibliografía

1. *Li, M.-X., Gui, H.-S., Kwan, J. S. H., Bao, S.-Y., & Sham, P. C. (2012). A comprehensive framework for prioritizing variants in exome sequencing studies of Mendelian diseases. Nucleic acids research, 40(7), e53. doi:10.1093/nar/gkr1257*

2. *McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., … DePristo, M. a. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297–1303. doi:10.1101/gr.107524.110.20*
