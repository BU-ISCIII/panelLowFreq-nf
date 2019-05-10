|Column|Name|
| --- | --- |
|Chrom|chromosome name|
|Position|position (1-based)|
|Ref|reference allele at this position|
|Var|variant allele observed|
|PoolCall|Cross-sample call using all data (Cons:Cov:Reads1:Reads2:Freq:P-value)|
||Cons – consensus genotype in IUPAC format|
||Cov - total depth of coverage|
||Reads1 - number of reads supporting reference|
||Reads2 - number of reads supporting variant|
||P-value - FET p-value of observed reads vs expected non-variant|
|StrandFilt|Information to look for strand bias using all reads, format R1+:R1-:R2+:R2-:pval|
||R1+ = reference supporting reads on forward strand|
||R1- = reference supporting reads on reverse strand|
||R2+ = variant supporting reads on forward strand|
||R2- = variant supporting reads on reverse strand|
||pval = FET p-value for strand distribution, R1 versus R2|
|SamplesRef|Number of samples called reference (wildtype)|
|SamplesHet|Number of samples called heterozygous-variant|
|SamplesHom|Number of samples called homozygous-variant|
|SamplesNC|Number of samples not covered / not called|
|SampleCalls|The calls for each sample in the mpileup, space-delimited|
