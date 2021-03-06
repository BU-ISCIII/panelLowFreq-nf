|Column|Meaning|
| --- | --- |
|Chromosome|chromosome number|
|StartPosition|Human genome reference position|
|ReferenceAlternativeAllele|reference/alternative allele|
|rsID|SNP rs ID|
|MostImportantFeatureGene|Gene Symbol|
|MostImportantGeneFeature|Gene feature {missense,intronic, ncRNA, etc}|
|RefGeneFeatures|Gene Features {codons,transcripts,etc}|
|SLR|Sitewise Likelihood-ratio (SLR) test statistic for testing natural selection on codons. A negative value indicates negative selection, and a positive value indicates positive selection. Larger magnitude of the value suggests stronger evidence.|
|SIFT_score|SIFT uses the 'Sorting Tolerant From Intolerant' (SIFT) algorithm to predict whether a single amino acid substitution affects protein function or not, based on the assumption that important amino acids in a protein sequence should be conserved throughout evolution and substitutions at highly conserved sites are expected to affect protein function.A small scoreindicates a high chance for a substitutionto damage the protein function.|
|Polyphen2_HDIV_score|"Polyphen2 score based on HumDiv, i.e. hdiv_prob. The score ranges from 0 to 1, and the corresponding prediction is ""probably damaging"" if it is in [0.957,1]; ""possibly damaging"" if it is in [0.453,0.956]; ""benign"" if it is in [0,0.452]. Score cutoff for binary classification is 0.5, i.e. the prediction is ""neutral"" if the score is smaller than 0.5 and ""deleterious"" if the score is larger than 0.5. Multiple entries separated by "";"""|
|Polyphen2_HVAR_score|Polyphen2 predicts the possible impact of an amino acid substitution on the structure and function of a human protein using straightforward physical and comparative considerations by an iterative greedy algorithm. In the present study, we use the original scores generated by the HumVar (instead ofHumDiv) trained model as it is preferred for the diagnosis of Mendelian diseases. The scores range from 0 to 1. A substitution with larger score has a higher possibility to damage the protein function.|
|LRT_score|LRT employed a likelihood ratio test to assess variant deleteriousnessbased on a comparative genomics data set of 32 vertebrate species. The identified deleterious mutations could disrupt highly conserved amino acids within protein-coding sequences, which are likely to be unconditionally deleterious.The scores range from 0 to 1. A larger score indicates a larger deleterious effect.|
|MutationTaster_score|MutationTaster assesses the impact of the disease-causing potential of a sequence variant by a naive Bayes classifier using multiple resources such as evolutionary conservation, splice-site changes, loss of protein features and changes that might affect mRNA level. The scores range from 0 to 1. The larger score suggests a higher probability to cause a human disease.|
|MutationAssessor_score|"MutationAssessor ""functional impact of a variant : predicted functional (high, medium), predicted non-functional (low, neutral)"" Please refer to Reva et al. Nucl. Acids Res. (2011) 39(17):e118 for details"|
|FATHMM_score|"FATHMM default score (weighted for human inherited-disease mutations with Disease Ontology); If a score is smaller than -1.5 the corresponding NS is predicted as ""D(AMAGING)""; otherwise it is predicted as ""T(OLERATED)"". If there's more than one scores associated with the same NS due to isoforms, the smallest score (most damaging) was used. Please refer to Shihab et al Hum. Mut. (2013) 34(1):57-65 for details"|
|VEST3|VEST 3.0 score. Score ranges from 0 to 1. The larger the score the more likely the mutation may cause functional change. In case there are multiple scores for the same variant, the largest score (most damaging) is presented. Please refer to Carter et al., (2013) BMC Genomics. 14(3) 1-16 for details. Please note this score is free for non-commercial use. For more details please refer to http://wiki.chasmsoftware.org/index.php/SoftwareLicense. Commercial users should contact the Johns Hopkins Technology Transfer office.|
|CADD_score|Combined Annotation Dependent Depletion (CADD) score for funtional prediction of a SNP. Please refer to Kircher et al. (2014) Nature Genetics 46(3):310-5  for details. The larger the score the more likely the SNP has damaging effect.|
|GERP++_NR|Neutral rate|
|GERP++_RS|RS score, the larger the score, the more conserved the site|
|phyloP|PhyloP estimates the evolutional conservation at each variant from multiple alignments of placental mammal genomes to the human genome based on a phylogenetic hidden Markov model.|
|29way_logOdds|SiPhy score based on 29 mammals genomes. The larger the score, the more conserved the site.|
|LRT_Omega|Estimated nonsynonymous-to-synonymous-rate ratio ( reported by LRT)|
|AffectedRefHomGtyNum|Number of affected individuals with reference homozygote at this variant;|
|AffectedHetGtyNum|Number of affected individuals with heterozygote at this variant;|
|AffectedAltHomGtyNum|Number of affected individuals with non-ref homozygote;|
|UnaffectedRefHomGtyNum|Number of unaffected individuals with reference homozygote at this variant;|
|UnaffectedHetGtyNum|Number of unaffected individuals with heterozygote at this variant;|
|UnaffectedAltHomGtyNum|Number of unaffected individuals with non-ref homozygote;|
|DenovoMutationEvent|"n the main output file, there is a column named DenovoMutationEvent to record the genotypes of a child and his or her parents. Example: N140_0:0/1:46,59&N140_1:0/0:57,0&N140_2:0/0:68,0. The child N140_0 has genotype 0/1 with 46 and 59 reads carrying reference alleles and alternative alleles respectively. The father N140_1 and mother N140_2 are homozygous 0/0."|
|UniProtFeatureForRefGene|Annotate a variant of coding gene using the UniProt protein annotations.|
|GeneDescription|Gene description|
|Pseudogenes|Pseudogenes listed in http://tables.pseudogene.org/set.py?id=Human61|
|DiseaseName(s)MIMid|"Disorder, <disorder MIM no.> (<phene mapping key>) Phenotype mapping method <phene mapping key>:1 - the disorder is placed on the map based on its association with a gene, but the underlying defect is not known.2 - the disorder has been placed on the map by linkage; no mutation has been found. 3 - the molecular basis for the disorder is known; a mutation has been found in the gene.4 - a contiguous gene deletion or duplication syndrome, multiple genes are deleted or duplicated causing the phenotype."|
|GeneMIMid|GeneMIMid : Gene/locus MIM no.|
|SIFT_pred|SIFT prediction filter|
|Polyphen2_HDIV_pred|"Polyphen2 prediction based on HumDiv, ""D"" (""porobably damaging""), ""P"" (""possibly damaging"") and ""B"" (""benign""). Multiple entries separated by "";"""|
|Polyphen2_HVAR_pred|"Polyphen2 prediction based on HumVar, ""D"" (""porobably damaging""), ""P"" (""possibly damaging"") and ""B"" (""benign""). Multiple entries separated by "";"""|
|LRT_pred|Classification using LRT (D = deleterious, N = neutral, or U = unknown)|
|MutationTaster_pred|Classification using MutationTaster (A = disease_causing_automatic, D = disease_causing, N = polymorphism, or P = polymorphism_automatic)|
|MutationAssessor_pred|"MutationAssessor ""functional impact of a variant : predicted functional (high, medium), predicted non-functional (low, neutral)"""|
|FATHMM_pred|FATHMM prediction filter.|
|DiseaseCausalProb_ExoVarTrainedModel|Conditional probability of being Mendelian disease-causing given the above prediction scores under a logistic regression model trained by our dataset ExoVar.|
|IsRareDiseaseCausal_ExoVarTrainedModel|Classification using the logistic regression model (Y = disease-causing or N = neutral)|
|BestCombinedTools:OptimalCutoff:TP:TN|The subset of original prediction tools (out of the 13 tools) used for the combined prediction by our Logistic Regression model which have the largest posterior probability among all possible combinatorial subsets: the cutoff leads to the maximal Matthews correlation coefficient (MCC): the corresponding true positive and true negative at the maximal MCC.|
|TFBSconsSite[tfbsName:rawScore:zScore]|Conserved TFBSs in the UCSC genome browser|
|vistaEnhancer[enhancerName:positive/negative]|Known enhancers in the VISTA enhancer browser|
|PubMedIDIdeogram|PubMed ID of articles in which the term and the cytogeneic position of the variant are co-mentioned|
|PubMedIDGene|PubMed ID of articles in which the term and the gene containing the variant are co-mentioned|
