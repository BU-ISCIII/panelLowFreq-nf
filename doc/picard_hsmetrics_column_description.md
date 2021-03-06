|BAIT_SET|The name of the bait set used in the hybrid selection.|
| --- | --- |
|GENOME_SIZE|The number of bases in the reference genome used for alignment.|
|BAIT_TERRITORY|The number of bases which are localized to one or more baits.|
|TARGET_TERRITORY|The unique number of target bases in the experiment, where the target sequence is usually exons etc.|
|BAIT_DESIGN_EFFICIENCY|The ratio of TARGET_TERRITORY/BAIT_TERRITORY. A value of 1 indicates a perfect design efficiency, while a valud of 0.5 indicates that half of bases within the bait region are not within the target region.|
|TOTAL_READS|The total number of reads in the SAM or BAM file examined.|
|PF_READS|The total number of reads that pass the vendor's filter.|
|PF_UNIQUE_READS|The number of PF reads that are not marked as duplicates.|
|PCT_PF_READS|The fraction of reads passing the vendor's filter, PF_READS/TOTAL_READS.|
|PCT_PF_UQ_READS|The fraction of PF_UNIQUE_READS from the TOTAL_READS, PF_UNIQUE_READS/TOTAL_READS.|
|PF_UQ_READS_ALIGNED|The number of PF_UNIQUE_READS that aligned to the reference genome with a mapping score > 0.|
|PCT_PF_UQ_READS_ALIGNED|The fraction of PF_UQ_READS_ALIGNED from the total number of PF reads.|
|PF_BASES_ALIGNED|The number of PF unique bases that are aligned to the reference genome with mapping scores > 0.|
|PF_UQ_BASES_ALIGNED|The number of bases in the PF_UQ_READS_ALIGNED reads. Accounts for clipping and gaps.|
|ON_BAIT_BASES|The number of PF_BASES_ALIGNED that are mapped to the baited regions of the genome.|
|NEAR_BAIT_BASES|The number of PF_BASES_ALIGNED that are mapped to within a fixed interval containing a baited region, but not within the baited section per se.|
|OFF_BAIT_BASES|The number of PF_BASES_ALIGNED that are mapped away from any baited region.|
|ON_TARGET_BASES|The number of PF_BASES_ALIGNED that are mapped to a targeted region of the genome.|
|PCT_SELECTED_BASES|The fraction of PF_BASES_ALIGNED located on or near a baited region (ON_BAIT_BASES + NEAR_BAIT_BASES)/PF_BASES_ALIGNED.|
|PCT_OFF_BAIT|The fraction of PF_BASES_ALIGNED that are mapped away from any baited region, OFF_BAIT_BASES/PF_BASES_ALIGNED.|
|ON_BAIT_VS_SELECTED|The fraction of bases on or near baits that are covered by baits, ON_BAIT_BASES/(ON_BAIT_BASES + NEAR_BAIT_BASES).|
|MEAN_BAIT_COVERAGE|The mean coverage of all baits in the experiment.|
|MEAN_TARGET_COVERAGE|The mean coverage of a target region.|
|MEDIAN_TARGET_COVERAGE|The median coverage of a target region.|
|MAX_TARGET_COVERAGE|The maximum coverage of reads that mapped to target regions of an experiment.|
|PCT_USABLE_BASES_ON_BAIT|The number of aligned, de-duped, on-bait bases out of the PF bases available.|
|PCT_USABLE_BASES_ON_TARGET|The number of aligned, de-duped, on-target bases out of all of the PF bases available.|
|FOLD_ENRICHMENT|The fold by which the baited region has been amplified above genomic background.|
|ZERO_CVG_TARGETS_PCT|The fraction of targets that did not reach coverage=1 over any base.|
|PCT_EXC_DUPE|The fraction of aligned bases that were filtered out because they were in reads marked as duplicates.|
|PCT_EXC_MAPQ|The fraction of aligned bases that were filtered out because they were in reads with low mapping quality.|
|PCT_EXC_BASEQ|The fraction of aligned bases that were filtered out because they were of low base quality.|
|PCT_EXC_OVERLAP|The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads.|
|PCT_EXC_OFF_TARGET|The fraction of aligned bases that were filtered out because they did not align over a target base.|
|FOLD_80_BASE_PENALTY|The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to the mean coverage level in those targets.|
|PCT_TARGET_BASES_1X|The fraction of all target bases achieving 1X or greater coverage.|
|PCT_TARGET_BASES_2X|The fraction of all target bases achieving 2X or greater coverage.|
|PCT_TARGET_BASES_10X|The fraction of all target bases achieving 10X or greater coverage.|
|PCT_TARGET_BASES_20X|The fraction of all target bases achieving 20X or greater coverage.|
|PCT_TARGET_BASES_30X|The fraction of all target bases achieving 30X or greater coverage.|
|PCT_TARGET_BASES_40X|The fraction of all target bases achieving 40X or greater coverage.|
|PCT_TARGET_BASES_50X|The fraction of all target bases achieving 50X or greater coverage.|
|PCT_TARGET_BASES_100X|The fraction of all target bases achieving 100X or greater coverage.|
|HS_LIBRARY_SIZE|The estimated number of unique molecules in the selected part of the library.|
|HS_PENALTY_10X|The "hybrid selection penalty" incurred to get 80% of target bases to 10X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 10X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 10 * HS_PENALTY_10X.|
|HS_PENALTY_20X|The "hybrid selection penalty" incurred to get 80% of target bases to 20X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 20X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 20 * HS_PENALTY_20X.|
|HS_PENALTY_30X|The "hybrid selection penalty" incurred to get 80% of target bases to 30X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 30X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 30 * HS_PENALTY_30X.|
|HS_PENALTY_40X|The "hybrid selection penalty" incurred to get 80% of target bases to 40X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 40X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 40 * HS_PENALTY_40X.|
|HS_PENALTY_50X|The "hybrid selection penalty" incurred to get 80% of target bases to 50X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 50X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 50 * HS_PENALTY_50X.|
|HS_PENALTY_100X|The "hybrid selection penalty" incurred to get 80% of target bases to 100X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 100X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 100 * HS_PENALTY_100X.|
|AT_DROPOUT|A measure of how undercovered <= 50% GC regions are relative to the mean. For each GC bin [0..50] we calculate a = % of target territory, and b = % of aligned reads aligned to these targets. AT DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total reads that should have mapped to GC<=50% regions mapped elsewhere.|
|GC_DROPOUT|A measure of how undercovered >= 50% GC regions are relative to the mean. For each GC bin [50..100] we calculate a = % of target territory, and b = % of aligned reads aligned to these targets. GC DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total reads that should have mapped to GC>=50% regions mapped elsewhere.|
|HET_SNP_SENSITIVITY|The theoretical HET SNP sensitivity.|
|HET_SNP_Q|The Phred Scaled Q Score of the theoretical HET SNP sensitivity.|
