# Integration_Site_Analysis
Snakemake pipeline to analyze the integration patterns of the synthetic transposon Sleeping Beauty in relation to several genomic features. The pipeline utilizes several custom python scripts and applies a Hidden Markov Model. Due to the sensitivity of the research, only the pipeline and scripts I have written are shown here.

Snakemake was chosen for speed, reproducibility, and ease-of-use for users. Parallelization and multi-threading ensures that the pipeline remains fast even as the number of inputs increaess. Reproducibility is ensured through a log file that tracks changes to the pipeline and through a yaml file that can be exported. Since snakemake is pythonic, it is relatively easy for users to understand. 

FastQC v0.11.9 provided initial quality control metrics, such as the presence of overrepresented sequences indicative of adapter contamination, distribution of read lengths, and more. Using Cutadapt v3.4, the SB sequences from the 5’ end of the read and any Illumina primers from the 3’ end were removed. For its computational speed, Bowtie2 v2.2.5 was utilized to the align the trimmed reads to the Mus musculus genome assembly GRCm38 (mm10). Due to the multiple PCR steps during library preparation, the frequency of nucleotide misincorporation was high enough to affect the removal of the adapters. To combat this, unmapped reads were trimmed again with less stringent criteria and mapped again. Strings of A’s and T’s were also removed. The mapped reads from each trimming-mapping step were combined for downstream analyses irrespective of mapping quality. The omission of a filter by mapping quality permitted the incorporation of multi-mapping reads increasing the sensitivity of the pipeline at the risk of false positives.

To further increase the accuracy, reads mapping to the same position were clustered together employing a custom python script and a consensus sequence was generated. Misalignment introduced by improper trimming were accounted for by clustering reads that were a few bases away. Instead of alignment score, the read with highest number of occurrence within the cluster was taken to have the correct starting position and insert length. The user can modify the thresholds for the read proximity for clustering and for the minimum number of occurrences needed to be reported. Sorting of the reads by chromosome position was performed by Samtools v1.12. Because it is known that SB integrates with a preference in TA dinucleotide repeats, an additional script was written to attempt to correct the position for reads that do not begin with TA. The script added or subtracted bases from the 5’ end until it locates a TA with a preference for a new read length closer to the theoretical trimmed sequencing insert length. Sequencing errors were suppressed by taking a consensus of the sequences using Bedtools v2.30.0. The consensus reads were filtered by length to only retain those close to the theoretical trimmed sequencing insert length. The genomic coordinates of the SB landing sites were stored in a Browser Extensible Data file as the “final” output of the pipeline. Additional steps can be added to the snakemake pipeline depending on the desired analysis.
