# ClusteringForPhylobook
This is a repository of R code that clusters nucleotide sequence data using the GAP procedure and K-mediods clustering.  The input is a folder containing fastA files for all the data within a project.  The output is a collection of files containing the sequence-to-cluster-ID for each clustering run (e.g. the GAP procedure, and various K's in K-mediod clustering).
This R code is designed to run from within R or R studio.  It expects that the working directory will be set to a folder that contains only fasta files.
It will process each fasta file within the directory to generate cluster assignments using the GAP procedure (as publicly available at - https://github.com/vrbiki/GapProcedure) or K-mediods clustering.  The output is a collection of csv files of format "sequence_name, clusterID" - one file for each type of clustering (e.g. GAP procedures and values of K from 2 to KMax).  KMax is adjustable by editing the code as needed.  The purpose of this program is to provide guidance to users of Phylobook during manual selection/editing of lineages.  See - https://github.com/MullinsLab/phylobook for more details.
KmedoidsGapsHandled is the latest release of the code.  Updates to this version include using a different method the calculat the distances.  In this method pairwise gaps are masked out.  Also an HIV specific substitution matrix is now used.  Finally, the code to estimate "kbest" e.g. the k at which kmediods is likely the best choice has been implemented.  The selection of k-best is simple in this code - it's the k after which new clusters contain only one member.
