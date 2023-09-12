# This is a script to run both the GAP procedure and K-mediods on DNA fasta files in a folder
# The output is a set of cluster files that are appropriate for use in Phylobook
# Each output file contains one set of clustering results per fasta file 
# The output file name convention is "inputfilename.cluster.kmNclusters" where
# N = the number of clusters
# I run this program in RStudio and it is assumed that the
# working directory is set to the directory containing the fasta files
# so make sure to set to the desired working directory first
#Updates to this version 8/31/2023
#  1) The distance measure has been switched to use dist.dna from the ape library.  This allows
#     the tool to better handle gaps (e.g. gaps are masked out in a pairwise fashion prior to calculating distance).
#     Thus the tool no longer filters out sequences with long gaps prior calculating distances
#  2) Kmax in kmedoid clustering is now set to 1 less than the number of sequences if there are less than 10 sequences
#  3) The code to calculate "kbest" is no longer commented out as phylobook properly handles the names.  Note "kbest" is a 
#     very simple calculation as it is chosen when increasing K results in singleton clusters (as opposed to BIC or some other
#     more common method).
# Known issues: Code will hang on fasta files with only one sequence and will throw errors on fasta files with <=3 sequences.

# load libraries
library(GapProcedure)
library(adegenet)
library(cluster)
library(tidyverse)
library(Biostrings)
library(ape)
library(seqinr)

#Read List of files to be operated on
all_files <-dir()
numfiles <-length(all_files)
#For each FASTA file
#Outermost loop
for (i in 1:numfiles) {
  #read in fasta and create dnabin
  seq_data <- fasta2DNAbin(all_files[i])
  names <- dimnames(seq_data)
  nams <- unlist(names)
  numberofsequences <- length(nams)
  print("Number of sequences")
  print(numberofsequences)
  GP_out <-GapProcedure(seq_data, submod = "aK80", distance.matrix = NULL,outlier.adj = 0.9)
  Gap_classification <- GP_out[["classification"]]
  temp = cbind(nams,Gap_classification)
  outputfilename <- paste(strsplit(all_files[i], ".fasta"),".cluster.GAPclusters", sep="")
  write.table( temp,file=outputfilename, sep=",",col.names=FALSE, row.names= FALSE, quote=FALSE)
  #dist <- ak80(seq_data)
  dist <- dist.dna(seq_data, model= "K80", variance = FALSE, gamma = FALSE, pairwise.deletion = TRUE)
  # select the staring and ending K for K-mediods
  kmin <- 2
  kmax <-10
  if(kmax >numberofsequences-1){kmax <- numberofsequences-1}
  numberofsingletons = c()
  
  # For each k in K mediods
  for (k in kmin:kmax) {
    kmout <- pam(dist,k , diss=TRUE)
    kmclusters <- kmout[["clustering"]]
    
    if (k<10){
      middle <- ".cluster.km0"}
    else{
      middle <- ".cluster.km"}
    
    outputfilename <- paste(strsplit(all_files[i], ".fasta"),middle,k,"clusters", sep="")
    temp = cbind(nams, kmclusters)
    write.table( temp,file=outputfilename, sep=",",col.names=FALSE, row.names= FALSE, quote=FALSE)
    
    
    #Creat an array to store the count of the number of sequences in each cluster
    numberincluster =c()
    
    #For each cluster in the current k mediods run count the number in each cluster
    for (k2 in 1:k){
      numberincluster[k2] = (sum(kmclusters == k2))
    }
    
    print (numberincluster)
    #Count the number of clusters that are singletons	in a given k mediods run  
    #number of singletons is an array that contains the number of clusters that 
    #were a singleton for each kmediods	run. This is used to estimate the best k with 
    #the assumption that  k is best when higher k generates mostly singletons
    numberofsingletons[k-1] = (sum(numberincluster ==1))
    
  }
  # end of looping through values of k
  
  print (numberofsingletons)
  # Now look for m where the number of singletons is increasing in two consecutive k mediods runs
  d0 <- 0
  
  # Loop to process n	   
  for (m in 2:length(numberofsingletons)) {
    d <- numberofsingletons[m]-numberofsingletons[m-1]
    print (paste(d0,d))
    if((d0 > 0) & (d>0)) {
      bestk = m-1
      print (paste("bestK = ", bestk))
      break
    }
    d0 <- d
  }
  # Now rename the clustering file with the "best k" to make it first in the pulldown list
  # eg add an extra 0 to the clustername 	   
  if(bestk<10){
    middle <- ".cluster.km0"
    middle2 <- ".cluster.km00"}
  else{
    middle <- ".cluster.km"
    middle2 < ".cluster.km00"}
  filenametochange <- paste(strsplit(all_files[i], ".fasta"),middle,bestk,"clusters", sep="")
  newfilename <- paste(strsplit(all_files[i], ".fasta"),middle2,bestk,"clusters", sep="")
  file.rename(filenametochange, newfilename)
}