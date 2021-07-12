install.packages("BiocManager",repos = "http://cran.us.r-project.org")

BiocManager::install("TFBSTools")
BiocManager::install("Biostrings")
BiocManager::install("topGO")
BiocManager::install("org.At.tair.db")
library(TFBSTools)
library(Biostrings)
library(topGO)
library("org.At.tair.db")

setwd("../")
in_matrix <- readJASPARMatrix("./promoter_annotation/Transcription_factor_weight_matrix_Arabidopsis_thaliana.txt",matrixClass="PWM")

fasta_files <- readDNAStringSet("./promoters/1kb_promoters.fa",format="fasta")


for(i in 1:length(fasta_files)) {
	site_set_list <- searchSeq(in_matrix,as.character(fasta_files[i],use.names=TRUE),seqname=names(fasta_files[i]),min.score="80%",strand = "+")
	file_path = paste("./promoter_motifs/",names(fasta_files[i]),".gff")
	write.table(writeGFF3(site_set_list),file = file_path,sep = "\t", quote = FALSE)
}

#convert each gene in promoter database into a list of it's associated go terms
possible_genes <- read.table("tf_in_matrix.txt")
pos_genes <- factor(possible_genes[,2])
names(pos_genes) <- possible_genes[,1]

GOdata <- new("topGOdata",ontology = "BP",allGenes = pos_genes,geneSel = function(x)x,annot = annFUN.org,mapping = "org.At.tair.db")
allGO <- genesInTerm(GOdata)
sapply(names(allGO),function(x) paste(x,paste(allGO[[x]],collapse="\t")))
test <- sapply(names(allGO),function(x) paste(x,paste(allGO[[x]],collapse=" ")))
lapply(test,write,"go-to-gene.map",append=TRUE)