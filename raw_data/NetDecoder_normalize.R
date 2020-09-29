#To use this script, basically copy and paste the code below in the R command line.
#Please, remember to adjust the paths to point to the correct locations of your files.
#You need to have the affy package installed. If you don't, type in the R command line:
#source('http://bioconductor.org/biocLite.R')
#biocLite('affy')

setwd('~/NetDecoder_Example/raw_data/CEL_files') #path to your .CEL files
source('../NetDecoder_utils.R')
library(affy)

prefix <- '../expDat'
files <- list.files('.', pattern='CEL')
raw.data <- ReadAffy(verbose = FALSE, filenames = files)
data.rma.norm <- rma(raw.data)
expRaw <- exprs(data.rma.norm)
fname <- paste(prefix, 'raw.R', sep='_')
save(expRaw, file=fname) #comment this line if you do not want to save RMA normalized data

platform <- "hgu133plus2" #change this if your microarray data was profiled in another platform
geneTab<-Norm_loadAnnotation(platform);
rownames(geneTab)<-as.vector(geneTab$probe_id)
ggs<-intersect(rownames(expRaw), as.vector(geneTab$probe_id));
expRaw2<-expRaw[ggs,];
gTab<-geneTab[ggs,];
expClean<-Norm_cleanExp(expRaw2, gTab, "symbol");
x <- sub('.cel', '', colnames(expClean), ignore.case=TRUE)
colnames(expClean) <- gsub("-", "_", x)
fname <- paste(prefix, 'clean.R', sep='_')
save(expClean, file=fname) 

#expClean is a clean expression matrix with gene symbol annotations
#this is the gene expression matrix that you have to provide to the script NetDecoder_Create_EWN.R
#where EWN means Edge-Weighted Networks.