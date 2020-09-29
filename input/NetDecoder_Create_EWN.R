
#Once you have your data normalized using the script NetDecoder_normalize.R,
#you are ready to create the edge-weighted interactomes required 
#as input for running NetDecoder. To that end, copy and paste the code below in the R command line and 

########## Code to create edge-weigthed interactomes from an iRefIndex interaction network
library(igraph)
coExpressionNetwork <- function(networkDat, expDat, genesPPI){
  aux <- data.frame()
  count <- 1
  for(edges in rownames(networkDat)){
    edge <- strsplit(edges, split = '_')[[1]]
    proteinA <- edge[1]
    proteinB <- edge[2]
    if((proteinA %in% genesPPI) & (proteinB %in% genesPPI)){
      cor <- cor.test(expDat[proteinA,], expDat[proteinB,])
      aux[edges,'proteinA'] <- proteinA
      aux[edges,'proteinB'] <- proteinB
      aux[edges, 'abs_cor'] <- abs(cor$estimate)
      aux[edges,'cor'] <- cor$estimate
      aux[edges, 'pvalue'] <- cor$p.value
      
      if((count %% 30000) == 0){cat(count, '\n')}
      count <- count + 1
    }
  }
  aux
}

path<-"~/NetDecoder_Example/input/"
setwd(path)
# Here, you can change "expBreastCancer_15Set2014.R" by the name you provided for the clean normalized expression
# matrix (expDat_clean.R, for example) you created for your data using the script NetDecoder_normalize.R.
# You will need to copy expDat_clean.R from the raw_data folder to the input folder to run this properly.

expQuery <- get(load("expBreastCancer_15Set2014.R")); #change the filename according to the filename you provided for the expClean object above

#You can see a .csv file corresponding to "stBreastCancer.R" in the folder input from the NetDecoder_Example folder.
#You can also load the corresponding .R file using the command below.
#stQuery <- get(load("stBreastCancer.R"))
stQuery <- read.csv("stBreastCancer.csv")
rownames(stQuery) <- stQuery$sample_id
network_data <- get(load("9606.mitab.04072015.R"))

network_table <- data.frame(get.edgelist(network_data$network))
rownames(network_table) <- paste(network_table$X1, network_table$X2 , sep = "_" )
genesPPI <- intersect(rownames(expQuery), network_data$all_genes)
expQuery <- expQuery[genesPPI,]

######################
## Generate co-expression networks for each disease state
for(condition in unique(stQuery$conditions)){
  cat(condition, '\n')
  sampTmp <- stQuery[which(stQuery$conditions == condition),]
  expTmp <- expQuery[,rownames(sampTmp)]
  cat('dimensions: ', dim(expTmp), '\n')
  co_exprNet <- coExpressionNetwork(network_table, expTmp, genesPPI)
  co_exprNet <- co_exprNet[complete.cases(co_exprNet),]
  filename <- paste('co_expression_network_breast_cancer_', condition, '_', Sys.Date(), '.txt', sep='')
  write.table(co_exprNet, file=filename, sep='\t', row.names=FALSE, col.names = FALSE, quote=FALSE)
}
