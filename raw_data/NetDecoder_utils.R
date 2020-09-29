#code from CellNet:Network biology applied to stem cell engineering, Cahan P, et al, 2014, Cell.
#In the RStudio environment, copy and paste the functions Norm_cleanExp and Norm_loadAnnotation in the R command line.
#Then, go the line starting with ## START HERE -> .CEL files processing
Norm_cleanExp<-function
### returns the gene averaged matrix
(expDat,
 ### exp matrix
 geneTab,
 ### gene annotation
 nameCol="symbol"
 ### gene ann table column name to average over
){
  if(!is.matrix(expDat)){
    expDat<-as.matrix(expDat);
  }
  # Make sure the ann table and expDat have same rows:
  ## altered 05-03-13
  ## rownames(geneTab)<-as.vector(geneTab$probe_id);
  ## sameProbes<-intersect(rownames(expDat), as.vector(geneTab$probe_id));
  rownames(geneTab)<-as.character(geneTab$probe_id);
  sameProbes<-intersect(rownames(expDat), as.character(geneTab$probe_id));
  expDat<-expDat[sameProbes,];
  ## cat(length(sameProbes),"\n");
  geneTab<-geneTab[sameProbes,];
  ## cat(length(sameProbes),"\n");
  eids<-unique(as.vector(geneTab[,nameCol]));
  uSymbols<-vector(length=length(eids));
  ## cat(length(eids),"\n");
  ans<-matrix(nrow=length(eids), ncol=ncol(expDat));
  for(i in seq(length(eids))){
    eid<-eids[i];
    #cat(".");
    xi <- which( geneTab[,nameCol]==eid );
    ## desProbes <- as.vector(geneTab[xi,]$probe_id);
    desProbes <- as.character(geneTab[xi,]$probe_id);
    if(length(xi)>1){
      ans[i,]<- apply(expDat[desProbes,], 2, mean);
    }
    else{
      ans[i,]<-expDat[desProbes,];
    }
    uSymbols[i]<-as.vector(geneTab[ xi[1] ,nameCol]);
  }
  rownames(ans)<-uSymbols;
  colnames(ans)<-colnames(expDat);
  ans;
  ### gene averaged matrix
}

Norm_loadAnnotation<-function
### load the specified gene annotation table
(pName
 ### platform name
){
  # Need to add code to specifically download and install proper library versions.
  geneTab<-'';
  if(pName=="mogene10stv1"){
    pName<-"mogene10sttranscriptcluster";
  }
  if(pName=="hugene10stv1"){
    pName<-"hugene10sttranscriptcluster";
  }
  libName<-paste(pName, ".db",sep='');
  ## cat("Loading ", libName,"\n");
  x<-require(libName, character.only=TRUE);
  if(FALSE){
    if(!x){
      source("http://bioconductor.org/biocLite.R")
      biocLite(libName);
      x<-require(libName, character.only=TRUE);
    }
  }
  if(!x){
    cat(".loadAnnotation\tUnable to install ",libName,"\n")
    return;
  }
  else{
    entrezCmd<-paste(pName, "ENTREZID", sep='');
    idType<-'entrezgeneid';
    probes<-eval(parse(text=entrezCmd));
    probeTable<-links(probes);
    symbolsCmd<-paste(pName,"SYMBOL", sep='');
    symbols<-eval(parse(text=symbolsCmd));
    symbols<-links(symbols);
    geneTab<-merge(probeTable, symbols);
  }
  geneTab;
  ### table of probe_id, entrez id and gene symbols
}