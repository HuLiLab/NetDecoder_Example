# NetDecoder_Example

Please refer to the Netdecoder paper <https://doi.org/10.1093/nar/gkw166> for details about the algorithms and figures.

## Prerepuisites
NetDecoder was developed and tested using Java version 1.8. The Oracle JDK, not the open JDK, was used. Therefore, werecommend to have Oracle JDK version 1.8 or higher installed. The R version used was 3.1.1 and the required R packagesare: gplots, ggplot2, grid, reshape, reshape2, plyr, RcolorBrewer, igraph and Vennerable. To install Vennerable, please usethe following command in the R command line:
```R
install.packages("Vennerable", repos="http://R-Forge.R-project.org")
```
The other packages are available through Bioconductor and can be installed using the following commands in the R command line:
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
```

## Preparing input files
### Processing of raw microarray data
You may need to install additional R packages depending on the microarray platform you are using. To run our example,you need to have installed the package
`affy`, `hgu133plus2.db` and `org.Hs.eg.db` .

Processing of raw microarray data can beachieved using the script NetDecoder_normalize.R inside the folder `raw_data` of the `NetDecoder_Example` directory. 
As anexample, we provide 5 .CEL fi les inside the folder `CEL_files` and used the RMA normalization method to process these.CEL files, 
which are from the hgu133plus2 microarray platform. Then, probeset ids mapping to the same gene symbols areaveraged and a clean normalized gene expression matrix generated and saved in the
`raw_data` directory. This procedure alsoshould work fi ne for the HuGene microarray platform. A description regarding other microarray platforms and organisms is beyond the scope of this tutorial. 
If you have your own .CEL files, you can put them in the `CEL_files` directory and run the script `NetDecoder_normalize.R` as you just did with our example .CEL files. 
You can also download microarray data from GEO and apply our script to such datasets. The script
`NetDecoder_normalize.R` is further commented to help users to run it properly.

### Processing of RNA-seq data
Netdecoder accepts normalized gene expression values from RNA-seq, such as RPKM and DESeq2-normalized values.

## Creating edge-weighted networks
Example files and a .R script to generate the input fi les from your data are provided in the folder `NetDecoder_Example/input`. 
To create the edge-weighted interactome used as input to NetDecoder, an interaction network and a normalized gene expression matrix (microarray, RNA-seq) are required. 
The script `NetDecoder_Create_EWN.R`, where EWN means Edge-Weighted Networks will help you to create the edge-weighted interactomes given an expression matrix.

First, load the file `NetDecoder_Create_EWN.R` into Rstudio and if you saved our NetDecoder example in your home directory, you are all set. Otherwise, you need to change the following line to reflect your local directory structure:
```R
path<-"~/NetDecoder_Example/input/"
```
This script takes as input three files: a `NxM` gene expression matrix (expBreastCancer_15Set2014.R) 
where `N` is the number of genes and `M` the number of samples, a sample annotation table (stBreastCancer.csv)
and a .R object containing the network (as an igraph object) and a list of all genes in the network (9606.mitab.04072015.R). 
Importantly, when creating your own sample annotation table, make sure to include a column called `conditions`.
This column contains the annotation for each sample such as normal_tissue, ERnegative or ERpositive. 
An example table is also provided inside the folder input (stBreastCancer.csv). 
Then, run this script (`source` function) a co-expression network for each phenotype (control, ER-negative and ER-positive) 
will be generated. This co-expression network consists of the absolute value of the Pearson-correlation coefficient (PCC) 
of each interaction reported in the protein interaction network derived from iRefIndex.

This computation may take several hours (we estimate about 3-5 hours, depending on the number of phenotypes/conditions studied)
because it requires to iterate over all interactome edges, compute a PCC for each reported interaction and assign the value to the respective edge.
In our example, once it is done, it will generate three files in the directory specified by the path above, one for each phenotype under analysis. 
You can move these files to the folder `breast_cancer` such that you can perform the analysis using the newly generated edge-weighted interactomes. 
In our previous example, a co-expression network for the control transcriptome was called
`co_expression_network_breast_cancer_control_2015-07-09.txt`
,inside the folder `breast_cancer`. The other input parameter is a gene list to be used as sources. 
We leave up to the users which gene list to use such as differentially expressed genes or any other set of genes of interest. 
Example gene lists can be foundin the folder breast_cancer, for example, the file
`sig_genes_control_ERnegative_2015-07-15.txt`
contains the genes identified as preferentially expressed in ER-negative breast cancer when compared to normal breast tissues. 
It is a .txt file with a column header, such as “genes” in our example gene lists. 
You may want to change the following line to reflect the phenotype you are interested in to study using NetDecoder, changing
`breast_cancer` to a meaningful name associated to your experiment or condition.
```R
filename <- paste('co_expression_network_breast_cancer_', condition, '_', Sys.Date(), '.txt', sep='')
```
Once you have your edge-weighted interaction networks for each phenotype of interest as well as your gene lists, 
you can change the scripts provided in the `A reproducible example` section to point to your input files. 
Importantly, NetDecoder always compute at least two subnetworks from the same set of source genes. In our breast cancer example, 
we used the genes differentially expressed between control and ER-negative breast cancer as sources for both the control and 
ER-negative edge-weighted interactomes. This is requested such that NetDecoder can evaluate path differences as a 
function of the transcriptome of each phenotype.

## Analysis: A reproducible example
NetDecoder takes as input an edge-weighted interactome and a gene set to be used as sources. Transcription-related genes 
are defined as targets (sinks) by default. A reproducible example is provided in this repository, 
which includes the iRefIndex network used in our paper. Download this repository into your home directory, or any other directory you want, and go to the
`NetDecoder_Example` folder.
```bash
git clone https://github.com/HuLiLab/NetDecoder_Example/
cd NetDecoder_Example
```

The first step is to compute the subnetworks associated to each phenotype of interest. We will use ER-negative and ER-positive breast cancers as an example 
with the respective normal tissue transcriptome profile. In this scenario, 4 subnetworks will the generated: two for the normal tissue transcriptome data (two different sets of gene signatures are used) 
and another two for ER-negative and ER-positive breast cancer transcriptomes. To run this step, go to the folder
`networks` and edit the script `NetDecoder_CreateNetworks.sh`. If you saved the example files in a directory other than your home directory, please change the paths accordingly 
to point to the correct locations of the files in this tutorial. For example, the lines bellow should be updated in the scripts.
```bash
LIB_DIR=~/NetDecoder_Example/lib
FILES_DIR=~/NetDecoder_Example/lib
OUT_DIR=~/NetDecoder_Example/breast_cancer/networks/
BASE_DIR=~/NetDecoder_Example/breast_cancer/networks/
```
Then, open a terminal and go to the folder `NetDecoder_Example/breast_cancer/networks/` , which contains the script
`NetDecoder_CreateNetworks.sh`. With the paths pointing to the correct locations, run the script in the terminal:
```bash
./NetDecoder_CreateNetworks.sh
```
This process could take a while, usually about 3-4 hours. The phenotype-specific subnetworks will be saved in the current directory `networks`
and will be are used as input for our downstream analysis to find interactions (edges), network routers, key targets and high impact genes, 
all of which are predicted to play a role in breast cancer pathogenesis.

In the next step, we are going to perform a typical NetDecoder analysis as reported in the Figures 2, 3 and 4 in the Netdecoder paper. 
If you saved the NetDecoder_Example in your home directory no additional changes should be required. Otherwise, change the paths in the script
`NetDecoder_Analysis.sh` to point to the correct locations for the input files. To run this analysis, enter directory `breast_cancer/analysis`,
and run 
```bash
./NetDecoder_Analysis.sh
```
This will generate several files and plots in two different folders, corresponding to ER-negative and ER-positive breastcancers. 
The first ones that you might be interested to look at are the plots in .pdf format. For example, we may look at the interactions with higher flows 
in ER-negative or ER-positive and these might provide clues about interactions that might be important for further experimental validation. 
Hereafter, we are going to call these networks as control and disease-subnetworks for ER-negative and ER-positive breast cancers. 

In addition to these results, this script also generates several important files, many of them contain the raw data from where the plots are generated. 
For instance, the files ending with `_keyEdges` contain the plot and the data used to generate the barchart for gene interactions shown in Figure 1a.
The file ending with `_Jaccard_ERnegative` contains the raw data togenerate the heatmap in Figure 1b and the files ending with
`_Network routers.pdf` (raw data in the file ending with `_HIDDEN.txt`) and `_Key targets`(raw data ending with `_SINKS.txt`) 
contains the Network routers and Key targets in Figure 1c. Gml files containing `_NR_` and `_KT_` as part of the filenames 
identify the network motifs for the respective genes in `_Control` and `_Disease` subnetworks. Bar charts ending with
`edge_flows_Control` or `edge_flows_Disease` show histograms for edge flow values in the respective network motif. 
For an example of these network motifs, please see Figure2. In addition, there are several temporary files that are used to 
generate the network motifs and in general, you do not needto worry about them. These files contain `_flowDifference_X`
or `_totalFlow_X` as part of the filename, where “X” is either `Control` or `Disease`. These files can be used to give an idea 
about the values behind the color coding used in the network motifs. Last, this script also generates a .gml file for the full network
`_FULL_` , the edge-centered networks `_EDGE_CENTERED_SUBNET_` and the prioritized subnetworks `_PRIORITIZED_NETWORK_`

By default, we also generate differential networks by comparing paths between two phenotype-specific subnetworks and generating
a subnetwork from paths that are specific to each phenotype. Therefore, a differential network is generated for phenotype 1 and phenotype 2 
with the overlap between paths in both subnetworks. These files are identified by `_DFN_` in the filename.

This analysis also computed the impact scores for all genes in both phenotype 1 and phenotype 2 subnetworks and found the genes called
high impact genes in our paper. A heatmap containing the genes with the highest or lowest impact scoresis generated in the folder
`analysis`.

Among the .gml files previously mentioned for network routers and key targets motifs, the folder ER-negative and ER-positive also 
contain network motifs for the high impact genes. High impact genes are identified by a `_IP_` in their filenames. 
The .pdf files contain histograms for edge flow values as described for network routers and key targets above.

**That's it! You successfully ran a NetDecoder analysis and have plenty of data to figure out what is happening in your experiment.**

## NetDecoder parameters
### Parameters for required input files
No changes are needed in these parameters. NetDecoder uses these files internally.
`-SYMBOL`: Gene association file. Used to perform mappings between genes and gene ontology terms.
`-GO`: The Gene Ontology database used by NetDecoder.

### Parameters to generate subnetworks and perform downstream analysis
`-gen`: It is used in the script `NetDecoder_CreateNetworks.sh`. It is used to generate the subnetworks connecting sources to sinks.

`-E`: It is used to run the analysis to obtain the edges, network routers and key targets described in our paper, Figure 2. 
      It also generate .gml files for edge-centered networks and network motifs for network routers and key targets. These .gml files can be
      imported into commonly used network visualization softwares, such as Cytoscape or Gephi. Plots for the distribution of edge flows are also provided.

`-C`: It is used to perform the analysis of high impact genes, as in the Figure 3 of our paper. Also generate network motifs as.gml files that could be imported into Cytoscape or Gephi.

`-cCCS`: It is used for plotting a heatmap containing the high impact genes, as in the Figure 3 of our paper.

`-ncp`: file containing the paths associated to the phenotype 1 subnetwork, such as control.

`-ndp`: file containing the paths associated to the phenotype 2 subnetwork, such as ER-negative.

`-control`: a string specifying the phenotype 1 state, such as control.

`-condition`: a string specifying the phenotype 2 state, such as ER-negative.

`-corThreshold`: this parameter selects all edges with fl ow higher than corThreshold. The default value is 0.5.

`-ratioThreshold`: this parameter is used to further select edges with high fl ow differences between phenotype 2 andphenotype 1 subnetworks. The default value is 5.

`-top`: it is used to select the top edges contributing the most to better distinguish phenotype 1 and phenotype 2 subnetworks.
      It is an additional filtering step to prioritize edges with higher flow values under phenotype 2 (often disease) state. 
      It also used to select the top intermediary and target genes with the highest total node flow differences (network routers and keytargets, respectively) 
      between phenotype 2 and phenotype 1 subnetworks. Importantly, the parameters

`-corThreshold` , `-ratioThreshold` and `-top` are used in combination to find out edges with higher flow differences between phenotype 2 andphenotype 1 subnetworks.

`-overlap`: when this parameter is provided, only edges shared between both subnetworks will be included in the analysis. 
      Otherwise, all edges in both subnetworks will be included. The NetDecoder example above illustrates a situation where `-overlap` is not provided.

`-g`: it is used to specify an input gene list. In `NetDecoder_CreateNetworks.sh`, `-g` is used to indicate the list of genes to beused as sources. 
      In `NetDecoder_Analysis.sh` or in `NetDecoder_Edges.sh` .It is used to provide a gene list from which theinteracting partners for each gene will be plotted using an adjacency matrix representation. 
      If a given gene is present in bothphenotype 1 and phenotype 2 subnetworks, the establishment of new interactions can observed, as illustrated in Figure 4.

`-out`: the output directory where the results will be saved.

`-f`: filename for the output


