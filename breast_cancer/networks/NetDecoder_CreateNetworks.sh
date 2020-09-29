#!/bin/sh
#
# run NetDecoder, run!
#

LIB_DIR=~/NetDecoder_Example/lib
FILES_DIR=~/NetDecoder_Example/lib
LIB=$LIB_DIR/NetDecoder.jar:$LIB_DIR/commons-cli-1.2.jar
SYMBOL=$FILES_DIR/gene_associationgoa_ref_HUMAN_25Jun2014.txt
GO=$FILES_DIR/gene_ontology_ext_FULL.obo
OUT_DIR=~/NetDecoder_Example/breast_cancer/networks/

##ERnegative
#ER- sources with control co-expression network
geneList=~/NetDecoder_Example/breast_cancer/sig_genes_control_ERnegative_2015-07-15.txt

PPI=~/NetDecoder_Example/breast_cancer/co_expression_network_breast_cancer_control_2015-07-09.txt
state=control
fCompositeNetwork=ControlNetwork_ERnegative_Sources		#output filename for context scores
echo "1) Computing subnetworks for $state state"
java -cp $LIB "netdecoder.NetDecoder" -PPI $PPI -SYMBOL $SYMBOL -GO $GO -g $geneList -out $OUT_DIR -gen -d $state -f $fCompositeNetwork

#ER- sources with ER- co-expression network
PPI=~/NetDecoder_Example/breast_cancer/co_expression_network_breast_cancer_ERnegative_2015-07-09.txt
state=ERnegative
fCompositeNetwork=ERnegativeNetwork_ERnegative_Sources		#output filename for context scores
echo "2) Computing subnetworks for $state state"
java -cp $LIB "netdecoder.NetDecoder" -PPI $PPI -SYMBOL $SYMBOL -GO $GO -g $geneList -out $OUT_DIR -gen -d $state -f $fCompositeNetwork

##ERpositive
#ER+ sources with control co-expression network
geneList=~/NetDecoder_Example/breast_cancer/sig_genes_control_ERpositive_2015-07-15.txt

PPI=~/NetDecoder_Example/breast_cancer/co_expression_network_breast_cancer_control_2015-07-09.txt
state=control
fCompositeNetwork=ControlNetwork_ERpositive_Sources		#output filename for context scores
echo "1) Computing subnetworks for $state state"
java -cp $LIB "netdecoder.NetDecoder" -PPI $PPI -SYMBOL $SYMBOL -GO $GO -g $geneList -out $OUT_DIR -gen -d $state -f $fCompositeNetwork

#ER+ sources with ER+ co-expression network
PPI=~/NetDecoder_Example/breast_cancer/co_expression_network_breast_cancer_ERpositive_2015-07-09.txt
state=ERpositive
fCompositeNetwork=ERpositiveNetwork_ERpositive_Sources		#output filename for context scores
echo "2) Computing subnetworks for $state state"
java -cp $LIB "netdecoder.NetDecoder" -PPI $PPI -SYMBOL $SYMBOL -GO $GO -g $geneList -out $OUT_DIR -gen -d $state -f $fCompositeNetwork

