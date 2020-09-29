#!/bin/sh
#
# run NetDecoder, run!
#

LIB_DIR=~/NetDecoder_Example/lib
FILES_DIR=~/NetDecoder_Example/lib
LIB=$LIB_DIR/NetDecoder.jar:$LIB_DIR/commons-cli-1.2.jar
SYMBOL=$FILES_DIR/gene_associationgoa_ref_HUMAN_25Jun2014.txt
GO=$FILES_DIR/gene_ontology_ext_FULL.obo

OUT_DIR=~/NetDecoder_Example/breast_cancer/analysis/
BASE_DIR=~/NetDecoder_Example/breast_cancer/networks/

###ERnegative
netpath1=$BASE_DIR/ControlNetwork_ERnegative_Sources_paths.ser
netpath2=$BASE_DIR/ERnegativeNetwork_ERnegative_Sources_paths.ser
ann1=control
ann2=ERnegative
ratioThreshold=5
corThreshold=0.5
top=10
geneList=none

java -cp $LIB "netdecoder.NetDecoder" -SYMBOL $SYMBOL -GO $GO -out $OUT_DIR -E -control $ann1 -condition $ann2 -ncp $netpath1 -ndp $netpath2 -corThreshold $corThreshold -ratioThreshold $ratioThreshold -top $top -g $geneList -overlap -f ERnegative

###ERpositive
netpath1=$BASE_DIR/ControlNetwork_ERpositive_Sources_paths.ser
netpath2=$BASE_DIR/ERpositiveNetwork_ERpositive_Sources_paths.ser
ann1=control
ann2=ERpositive
ratioThreshold=5
corThreshold=0.5
top=10

java -cp $LIB "netdecoder.NetDecoder" -SYMBOL $SYMBOL -GO $GO -out $OUT_DIR -E -control $ann1 -condition $ann2 -ncp $netpath1 -ndp $netpath2 -corThreshold $corThreshold -ratioThreshold $ratioThreshold -top $top -g $geneList -overlap -f ERpositive

