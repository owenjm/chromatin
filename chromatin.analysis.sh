#!/bin/bash

STATE=20
NAME=analysis
ITER=1000
INIT=10
DATE=`date +%d-%m-%Y`
SCRIPT=chromatin.all.states.hmms.r

# parse the options
while getopts 's:n:i:j:k:g:h:' opt ; do
  case $opt in
    s) STATE=$OPTARG ;;
    n) NAME=$OPTARG ;;
    i) ITER=$OPTARG ;;
    j) INIT=$OPTARG ;;
	k) SCRIPT=$OPTARG ;;
    g) GFFPATH=$OPTARG ;;
    h) echo "options:
	-s [state, DEFAULT=$STATE]
	-n [name]
	-i [iter, DEFAULT=$ITER]
	-j [init iterations, DEFAULT=$INIT]
	-k [script, DEFAULT=$SCRIPT]
	-g [gff_path]
	[-help]"; exit 1 ;;
  esac
done

echo "state $STATE
name: $NAME
iter: $ITER
init: $INIT
path: $GFFPATH";

DIR=analysis.$NAME.$DATE.iter.$ITER.init.$INIT.st.$STATE;
mkdir $DIR;
cd $DIR;

Rscript --max-ppsize=500000 $SCRIPT\
    --name=$NAME \
    --st.start=$STATE \
    --st.end=$STATE  \
    --gff.path=$GFFPATH \
    --rhmm.iter=$ITER \
    --rhmm.iter.init=$INIT;


