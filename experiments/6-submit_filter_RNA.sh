#!/bin/bash

THREADS=32
######################### SLURM PARAMETERS #########################
SLURM_ACC="your_account"
SLURM_QOS="your_qos"
SLURM_NODE_TYPE="your_node"
LOGS_DIR="./logs/filter/"
mkdir -p $LOGS_DIR

DNA_DIR=$(pwd)'/out_DNA/'
RNA_DIR="/path/to/RNA/output/"
REF_DIR=$(pwd)'/../../data/References/'
SAMPLES_LIST=$(pwd)'/data/samples.txt'
RNA_COV_THRESHOLD="30"
DNA_COV_THRESHOLD="30"
EDIT_THRESHOLD="0.03"
OUT_RNA_DIR=$(pwd)"/out_RNA_filtered/"
OUT_DNA_DIR=$(pwd)"/out_DNA_filtered/"

REF_LIST=${REF_DIR}/ref_list.txt

while IFS=, read -r BACTERIA_NAME REF_NAME; do
    sbatch -J filterSites --account=$SLURM_ACC \
                --qos=$SLURM_QOS \
                -p $SLURM_NODE_TYPE \
                -N 1 \
                --cpus-per-task $THREADS \
                -o "${LOGS_DIR}/stdout-${BACTERIA_NAME}.txt" \
                -e "${LOGS_DIR}/stderr-${BACTERIA_NAME}.txt" \
              ./submit_filter_one.sh $DNA_DIR $RNA_DIR $REF_DIR $REF_NAME $BACTERIA_NAME $SAMPLES_LIST $RNA_COV_THRESHOLD $DNA_COV_THRESHOLD $EDIT_THRESHOLD $OUT_RNA_DIR $OUT_DNA_DIR
done < $REF_LIST


