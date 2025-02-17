#!/bin/bash

module load singularity-3.5.3

SAMPLE_TYPE="RNA"
mkdir -p ./logs
LOGS_DIR=./logs/${SAMPLE_TYPE}/
mkdir -p $LOGS_DIR
THREADS=32

######################### SLURM PARAMETERS #########################
SLURM_ACC="your_account"
SLURM_QOS="your_qos"
SLURM_NODE_TYPE="your_node"
####################################################################

################################ INPUT ##############################
SAMPLES_DIR="/path/to/samples/dir"
REF_DIR="../../data/References/"
REF_LIST=${REF_DIR}/ref_list.txt
OUT_DIR="/path/to/output/dir"
mkdir -p $OUT_DIR

SAMPLES_DIR=${SAMPLES_DIR}/${SAMPLE_TYPE}/
TRIMMED_DIR=${OUT_DIR}/trimmed_RNA/
mkdir -p $TRIMMED_DIR

sbatch -J metaEdit -a 0-748 \
                --account=$SLURM_ACC \
                --qos=$SLURM_QOS \
                -p $SLURM_NODE_TYPE \
                -N 1 \
                --cpus-per-task $THREADS \
                -o "${LOGS_DIR}/stdout-%a.txt" \
                -e "${LOGS_DIR}/stderr-%a.txt" \
                ./metaEdit_one.sh $SAMPLES_DIR $REF_DIR $REF_LIST $TRIMMED_DIR $OUT_DIR $THREADS


