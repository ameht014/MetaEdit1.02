#!/bin/bash

module load singularity-3.5.3

THREADS=32

SAMPLE_TYPE="DNA"
mkdir -p ./logs
LOGS_DIR=./logs/${SAMPLE_TYPE}/
mkdir -p $LOGS_DIR

######################### SLURM PARAMETERS #########################
SLURM_ACC="iacc_giri"
SLURM_QOS="pq_giri"
SLURM_NODE_TYPE="investor"
####################################################################

################################ INPUT ##############################
SAMPLES_DIR="/home/vsteb002/samples/iHMP/IBD"

REF_DIR="../../data/References/"
REF_LIST=${REF_DIR}/ref_list.txt
OUT_DIR="./out_DNA/"
mkdir -p $OUT_DIR

SAMPLES_DIR=${SAMPLES_DIR}/${SAMPLE_TYPE}/
TRIMMED_DIR=${SAMPLES_DIR}/trimmed/
mkdir -p $TRIMMED_DIR

#####################################################################

#sbatch -J metaEdit -a 0-400 \
#                --account=$SLURM_ACC \
#                --qos=$SLURM_QOS \
#                -p $SLURM_NODE_TYPE \
#                -N 1 \
#                -n $THREADS \
#                -o "${LOGS_DIR}/stdout-%a.txt" \
#                -e "${LOGS_DIR}/stderr-%a.txt" \
#                ./metaEdit_one.sh $SAMPLES_DIR $REF_DIR $REF_LIST $TRIMMED_DIR $OUT_DIR $THREADS

sbatch -J metaEdit -a 400-748 \
                --account=$SLURM_ACC \
                --qos=$SLURM_QOS \
                -p $SLURM_NODE_TYPE \
                -N 1 \
                --cpus-per-task $THREADS \
                -o "${LOGS_DIR}/stdout-%a.txt" \
                -e "${LOGS_DIR}/stderr-%a.txt" \
                ./metaEdit_one.sh $SAMPLES_DIR $REF_DIR $REF_LIST $TRIMMED_DIR $OUT_DIR $THREADS
