#!/bin/bash

IFS=$'\n' read -d '' -r -a ALL_SAMPLES < ./data/samples.txt
SAMPLE=${ALL_SAMPLES[${SLURM_ARRAY_TASK_ID}]}

# $SAMPLES_DIR $REF_DIR $REF_LIST $OUT_DIR

SAMPLES_DIR=$1

REF_DIR=$2
REF_LIST=$3
TRIMMED_DIR=$4
OUT_DIR=$5
THREADS=$6


mkdir -p $OUT_DIR

SAMPLES_DIR=${SAMPLES_DIR}/${SAMPLE}/

#if [ ! -d $SAMPLES_DIR ]; then
#  tar -xvf ${SAMPLES_DIR}/${SAMPLE}.tar

#SAMPLES_DIR=${SAMPLES_DIR}/${SAMPLE_TYPE}/
TRIMMED_DIR="${TRIMMED_DIR}/${SAMPLE}"
mkdir -p $TRIMMED_DIR
# Try different affixes and see if they exist
AFFIX1="_1.fastq.gz"
AFFIX2='_2.fastq.gz'

if [ ! -f "${SAMPLES_DIR}/${SAMPLE}${AFFIX1}" ]; then
  AFFIX1="_R1.fastq.gz"
  AFFIX2='_R2.fastq.gz'
  if [ ! -f "${SAMPLES_DIR}/${SAMPLE}${AFFIX1}" ]; then
    AFFIX1="_1_sequence.txt.bz2"
    AFFIX2='_2_sequence.txt.bz2'
  fi
fi

mkdir -p $OUT_DIR
AIMAP_PATH="../../src/third_party/aimap/"
chmod +x ${AIMAP_PATH}/bin/Adenosine_to_inosine.py
while IFS=, read -r DIR REF_NAME; do

    REF_DIR_I="${REF_DIR}/${DIR}/"
    echo $REF_DIR_I $REF_NAME
    OUT_DIR_I=${OUT_DIR}/${SAMPLE}.${REF_NAME}

    singularity exec --bind $SAMPLES_DIR ../../env/metaEdit.sif python3 ${AIMAP_PATH}/bin/Adenosine_to_inosine.py  \
            -g ${REF_DIR_I}/${REF_NAME}.fna \
            -tr ${TRIMMED_DIR} \
            -l 18 \
            -a ${REF_DIR_I}/${REF_NAME}.gff  \
            -t $THREADS \
            -o $OUT_DIR_I  \
            --outfile_name ${SAMPLE}.${REF_NAME} \
            -m paired \
             -f1  ${SAMPLES_DIR}/${SAMPLE}${AFFIX1} \
             -f2  ${SAMPLES_DIR}/${SAMPLE}${AFFIX2}

    rm -f ${OUT_DIR_I}/*.sam
    rm -f ${OUT_DIR_I}/*.bam
  done < $REF_LIST

rm -r $TRIMMED_DIR

#/bin/singularity exec ./../env/metaEdit.sif cutadapt --help

