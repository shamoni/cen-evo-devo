#!/bin/bash

CHIP=$1
INPUT=$2
OUTPREFIX=$3

echo -e "Sample:"$CHIP "Control:"$INPUT

module load macs2

macs2 callpeak -t $CHIP -c $INPUT -f BAM --outdir ${OUTPREFIX} -n ${OUTPREFIX} -B --nomodel --extsize 165 --keep-dup all 2> ${OUTPREFIX}.stdout

macs2 bdgcmp -t ${OUTPREFIX}/${OUTPREFIX}_treat_pileup.bdg -c ${OUTPREFIX}/${OUTPREFIX}_control_lambda.bdg -m FE -o ${OUTPREFIX}_FE.bdg

LC_COLLATE=C sort -k1,1 -k2,2n ${OUTPREFIX}_FE.bdg > ${OUTPREFIX}_FE_sorted.bdg

~/bin/bedGraphToBigWig ${OUTPREFIX}_FE_sorted.bdg T10.chrom.sizes ${OUTPREFIX}_FE.bw

rm ${OUTPREFIX}_FE.bdg

mv ${OUTPREFIX}_FE_sorted.bdg ${OUTPREFIX}_FE.bdg

mv ${OUTPREFIX}* ${OUTPREFIX}/
