#!/bin/bash
#

# Description: run DADA2 pipeline on a MiSeq run data
# Instructions: go to the Analysis folder and run:
# 01_run_DADA2_16S_V8.3.sh $POOLING $DATAFOLDER $RESOURCES

HOMEFOLDER=$(pwd)
POOLING=${1:-pseudo} # pseudo, TRUE, FALSE
DATAFOLDER=${2:-"../../Data/16S/noPrimers/"}
RESOURCES=${3:-"~/Resources/"} # RESOURCES should be referred to from $DATAFOLDER/DADA2_pseudo!
SCRIPTS=${HOMEFOLDER}
GENE=16S
MOCKID="" # change to "" if no mock sample was used
VERSION=V8.4

if [ -d $DATAFOLDER/DADA2_${POOLING} ]
then
 rm -rf $DATAFOLDER/DADA2_${POOLING}
fi
mkdir $DATAFOLDER/DADA2_${POOLING}

LOG=01_run_DADA2_${GENE}_${POOLING}_${VERSION}.log
touch ${HOMEFOLDER}/${LOG}

echo -e "Run DADA2 ${VERSION} on ${GENE} dataset with ${POOLING} option  \n" >${HOMEFOLDER}/${LOG}
echo -e "Will process reads from the following folder:"  >> ${HOMEFOLDER}/${LOG}
echo -e $DATAFOLDER >> ${HOMEFOLDER}/${LOG}

echo -e "Starting R script \n" >> ${HOMEFOLDER}/${LOG}
cd $DATAFOLDER/DADA2_${POOLING}
/usr/bin/time -v Rscript --vanilla ${SCRIPTS}/01_DADA2_16S_${VERSION}.R $POOLING ../ $RESOURCES $MOCKID >> ${HOMEFOLDER}/${LOG} 2>&1
cd ${HOMEFOLDER}
mv ${DATAFOLDER}/DADA2_${POOLING} ${HOMEFOLDER}
rm -r ${DATAFOLDER}/filtered

