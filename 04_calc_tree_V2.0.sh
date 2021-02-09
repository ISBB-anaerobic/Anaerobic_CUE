#!/bin/bash -   
#title          :05_calc_tree_V2.0.sh
#description    :Align DADA2 output sequences using mafft and build a an ML tree using IQ-TREE
#author         :Roey Angel
#date           :20201214
#version        :2.0    
#usage          :After running DADA2 pipeline, run from the Analysis folder: 05_calc_FastTree.sh $full_path_to_rep_seqs.fa ./05_calc_fasttree.sh
#notes          :  $rep_seqs.fa is the full path (absolute or relative) to a fasta file containing OTU representatives or ASVs     
#bash_version   :5.0.17(1)-release
#============================================================================


## TODO


eval "$(conda shell.bash hook)" # this is needed to be able to use conda activate

HOMEFOLDER=`pwd`
DBPATH="/data/Resources/SILVA/"
DB="silva_refnr_132.arb"


# Start
logFile="05_calc_tree.log"

if [ -z ${1+x} ]
    then 
        echo "No fasta file provided. Exiting" > ${HOMEFOLDER}/${logFile}
        exit 1 >> ${HOMEFOLDER}/${logFile}
    else 
        inputPath=$1
fi
DATAFOLDER="${inputPath%/*}"
INPUTSEQS=$(basename $inputPath)

echo -e "05_calc_tree.sh ${INPUTSEQS}" > ${HOMEFOLDER}/${logFile}
echo -e "Will process DADA2 sequences from the following file:" >> ${HOMEFOLDER}/${logFile}
echo -e $(ls ${inputPath}) >> ${HOMEFOLDER}/${logFile}

if [ -d ${DATAFOLDER}/Tree ]
then
 rm -rf ${DATAFOLDER}/Tree
fi
mkdir ${DATAFOLDER}/Tree

cp ${inputPath} ${DATAFOLDER}/Tree
cd ${DATAFOLDER}/Tree

# Align the sequences against a DB
echo -e "Align the sequences against a DB" >> ${HOMEFOLDER}/${logFile}
#mafft --add  --thread 50 ${DBPATH}${DB} > ${INPUTSEQS%.fa*}.align &>> ${HOMEFOLDER}/${logFile}
conda activate sina
/usr/bin/time -v sina -i ${INPUTSEQS} --db ${DBPATH}${DB} -o ${INPUTSEQS%.fa*}.align &>> ${HOMEFOLDER}/${logFile}
sed -i 's/;size=[0-9]\+;//g' ${INPUTSEQS%.fa*}.align # remove ;size=[0-9]\+; added by sina
conda deactivate
mothur "# filter.seqs(fasta=${INPUTSEQS%.fa}.align, vertical=T, processors=20)"
cat mothur.* >> ${HOMEFOLDER}/${logFile}
rm mothur.*
mv ${INPUTSEQS%.fa}.filter.fasta ${INPUTSEQS%.fa}.filtered.align

# Make a quick tree using FastTree
echo -e "\n Make a quick tree using FastTree" >> ${HOMEFOLDER}/${logFile}
/usr/bin/time -v FastTree -gtr -nt < ${INPUTSEQS%.fa*}.filtered.align > ${INPUTSEQS%.fa*}.FT.tree | tee -a ${HOMEFOLDER}/${logFile}
sed -i '/(Seq_/!p' ${HOMEFOLDER}/${logFile} # becuase the tree also gets printed out

# Calculate an iqtree using fast and normal bootstrapping (ultrafast bs doesn't work with fast)
echo -e "\n Calculate an iqtree using fast and normal bootstrapping" >> ${HOMEFOLDER}/${logFile}
#/usr/bin/time -v  iqtree -s ${INPUTSEQS%.fa*}.align -m GTR+R10 -alrt 1000 -bb 1000 -nt AUTO >> ${HOMEFOLDER}/${logFile}
#/usr/bin/time -v iqtree -s ${INPUTSEQS%.fa*}.filtered.align -fast -mset GTR -mrate R -cmin 5 -cmax 15 -b 100 -nt AUTO >> ${HOMEFOLDER}/${logFile}
#/usr/bin/time -v mpirun -np 20 iqtree-mpi -nt 1 -s ${INPUTSEQS%.fa*}.filtered.align -fast >> ${HOMEFOLDER}/${logFile}
/usr/bin/time -v iqtree -nt AUTO -s ${INPUTSEQS%.fa*}.filtered.align -fast >> ${HOMEFOLDER}/${logFile}

rm ${INPUTSEQS} # rm copied seq file
zip -m ${INPUTSEQS%.fa*}.align.zip ${INPUTSEQS%.fa*}.align # zip large alignment file

