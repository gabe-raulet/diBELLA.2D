#!/bin/bash

DATASET=~/Desktop/CONTIGGING/read-datasets/ecoli_hifi_29x.numbered.fasta
EXE=~/Desktop/CONTIGGING/diBELLA.2D/build_release/dibella

RESULT=dibella.grcontig.ecoli.n1.L20.U30.D0.7.K31.X15
rm -rf $RESULT
mkdir $RESULT

NUM_READS=`grep '>' $DATASET | wc -l`

export OMP_NUM_THREADS=1

cd $RESULT
mpirun -np 1 $EXE -i $DATASET -k 31 --idxmap $RESULT.idxmap -c $NUM_READS --alph dna -s 1 -O 100000 --afreq 100000 --xa 15 --af $RESULT.af
