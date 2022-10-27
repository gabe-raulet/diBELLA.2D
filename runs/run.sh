#!/bin/bash

mpirun -np 4 ../build_release/elba -i chimer_reads.fa --idxmap idmap -c 600 -k 31 --xa 15 -s 1 -O 100000 --afreq 100000 --pb --tip 5 2>&1 | tee elba.out

rm bridges.txt
rm contig-assignment.txt
rm elba_rank_*_log.txt
rm readNameMap_*
rm small_to_large_contig_map.txt
rm {string,overlap}.mtx

grep "^AA" elba.out | cut -f 2-5 | sort -V > seeds.tsv
