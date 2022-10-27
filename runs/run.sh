#!/bin/bash
#
mpirun -np 4 ../build_release/elba -i chimer_reads.fa --idxmap idmap -c 600 -k 31 --xa 15 -s 1 -O 100000 --afreq 100000 --pb --tip 5 2>&1 | tee elba.out
