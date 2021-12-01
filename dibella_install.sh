#!/bin/bash

LOWER_KMER_FREQ=$1
UPPER_KMER_FREQ=$2
DELTACHERNOFF=$3
FUZZ=$4

cd /Users/gabrielraulet/Desktop/CONTIGGING/diBELLA.2D

export COMBBLAS_HOME=$PWD
export BLOOM_HOME=$PWD/src/libbloom/
export SEQAN_HOME=$PWD/seqan
export MYPROJDIR=$PWD

CURDIR=$PWD

rm -rf build_release
mkdir -p build_release
cd build_release

export CC=gcc-11
export CXX=g++-11

cmake .. -DLOWER_KMER_FREQ=$LOWER_KMER_FREQ -DUPPER_KMER_FREQ=$UPPER_KMER_FREQ -DDELTACHERNOFF=$DELTACHERNOFF -DFUZZ=$FUZZ
make -j4
cd $CURDIR
