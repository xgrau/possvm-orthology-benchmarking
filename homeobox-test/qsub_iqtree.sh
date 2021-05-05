#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q long-sl7
#$ -l virtual_free=10G,h_rt=120:00:00
#$ -o tmp/
#$ -e tmp/

# input
i=$1 # ALIGNED FASTA
c=$2 # NUM CPUS

function do_phy {

	local ilt=$1
	local tre=$3
	local c=$2

	# iqtree
	/users/asebe/xgraubove/Programes/iqtree-2.1.0-Linux/bin/iqtree2 \
		-s $ilt \
		-m TEST \
		-mset LG,WAG,JTT \
		-nt AUTO \
		-ntmax $c \
		-bb 1000 \
		-pre $tre \
		-nm 10000 \
		-nstop 200 \
		-cptime 1800

}

# do alignments and trees
do_phy $i $c ${i%%.fasta}.iqtree
