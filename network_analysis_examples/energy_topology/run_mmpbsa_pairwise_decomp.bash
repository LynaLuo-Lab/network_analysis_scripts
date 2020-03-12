#!/bin/bash

export AMBERHOME=/Users/wbotellosmith/anaconda2/envs/py3

sysName=$1

trajExt='dcd'
topExt='parm7'

echo "running: MMPBSA.py -O -i mmpbsa.pairwise_decomp.in -o ${sysName}.mmpbsa.pairwise_decomp.out -cp ${sysName}.${topExt} \"
echo "-do ${sysName}.mmpbsa.pairwise_decomp.decomp_summary.csv -eo ${sysName}.mmpbsa.pairwise_decomp.energy_series.csv \"
echo "-deo ${sysName}.mmpbsa.pairwise_decomp.decomp_energy_series.csv -y ${sysName}.${trajExt} "

MMPBSA.py -O -i mmpbsa.pairwise_decomp.in -o ${sysName}.mmpbsa.pairwise_decomp.out -cp ${sysName}.${topExt} \
	-do ${sysName}.mmpbsa.pairwise_decomp.decomp_summary.csv -eo ${sysName}.mmpbsa.pairwise_decomp.energy_series.csv \
	-deo ${sysName}.mmpbsa.pairwise_decomp.decomp_energy_series.csv -y ${sysName}.${trajExt}
echo "done"
