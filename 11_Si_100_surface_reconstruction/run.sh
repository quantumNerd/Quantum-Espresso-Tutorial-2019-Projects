#!/bin/sh
mkdir /scratch/anonymous
mkdir /scratch/anonymous/outdir_Si100

mpirun -np 12 ../qe-6.1_intel/bin/pw.x -inp si100.2_slab_relax.in > si100.2_slab_relax.out

rm -r /scratch/anonymous
