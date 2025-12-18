#!/bin/bash
#
#SBATCH --job-name=L20FHmixD1D
#SBATCH --comment=""
#SBATCH --time=6-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hannah.lange@physik.uni-muenchen.de
#SBATCH --chdir=
#SBATCH --output=slurm.%j.%N.out
#SBATCH --constraint=x86-64-v3
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=99GB


source /project/th-scratch/s/Sebastian.Paeckel/init_syten.sh

python3 setup.py -l 41 -t 1.0 -timp 1.0 -Delta -4.0 -U 0.0 -Uimp 8.0 -N 29 -dmrg 1

