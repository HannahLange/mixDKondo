#!/bin/bash
#
#SBATCH --job-name=TL20tJKondo
#SBATCH --comment=""
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hannah.lange@physik.uni-muenchen.de
#SBATCH --chdir=
#SBATCH --output=slurm.%j.%N.out
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=100GB

source /project/th-scratch/s/Sebastian.Paeckel/init_syten.sh

python3 setup.py

