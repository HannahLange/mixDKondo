#!/bin/bash
#
#SBATCH --job-name=TL20tJKondo
#SBATCH --comment=""
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hannah.lange@physik.uni-muenchen.de
#SBATCH --chdir=
#SBATCH --output=slurm.%j.%N.out
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=100GB

source /project/th-scratch/s/Sebastian.Paeckel/init_syten.sh

python3 setup.py -l 20 -t1 1.0 -t2 0.0 -J1 0.5 -J2 0.0 -Jperp 1.5 -N 18 -dmrg 0 -temp 1 -imps 3 4
python3 setup.py -l 20 -t1 1.0 -t2 0.0 -J1 0.5 -J2 0.0 -Jperp 1.5 -N 18 -dmrg 0 -temp 1 -imps 3 5
python3 setup.py -l 20 -t1 1.0 -t2 0.0 -J1 0.5 -J2 0.0 -Jperp 1.5 -N 18 -dmrg 0 -temp 1 -imps 3 6
python3 setup.py -l 20 -t1 1.0 -t2 0.0 -J1 0.5 -J2 0.0 -Jperp 1.5 -N 18 -dmrg 0 -temp 1 -imps 3 7
python3 setup.py -l 20 -t1 1.0 -t2 0.0 -J1 0.5 -J2 0.0 -Jperp 1.5 -N 18 -dmrg 0 -temp 1 -imps 3 8
python3 setup.py -l 20 -t1 1.0 -t2 0.0 -J1 0.5 -J2 0.0 -Jperp 1.5 -N 18 -dmrg 0 -temp 1 -imps 3 9
python3 setup.py -l 20 -t1 1.0 -t2 0.0 -J1 0.5 -J2 0.0 -Jperp 1.5 -N 18 -dmrg 0 -temp 1 -imps 3 10
python3 setup.py -l 20 -t1 1.0 -t2 0.0 -J1 0.5 -J2 0.0 -Jperp 1.5 -N 18 -dmrg 0 -temp 1 -imps 3 11
python3 setup.py -l 20 -t1 1.0 -t2 0.0 -J1 0.5 -J2 0.0 -Jperp 1.5 -N 18 -dmrg 0 -temp 1 -imps 3 12
python3 setup.py -l 20 -t1 1.0 -t2 0.0 -J1 0.5 -J2 0.0 -Jperp 1.5 -N 18 -dmrg 0 -temp 1 -imps 3 13
python3 setup.py -l 20 -t1 1.0 -t2 0.0 -J1 0.5 -J2 0.0 -Jperp 1.5 -N 18 -dmrg 0 -temp 1 -imps 3 14
python3 setup.py -l 20 -t1 1.0 -t2 0.0 -J1 0.5 -J2 0.0 -Jperp 1.5 -N 18 -dmrg 0 -temp 1 -imps 3 15
python3 setup.py -l 20 -t1 1.0 -t2 0.0 -J1 0.5 -J2 0.0 -Jperp 1.5 -N 18 -dmrg 0 -temp 1 -imps 3 16
python3 setup.py -l 20 -t1 1.0 -t2 0.0 -J1 0.5 -J2 0.0 -Jperp 1.5 -N 18 -dmrg 0 -temp 1 -imps 3 17

