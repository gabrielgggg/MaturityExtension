#!/bin/bash
#SBATCH -p standard
#SBATCH -N 3 
#SBATCH --ntasks-per-node=1
#SBATCH -c 24 
#SBATCH -t 0-03:00:00
#SBATCH -J jmp_arrears
#SBATCH -o out.txt
#SBATCH -e err.txt
#SBATCH --mem=64gb
#SBATCH 

module load circ slurm intel openmpi/2.0.1/b2 matlab

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

echo "$SLURM_JOB_NUM_NODES nodes with $SLURM_CPUS_ON_NODE processors each." > log.txt
date >> log.txt

make clean &>> log.txt
make -j &>> log.txt

srun bin/jmp &>> log.txt

cd results/
srun -n 1 matlab -nojvm -nodesktop -r "toMatlab; loadSim; loadEvents; exit" &>> ../log.txt

date >> ../log.txt
echo "End." &>> ../log.txt
