#!/bin/bash -login
#SBATCH --mem=64GB
#SBATCH --job-name=prefetch_fasterqdump
#SBATCH --output=%x-%j.out
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32

# print start time
start=`date +%s`

# load the SRA-toolkit module
module purge
module load SRA-Toolkit/3.0.10-gompi-2023a

# set the dir you want the data to go in
RAWDATA=$SCRATCH/BCC105_scRNAseq_training/raw_data/Lukassen_testes

# save the sra file to the appropriate data directory
prefetch SRR6129050 -O $RAWDATA 
prefetch SRR6129051 -O $RAWDATA

# extract the fastq files from the sra file
cd $RAWDATA
fasterq-dump SRR6129050 -e 32
fasterq-dump SRR6129051 -e 32

# print end time
echo "Finished"
end=`date +%s`
runtime=$((end-start))
echo execution time was `expr $end - $start` 