#!/bin/bash -login
#SBATCH --mem=64GB
#SBATCH --job-name=kb_count
#SBATCH --output=%x-%j.out
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32

# print start time
start=`date +%s`

# load the conda module 
module purge
module load Conda/3

# make count matrix
kb count \
  -o /mnt/research/bioinformaticsCore/projects/petroffm/BCC105_scRNAseq_training/data/Lukassen_testes/SRR6129051_mouse2 \
  -t 32 \
  -i /mnt/research/bioinformaticsCore/projects/petroffm/BCC105_scRNAseq_training/data/kallisto_mouse.idx  \
  -g /mnt/research/bioinformaticsCore/projects/petroffm/BCC105_scRNAseq_training/data/kallisto_mouse_t2g.txt \
  -x 10XV2 \
  --filter bustools \
  $SCRATCH/BCC105_scRNAseq_training/raw_data/Lukassen_testes/SRR6129051_1.fastq \
  $SCRATCH/BCC105_scRNAseq_training/raw_data/Lukassen_testes/SRR6129051_2.fastq 

# print end time
echo "Finished"
end=`date +%s`
runtime=$((end-start))
echo execution time was `expr $end - $start` 