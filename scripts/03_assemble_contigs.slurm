#!/bin/bash
#SBATCH --time=3-00:00:00 # days-hours:minutes:seconds (or just minutes)
#SBATCH --nodes=1      
#SBATCH -o slurm-%j.out-%N 
#SBATCH -e slurm-%j.err-%N 
#SBATCH --ntasks=20
#SBATCH --account=uvu     
#SBATCH --partition=lonepeak
#SBATCH --mail-type=END,FAIL                                                     
#SBATCH --mail-user=zahn.geoff@gmail.com                  	
#SBATCH --job-name=metaspades
#SBATCH --mem=250000

#

# load appropriate modules

module load spades

# run spades separately on each sample

# change file pattern matching to run on cleaned/normalized data

for file1 in /scratch/kingspeak/serial/u6033249/Mangrove_Metagenome/RawData/*_normalized_1.fastq.gz
do 
file2=${file1/_1.fastq.gz/_2.fastq.gz} 
out=${file1%%.fastq.gz}_spadesout 
metaspades.py -1 $file1 -2 $file2 -o $out -t 20 -m 299
done
