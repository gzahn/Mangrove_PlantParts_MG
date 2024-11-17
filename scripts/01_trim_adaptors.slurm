#!/bin/bash
#SBATCH --time=300 # days-hours:minutes:seconds (or just minutes)
#SBATCH --nodes=1      
#SBATCH -o slurm-%j.out-%N 
#SBATCH -e slurm-%j.err-%N 
#SBATCH --ntasks=20
#SBATCH --account=uvu     
#SBATCH --partition=kingspeak-shared
#SBATCH --mail-type=END,FAIL,START                                                     
#SBATCH --mail-user=zahn.geoff@gmail.com                  	
#SBATCH --job-name=trim_adaptors

#


# load bbtools
module load bbtools

# while loop to parse file names for bbtools

while read fn;
do
fn2=${fn/_1.fastq.gz/_2.fastq.gz}
basenames=$(echo $(basename $fn))
base=${basenames/_1.fastq.gz/}
outpath=${fn%/*}
full_out=${outpath}/${base}_trimmed_1.fastq.gz
full_out2=${full_out/_trimmed_1.fastq.gz/_trimmed_2.fastq.gz}

# run bbduk.sh to trim adaptors	(and merge???)
bbduk.sh in=$fn in2=$fn2 out=$full_out out2=$full_out2 ktrim=r k=23 mink=11 hdist=1 ref=adapters  

done < /uufs/chpc.utah.edu/common/home/u6033249/Mangrove_Metagenome/abs_filepaths

#abs_filepaths is just a list of abs paths to fwd reads

# run took ~3 hours to complete
