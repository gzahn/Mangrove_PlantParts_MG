#!/bin/bash
#SBATCH --time=3-00:00:00 # days-hours:minutes:seconds (or just minutes)
#SBATCH --nodes=1      
#SBATCH -o slurm-%j.out-%N 
#SBATCH -e slurm-%j.err-%N 
#SBATCH --ntasks=20
#SBATCH --account=uvu     
#SBATCH --partition=kingspeak-shared
#SBATCH --mail-type=ALL                                                     
#SBATCH --mail-user=zahn.geoff@gmail.com                  	
#SBATCH --job-name=normalize_all

#


# load bbtools
module load bbtools

# while loop to parse file names for bbtools

while read fn;
do
basenames=$(echo $(basename $fn))
base=${basenames/_1.fastq.gz/}
outpath=${fn%/*}
trimmed=${outpath}/${base}_trimmed_1.fastq.gz #cleaned reads from adaptor trimming step
trimmed2=${trimmed/_trimmed_1/_trimmed_2} 
normalized=${outpath}/${base}_normalized_1.fastq.gz # new filename for normalized kmers reads
normalized2=${normalized/normalized_1/normalized_2}

#Normalize
bbnorm.sh in=$trimmed in2=$trimmed2 out=$normalized out2=$normalized2 target=100 min=5
        

done < /uufs/chpc.utah.edu/common/home/u6033249/Mangrove_Metagenome/abs_filepaths

#abs_filepaths is just a list of abs paths to fwd reads


# took ~24 hours to run on all samples
