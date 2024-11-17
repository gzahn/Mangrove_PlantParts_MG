
#!/bin/bash
#SBATCH --time=14-00:00:00 # days-hours:minutes:seconds (or just minutes)
#SBATCH --nodes=1
#SBATCH -o /scratch/general/vast/Zahn/Mangrove_Annotations/pfam/slurm-%j.out-%N
#SBATCH -e /scratch/general/vast/Zahn/Mangrove_Annotations/pfam/slurm-%j.err-%N
#SBATCH --ntasks=64
#SBATCH --account=uvucos-np
#SBATCH --partition=uvucos-np
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zahn.geoff@gmail.com
#SBATCH --job-name=KOfamscan




# map raw reads to gene fastas | count hits
module load bbtools
RAW_DIR=./data/raw
REF_DIR=./data/fna


cd $REF_DIR

for REF_GENES in ${REF_DIR}*contigs.predsResults.codon.fas
do

SAMPLEID=$(basename ${REF_GENES/_contigs.predsResults.codon.fas/})
INPUT1=${RAW_DIR}${SAMPLEID}_1.fastq.gz
INPUT2=${RAW_DIR}${SAMPLEID}_2.fastq.gz

bbsplit.sh in=$INPUT1 in2=$INPUT2 ref=$REF_GENES scafstats=${SAMPLEID}_scafstats.txt refstats=${SAMPLEID}_refstats.txt ambiguous=best threads=64

done
