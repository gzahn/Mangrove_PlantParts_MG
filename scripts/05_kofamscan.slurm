#!/bin/bash
#SBATCH --time=14-00:00:00 # days-hours:minutes:seconds (or just minutes)
#SBATCH --nodes=1
#SBATCH -o /scratch/general/vast/Zahn/Mangrove_Annotations/KEGG_Profiles/slurm-%j.out-%N
#SBATCH -e /scratch/general/vast/Zahn/Mangrove_Annotations/KEGG_Profiles/slurm-%j.err-%N
#SBATCH --ntasks=64
#SBATCH --account=uvucos-np
#SBATCH --partition=uvucos-np
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zahn.geoff@gmail.com
#SBATCH --job-name=KOfamscan


module load hmmer3/3.2.1 zlib/1.2.13 gcc/13.1.0 python/3.12.4
exec_annotation="/scratch/general/vast/Zahn/Mangrove_Annotations/ko/kofam_scan-1.3.0/exec_annotation"
working_dir=/scratch/general/vast/Zahn/Mangrove_Annotations/ko/kofam_scan-1.3.0
file_dir=/scratch/general/vast/Zahn/Mangrove_Annotations/metaeuk
outdir=/scratch/general/vast/Zahn/Mangrove_Annotations/KEGG_Profiles

cd $working_dir

for i in $file_dir/*_contigs.predsResults.fas;
do
fna=$i
out=${i/_contigs.predsResults.fas/_kofamscan.txt}
echo $out

./exec_annotation -f detail-tsv --ko-list ./ko_list -p ../profiles -o $out --cpu 64 $fna

done
