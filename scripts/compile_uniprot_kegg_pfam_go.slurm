#!/bin/bash
#SBATCH --time=14-00:00:00 # days-hours:minutes:seconds (or just minutes)
#SBATCH --nodes=1
#SBATCH -o /scratch/general/vast/Zahn/Mangrove_Annotations/KEGG_Profiles/slurm-%j.out-%N
#SBATCH -e /scratch/general/vast/Zahn/Mangrove_Annotations/KEGG_Profiles/slurm-%j.err-%N
#SBATCH --ntasks=96
#SBATCH --account=uvucos-np
#SBATCH --partition=uvucos-np
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zahn.geoff@gmail.com
#SBATCH --job-name=build_annot_table


# double-check directory
file_dir=/scratch/general/vast/Zahn/Mangrove_Annotations/metaeuk
cd $file_dir

# for each sample file, add info based on UniProt ID assigned my metaeuk

#for FILE in ${file_dir}/*_kofamscan.txt
# restart, using list of files (change to for-loop for finished codebase)
while read FILE
do

sampleid=$(basename ${FILE/_kofamscan.txt/})
echo $sampleid

# get unique list of uniprot query ids
cat $FILE | tail -n +3 | cut -f2 | cut -d "|" -f1 | sort -u > uniprot_list
echo "${sampleid} uniprot_list completed"


# find only first match of those in kofamscan results
while read id; do grep -m 1 "$id" $FILE; done < uniprot_list > best_ko_table
echo "${sampleid} best_ko_table completed"

# find only first match of those in pfam results
while read id; do grep -m 1 "$id" ../pfam/uniprotID_to_pfam_edited.tsv; done < uniprot_list | cut -f2 | cut -d " " -f5 | cut -d ";" -f1 > best_ko_pfamid
while read id; do grep -m 1 "$id" ../pfam/uniprotID_to_pfam_edited.tsv; done < uniprot_list | cut -f2 | cut -d " " -f6 | cut -d ";" -f1 > best_ko_pfamdesc
echo "${sampleid} first pfam matches completed"


# and the GO: assignment as well
while read id; do grep -m 1 "$id" ../pfam/uniprotID_to_pfam.tsv; done < uniprot_list  | cut -f4 | cut -d ";" -f 2 | cut -d " " -f2 > best_ko_goterm
echo "${sampleid} GO term lookup completed"

# cut up and rebuild into convenient format
cat best_ko_table | cut -f2 | cut -d "|" -f1 > best_ko_upid
cat best_ko_table | cut -f2 | cut -d "|" -f2 > best_ko_locus
cat best_ko_table | cut -f3 > best_ko_koid
cat best_ko_table | cut -f7 > best_ko_kodesc

# final table
cat <(echo "predicted_gene  uniprot_id  ko_id   ko_description  pfam_id pfam_description    go_term") <(paste best_ko_locus best_ko_upid best_ko_koid best_ko_kodesc best_ko_pfamid best_ko_pfamdesc best_ko_goterm) > ${file_dir}/${sampleid}_annotation_table.tsv

rm uniprot_list best_ko_table best_ko_pfamid best_ko_pfamdesc best_ko_goterm best_ko_upid best_ko_locus best_ko_koid best_ko_kodesc

echo "${sampleid}_annotation_table.tsv has been built"

done < inputfiles
