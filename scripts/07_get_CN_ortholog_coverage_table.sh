#!/bin/bash
set -ueo pipefail

# start from scripts directory

#############################
# Setup variables and paths #
#############################

pathways_dir="../data/C_Pathways"
annot_dir="../data/annotations"
C_paths=${pathways_dir}/paths.txt


echo "step1"

#############################################
# Get associated genes from C pathways list #
#############################################

while read pathway
do
# download with wget
#https://www.genome.jp/dbget-bin/get_linkdb?-t+orthology+path:map00910

if ! [ -f ${pathways_dir}/${pathway}.paths ]; then
wget $(echo https://www.genome.jp/dbget-bin/get_linkdb?-t+orthology+path:map$pathway) -O ${pathways_dir}/${pathway}.paths
fi
done < $C_paths

echo "step2"


###########################################
# Pull out C genes from annotation tables #
###########################################

# make list of all KEGG genes from the listed C pathways
cat ${pathways_dir}/*.paths | cut -d ">" -f2 | cut -d "<" -f1 | grep "^K" > ALL_KO_GENES

echo "step3"

for sample in ${annot_dir}/*[0,1,2,3,4,5,6,7,8,9]_annotation_table.tsv; 
do 
SAMPLEID=$(basename ${sample/_annotation_table.tsv//}); 
echo "Processing $SAMPLEID"
grep -Fwf ALL_KO_GENES $sample > ${annot_dir}/${SAMPLEID}_C_genes.txt; 

echo "step4"

# look at new C_genes.txt file
# find the KEGG ID
# use it to find the total number of assigned reads in _scafstats.txt
# add that value to the _C_genes.txt file


#   WHERE DID I GET SCAFSTATS FILES!?????


# build lookup terms
paste -d "|" <(cat ${annot_dir}/${SAMPLEID}_C_genes.txt | cut -f2) <(cat ${annot_dir}/${SAMPLEID}_C_genes.txt | cut -f1) > ${annot_dir}/${SAMPLEID}_lookup_terms.txt
grep -f ${annot_dir}/${SAMPLEID}_lookup_terms.txt ${annot_dir}/${SAMPLEID}_scafstats.txt > ${annot_dir}/${SAMPLEID}_C_gene_scafstats.txt

echo "step5"

# build simple lookup table
cat ${annot_dir}/${SAMPLEID}_C_gene_scafstats.txt | cut -f 1 | cut -d "|" -f1 > a
cat ${annot_dir}/${SAMPLEID}_C_gene_scafstats.txt | cut -f 8 > b
paste a b | sort -k1 > ${annot_dir}/${SAMPLEID}_C_gene_readcounts.txt # sorted into same order as _C_genes.txt

echo "step6"

# add to _C_genes.txt
cat <(echo "predicted_gene  uniprot_id  ko_id   ko_description  pfam_id pfam_description    go_term   uniprot_id_2    assigned_reads") <(paste ${annot_dir}/${SAMPLEID}_C_genes.txt ${annot_dir}/${SAMPLEID}_C_gene_readcounts.txt) > ${annot_dir}/${SAMPLEID}_C_gene_annotation_table.tsv

echo "step7"

# cleanup
rm a b ${annot_dir}/${SAMPLEID}_C_genes.txt ${annot_dir}/${SAMPLEID}_C_gene_readcounts.txt ${annot_dir}/${SAMPLEID}_lookup_terms.txt ${annot_dir}/${SAMPLEID}_C_gene_scafstats.txt

echo "step8"

done

# final cleanup (outside of loop)
rm ALL_KO_GENES

