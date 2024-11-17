# Assign variables
metaeuk=/uufs/chpc.utah.edu/common/home/u6033249/programs/metaeuk_new/metaeuk/bin/metaeuk
mmseqs=/uufs/chpc.utah.edu/common/home/u6033249/programs/mmseqs/bin/mmseqs
metaeuk_db=/scratch/general/vast/Zahn/Databases/swissprot


ass_dir=/scratch/general/vast/Zahn/Mangrove_Transplant_MG/assemblies
annot_dir=/scratch/general/vast/Zahn/Mangrove_Transplant_MG/annotations

for i in ${ass_dir}/*.fasta.gz
do 
        echo $i
        contigs=$i
        sample_id=$(basename $i | cut -d "_" -f1)

# metauk command
        $metaeuk easy-predict $contigs $metaeuk_db ${annot_dir}/${sample_id}_predsResults ${annot_dir}/${sample_id}_tempFolder

done


