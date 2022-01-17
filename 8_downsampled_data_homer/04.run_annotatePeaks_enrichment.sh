
#!/bin/bash

# Running findMotifsGenome HOMER script
# Finding enriched motifs in DARs previously reported.

source ~/.bashrc

cd ../final_dars/

#for file in $(ls -d *B_cell__DARs.bed.tmp)
for file in $(ls -d *B_cell__DARs.TOP2K.bed.tmp)
do
	echo $file
	query="${file}"

	# Generating background file (merging all the techniques and all the cell types)
	cat libds__merged__B_cell__DARs.bed.tmp libds__merged__CD4+_T_cell__DARs.bed.tmp libds__merged__CD14+_monocyte__DARs.bed.tmp libds__merged__CD16+_monocyte__DARs.bed.tmp libds__merged__Cytotoxic_T_cell__DARs.bed.tmp libds__merged__Natural_killer_cell__DARs.bed.tmp libds__merged__Dendritic_cell__DARs.bed.tmp > background_libds__merged.tmp 

	awk 'BEGIN {OFS="\t"}; {print NR,$2,$3,$4,$5}' background_libds__merged.tmp >> background_libds__merged.txt
	rm background_libds__merged.tmp 

	name_output=${file/.txt.tmp/}

	mkdir "scATAC_Benchmarking/HOMER_results/enrichment_motifs_2k/${name_output}"
	output_enr="scATAC_Benchmarking/HOMER_results/enrichment_motifs_2k/${name_output}/"

perl HOMER/bin/findMotifsGenome.pl $query hg38 $output_enr -size given -bg background_libds__merged.txt
done