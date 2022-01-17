
#!/bin/bash

source ~/.bashrc

# Running AnnotatePeaks HOMER script
# Gene ontology analysis of associated genes
# Statistics associated to the peak annotation

cd ../final_dars/

#for file in $(ls -d *B_cell__DARs.bed.tmp)
for file in $(ls -d *B_cell__DARs.TOP2K.bed.tmp)

do
	echo $file
	query="${file}"

	name_output=${file/.bed.tmp/}

	output="scATAC_Benchmarking/HOMER_results/stats/${name_output}_stats.txt"

	mkdir "scATAC_Benchmarking/HOMER_results/GO/${name_output}"
	output_go="scATAC_Benchmarking/HOMER_results/GO/${name_output}/"

perl HOMER/bin/annotatePeaks.pl $query hg38 -annStats $output -go $output_go
done