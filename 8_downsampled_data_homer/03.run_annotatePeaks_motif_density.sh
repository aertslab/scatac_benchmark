
#!/bin/bash

# Running annotatePeaks HOMER script
# Targetted analysis: Motif distribution profile 
# The motifs: classical motifs detected in the B cell lineage. 

source ~/.bashrc

cd ../final_dars/

homer_motifs="HOMER/motifs/"

for file in $(ls -d *B_cell__DARs.TOP2K.bed)
do
	echo $file
	query="${file}"

	name_output=${file/.txt.tmp/}
	output="scATAC_Benchmarking/HOMER_results/density_motifs_2k/${name_output}_motif_density.txt"

perl HOMER/bin/annotatePeaks.pl $query hg38 -size 1000 -hist 10 -m ${homer_motifs}oct11.motif ${homer_motifs}irf1.motif ${homer_motifs}pax5.motif ${homer_motifs}ebf.motif ${homer_motifs}spib.motif ${homer_motifs}pu1-irf.motif ${homer_motifs}oct2.motif ${homer_motifs}e2a-near-pu1.motif > $output 

done
