
#!/bin/bash

# Adapt the input file to a HOMER format file

# ------------
# HOMER format
# ------------
# Column1: Unique Peak ID
# Column2: chromosome
# Column3: starting position
# Column4: ending position
# Column5: Strand (+/- or 0/1, where 0="+", 1="-")

# We added manually the header to the Hydrop samples and libds__libdsmerged
# header -> chrom	chromStart	chromEnd	Contrast	Log2FC	strand	Adjusted_pval	


cd ../final_dars/

for file in $(ls -d *B_cell__DARs.bed)
#for file in $(ls -d *B_cell__DARs.TOP2K.bed)
#for file in $(ls -d libds__merged__*)

do
    echo $file	
	awk 'BEGIN {OFS="\t"} {{print NR,$1,$2,$3,"+"}}' $file >> "$file.tmp"

done

