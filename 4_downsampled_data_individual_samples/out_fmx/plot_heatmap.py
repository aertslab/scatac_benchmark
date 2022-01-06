#!/usr/bin/env python3

import argparse
import pandas as pd
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
# sns.set(color_codes=True)

parser = argparse.ArgumentParser(description='')
parser.add_argument(
    "--input",
    type=str,
    help='Input table (non-reference discordance rate from vcf-compare).'
)

parser.add_argument(
    "--output",
    type=str,
    help='Output plot (pdf).'
)
args = parser.parse_args()


#f_nrd = 'genotype_concordance_nrd.txt'
#f_nvar = 'genotype_concordance_nvar.txt'

nrd = pd.read_csv(args.input, sep='\t', header=0, index_col=0)
#nvar = pd.read_csv(f_nvar, sep='\t', header=0, index_col=0)

g = sns.clustermap(nrd, annot=True,  square=True,  linecolor='gray',
                   yticklabels=True, xticklabels=True, 
                   vmin=0.0, vmax=nrd.max().max(), # row_colors=colormap,
                   cmap=sns.color_palette('Reds'),
                   #cmap='vlag',
                   cbar_kws={'label': 'Non-ref. discord rate (%)'},
                   figsize=(8,8) )
g.cax.set_visible(True)
g.savefig(args.output, dpi=300, bbox_inches = "tight")






