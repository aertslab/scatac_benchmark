


singularity exec --cleanenv -H $PWD -B /lustre1,/staging /staging/leuven/res_00001/software/vsn_containers/vibsinglecellnf-scanpy-0.5.2.img ipython

import scanpy as sc
import pandas as pd

f_loom = '/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/Sanger_multiome/out/loom/Sanger_multiome_RNA.SCope_output.loom'
adata = sc.read_loom(f_loom, validate=False)

f_adata = '/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/Sanger_multiome/out/data/Sanger_multiome_RNA.HARMONY.h5ad'
adata2 = sc.read_h5ad(f_adata)


exprmat = pd.DataFrame(adata.X.todense(), index=adata.obs_names, columns=adata.var_names)
exprmat.to_csv('/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/Sanger_multiome_RNA.SCope_output.tsv.gz',
               sep='\t')



pd.DataFrame(adata2.obsm['X_tsne'], index=adata.obs_names, columns=['X','Y']).to_csv(
    '/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/Sanger_multiome_RNA.X_tsne.tsv.gz', sep='\t')

pd.DataFrame(adata2.obsm['X_umap'], index=adata.obs_names, columns=['X','Y']).to_csv(
    '/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/Sanger_multiome_RNA.X_umap.tsv.gz', sep='\t')






################################################################################
################################################################################

singularity exec -B /lustre1,/staging /staging/leuven/stg_00002/lcb/cflerin/containers/monocle-garnett.sif R

library(garnett)

# load classifier:
#load("/staging/leuven/stg_00002/lcb/cflerin/resources/garnett/hsPBMC")

hsPBMC = readRDS("/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/garnett/hsPBMC_20190911.RDS")

countMat = read.delim('/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/Sanger_multiome_RNA.SCope_output.tsv.gz',
                       sep="\t", row.names=1 )

.embd = '/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/'

l = list(
    sctsne = 'Sanger_multiome_RNA.X_tsne.tsv.gz',
    scumap = 'Sanger_multiome_RNA.X_umap.tsv.gz'
)
emb = lapply( l, function(x) read.delim(paste0(.embd,x), sep="\t", stringsAsFactors=FALSE ))


pd = new( "AnnotatedDataFrame",
        data = data.frame(
            sampleNames = rownames(countMat) ,
            #CellType = sapply(strsplit(rownames(countMat),"_"),function(x) paste(x[1:(length(x)-1)],collapse="_") ) ,
            row.names=rownames(countMat) ,
            # scanpy
            scumap_X = emb$scumap[,2] ,
            scumap_Y = emb$scumap[,3] ,
            sctsne_X = emb$sctsne[,2] ,
            sctsne_Y = emb$sctsne[,3] ,
            stringsAsFactors=FALSE
            )
        )

fd = new( "AnnotatedDataFrame",
        data = data.frame(
            gene_short_name = colnames(countMat) ,
            row.names=colnames(countMat) ,
            stringsAsFactors=FALSE
            )
        )

pbmc <- newCellDataSet(
        as( t(as.matrix(countMat)), "sparseMatrix" ) ,
        phenoData = pd,
        featureData = fd ,
        expressionFamily = negbinomial.size()
        )

# generate size factors for normalization later
pbmc <- estimateSizeFactors(pbmc)

# classify
library(org.Hs.eg.db)
pbmc <- classify_cells(pbmc, hsPBMC,
        db = org.Hs.eg.db,
        cluster_extend = TRUE,
        cds_gene_id_type = "SYMBOL")


head(pData(pbmc))

table( pData(pbmc)$garnett_cluster )

table(pData(pbmc)$cell_type)

table(pData(pbmc)$cluster_ext_type)


table(pData(pbmc)[,c('cell_type','cluster_ext_type')])


################################################################################
marker_file_path <- system.file("extdata", "pbmc_bad_markers.txt",
                                package = "garnett")
marker_check <- check_markers(pbmc, marker_file_path,
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")

plot_markers(marker_check)

################################################################################

.wd = '/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/garnett/'

write.table( pData(pbmc), file=paste0(.wd,"pbmc_garnett_results_full.txt"), sep="\t", quote=FALSE )

cellannot=pData(pbmc)[,c('sampleNames','cell_type','cluster_ext_type')]
colnames(cellannot) = c('composite_sample_id', 'cell_type','cell_type_extended')

cellannot$barcode = sapply(strsplit(cellannot$composite_sample_id,"___"), '[', 1)
cellannot$sample_id = sapply(strsplit(cellannot$composite_sample_id,"___"), '[', 2)
write.table(cellannot, file=paste0(.wd,"pbmc_garnett_results__cell_annotations.txt"), sep="\t", quote=FALSE, row.names=FALSE )














################################################################################

devtool::install_github("aertslab/SCopeLoomR", lib="~/R")



library(SCopeLoomR)


loom_path <- file.path(system.file('extdata', package='SCopeLoomR'), "example.loom")
loom <- open_loom(loom_path, mode="r+")
dgem <- get_dgem(loom)
close_loom(loom)

dgem[1:5,1:5]



