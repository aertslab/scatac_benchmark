

singularity run --cleanenv -H $PWD -B /lustre1,/staging,${VSC_SCRATCH}/tmp:/tmp /staging/leuven/stg_00002/lcb/cflerin/containers/cflerin-seurat-4.0.3-plus.sif R

################################################################################

library(Seurat)
library(SeuratDisk)
library(Signac)
library(ggplot2)



pbmc.rna = readRDS('/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/seurat/reference/pbmc_ssc_mat__integrated.rds')


annotations = readRDS("granges_annotation.rds")
pred_thr = 0.7


sample_ids = c(
    'Broad_1', 'Broad_2', 
    'Broad_mito_1', 'Broad_mito_2', 'CNAG_1', 'CNAG_2', 's3atac', 'Sanger_1', 'Sanger_2', 'Stanford_1', 'Stanford_2', 'VIB_1', 'VIB_2')

f_fragdir = '/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/data_freeze_Jun2021/atac_preprocess/multiplet_tagged/fragments/'


# with screen regions
#f_loomdir = '/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/data_freeze_Jun2021/atac_qc_multiplet_merged/jupyter/cell_region_loom__screen/'
#f_outdir = 'predictions__screen/'


# with consensus regions
f_loomdir = '/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/data_freeze_Jun2021/atac_qc_multiplet_merged/jupyter/cell_region_loom__consensus/'

f_outdir = 'predictions__consensus/'

################################################################################
################################################################################
for(sample_id in sample_ids) {

print(sample_id)
f_loom = paste0(f_loomdir, sample_id, 'cell_region-all.loom')
f_frag = paste0(f_fragdir, sample_id, '.sinto.mm.fragments.tsv.gz')

### get data from loom:
atacloomcon <- Connect(filename = f_loom, mode = "r")
atacloomcon
atac_tmp <- as.Seurat(atacloomcon, assay='ATAC')
atacloomcon$close_all()

# correctly parse regions (default delims are '-','-')
# create chromatin assay
chromatinassay = CreateChromatinAssay(
    counts=GetAssayData(atac_tmp, slot = "counts", assay='ATAC'),
    genome='hg38',
    fragments = f_frag,
    #ranges=regions,
    sep=c(':','-'),
    validate.fragments=FALSE
    )

# findOverlaps(sort(chromatinassay@ranges))


atac <- CreateSeuratObject(counts = chromatinassay, assay='ATAC')
Annotation(atac) <- annotations

# We exclude the first dimension as this is typically correlated with sequencing depth
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


##################################################
# Identify RNA-ATAC anchors

# quantify gene activity
gene.activities <- GeneActivity(atac, features = VariableFeatures(pbmc.rna))

# add gene activities as a new assay
atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac)
atac <- ScaleData(atac, features = rownames(atac))

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = pbmc.rna,
    query = atac,
    features = VariableFeatures(object = pbmc.rna),
    reference.assay = "RNA",
    query.assay = "ACTIVITY",
    reduction = "cca")


# predict celltype
celltype.predictions <- TransferData(
    anchorset = transfer.anchors,
    refdata = pbmc.rna$CellType,
    weight.reduction = atac[["lsi"]],
    dims = 2:30)

#pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)

md = celltype.predictions

tmp = data.frame(
      composite_sample_id = paste0(rownames(md),'-',sample_id),
      barcode = rownames(md),
      sample_id = sample_id,
      cell_type = md$predicted.id,
      cell_type_pred_score = md$prediction.score.max
      )
tmp$cell_type_hiconf_70 = tmp$cell_type
tmp$cell_type_hiconf_70[tmp$cell_type_pred_score<pred_thr] = 'Unknown'


write.table(tmp,
            file=paste0(f_outdir,sample_id,'__cell_type_seurat.txt'),
            sep='\t', row.names=FALSE, quote=FALSE
            )

}

#table(tmp$cell_type)
#table(tmp$cell_type_hiconf_70)
