

module load Seurat

##################################################

singularity run -B /lustre1,/staging /staging/leuven/stg_00002/lcb/cflerin/containers/satijalab-seurat-4.0.3.sif R

singularity run -B /lustre1,/staging,/usr /staging/leuven/stg_00002/lcb/cflerin/containers/satijalab-seurat-4.0.3.sif R

singularity run --cleanenv -H $PWD -B /lustre1,/staging,${VSC_SCRATCH}/tmp:/tmp /staging/leuven/stg_00002/lcb/cflerin/containers/cflerin-seurat-4.0.3-plus.sif R


################################################################################
# https://satijalab.org/seurat/articles/integration_mapping.html


library(Seurat)
library(ggplot2)
#library(sctransform)

################################################################################
### reference data:
f_ref = '/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/seurat/reference/pbmc_ssc_mat.rds'
f_meta = '/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/seurat/reference/pbmc_ssc_metadata.rds'
pbmc.data <- readRDS(f_ref)
pbmc.metadata <- readRDS(f_meta)

pbmc <- CreateSeuratObject(counts = pbmc.data, meta.data = pbmc.metadata)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200)

pbmc.list <- SplitObject(pbmc, split.by = "Method")
#pbmc.list$bench = qdata

for (i in names(pbmc.list)) {
    #pbmc.list[[i]] <- SCTransform(pbmc.list[[i]], verbose = FALSE)
    pbmc.list[[i]] <- NormalizeData(pbmc.list[[i]])
    pbmc.list[[i]] = FindVariableFeatures(pbmc.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

##################################################
### integrate reference:

reference.list <- pbmc.list[
    c("10x Chromium (v2) A", "10x Chromium (v2) B", "10x Chromium (v3)", "10x Chromium (v2)")
]
#reference.list$query = qdata
#pbmc.ref = pbmc.list$`10x Chromium (v3)`
#pbmc.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

# integrate all available assays in the reference:
pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, dims = 1:30)

pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, dims = 1:30)
DefaultAssay(pbmc.integrated) <- "integrated"

pbmc.integrated <- ScaleData(pbmc.integrated, verbose = FALSE)
pbmc.integrated <- RunPCA(pbmc.integrated, npcs = 30, verbose = FALSE)
pbmc.integrated <- RunUMAP(pbmc.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)


table(pbmc.integrated@meta.data$CellType)

table(pbmc.integrated@meta.data$CellType)/length(pbmc.integrated@meta.data$CellType)

# remove 'Megakaryocyte' and 'Plasmacytoid dendritic cell'

pbmc.integrated = subset(pbmc.integrated, subset = CellType!="Unassigned")
pbmc.integrated = subset(pbmc.integrated, subset = CellType!="Plasmacytoid dendritic cell")
pbmc.integrated = subset(pbmc.integrated, subset = CellType!="Megakaryocyte")

table(pbmc.integrated@meta.data$CellType)

saveRDS(pbmc.integrated, file='/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/seurat/reference/pbmc_ssc_mat__integrated.rds')



##################################################
### query data:

f_gex = list(
Sanger_1 = '/staging/leuven/stg_00002/lcb/lcb_projects/BAP/data_share/data_processed/cellranger_arc_2.0.0/Sanger_1/outs/filtered_feature_bc_matrix.h5',
Sanger_2 = '/staging/leuven/stg_00002/lcb/lcb_projects/BAP/data_share/data_processed/cellranger_arc_2.0.0/Sanger_2/outs/filtered_feature_bc_matrix.h5',
PBMCmix_1 = '/staging/leuven/stg_00002/lcb/lcb_projects/BAP/data_unsorted/scRNA/SCGRES_13/cellranger/PBMCmix_1/outs/filtered_feature_bc_matrix.h5',
PBMCmix_2 = '/staging/leuven/stg_00002/lcb/lcb_projects/BAP/data_unsorted/scRNA/SCGRES_13/cellranger/PBMCmix_2/outs/filtered_feature_bc_matrix.h5',
PBMCmix_3 = '/staging/leuven/stg_00002/lcb/lcb_projects/BAP/data_unsorted/scRNA/SCGRES_13/cellranger/PBMCmix_3/outs/filtered_feature_bc_matrix.h5'
)

query.list = list()
for(i in 1:length(f_gex)) {
    tmp_h5 = Read10X_h5(f_gex[[i]])
    if(is.list(tmp_h5)) {
        query.list[[i]] = CreateSeuratObject(tmp_h5$`Gene Expression`)
    } else {
        query.list[[i]] = CreateSeuratObject(tmp_h5)
    }
    query.list[[i]] = NormalizeData(query.list[[i]])
    query.list[[i]] = FindVariableFeatures(query.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

names(query.list) = names(f_gex)

##################################################
# transfer


for(i in 1:length(query.list)) {
    transfer.anchors <- FindTransferAnchors(reference=pbmc.integrated, query=query.list[[i]], dims=1:30, reference.reduction="pca")
    predictions <- TransferData(anchorset=transfer.anchors, refdata=pbmc.integrated$CellType, dims=1:30)
    query.list[[i]] = AddMetaData(query.list[[i]], metadata = predictions)
}

sapply(query.list, function(x) sum(x$prediction.score.max>0.5))
sapply(query.list, function(x) sum(x$prediction.score.max>0.5)/length(x$prediction.score.max))
sapply(query.list, function(x) ncol(x))

3306 Sanger_1
4079 Sanger_2

get_difference_to_next_prediction_score = function(x) {
    y = x[,grep('prediction.score',colnames(x))]
    xcols = grep('prediction.score',colnames(x))
    xcols = xcols[ 1:(length(xcols)-1) ]
    pred_score_next = numeric(nrow(x))
    for(i in 1:nrow(x)) {
        pred_score_next[i] = sort(x$prediction.score.max[i] - as.numeric(x[i,xcols]))[2]
    }
    return(pred_score_next)
}

pred_thr = 0.7
diff_thr = 0.1
for(i in 1:length(query.list)) {
    query.list[[i]]@meta.data$diff_to_next_pred_score = 
        get_difference_to_next_prediction_score(query.list[[i]]@meta.data)
    pf = (query.list[[i]]$prediction.score.max>pred_thr) & (query.list[[i]]$diff_to_next_pred_score>diff_thr)
    cat(names(query.list)[i],": ", 
        length(query.list[[i]]$prediction.score.max), " | ",
        sum(query.list[[i]]$prediction.score.max>=pred_thr), "> ",pred_thr," | ",
        sum(query.list[[i]]$diff_to_next_pred_score>=diff_thr), "> ",diff_thr," | both:",
        sum(pf), " ", sum(pf)/length(pf), "\n")
}

plot(x$prediction.score.max, pred_score_next)

'predicted.id'
'prediction.score.max'


cell.annot = list()
for(i in 1:length(query.list)) {
    md = query.list[[i]]@meta.data

    tmp = data.frame(
          composite_sample_id = paste0(rownames(md),'___',names(query.list)[i]),
          barcode = rownames(md),
          sample_id = names(query.list)[i],
          cell_type = md$predicted.id,
          cell_type_pred_score = md$prediction.score.max
          )
    tmp$cell_type_hiconf_70 = tmp$cell_type
    tmp$cell_type_hiconf_70[tmp$cell_type_pred_score<pred_thr] = 'Unknown'

    cell.annot[[i]] = tmp
}

# Sanger MO samples:
write.table(
            rbind(cell.annot[[1]],cell.annot[[2]]),
            file='10xMO_cell_type_seurat__filtcelltypes_fullintegratedref.txt',
            sep='\t', row.names=FALSE, quote=FALSE
            )

# CNAG RNA samples:
write.table(
            do.call('rbind',cell.annot[3:5]),
            file='CNAG_RNA_cell_type_seurat__filtcelltypes_fullintegratedref.txt',
            sep='\t', row.names=FALSE, quote=FALSE
            )
    

################################################################################
################################################################################

sapply(cell.annot, function(x) sum(x$cell_type_hiconf=='Unknown'))


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


##################################################
### query data:

countMat = read.delim('/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/Sanger_multiome_RNA.SCope_output.tsv.gz', sep="\t", row.names=1 )

# input matrix is rows=genes, columns=cells
qdata = CreateSeuratObject(counts=t(countMat))

qdata <- NormalizeData(qdata)
qdata = FindVariableFeatures(qdata, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

#qdata <- ScaleData(qdata)
#qdata = RunPCA(qdata, npcs = 30, verbose = FALSE)
#qdata <- RunUMAP(qdata, reduction = "pca", dims = 1:30, verbose = FALSE)

#f_rna1 = '/staging/leuven/stg_00002/lcb/lcb_projects/BAP/data_share/data_processed/cellranger_arc_2.0.0/Sanger_1/outs/filtered_feature_bc_matrix.h5'
#q1 = Read10X_h5(f_ref)

##################################################



##################################################

#pbmc.query = pbmc.list[['bench']]
#pbmc.query = pbmc.list$inDrops
pbmc.query = qdata

pbmc.transfer.anchors <- FindTransferAnchors(reference=pbmc.integrated, query=pbmc.query, dims=1:30, reference.reduction="pca")

predictions <- TransferData(anchorset=pbmc.transfer.anchors, refdata=pbmc.integrated$CellType, dims=1:30)


qdata = AddMetaData(qdata, metadata = predictions)

'predicted.id'
'prediction.score.max'


################################################################################
library(SeuratData)
InstallData("panc8")
data("panc8")
pancreas.list <- SplitObject(panc8, split.by = "tech")
pancreas.list <- pancreas.list[c("celseq", "celseq2", "fluidigmc1", "smartseq2")]

pancreas.query <- pancreas.list[["fluidigmc1"]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query,
    dims = 1:30, reference.reduction = "pca")

##################################################

################################################################################
# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)





### reference data:
f_ref = '/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/seruat/reference/pbmc_ssc_mat.rds'
f_meta = '/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/seruat/reference/pbmc_ssc_metadata.rds'
pbmc.data <- readRDS(f_ref)
pbmc.metadata <- readRDS(f_meta)

pbmc <- CreateSeuratObject(counts = pbmc.data, meta.data = pbmc.metadata)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200)

pbmc.list <- SplitObject(pbmc, split.by = "Method")
#pbmc.list$bench = qdata

for (i in names(pbmc.list)) {
    pbmc.list[[i]] <- SCTransform(pbmc.list[[i]], verbose = FALSE)
    #pbmc.list[[i]] <- NormalizeData(pbmc.list[[i]])
    #pbmc.list[[i]] = FindVariableFeatures(pbmc.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}


### query data:
countMat = read.delim('/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/Sanger_multiome_RNA.SCope_output.tsv.gz',
                       sep="\t", row.names=1 )
qdata = CreateSeuratObject(counts=countMat)


# normalize (apply similar step as reference)
qdata <- SCTransform(qdata, verbose = FALSE)
#qdata = RunPCA(qdata, npcs = 50, verbose = FALSE)

reference = pbmc.list$`10x Chromium (v3)`
reference = RunPCA(reference, npcs = 50, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = reference,
  query = qdata,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)









# example
pancreas.query <- pancreas.list[["fluidigmc1"]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query,
    dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$celltype,
    dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)











##################################################

### reference data:
f_ref = '/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/seruat/reference/pbmc_ssc_mat.rds'
f_meta = '/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis2/scRNA/cell_type_classification/seruat/reference/pbmc_ssc_metadata.rds'
pbmc.data <- readRDS(f_ref)
pbmc.metadata <- readRDS(f_meta)

pbmc <- CreateSeuratObject(counts = pbmc.data, meta.data = pbmc.metadata)

pbmc <- subset(pbmc, subset = nFeature_RNA > 200)
pbmc.list <- SplitObject(pbmc, split.by = "Method")
for (i in names(pbmc.list)) {
    pbmc.list[[i]] <- SCTransform(pbmc.list[[i]], verbose = FALSE)
}
pbmc.features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 3000)
pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = pbmc.features)

# This command returns dataset 5.  We can also specify multiple refs. (i.e. c(5,6))
reference_dataset <- which(names(pbmc.list) == "10x Chromium (v3)")

pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "SCT",
    anchor.features = pbmc.features, reference = reference_dataset)
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT")

pbmc.integrated <- RunPCA(object = pbmc.integrated, verbose = FALSE)
pbmc.integrated <- RunUMAP(object = pbmc.integrated, dims = 1:30)




plots <- DimPlot(pbmc.integrated, group.by = c("Method", "CellType"), combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, 
    byrow = TRUE, override.aes = list(size = 3))))
CombinePlots(plots)








