library(RCTD)
library(Seurat)
library(data.table)
library(magrittr)
#' convert data.table to matrix assumes first column is rownames
#'
#' @param dt  data.table 
#'
#' @return matrix
#' @export
#'
#' @examples
dt2mat <- function(dt) {
  dt[,-1, with=F] %>% as.matrix %>% 
    set_rownames(unlist(dt[,1]) )
}

preprocessing.geomx.seurat = function (sco, num.dim = 50, nfeatures=500) {
  require(Seurat)
    sco %<>% Seurat::NormalizeData(verbose = FALSE) %>% FindVariableFeatures(selection.method = "vst", 
        nfeatures = nfeatures) %>% ScaleData(verbose = FALSE)
    sco %<>% RunPCA(pc.genes = sco@var.features, npcs = num.dim, 
        verbose = FALSE) %>% RunUMAP(reduction = "pca", dims = seq(num.dim)) %>% 
        RunTSNE(reduction = "pca", dims = seq(num.dim))  
        sco
}

transfer.label = function(reference, query, refdata = "seurat_clusters") {
  require(Seurat)
  refdata = reference@meta.data[,refdata]

  anchors <- FindTransferAnchors(reference = reference, query = query,
    reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors, refdata = refdata,
    dims = 1:30)
  predictions
}


get.spatial.obj = function(sco.spatial){
  count.curr = sco.spatial@assays$RNA@counts %>%  round
  coords.curr = sco.spatial@reductions$umap@cell.embeddings %>% as.data.frame
  nUMI <- colSums(count.curr) %>%round# In this case, total counts per pixel is nUMI

### Create SpatialRNA object
  puck <- SpatialRNA(coords.curr, count.curr, nUMI)
  puck
}

#' Get RTCD result
#'
#' <full description>
#'
#' @param sco reference Seurat object
#' @param annotation name of meta.data in sco that should be used as cell annotation
#' @param resultsdir dir to store result 
#' @param width=16 width and height to store matrix
#' @param height=6 
#'
#' @export
#' @return

get.RTCD.results = function(sco, rctd.spatial, annotation, resultsdir, tumor_or_immune, width=16, height=6){
  dir.create(resultsdir)
## create reference
  counts = sco@assays$RNA@counts %>% round
  cell_type = sco@meta.data[[annotation]] %>% as.factor %>% 
  set_names(colnames(counts))
   levels(cell_type) %<>% gsub("/", ., replacement="_")

  p = DimPlot(sco, group.by=annotation, label=T) + NoLegend()
  ggsave(sprintf("%s/dimplot.pdf", resultsdir))

  cell_type.sel = which(table(cell_type) > 25) %>% names 
  sel.inx = which(cell_type %in%  cell_type.sel)
  cell_type.curr =  cell_type[sel.inx] 
  cell_type.curr = factor(cell_type.curr, levels=unique(cell_type.curr))
  counts.curr = counts[,sel.inx]
  reference <- Reference(counts=counts.curr, cell_types=cell_type.curr)
## run RTCD
  myRCTD <- create.RCTD(rctd.spatial, reference, max_cores = 10)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  results <- myRCTD@results
  norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
  cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
  spatialRNA <- myRCTD@spatialRNA

## plots 
try({
avi.dt = data.table(tumor=tumor_or_immune%>% as.factor, mal.weights=norm_weights[,"Malignant"])
avi.dt = avi.dt[!is.na(tumor)]
library(ggpubr)
p = ggboxplot(avi.dt, x="tumor", y = "mal.weights", add = "dotplot")
filename = sprintf("%s/tumor_weights.pdf", resultsdir)
ggsave(file=filename, p)
})
  library(corrplot)
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF",  "#FFFFFF", "#77AADD", "#4477AA") %>% rev)
  ord=tumor_or_immune %>% order(na.last="FALSE")
  label[1] = "control"
  weights.ord = as.matrix(norm_weights)[ord,] %>% t %>% 
  set_colnames(label)
  pdf(sprintf("%s/weights_matrix.pdf", resultsdir), width=width, height=height)
  corrplot(weights.ord, method="color", col=col(200),
   is.corr = F,
   tl.col="black", tl.srt=90, #Text label color and rotation
   tl.cex = 0.8,
   diag=T
   )
  dev.off()

  weights.ord = as.matrix(results$weights)[ord,] %>% t %>% 
  set_colnames(label)
  pdf(sprintf("%s/nonnorm_weights_matrix.pdf", resultsdir), width=width, height=height)
  corrplot(weights.ord, method="color", col=col(200),
   is.corr = F,
   tl.col="black", tl.srt=90, #Text label color and rotation
   tl.cex = 0.8,
   diag=T
   )
  dev.off()
  myRCTD
}


get.sco.label = function(scRNA.file) {
  sco = readRDS(scRNA.file)
  if(is.list(sco)) sco = sco[[1]]
  sco %<>% NormalizeData %>% 
FindVariableFeatures(selection.method = "vst", 
        nfeatures = 3000, verbose = FALSE)
xx = transfer.label(reference=covid.obj, query=sco, refdata="celltype")
sco$predicted.id = xx$predicted.id
sco  
}