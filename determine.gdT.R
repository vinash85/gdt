setwd("~/project/gdt")
gs.enrichment.score = function(exp.count, genes.set){
	colSums(exp.count[rownames(exp.count) %in% genes.set, ])/colSums(exp.count)
}
norm.gs.enrichment.score = function(gs.score){
	(gs.score - min(gs.score))/(max(gs.score) - min(gs.score))
}

gdT.gs1 = c("CD3D", "CD3E", "TRDC", "TRGC1", "TRGC2")
gdT.gs2= c("CD8A", "CD8B")

trg.genes = grep("^TRG[V,D,J, C]", rownames(pbmc.combined), value=T)
trd.genes = grep("^TRD[V,D,J,C]", rownames(pbmc.combined), value=T)
get.gdT.score  = function(sco) {
	gdT.gs1 = c("CD3D", "CD3E", "TRDC", "TRGC1", "TRGC2")
	gdT.gs2= c("CD8A", "CD8B")
	sco$gdt.gs1  = gs.enrichment.score(sco@assays$RNA@counts, gdT.gs1) 

	sco$gdt.gs2 = (-1 * gs.enrichment.score(sco@assays$RNA@counts, gdT.gs2))

	sco$gdt = norm.gs.enrichment.score(sco$gdt.gs1) * norm.gs.enrichment.score(sco$gdt.gs2)
	sco
}


library(Seurat)
library(data.table)
library(magrittr)
library(Matrix)
exp.mat  = fread("~/liulab_home/data/single_cell/GSE128223/data/GSE128223_3donors_d1_d2_v2.tsv.gz")
row.name = exp.mat$V1
exp.mat = exp.mat[,V1:=NULL] %>%  as.matrix %>%
set_rownames(row.name)
gdt.pnas <- CreateSeuratObject(counts = exp.mat, project = "PNAS")
gdt.pnas %<>% preprocessing.std.seurat

norm.exp = exp(gdt.pnas@assays$RNA@data) -1 
gs1.score = gs.enrichment.score(norm.exp, gdT.gs1) %>%
norm.gs.enrichment.score

gs2.score = (-1 * gs.enrichment.score(norm.exp, gdT.gs2)) %>%
norm.gs.enrichment.score



## pbmc datasets analysis

library(Matrix)
library(Seurat)
library(data.table)
library(magrittr)
library(cowplot)
library(harmony)
library(uwot)
library(parallel)
library(ggplot2)

pbmc4k.data <- Read10X(data.dir = "/homes/cwang/projects/DATA/SCRNAseq/Data_normal/10XHCA_PBMC/Data/pbmc4k/outs/filtered_gene_bc_matrices/GRCh38/")
pbmc4k <- CreateSeuratObject(counts = pbmc4k.data, project = "PBMC4K")


pbmc8k.data <- Read10X(data.dir = "/homes/cwang/projects/DATA/SCRNAseq/Data_normal/10XHCA_PBMC/Data/pbmc8k/outs/filtered_gene_bc_matrices/GRCh38/")
pbmc8k <- CreateSeuratObject(counts = pbmc8k.data, project = "PBMC8K")




pbmc.combined <- merge(pbmc4k, y = pbmc8k, add.cell.ids = c("4K", "8K"), project = "PBMC12K")

pbmc.combined$gdt.gs1  = gs.enrichment.score(exp(pbmc.combined@assays$RNA@data) -1, gdT.gs1) 

pbmc.combined$gdt.gs2 = (-1 * gs.enrichment.score(exp(pbmc.combined@assays$RNA@data)-1, gdT.gs2))

pbmc.combined$gdt = norm.gs.enrichment.score(pbmc.combined$gdt.gs1) * norm.gs.enrichment.score(pbmc.combined$gdt.gs2)

avg.exp <- function(sco, pattern){
 tryCatch({
aa = exp(sco@assays$RNA@data[grep(pattern, rownames(sco)),, drop=F]) - 1  
 if(ncol(aa) >1) aa= colSums(aa)
 	aa
}, error = function(e) NULL) 
} 
for (gd in c("TRG", "TRD", "TRA", "TRB")){
	for (type in c("V", "D", "J", "C")){
		label = paste0(gd,type)
		out = avg.exp(pbmc.combined, label)
		if(!is.null(out))
			pbmc.combined %<>% AddMetaData(metadata=out,
				col.name =label)
	}
}

pbmc.combined  %<>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)%>%
    RunUMAP(reduction = "pca", dims = 1:20) %>% 
    FindNeighbors(reduction = "pca", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()

pbmc.combined  %<>%  RunTSNE(reduction = "pca", dims = 1:20)


   pbmc.combined  %<>%  FindNeighbors(reduction = "tsne",  dims = 1:2) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()

FeaturePlot(pbmc.combined, reduction = "tsne", features= c("gdt", "gdt.gs1", "gdt.gs2")) %>% 
ggsave("results/.figs/pbmc.combined.tsne.gdT.gs.score.pdf", plot=.)

FeaturePlot(pbmc.combined, reduction = "tsne", features= grep("^TR", colnames(pbmc.combined@meta.data), value=T)) %>% 
ggsave("results/.figs/pbmc.combined.tsne.TRXs.pdf", plot=.,width=16, height=10)

FeaturePlot(pbmc.combined, features = c("gdt", "TRDC", "TRGC1", "TRGC2","CD3E" ,  "CD8A", "GNLY", "NKG7", "PRF1", "GZMB", "GZMK", "GZMH")) %>% 
ggsave("results/.figs/pbmc.combined.umap.features.pdf", plot=., width=16, height=10)

FeaturePlot(pbmc.combined, reduction = "tsne", features= "gdt") %>% 
ggsave("results/.figs/pbmc.combined.tsne.gdT.gs.score.pdf", plot=.)

FeaturePlot(pbmc.combined, reduction = "tsne",features = c("gdt", "TRDC", "TRGC1", "TRGC2","CD3E" ,  "CD8A", "GNLY", "NKG7", "PRF1", "GZMB", "GZMK", "GZMH")) %>% 
ggsave("results/.figs/pbmc.combined.tsne.features.pdf", plot=., width=16, height=10)
 
 p =DimPlot(object = pbmc.combined, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend() 
 ggsave("results/.figs/pbmc.combined.tsne.cluster.pdf", plot=p)

pbmc.combined$marked.cluster = ifelse(pbmc.combined$seurat_clusters %in% c(12,23,11, 32, 19 , 13, 7), pbmc.combined$seurat_clusters, 0)
 p = ggplot(pbmc.combined@meta.data, aes(gdt.gs1, gdt.gs2)) + geom_point() + facet_wrap(~marked.cluster)
  ggsave("results/.figs/pbmc.combined.gs12.pdf", plot=p)





pbmc.tgd.combined <- merge(pbmc.combined[rownames(pbmc.combined) %in% rownames(gdt.pnas),], y = gdt.pnas, add.cell.ids = c("12k", "pnas"), project = "PBMCTGD")

pbmc.tgd.combined$gdt.gs1  = gs.enrichment.score(pbmc.tgd.combined@assays$RNA@counts, gdT.gs1) 

pbmc.tgd.combined$gdt.gs2 = (-1 * gs.enrichment.score(pbmc.tgd.combined@assays$RNA@counts, gdT.gs2))

pbmc.tgd.combined$gdt = norm.gs.enrichment.score(pbmc.tgd.combined$gdt.gs1) * norm.gs.enrichment.score(pbmc.tgd.combined$gdt.gs2)


pbmc.tgd.combined  %<>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)%>%
    RunUMAP(reduction = "pca", dims = 1:20) %>% 
    FindNeighbors(reduction = "pca", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()

FeaturePlot(pbmc.tgd.combined, reduction = "umap", features= "gdt") %>% 
ggsave("results/.figs/pbmc.tgd.combined.umap.gdT.gs.score.pdf", plot=.)

FeaturePlot(pbmc.tgd.combined, features = c("gdt", "TRDC", "TRGC1", "TRGC2","CD3E" ,  "CD8A", "GNLY", "NKG7", "PRF1", "GZMB", "GZMK", "GZMH")) %>% 
ggsave("results/.figs/pbmc.tgd.combined.umap.features.pdf", plot=., width=16, height=10)









##

## Shipp's data 

lee.sco = readRDS("~/liulab_home/data/single_cell/Lee_data/lee.seurat.cd3.RDS")
lee.sco %<>% get.gdT.score
FeaturePlot(lee.sco, reduction = "umap", features= "gdt") %>% 
ggsave("results/.figs/lee.sco.umap.gdT.gs.score.pdf", plot=.)

FeaturePlot(lee.sco, features = c("gdt", "TRDC", "TRGC1", "TRGC2","CD3E" ,  "CD8A", "GNLY", "NKG7", "PRF1", "GZMB", "GZMK", "GZMH")) %>% 
ggsave("results/.figs/lee.sco.umap.features.pdf", plot=., width=16, height=10)


##GSE145281 

GSE145281 = readRDS(file="~/liulab_home/data/single_cell/GSE145281/seurat.RDS") 




## creating supervised clusters  
pbmc.combined.sel = pbmc.combined[,pbmc.combined$seurat_clusters %in% c(12, 23, 11, 32, 19 , 13, 7, 24, 2, 15)]
pbmc.combined.sel$is.tgd = ifelse(pbmc.combined.sel$seurat_clusters %in% c(12,23), 1, 0)
pbmc.combined.sel %>% FindVariableFeatures(method="vst")
genes.sel = intersect(VariableFeatures(pbmc.combined.sel), VariableFeatures(lee.sco)) %>%
intersect(., rownames(GSE145281))
dataset = pbmc.combined.sel@assays$RNA@data[genes.sel,] %>% t %>%
as.data.table %>% 
.[,type:="PBMC"] %>%
.[,is.tgd:=pbmc.combined.sel$is.tgd] %>%
.[,c("type", genes.sel, "is.tgd"), with=F]
 source("~/liulab_home/softwares/avinash/R/source.deepImmune.R")

output.dir = "/liulab/asahu/projects/icb/data/gdt/pbmc.gdt2_t_nk"
write.dataset(output.dir =output.dir, dataset = dataset, sample.name = colnames(pbmc.combined.sel))

file.copy("/liulab/asahu/projects/icb/data/gdt/pbmc.gdt_t_nk/params.json", output.dir)
file.copy("/liulab/asahu/projects/icb/data/gdt/pbmc.gdt_t_nk/datasets_tsne_list.txt", output.dir)
file.copy("/liulab/asahu/projects/icb/data/gdt/pbmc.gdt_t_nk/datasets_test_list.txt", output.dir)

## Shell script 
# python train.py  --data_dir  ~/project/deeplearning/icb/data/gdt/pbmc.gdt2_t_nk/datasets_tsne_list.txt --model_dir ~/project/deeplearning/icb/data/gdt/pbmc.gdt2_t_nk/. 
## 


output.dir1 = "/liulab/asahu/projects/icb/data/gdt/pbmc.gdt2_t_nk/GSE145281"
dataset = GSE145281@assays$RNA@data[genes.sel,] %>% t %>%
as.data.table %>% 
.[,type:="PBMC"] %>%
.[,is.tgd:=0] %>%
.[,c("type", genes.sel, "is.tgd"), with=F]
write.dataset(output.dir = output.dir1, dataset = dataset, sample.name = colnames(GSE145281), training=F)

output.dir1 = "/liulab/asahu/projects/icb/data/gdt/pbmc.gdt2_t_nk/lee.sco"
dataset = lee.sco@assays$RNA@data[genes.sel,] %>% t %>%
as.data.table %>% 
.[,type:="PBMC"] %>%
.[,is.tgd:=0] %>%
.[,c("type", genes.sel, "is.tgd"), with=F]
write.dataset(output.dir = output.dir1, dataset = dataset, sample.name = colnames(lee.sco), training=F)

# python evaluate.py  --data_dir  ~/project/deeplearning/icb/data/gdt/pbmc.gdt2_t_nk/datasets_test_list.txt --model_dir ~/project/deeplearning/icb/data/gdt/pbmc.gdt2_t_nk/tensorboardLog/20200606-173423/. --restore_file  ~/project/deeplearning/icb/data/gdt/pbmc.gdt2_t_nk/tensorboardLog/20200606-173423/epoch-34.pth.tar --output_dir  ~/project/deeplearning/icb/data/gdt/pbmc.gdt2_t_nk/tensorboardLog/20200606-173423/outputs/



## read dataset 
gdt.out = fread("/liulab/asahu/projects/icb/data/gdt/pbmc.gdt2_t_nk/tensorboardLog/20200606-173423/outputs/val_1_prediction.csv")
# samp.names = fread("/liulab/asahu/projects/icb/data/gdt/pbmc.gdt2_t_nk/GSE145281/samples_name.txt")

GSE145281$gdt.di = gdt.out$is.tgd.output

FeaturePlot(GSE145281, features = c("gdt.di", "TRDC", "TRGC1", "TRGC2","CD3E" ,  "CD8A", "GNLY", "NKG7", "PRF1", "GZMB", "GZMK", "GZMH"), reduction="umap") %>% 
ggsave("results/.figs/GSE145281.umap.features.pdf", plot=., width=16, height=10)


## read dataset 

gdt.out = fread("/liulab/asahu/projects/icb/data/gdt/pbmc.gdt2_t_nk/tensorboardLog/20200606-173423/outputs/val_2_prediction.csv")

lee.sco$gdt.di = gdt.out$is.tgd.output

FeaturePlot(lee.sco, features = c("gdt.di", "gdt", "TRDC", "TRGC1", "TRGC2","CD3E" , "CD3D",  "CD8A", "GNLY", "NKG7", "PRF1", "GZMB", "GZMK", "GZMH"), reduction="umap") %>% 
ggsave("results/.figs/lee.sco.umap.features.pdf", plot=., width=16, height=10)


avg.exp <- function(sco, pattern){
 tryCatch({
aa = exp(sco@assays$RNA@data[grep(pattern, rownames(sco)),, drop=F]) - 1  
 if(ncol(aa) >1) aa= colSums(aa)
 	aa
}, error = function(e) NULL) 
} 
for (gd in c("TRG", "TRD", "TRA", "TRB")){
	for (type in c("V", "D", "J", "C")){
		label = paste0(gd,type)
		out = avg.exp(lee.sco, label)
		if(!is.null(out))
			lee.sco %<>% AddMetaData(metadata=out,
				col.name =label)
	}
}


FeaturePlot(lee.sco, reduction = "umap", features= sort(grep("^TR", colnames(lee.sco@meta.data), value=T))) %>% 
ggsave("results/.figs/lee.sco.umap.TRXs.pdf", plot=.,width=16, height=10)

sort(grep("TRG[V,C]", rownames(lee.sco), value=T))

FeaturePlot(lee.sco, reduction = "umap", features= sort(grep("TRG[V,C]", rownames(lee.sco), value=T))) %>% 
ggsave("results/.figs/lee.sco.umap.TRG_VC.pdf", plot=.,width=16, height=10)
