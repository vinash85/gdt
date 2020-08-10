##  Analysis of surv genes  
### compare OS PFI genes  
library(ggrepel)
tcga.surv.genes = readRDS("~/liulab_home/data/tcga/tcga.cox.genes.Rds")
avi.dt = tcga.surv.genes$OS
avi.dt$PFI.z = tcga.surv.genes$PFI$z
avi.dt$PFI.P = tcga.surv.genes$PFI$P
avi.dt$OS.z = tcga.surv.genes$OS$z
avi.dt$OS.P = tcga.surv.genes$OS$P

aa = fit.loess(avi.dt$PFI.z, avi.dt$OS.z)
inx.pos = (avi.dt$OS.z - aa$fitted) %>% order(decreasing=T) %>% 
  .[1:40]
inx.neg = (aa$fitted- avi.dt$OS.z ) %>% order(decreasing=T) %>% 
  .[1:30]

avi.dt.sub = avi.dt[unique(c(
  order(avi.dt$PFI.z,decreasing=F)[1:30],
  order(avi.dt$OS.z,decreasing=F)[1:30],
  order(avi.dt$PFI.z,decreasing=T)[1:40],
  order(avi.dt$OS.z,decreasing=T)[1:40], 
  inx.pos,
  inx.neg
))]
p.sc.volcano = ggplot(data=avi.dt, aes(x=OS.z, y=PFI.z))  + 
  geom_point(alpha=0.5) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
  geom_text_repel(
    data = avi.dt.sub,
    aes(y = PFI.z, x = OS.z, label = genes),    
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  ylab("PFI stat")  + 
  xlab("OS stat") + 
  theme_bw()

### enrichments
library(clusterProfiler)
gene.eg.dt = clusterProfiler::bitr(unique(avi.dt$genes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% 
  data.table
gene.eg.dt = gene.eg.dt[match(avi.dt$genes, SYMBOL)]
kk <- clusterProfiler::enrichKEGG(gene = gene.eg.dt[avi.dt$OS.z< -4,]$ENTREZID,
                                    universe= gene.eg.dt$ENTREZID,
                                  minGSSize = 3,
                                  organism     = 'hsa',
                                  pvalueCutoff = 0.05)

p = clusterProfiler::dotplot(kk, showCategory=30)
ggsave(filename = "/liulab/asahu/projects/results/.figs/tcga.os.pos.genes.pathway.pdf",p)
kk <- clusterProfiler::enrichGO(gene = gene.eg.dt[avi.dt$OS.z< -4,]$ENTREZID,
                                universe= gene.eg.dt$ENTREZID,
                                minGSSize = 3,
                                OrgDb    = org.Hs.eg.db,
                                ont      = "BP",
                                pvalueCutoff = 0.05)
p = clusterProfiler::dotplot(kk, showCategory=30)


kk2 <- clusterProfiler::enrichGO(gene = gene.eg.dt[avi.dt$OS.z > 4,]$ENTREZID,
                                universe= gene.eg.dt$ENTREZID,
                                minGSSize = 3,
                                OrgDb    = org.Hs.eg.db,
                                ont      = "BP",
                                pvalueCutoff = 0.05)
p = clusterProfiler::dotplot(kk2, showCategory=30)