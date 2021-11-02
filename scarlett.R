library(tidyr)
crispr = data.table(readxl::read_excel("../../data/Murad_Screen_Data.xlsx", sheet = 1))
crispr.murad.hits = data.table(readxl::read_excel("../../data/2020-0817.Murad_Top_Hits.xlsx", sheet = 1))
crispr.new =crispr[!(id %in% crispr.murad.hits$id)]
genes.curr = intersect(cor.dt$gene, crispr$id)
comb.dt = cbind(cor.dt[match(genes.curr,gene)], crispr[match(genes.curr, id)])
comb.dt[,crispr.p:=ifelse(`pos|p-value` < 0.05, `pos|p-value`, `neg|p-value`)]
comb.dt[crispr.p < 1E-4 & p < 0.05]

kai.screen = fread("~/liulab_home/data/zzeng/29301958_KaiWWucherpfennig_Science_2018/B16F10_Pmel-1.IFNg//rra.gene_summary.txt")
common.genes = intersect(crispr$id %>% toupper(), kai.screen$id %>% toupper())

kai.screen[,p:=ifelse(`neg|lfc`> 0, `pos|p-value`, `neg|p-value` )][,direction.p:=-sign(`neg|lfc`)*log10(p)]
crispr[,p:=ifelse(`neg|lfc`> 0, `pos|p-value`, `neg|p-value` )][,direction.p:=-sign(`neg|lfc`)*log10(p)]
kai.screen.match =kai.screen[match(common.genes, id %>% toupper)]
crispr.match =crispr[match(common.genes, id %>% toupper)]
crispr.match$kai = kai.screen.match$direction.p
crispr.match.sub = crispr.match[direction.p< -7 | kai < -4]
p.com = ggplot(data=crispr.match, aes(x=direction.p, y=kai))  + 
  geom_point(alpha=0.5) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
  geom_text_repel(
    data = crispr.match.sub,
    aes(x=direction.p, y=kai, label = id),    
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  ylab("Kai")  + 
  xlab("Murad") + 
  theme_bw()

top.thrs = 500
sel = which(
  (kai.screen.match$`neg|rank` <= top.thrs) &
    (crispr.match$`neg|rank` <= top.thrs))

gene.eg.dt = clusterProfiler::bitr(crispr.match[sel]$id, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

kk <- clusterProfiler::enrichKEGG(gene = gene.eg.dt$ENTREZID,
                                  organism     = 'hsa',
                                  minGSSize = 2,
                                  pvalueCutoff = 0.05)
p = clusterProfiler::dotplot(kk, showCategory=30)
kk1 <- clusterProfiler::enrichGO(gene = gene.eg.dt$ENTREZID,
                                 OrgDb    = org.Hs.eg.db,
                                 ont      = c("BP"),
                                 minGSSize = 2,
                                 pvalueCutoff = 0.05)
p1 = clusterProfiler::dotplot(kk1, showCategory=30)
kk2 <- clusterProfiler::enrichGO(gene = gene.eg.dt$ENTREZID,
                                 OrgDb    = org.Hs.eg.db,
                                 ont      = c("MF"),
                                 minGSSize = 2,
                                 pvalueCutoff = 0.05)
p2 = clusterProfiler::dotplot(kk2, showCategory=30)

sel = which(
  (kai.screen.match$`pos|rank` <= top.thrs) &
    (crispr.match$`pos|rank` <= top.thrs))
# plot pathway above

gene.eg.dt = clusterProfiler::bitr(crispr.match[sel]$id, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

kk <- clusterProfiler::enrichKEGG(gene = gene.eg.dt$ENTREZID,
                                  organism     = 'hsa',
                                  minGSSize = 2,
                                  pvalueCutoff = 0.05)
q = clusterProfiler::dotplot(kk, showCategory=30)
kk1 <- clusterProfiler::enrichGO(gene = gene.eg.dt$ENTREZID,
                                 OrgDb    = org.Hs.eg.db,
                                 ont      = c("BP"),
                                 minGSSize = 2,
                                 pvalueCutoff = 0.05)
q1 = clusterProfiler::dotplot(kk1, showCategory=30)
kk2 <- clusterProfiler::enrichGO(gene = gene.eg.dt$ENTREZID,
                                 OrgDb    = org.Hs.eg.db,
                                 ont      = c("MF"),
                                 minGSSize = 2,
                                 pvalueCutoff = 0.05)
q2 = clusterProfiler::dotplot(kk2, showCategory=30)


# highlight those are receptor ligand
# highligh those are receptor ligand of gd-T expressed complements.
complexes = fread("~/liulab_home/softwares/cellphonedb-data/data/complex_input.csv")
genesinp = fread("~/liulab_home/softwares/cellphonedb-data/data/gene_input.csv")
proteins = fread("~/liulab_home/softwares/cellphonedb-data/data/protein_input.csv")
interactions = fread("~/liulab_home/softwares/cellphonedb-data/data/interaction_input.csv")

cytokines = fread("/liulab/xmwang/oxphos_proj/loading_data/surface/cytokine.txt", header=F)
surface.genes = fread("/liulab/xmwang/oxphos_proj/loading_data/surface/ExpressionLigRec.txt", header=T)

crispr.match[,recLig:=ifelse(toupper(id) %in% toupper(genesinp$gene_name), 1, 0)]
crispr.match[,recLig:=ifelse(toupper(id) %in% toupper(cytokines$V2), 1, 0)]
crispr.match[,recLig:=ifelse(toupper(id) %in% toupper(surface.genes$ApprovedSymbol), 1, 0)]
crispr.match.sub1 = crispr.match[abs(direction.p) > 2 | abs(kai) > 3][recLig==1]
recLig.p = ggplot(data=crispr.match, aes(x=direction.p, y=kai))  + 
  geom_point(alpha=0.3, aes(color=recLig)) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
  geom_text_repel(
    data = crispr.match.sub1,
    aes(x=direction.p, y=kai, label = id),    
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  ylab("Kai")  + 
  xlab("Murad") + 
  theme_bw()





#### 
# To scarlett for gdT cells
all.recplig = proteins$protein_name %>% gsub("_HUMAN", ., replacement = "") %>% 
  c(.,surface.genes$ApprovedSymbol, cytokines$V2) %>% toupper %>% unique
aa = crispr[id %in% all.recplig][p<1E-2]
write.table(aa, file="../../results/murad/gdt.surface.protein.txt", quote = F, row.names = F, col.names = T)
memberane= fread("../../data/gene_lists/bloodatlas_membrane.tsv")
secreted= fread("../../data/gene_lists/bloodatlas_secreted.tsv")
membrane.secreted= fread("../../data/gene_lists/bloodatlas_membrane_and_secreted.tsv")
all.membrane.secreted =  c(memberane$Gene, secreted$Gene, membrane.secreted$Gene) %>% toupper %>% unique
bb = crispr[id %in% all.membrane.secreted][p<1E-3]
write.table(bb, file="../../results/murad/gdt.membrane.secreted.protein.txt", quote = F, row.names = F, col.names = T)
