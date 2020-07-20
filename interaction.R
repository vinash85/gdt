library(coxme)
library(data.table)
library(magrittr)
library(parallel)
library(survival)

extract_coxme_table <- function (mod){
    beta <- mod$coefficients #$fixed is not needed
    nvar <- length(beta)
    nfrail <- nrow(mod$var) - nvar
    se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
    z<- round(beta/se, 2)
    p<- signif(1 - pchisq((beta/se)^2, 1), 2)
    table=data.frame(cbind(beta,se,z,p))
    return(table)
}



eval.coxme = function(dat, dtx){
    tryCatch(
        {
            dtx$col = dat
            dtx = dtx[!is.na(col)]
            dtx = dtx[dataset %in% names(which(table(dtx$dataset) > 10))] 
            aa =  coxme(Surv(Survival, Event) ~ col + (1|dataset), dtx)
            extract_coxme_table(aa)
        },error=  function(e) rep(NA,4))
}



calculate.cox.interaction.pancancer = function(dtx){
    tryCatch(
    {
        dtx = dtx[!is.na(gene)]
        dtx = dtx[dataset %in% names(which(table(dtx$dataset) > 10))]
        dtx.summ =dtx[,.(var(gene)), by=dataset][V1 >1E-10]
        dtx = dtx[dataset %in% dtx.summ$dataset]
        aa =  coxme(Surv(Survival, Event) ~ gene + cell.type + gene*cell.type + (1|dataset), dtx)
        tt = extract_coxme_table(aa)
        tt["gene:cell.type",]
    },
    error=  function(e) rep(NA,4),
    warning = function(e) rep(NA,4)
    )

}




tcga.cor.unique = readRDS("~/liulab_home/data/tcga/tcga.singler.Rds")
tcga.dataset = fread("~/liulab_home/projects/icb/data/tcga/scrna.v4.allgenes.nopcs/dataset.txt")
tcga.sample = fread("~/liulab_home/projects/icb/data/tcga/scrna.v4.allgenes.nopcs/samples_name.txt")
tcga.surv.genes.dt = readRDS("~/liulab_home/data/tcga/tcga.cox.genes.Rds")

if.start.inx = which(colnames(tcga.dataset)== "B_cells_naive")
dataset.col = colnames(tcga.dataset)
surv.factors = tcga.dataset[,length(dataset.col) + seq(8)-8, with =F]
exp.tcga = tcga.dataset[,seq(2,if.start.inx-1), with=F]
exp.tcga = exp.tcga %>% as.matrix
run.seq = seq(2, if.start.inx-1)
genes = colnames(tcga.dataset)[run.seq]

## cell type important for response on its own. 


ii = 1 #OS only 

# for (ii in seq(4)) {
    sur.type = gsub(".time", colnames(surv.factors)[2*(ii-1)+1], replacement="")
    pfs.response.dt = surv.factors[,seq(2) + 2*(ii-1),with=F] %>%
    set_colnames(c("Survival", "Event")) %>%
    .[,dataset:=tcga.dataset$cancertype]
    pfs.response.genes.dt = mclapply(seq(nrow(tcga.cor.unique)), function(tt) eval.coxme(unlist(tcga.cor.unique[tt,pfs.response.dt), mc.cores=60) %>% 
    do.call(rbind, .) %>%
    set_colnames(c("estimate","se", "z", "P")) %>%
    data.table() %>%
    .[,cell.type:=rownames(tcga.cor.unique)]

    tcga.celltype.survival.out = pfs.response.genes.dt %>%
    .[,effect:=sign(estimate)*ifelse(abs(estimate)-se < 0, 0, abs(estimate)-se)]
# }

cell.type.sel = tcga.celltype.survival.out[P<0.001]$cell.type 

cell.types = rownames(tcga.cor.unique) %>% unique

dt.curr = surv.factors[,.(Survival=OS.time,Event=OS.filtered, dataset=tcga.dataset$cancertype)]
library(parallel) 
interaction.out = mclapply(cell.types, function(cc){ 
    dt.curr[,cell.type:=tcga.cor.unique[cc,]]
    out = apply(exp.tcga,2, function(exp.curr) {
        dt.curr[,gene:=exp.curr]%>%
        calculate.cox.interaction.pancancer
    }) %>% 
    do.call(rbind, .) %>%
    set_colnames(c("estimate","se", "z", "P")) 
}, mc.cores = 60) 

saveRDS(file="~/liulab_home/data/tcga/tcga.singler.gene.interaction.RDS", interaction.out)

interaction.mat = interaction.out %>% sapply(., function(tt) tt[,3]) %>% t

enrichment.p = WGCNA::corAndPvalue(t(interaction.mat), exp.gamma.delta, 
             use = "pairwise.complete.obs", 
             alternative = "two.sided",
             method="spearman")

