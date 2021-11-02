library(coxme)
library(data.table)
library(magrittr)
library(parallel)

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


tcga.dataset = fread("~/project/deeplearning/icb/data/tcga/scrna.v4.allgenes.nopcs/dataset.txt")

if.start.inx = which(colnames(tcga.dataset)== "B_cells_naive")
dataset.col = colnames(tcga.dataset)
# if.factors = tcga.dataset[,seq(if.start.inx, ncol(tcga.dataset)),with=F]
surv.factors = tcga.dataset[,length(dataset.col) + seq(8)-8, with =F]
exp.tcga = tcga.dataset[,seq(2,if.start.inx-1), with=F]


run.seq = seq(2, if.start.inx-1)
# run.seq = seq(2, 20) 
genes = colnames(tcga.dataset)[run.seq]
tcga.survival.out = list()

for (ii in seq(4)) {
    sur.type = gsub(".time", colnames(surv.factors)[2*(ii-1)+1], replacement="")
    pfs.response.dt = surv.factors[,seq(2) + 2*(ii-1),with=F] %>%
    set_colnames(c("Survival", "Event")) %>%
    .[,dataset:=tcga.dataset$cancertype]
    pfs.response.genes.dt = mclapply(run.seq, function(tt) eval.coxme(unlist(tcga.dataset[,tt, with=F]), dtx=pfs.response.dt), mc.cores=40) %>% 
    do.call(rbind, .) %>%
    set_colnames(c("estimate","se", "z", "P")) %>%
    data.table() %>%
    .[,genes:=genes]

    tcga.survival.out[[sur.type]]= pfs.response.genes.dt %>%
    # .[!is.na(estimate)] %>%
    # .[order(P)]  %>% 
    .[,effect:=sign(estimate)*ifelse(abs(estimate)-se < 0, 0, abs(estimate)-se)]
}

saveRDS(file = "~/liulab_home/data/tcga/tcga.cox.genes.Rds", tcga.survival.out)