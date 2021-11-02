setwd("~/liulab_home/data/tcga/TCGA_gdc/tcgabiolinks/")
library(TCGAbiolinks)

projects <- getGDCprojects()
projects.id = grep("TCGA",projects$project_id, value=T)

save.file ="TCGA_HTSeq_FPKM-UQ.rda"
if(!file.exists(save.file)){
query <- GDCquery(project = projects.id,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(query)

project.exp <- GDCprepare(query, 
                      save = TRUE, 
                      summarizedExperiment = TRUE, 
                      save.filename = save.file)
}else{
	load(save.file)
	project.exp = data
	rm(data);gc()

}
saveRDS(file="~/liulab_home/data/tcga/TCGA_gdc/tcgabiolinks/all.tcga.expression.RDS", project.exp)



## since not working separately. 

## protien interaction 

TCGAbiolinks:::getProjectSummary("CPTAC-2")
query <- GDCquery(project = "CPTAC-2",
                  data.category = "Transcriptome Profiling") 
GDCdownload(query)