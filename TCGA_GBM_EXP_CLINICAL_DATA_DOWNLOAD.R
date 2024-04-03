getwd()


# Load packages
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library(sesame)
library(DT)
library(GDCquery_clinic)

install.packages("glmnet")
install.packages("factoextra")
install.packages("FactoMineR")
install.packages("caret")
install.packages("survminer")
install.packages("gProfileR")
install.packages("genefilter")
BiocManager::install("genefilter")
BiocManager::install("sesameData")
BiocManager::install("SeSAMe")
install.packages("DT")



library(TCGAbiolinks)



# A total of 2.27 GB
query_exp_lgg <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)

GDCdownload(query_exp_lgg)
exp_lgg <- GDCprepare(
  query = query_exp_lgg
)

query_exp_gbm <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
GDCdownload(query_exp_gbm)
exp_gbm <- GDCprepare(
  query = query_exp_gbm
)

# The following clinical data is not available in GBM
missing_cols <- setdiff(colnames(colData(exp_lgg)),colnames(colData(exp_gbm)))
for(i in missing_cols){
  exp_lgg[[i]] <- NULL
}



library(SummarizedExperiment)

# Load object from TCGAWorkflowData package
# This object will be created in subsequent sections for enhanced clarity and understanding.


# get genes information
genes.info <- rowRanges(exp_gbm)
genes.info

# get sample information
sample.info <- colData(exp_gbm)
datatable(
  data = as.data.frame(sample.info), 
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)

# get indexed clinical patient data for GBM samples
gbm_clin <- GDCquery_clinic(
  project = "TCGA-GBM", 
  type = "Clinical"
)

# get indexed clinical patient data for LGG samples
lgg_clin <- GDCquery_clinic(
  project = "TCGA-LGG", 
  type = "Clinical"
)


gbm_res = getResults(query_exp_gbm) # make results as table
# head(lihc_res) # data of the first 6 patients.
colnames(gbm_res) # columns present in the table

lgg_res = getResults(query_exp_lgg) # make results as table
# head(lihc_res) # data of the first 6 patients.
colnames(lgg_res) # columns present in the table



head(gbm_res$sample_type) # first 6 types of tissue.
head(lgg_res$sample_type) # first 6 types of tissue.

summary(factor(gbm_res$sample_type)) # summary of distinct tissues types present in this study
summary(factor(lgg_res$sample_type)) # summary of distinct tissues types present in this study


query_exp_gbm = GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

query_exp_lgg = GDCquery(
  project = "TCGA-LGG",
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

GDCdownload(query = query_exp_gbm)
GDCdownload(query = query_exp_lgg)


tcga_data_gbm = GDCprepare(query_exp_gbm)
tcga_data_LGG = GDCprepare(query_exp_lgg)

dim(tcga_data_gbm)
dim(tcga_data_LGG)

class(tcga_data_gbm)
# Convert to dataframe
tcga_data_gbm_df <- as.data.frame(assay(tcga_data_gbm))
# Add rownames as a column
tcga_data_gbm_df$Gene <- rownames(tcga_data_gbm_df)
write.csv(tcga_data_gbm_df, "tcga_data_gbm.csv")                 # download exp df
write.csv(query_exp_gbm[[1]][[1]], "tcga_gbm_sampleinfo.csv")    # download sampleinfo df
write.csv(gbm_clin, "gbm_clinical.csv")


class(tcga_data_lgg)
# Convert to dataframe
tcga_data_lgg_df <- as.data.frame(assay(tcga_data_LGG))
# Add rownames as a column
tcga_data_lgg_df$Gene <- rownames(tcga_data_LGG_df)
write.csv(tcga_data_lgg_df, "tcga_data_lgg.csv")                # download exp df
write.csv(query_exp_lgg[[1]][[1]], "tcga_lgg_sampleinfo.csv")  # download sampleinfo df
write.csv(lgg_clin, "lgg_clinical.csv")






