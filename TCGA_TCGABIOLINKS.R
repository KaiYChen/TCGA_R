library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)

clinical <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")
TCGAbiolinks:::getProjectSummary("TCGA-COAD")

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "HTSeq - Counts",
                  legacy=FALSE)
##### Normalized Reads
# query <- GDCquery(project = "TCGA-COAD",
#                   data.category = "Gene expression",
#                   data.type = "Gene expression quantification",
#                   platform = "Illumina HiSeq", file.type  = "normalized_results", 
#                   experimental.strategy = "RNA-Seq",
#                   legacy = TRUE)

GDCdownload(query)
COAD_Rnaseq<- GDCprepare(query)
COAD_Matrix <- assay(COAD_Rnaseq)

#### miRNASeq
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification",
                  legacy=FALSE)
GDCdownload(query)
COAD_MiRnaseq<- GDCprepare(query)
COAD_MiRNA_Matrix <- COAD_MiRnaseq[,c(1,seq(2,dim(COAD_MiRnaseq)[2],3))]
