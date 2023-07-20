## Get Genelist from the paper: COVID-19 Transcriptomic Atlas: A Comprehensive Analysis of COVID-19 Related Transcriptomics Datasets

genelist="IFI6,OAS3,OAS1,CRTAM,TNFRSF9,CCR7,IFIT3,GBP1, CRTAM, OAS1, IFI44L, IFI6, OAS3, IFIT2, CXCL11, MX1, ISG15, IFIT1, CXCL10, OAS2, TNFSF18, TNFRSF9, LAG3, HLA-G, MS4A1, CD274, GZMK, NCR3, FCGR2B, IL10, MS4A1, HLA-DOB, CCR7, CRTAM,CRTAM, IFI44L, IFIT3, IFI6, OAS1, OAS3, NCR3, MX1, GBP1, IL10, IFIT2, TNFSF18, ISG15, IFIT1, CXCL11, CXCL10, TNFRSF9, HLA-DPA1, MS4A1, CRTAM, IFI27, CXCL11, HLA-DOB, IFI6, OAS3, ISG15, IFIT2, MX1, OAS1, IFIT3, IFIT1, TNFSF14, CRTAM, IL10, IDO1, GBP1, IFI44L, CXCL10, IL1A, IL23A,HSPA14, GBP5, PATL2, FHOD3, PHF11,ALPK3, MEGF6, H2AW"
genelist=gsub(" ", "", genelist, fixed = TRUE)
genelist_final=strsplit(genelist,",")[[1]]
genelist_final=unique(genelist_final)
rm(genelist)

#read level 5 data
setwd("/Users/taoxu/Desktop/COVID19PAPER")
library(cmapR)
my_ds_10_columns <- parse_gctx("Level5.gctx", cid=1:10)

#rm(my_ds_10_columns)

gene_info=read.csv("Geneinfo.csv", header=TRUE,sep=",")

## Pick those in genelist

gene_id_selected=gene_info$pr_gene_id[which(gene_info$pr_gene_symbol %in% genelist_final)]

## 
A=which(my_ds_10_columns@rid %in%gene_id_selected)
final_data= parse_gctx("Level5.gctx", rid=A)

cellinfo=read.csv("CellInfo.csv",header=TRUE)
InstInfo=read.csv("InstInfo.csv",header=TRUE)
PertInfo=read.csv("PertInfo.csv",header=TRUE)
PertMetric=read.csv("PertMetric.csv",header=TRUE)
SigInfo=read.csv("Siginfo.csv",header=TRUE)
SigMetric=read.csv("SigMetric.csv",header=TRUE)
## Pert ID: small molecule identifier

ruxo=which(final_data@cid %in% SigInfo[which(SigInfo$pert_iname=="ruxolitinib"),1])
riba=which(final_data@cid %in% SigInfo[which(SigInfo$pert_iname=="ribavirin"),1])
chlo=which(final_data@cid %in% SigInfo[which(SigInfo$pert_iname=="chloroquine"),1])

final_data=parse_gctx("Level5.gctx",rid=A,cid=c(ruxo,riba,chlo))





