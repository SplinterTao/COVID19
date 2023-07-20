library(openxlsx)
## Rank p-value of upregulated genes
gene_source=read.xlsx("Table_3.xlsx",colNames=TRUE)
upregulated=gene_source[,c(1,6)]
colnames(upregulated)=c("symbol","pval")
up=upregulated[-1,]
rownames(up)=1:nrow(up)
pval=up[,2]
pval=as.numeric(pval)
up[,2]=pval
up=up[1:475,]
up=up[order(up$pval),]


## Rank p-value of downregulated genes
gene_source=read.xlsx("Table_3.xlsx",colNames=TRUE)
downregulated=gene_source[,c(10,15)]
colnames(downregulated)=c("symbol","pval")
down=downregulated[-1,]
rownames(down)=1:nrow(down)
pval=down[,2]
pval=as.numeric(pval)
down[,2]=pval
down=down[order(down$pval),]

## Select genes with the least p-value
library(cmapR)

my_ds_10_columns <- parse_gctx("Level5.gctx", cid=1:10)
up_selected=up[which(up$pval<1e-3),]
down_selected=down[which(down$pval<2e-4),]
gene_info=read.csv("Geneinfo.csv", header=TRUE,sep=",")

gene_up_selected=gene_info$pr_gene_id[which(gene_info$pr_gene_symbol %in% up_selected[,1])]
gene_down_selected=gene_info$pr_gene_id[which(gene_info$pr_gene_symbol %in% down_selected[,1])]

gene_up_selected_name=gene_info$pr_gene_symbol[which(gene_info$pr_gene_symbol %in% up_selected[,1])]
gene_down_selected_name=gene_info$pr_gene_symbol[which(gene_info$pr_gene_symbol %in% down_selected[,1])]

up_data=which(my_ds_10_columns@rid %in%gene_up_selected)
down_data=which(my_ds_10_columns@rid %in%gene_down_selected)



final_data_up= parse_gctx("Level5.gctx", rid=up_data)
final_data_down= parse_gctx("Level5.gctx", rid=down_data)


## select data from
SigInfo=read.csv("Siginfo.csv",header=TRUE)
sort(table(SigInfo$pert_iname),decreasing=TRUE)
druglist=c("vorinostat","trichostatin-a","wortmannin","geldanamycin","sirolimus",
           "curcumin","estradiol","fulvestrant","tozasertib","sulforaphane","tamoxifen","raloxifene",
           "withaferin-a","dexamethasone","tretinoin","thioridazine","tanespimycin",
           "troglitazone","olaparib","parthenolide","panobinostat","resveratrol",
           "selumetinib","gemcitabine","doxorubicin","tert-butylhydroquinone","pifithrin-mu","manumycin-a",
           "veliparib","cyclosporin-a","gefitinib","emetine","forskolin","vemurafenib",
           "flutamide","niclosamide","toremifene","temozolomide","trifluoperazine","sildenafil",
           "bucladesine","entinostat","daunorubicin","menadione","tacedinaline")
druglist_main=c("chloroquine","ribavirin","ruxolitinib","ritonavir")

default_up=which(final_data_up@cid %in% SigInfo[which(SigInfo$pert_iname=="UnTrt"),1])
default_down=which(final_data_down@cid %in% SigInfo[which(SigInfo$pert_iname=="UnTrt"),1])
default_n=length(default_up)
default_up_data=parse_gctx("Level5.gctx",rid=up_data,cid=default_up)
default_down_data=parse_gctx("Level5.gctx",rid=down_data,cid=default_down)

for (i in druglist_main){
  temp_up=which(final_data_up@cid %in% SigInfo[which(SigInfo$pert_iname==i),1])
  temp_down=which(final_data_down@cid %in% SigInfo[which(SigInfo$pert_iname==i),1])
  temp_n=length(temp_up)
  temp_up_data=parse_gctx("Level5.gctx",rid=up_data,cid=temp_up)
  temp_down_data=parse_gctx("Level5.gctx",rid=down_data,cid=temp_down)
  data_full_temp_up=cbind(t(cbind(temp_up_data@mat,default_up_data@mat)),c(rep(1,temp_n),rep(0,default_n)))
  data_full_temp_down=cbind(t(cbind(temp_down_data@mat,default_down_data@mat)),c(rep(1,temp_n),rep(0,default_n)))
  write.csv(data_full_temp_up,paste("./alldata/COVID19Drug/",i,"up.csv",sep=""))
  write.csv(data_full_temp_down,paste("./alldata/COVID19Drug/",i,"down.csv",sep=""))
  print(i)
  print(ncol(temp_up_data@mat))
  
}

#### After running the python code please run this:
files=list.files("./Result", pattern="*.csv")
read_in=read.csv(paste("./Result/",file,sep=""))
read_in
for (file in files){
  read_in=read.csv(paste("./Result/",file,sep=""),header=TRUE,row.names=1)
  colnames(read_in)=sub("X","",colnames(read_in))
  for (i in 1:length(colnames(read_in))){
    geneid=colnames(read_in)[i]
    colnames(read_in)[i]=gene_info$pr_gene_symbol[which(gene_info$pr_gene_id==as.numeric(geneid))]
    }
    if (nrow(read_in)==ncol(read_in)){
    rownames(read_in)=colnames(read_in)
    write.csv(read_in,paste("./Result/ResultProcessed/Coefficient/",file,sep=""))
    }
    if (nrow(read_in)==1){
      write.csv(read_in,paste("./Result/ResultProcessed/Covariate/",file,sep=""))
    }
}

files=list.files("./Result/ResultProcessed/Covariate", pattern="*.csv")
for (file in files){
  read_in=read.csv(paste("./Result/ResultProcessed/Covariate/",file,sep=""),row.names=1)
  read_in_data=c(t(read_in[1,]))
  mygene=colnames(read_in)[order(abs(read_in_data),decreasing=TRUE)[1:5]]
  write.table(mygene,paste("./Result/ResultProcessed/Covariate/",file,".txt",sep=""))
}

for (i in druglist){
  print(i)
  temp_up=which(final_data_up@cid %in% SigInfo[which(SigInfo$pert_iname==i),1])
  temp_down=which(final_data_down@cid %in% SigInfo[which(SigInfo$pert_iname==i),1])
  temp_n=length(temp_up)
  temp_up_data=parse_gctx("Level5.gctx",rid=up_data,cid=temp_up)
  temp_down_data=parse_gctx("Level5.gctx",rid=down_data,cid=temp_down)
  data_full_temp_up=cbind(t(cbind(temp_up_data@mat,default_up_data@mat)),c(rep(1,temp_n),rep(0,default_n)))
  data_full_temp_down=cbind(t(cbind(temp_down_data@mat,default_down_data@mat)),c(rep(1,temp_n),rep(0,default_n)))
  write.csv(data_full_temp_up,paste("./alldata/DRUGS/",i,"up.csv",sep=""))
  write.csv(data_full_temp_down,paste("./alldata/DRUGS/",i,"down.csv",sep=""))
}
  
data_disease= read.table(file = 'GSE147507_RawReadCounts_Human.tsv', sep = '\t', header = TRUE, row.names=1)
data_disease_up=data_disease[which(rownames(data_disease)%in% gene_up_selected_name),]

data_disease_up_log=log(data_disease_up+1)


data_disease_down=data_disease[which(rownames(data_disease)%in% gene_down_selected_name),]
data_disease_down_log=log(data_disease_down+1)
write.csv(t(data_disease_down_log),"downdiseasenetworkdata.csv")
write.csv(t(data_disease_up_log),"updiseasenetworkdata.csv")

final_data_up_new= parse_gctx("GSE70138.gctx")


## Ritonavir 43+22
## chloroquine 42+127
## ruxolitinib ##129+84
## ribavirin ##42+17









