
datadir="" #provide directory in which data is depositied 
resultsdir="" #provide path to send the results 
RscriptsPath="" #provide path for R scripts 

#source auxiliary functions
source(paste(RscriptsPath,"auxiliary_scripts_onc_chRCC_paper .R",sep=""))

#########################################################
#Processing the GenomeStudio output to produce beta value matrix using minfi
######################################################

library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfiData)
library(sva)

baseDir2="" #Directory that contains the methylation array data 
targets2=read.metharray.sheet(baseDir2)
#this is concatenating long numbers, so I think I need to try do it manually, or else change the directory names

#Force R not to use exponential notation
options(scipen = 999)
setwd(baseDir2)
targets2 <- read.csv("Sample sheet.csv", skip = 10, header = T, stringsAsFactors = F)

library(plyr)
targets2$Sample_Group=revalue(as.factor(targets2$Sample_Group),c("1"="GroupA", "2"="GroupB", "3"="GroupC", "4"="GroupD"))

#Need to make Basename colname, which are the paths to Genome studio files
targets2=targets2[,c("Sample_Name","Sample_Well","Sample_Plate", "Sample_Group","Sentrix_Position","Sentrix_ID")]
colnames(targets2)=c("Sample_Name","Sample_Well","Sample_Plate", "Sample_Group","Array","Slide")

targets2$Basename=rep(NA,nrow(targets2))

targets2$Slide=revalue(as.factor(targets2$Sample_Group), c("GroupA"="200397540076", "GroupB"="200397540082","GroupC"="200397540095","GroupD"="200394970099"))
targets2$Basename=paste(baseDir2, targets2$Slide, "/", targets2$Slide, "_", targets2$Array ,sep="")

RGSet2 <- read.metharray.exp(targets = targets2)

phenoData2 <- pData(RGSet2)
#This is phenotype data

manifest2 <- getManifest(RGSet2)
#This is just a fairly useless manifest (array info)

head(getProbeInfo(manifest2))

MSet <- preprocessIllumina(RGSet2, bg.correct = TRUE, normalize = "controls")
qc <- getQC(MSet)
png(filename="~/QC.png2")
plotQC(qc)
dev.off()

#extract methylation and copy number data 
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
RSet <- updateObject(RSet)

#GenomicMethylSet to GenomicRatioSet
GRset <- mapToGenome(RSet)
snps <- getSnpInfo(GRset)
GRset <- addSnpInfo(GRset)
GRset <- dropLociWithSnps(GRset, snps=c("SBE","CpG"), maf=0)
GRset <- keepSeqlevels(GRset, paste0('chr', 1:22))

beta <- getBeta(GRset)
m <- getM(GRset)
CN <- getCN(GRset)

png(filename="~/DP2.png")
densityPlot(beta)
dev.off()

saveRDS(beta, "datadir/data_betas_kb.rds")
saveRDS(qc, "datadir/data_qc_kb.rds")

###################
#Investigate the proportions of hyermethylated and hypomethylated CpG sites in different TCGA RCC types
#Supplementary figure 1
######################

Barcodes_METcancer_KIRP=colnames(METcancer_KIRP)
Barcodes_METcancer_KIRC=colnames(METcancer_KIRC)
Barcodes_METcancer_KICH=colnames(METcancer_KICH)

dat=load("~/Documents/Projects/Oncocytoma/results/2017-08-30/MethylMix_TCGA_PanCanKidney_AllCpGs.RData")
MethylationStates=MethylMixResults$MethylationStates

#could do venn diagram to see overlap of CpG probes that are frequently hyper or hypo in each cancer type
MixtureStates=MethylMixResults$MixtureStates

HyperProbes=names(MixtureStates[unlist(lapply(1:length(MixtureStates), function(x) max(MixtureStates[[x]])>0.2))])
HypoProbes=names(MixtureStates[unlist(lapply(1:length(MixtureStates), function(x) min(MixtureStates[[x]])<(-0.2)))])
HyperAndHypoProbes=intersect(HyperProbes, HypoProbes)
JustHyperProbes=setdiff(HyperProbes, HyperAndHypoProbes)
JustHypoProbes=setdiff(HypoProbes, HyperAndHypoProbes)

MethylationStatesHyper=MethylationStates[HyperProbes,]
#40247 probes
MethylationStatesHypo=MethylationStates[HypoProbes,]
#1040 probes 

AbnormalProbes=c(JustHyperProbes, JustHypoProbes, HyperAndHypoProbes)
MethylationStatesAbnormal=MethylationStates[AbnormalProbes,]
#source(paste(RscriptsPath, "TCGA_enhancer_scripts.R", sep=""))

#Make version with both hyper and hypo

PropAbnormal=array(NA, c(nrow(MethylationStatesAbnormal),3))
rownames(PropAbnormal)=rownames(MethylationStatesAbnormal)
colnames(PropAbnormal)=c("pRCC/KIRP","ccRCC/KIRC","chRCC/KICH")

for(i in 1:nrow(MethylationStatesAbnormal)){
BarcodesAbnormal=names(which(MethylationStatesAbnormal[i,]!=0))
PropAbnormal[i,1]=length(intersect(BarcodesAbnormal, substr(Barcodes_METcancer_KIRP,1,12)))
PropAbnormal[i,2]=length(intersect(BarcodesAbnormal, substr(Barcodes_METcancer_KIRC,1,12)))
PropAbnormal[i,3]=length(intersect(BarcodesAbnormal, substr(Barcodes_METcancer_KICH,1,12)))
}

ns=c(length(Barcodes_METcancer_KIRP), length(Barcodes_METcancer_KIRC), length(Barcodes_METcancer_KICH))
#get the proportions
PAbAbnormal=sweep(PropAbnormal, MARGIN = 2, ns,"/")
PAbAbnormal=as.data.frame(PAbAbnormal)

#Make violin plots 

HyperMethGenesAll=list()
for(i in 1:ncol(MethylationStates)){
HyperMethGenesAll[[i]]=names(MethylationStates[,i][MethylationStates[,i]>0])
}
names(HyperMethGenesAll)=colnames(MethylationStates)

HypoMethGenesAll=list()
for(i in 1:ncol(MethylationStates)){
  HypoMethGenesAll[[i]]=names(MethylationStates[,i][MethylationStates[,i]<0])
}
names(HypoMethGenesAll)=colnames(MethylationStates)

#MethylationState_genes=rownames(MethylationStates)
Abnormal_MethylMix_state_N=array(NA, c(length(HyperMethGenesAll),2))
rownames(Abnormal_MethylMix_state_N)=names(HyperMethGenesAll)
colnames(Abnormal_MethylMix_state_N)=c("N_Hypermethylated","N_Hypomethylated")

for(i in 1:length(HyperMethGenesAll)){
  Abnormal_MethylMix_state_N[i,1]=length(HyperMethGenesAll[[i]])
  Abnormal_MethylMix_state_N[i,2]=length(HypoMethGenesAll[[i]])
}
Abnormal_MethylMix_state_N=as.data.frame(Abnormal_MethylMix_state_N)

Abnormal_MethylMix_state_N$cancer=rep(NA, nrow(Abnormal_MethylMix_state_N))

Abnormal_MethylMix_state_N[rownames(Abnormal_MethylMix_state_N) %in% gsub("-01","",Barcodes_METcancer_KIRP),"cancer"]="pRCC"
Abnormal_MethylMix_state_N[rownames(Abnormal_MethylMix_state_N) %in% gsub("-01","",Barcodes_METcancer_KIRC),"cancer"]="ccRCC"
Abnormal_MethylMix_state_N[rownames(Abnormal_MethylMix_state_N) %in% gsub("-01","",Barcodes_METcancer_KICH),"cancer"]="chRCC"

nhyper=Abnormal_MethylMix_state_N[,c("N_Hypermethylated","cancer")]
nhypo=Abnormal_MethylMix_state_N[,c("N_Hypomethylated","cancer")]

nhyper$MethylationState=rep("hyper",nrow(nhyper))
colnames(nhyper)[1]="nAbnormal"
nhypo$MethylationState=rep("hypo",nrow(nhypo))
colnames(nhypo)[1]="nAbnormal"

AbnormalStates=rbind(nhyper, nhypo)

file="~/Documents/Projects/oncocytoma/results/2017-08-30/figures/Violin_plot_Hypo_MethylMix_probes_number_hypo_hyper_per_patient"
png(file=paste(file,'.png',sep=''), units="in", width=5, height=5, res=600)
ggplot(AbnormalStates, aes(x=factor(MethylationState), y=nAbnormal, fill =MethylationState)) + geom_violin(trim=FALSE) + facet_grid(. ~ cancer)+ scale_fill_manual(values = c("red", "blue"), name = "Methylation state", labels = c("hyper" = "hypermethylated", "hypo" = "hypomethylated")) + theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y = black.bold.16.text, legend.text=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold")) + theme(strip.text.x = element_text(size = 14, colour = "black", face="bold")) 
dev.off()

###########################
#Get beta value data
###########################

beta=readRDS("datadir/data_betas_kb.rds")
#remove HNSCs
beta=beta[,colnames(beta) %in% rownames(Sample_file)]
#40 samples
#Just removing rows for which any data is missing, sample size to small to impute, and affects downstream analysis
#HNSCs
summary(apply(beta, 2, function(x) length(which(is.na(x)))))
#between 1 and 44 CpG probes missing per sample (mean 10)
NonMissingDataFeatures=which(apply(beta, 1, function(x) length(which(is.na(x))))==0)
#456318 (99.9% of CpG probes for which data was present for all features)
#195 CpG probes missing in any sample
#beta=beta[apply(beta, 1, function(x) length(which(is.na(x))))==0,]
beta=beta[NonMissingDataFeatures,]
beta_scale=scale(beta)
beta_Var=TCGA_GENERIC_GeneFiltering('MAD', beta_scale, 10)

########
#Set up annotation for 450k array
###########

library(plyr)
library(RColorBrewer)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library("gplots")
library("devtools")

#Get probe annotation for 450k array 
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#make a dataframe with annotation to be added to tables
ann450k2=ann450k[,c("UCSC_RefGene_Name", "chr","pos","UCSC_RefGene_Group", "Relation_to_Island")]
ann450k2=as.data.frame(ann450k2)
ann450k2$Gene=lapply(ann450k$UCSC_RefGene_Name, function(x) paste(unique(unlist(strsplit(x,";"))), collapse=", "))
ann450k2$UCSC_RefGene_Group=lapply(ann450k$UCSC_RefGene_Group, function(x) paste(unique(unlist(strsplit(x,";"))), collapse=", "))
ann450k2=ann450k2[,2:ncol(ann450k2)]

####################
#Make sample annotion file
#######################

n=60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

Sample_file=read.table("datadir/Sample_file.txt", sep="\t", header=T)
Sample_file$Sample_Group2=revalue(as.factor(Sample_file$Sample_Group), c("1"="200397540076", "2"="200397540082","3"="200397540095","4"="200394970099"))
rownames(Sample_file)=paste(Sample_file$Sample_Group2, Sample_file$Sentrix_Position, sep="_")
Sample_file=Sample_file[Sample_file$Sample_Type!="Head & neck tumor",]
Sample_file$PatientID=gsub("-.*","",Sample_file$Sample_Name)
Sample_file$Matching_label=factor(Sample_file$PatientID, labels=c(1:length(levels(factor(Sample_file$PatientID)))))
Sample_file$Sample_col=Sample_file$Sample_Type
levels(Sample_file$Sample_col)=col_vector[1:length(levels(Sample_file$Sample_col))]
Sample_file$Sample_Type=factor(Sample_file$Sample_Type)
ids=as.character(Sample_file[Sample_file$Sample_Type=="Normal kidney parenchyma","PatientID"])
summary(as.factor(Sample_file[Sample_file$PatientID %in% ids & Sample_file$Sample_Type!="Normal kidney parenchyma","Sample_Type"]))

Sample_file$Sample_Type2=revalue(as.factor(Sample_file$Sample_Type), c("Clear cell RCC"="ccRCC", "Hybrid oncocytic renal neoplasm"="Hybrid_RO_RN", "Hybrid oncocytic/Chromophobe type"="Hybrid_RO_chRCC_type", "Normal kidney parenchyma"="NKP", "Oncocytoma"="RO", "RCC-Chromophobe"="chRCC"))
Sample_file$PatientID=as.factor(Sample_file$PatientID)
Sample_file$PatientID2=Sample_file$PatientID

levels(Sample_file$PatientID2)=1:length(levels(as.factor(Sample_file$PatientID)))
Sample_file$Sample_ID2=paste(Sample_file$Sample_Type2, Sample_file$PatientID2, sep="_")

Sample_file=read.table("datadir/Sample_file.txt", sep="\t", header=T)

Sample_file$Sample_col=Sample_file$Sample_Type
Sample_file$Sample_col=revalue(as.factor(Sample_file$Sample_col), c("Clear cell RCC"="#E69F00","RCC-Chromophobe"="#56B4E9", "Normal kidney parenchyma"="#009E73", "Oncocytoma"="#FD61D1", "Hybrid oncocytic/Chromophobe type"="#A3A500", "Hybrid oncocytic renal neoplasm"="#D55E00"))
Sample_file$Sample_Type2=revalue(as.factor(Sample_file$Sample_Type), c("Clear cell RCC"="ccRCC", "RCC-Chromophobe"="chRCC", "Normal kidney parenchyma"="NKP", "Oncocytoma"="Oncocytoma", "Hybrid oncocytic/Chromophobe type"="Hybrid oncocytic/ChRCC", "Hybrid oncocytic renal neoplasm"="Hybrid oncocytic/renal"))
Sample_file$Sample_col=factor(Sample_file$Sample_col)
levels(Sample_file$Sample_Type2)=levels(Sample_file$Sample_Type2)[c(1,6,4,5,3,2)]
levels(Sample_file$Sample_col)=levels(Sample_file$Sample_col)[c(1,6,4,5,3,2)]

write.table(Sample_file, file="datadir/Sample_file.txt", sep="\t")

##########################################################
#Download and process TCGA 450k methylation data 
#################################################################

datadir="/srv/gevaertlab/data/TCGA/Kevin_Data/MethylationData450k/"

cancerType=c("KICH","KIRC","KIRP")
dat=load(paste(RscriptsPath,"BatchData.rda",sep=""))

for(i in 1:length(cancerType)){
  CancerSite=cancerType[i] 
  TargetDirectory=paste(datadir, CancerSite,"/", sep="")
  ifelse(!file.exists(TargetDirectory), dir.create(TargetDirectory),NA)
}

for(i in 1:length(cancerType)){
  CancerSite=toupper(cancerType[i])
  tag27kor450k="450"
  TargetDirectory=paste(datadir, CancerSite,"/", sep="")
  dataType="stddata"
  
  METdirectories=Download_CancerSite_DNAmethylation_KB(CancerSite,TargetDirectory,downloadData=TRUE, tag27kor450k)
  #Worked
}

#Process methylation data
library(data.table)

for(i in 1:length(cancerType)){
  CancerSite=cancerType[i]
  path=paste(datadir, CancerSite,"/", sep="")
  METdirectory=paste(path, "/", list.files(path), "/", sep="")
  ProcessedData=Preprocess_CancerSite_Methylation450k_KB(CancerSite,METdirectory)
  save(ProcessedData, file=paste(datadir, "MET_", CancerSite, "Processed.RData", sep=""))
}

####################################################################################################
#Run MethylMix on all three TCGA kidney cancers, versus TCGA kidney normal
####################################################################################################

#Restrict to CpG sites that are mean <0.2 in normal kidney parenchyma in TCGA
#load meth data

#Order KIRC, KIRP, KICH
dat=load("datadir/MET_KIRP_Processed.Rdata")
METcancer_KIRP=ProcessedData$MET_Data_Cancer
METnormal_KIRP=ProcessedData$MET_Data_Normal
rm(ProcessedData)
METcancer_KIRP=TCGA_GENERIC_CleanUpSampleNames(METcancer_KIRP, 12) 
#n=275
METnormal_KIRP=TCGA_GENERIC_CleanUpSampleNames(METnormal_KIRP, 12) 
#n=44

#KICH
dat=load("datadir/MET_KICH_Processed.Rdata")
METcancer_KICH=ProcessedData$MET_Data_Cancer
#No normal data for KICH
rm(ProcessedData)
METcancer_KICH=TCGA_GENERIC_CleanUpSampleNames(METcancer_KICH, 12) 
#n=65 

dat=load("datadir/ProcessedData_KIRC.RData")
METcancer_KIRC=ProcessedData$MET_Data_Cancer
METnormal_KIRC=ProcessedData$MET_Data_Normal
rm(ProcessedData)

METcancer_KIRC=TCGA_GENERIC_CleanUpSampleNames(METcancer_KIRC, 12) 
#n=315
METnormal_KIRC=TCGA_GENERIC_CleanUpSampleNames(METnormal_KIRC, 12) 
#n=160

#
OverlapProbes=intersect(rownames(METnormal_KIRP),rownames(METnormal_KIRC))
METnormal_KIRP=METnormal_KIRP[OverlapProbes,]
METnormal_KIRC=METnormal_KIRC[OverlapProbes,]

METnormal=cbind(METnormal_KIRP, METnormal_KIRC)

commoncpgs=intersect(
     intersect(rownames(METcancer_KIRP), intersect(rownames(METcancer_KIRC), rownames(METcancer_KICH))),
     rownames(METnormal))
commoncpgs=commoncpgs[-c(grep("rs", commoncpgs))]
#n=395143

#restrict meth to just these common variable CpGs
METcancer_KIRP=METcancer_KIRP[commoncpgs,]
METcancer_KIRC=METcancer_KIRC[commoncpgs,]
METcancer_KICH=METcancer_KICH[commoncpgs,]
METnormal_KIRP=METnormal_KIRP[commoncpgs,]
#METnormal_KIRP=scale(METnormal_KIRP)
METnormal_KIRC=METnormal_KIRC[commoncpgs,]
#METnormal_KIRC=scale(METnormal_KIRC)

#Need to specify which colnames represent cancer and normal
colnames(METcancer_KIRP)=paste(colnames(METcancer_KIRP),"-01", sep="")
colnames(METcancer_KIRC)=paste(colnames(METcancer_KIRC),"-01", sep="")
colnames(METcancer_KICH)=paste(colnames(METcancer_KICH),"-01", sep="")
colnames(METnormal_KIRP)=paste(colnames(METnormal_KIRP),"-11", sep="")
colnames(METnormal_KIRC)=paste(colnames(METnormal_KIRC),"-11", sep="")

Barcodes_METcancer_KIRP=colnames(METcancer_KIRP)
Barcodes_METcancer_KIRC=colnames(METcancer_KIRC)
Barcodes_METcancer_KICH=colnames(METcancer_KICH)
Barcodes_METnormal_KIRP=colnames(METnormal_KIRP)
Barcodes_METnormal_KIRC=colnames(METnormal_KIRC)
Barcodes_METnormal=c(Barcodes_METnormal_KIRC, Barcodes_METnormal_KIRP)

METall=cbind(
     METcancer_KIRP,
     METcancer_KIRC,
     METcancer_KICH,
     METnormal_KIRP,
     METnormal_KIRC)
#859 samples

#get CpGs with low mean methylation in both normal and oncocytoma tissue 
cgsNormLow=names(rowMeans(METnormal)[which(rowMeans(METnormal)<0.2)])
betaonc=beta[,colnames(beta) %in% BarcodesOncocytoma]
cgsOncLow=names(rowMeans(betaonc)[which(rowMeans(betaonc)<0.2)])
cgsOncNormLow=intersect(cgsNormLow, cgsOncLow)

METcancerAll=cbind(
     METcancer_KIRP,
     METcancer_KIRC,
     METcancer_KICH)

METcancerAll_cgsOncNormLow=METcancerAll[cgsOncNormLow,]
METcancerAll_cgsOncNormLow0.3=METcancerAll_cgsOncNormLow[names(rowMeans(METcancerAll_cgsOncNormLow)[which(rowMeans(METcancerAll_cgsOncNormLow)>0.3)]),]
#Running MethylMix on 95 CpGs that are <0.2 methylation in normal kidney parenchyma and oncocytoma

kidney_unmeth_norm_probes=rownames(METcancerAll_cgsOncNormLow0.3)
save(kidney_unmeth_norm_probes, file="resultsdir/kidney_unmeth_norm_probes.RData")
#Take this set of CpGs over to crosswood to 
#scp resultsdir/kidney_unmeth_norm_probes.RData kbren@crosswood.stanford.edu:/srv/gevaertlab/data/TCGA/Kevin_Data/

library(parallel)
library(doParallel)
library(RPMM)

#Ran MethylMix on CpGs that were unmethylated in onc and normal
MethylMixResults=MethylMix("univ.beta",METcancerAll_cgsOncNormLow0.3, METnormal, Parallel=T ,filter="none2")
names(MethylMixResults$MixtureStates)=MethylMixResults$MethylationDrivers
save(MethylMixResults, file="resultsdir/MethylMix_TCGA_PanCanKidney.RData")
#They are all retained as hypermethylated by MethylMix

dat=load("resultsdir/MethylMix_TCGA_PanCanKidney.RData")

##########################################
#Get TCGA kidney cancer processed methylation data
#This is just making a big dataframe with all of the TCGA kidney data, and a Sample_File to determine each class
###########################################

dat=load("datadir/MET_KIRP_Processed.Rdata")
METcancer_KIRP=ProcessedData$MET_Data_Cancer
METnormal_KIRP=ProcessedData$MET_Data_Normal
rm(ProcessedData)
METcancer_KIRP=TCGA_GENERIC_CleanUpSampleNames(METcancer_KIRP, 12) 
#n=275
METnormal_KIRP=TCGA_GENERIC_CleanUpSampleNames(METnormal_KIRP, 12) 
#n=44

#KICH
dat=load("datadir/MET_KICH_Processed.Rdata")
METcancer_KICH=ProcessedData$MET_Data_Cancer
#No normal data for KICH
rm(ProcessedData)
METcancer_KICH=TCGA_GENERIC_CleanUpSampleNames(METcancer_KICH, 12) 
#n=65 

#scp -r kbren@crosswood.stanford.edu:/srv/gevaertlab/data/TCGA/Kevin_Data/MethylationData450k/KIRC/gdac_20160715/ProcessedData_KIRC.RData datadir/
dat=load("datadir/ProcessedData_KIRC.RData")
METcancer_KIRC=ProcessedData$MET_Data_Cancer
METnormal_KIRC=ProcessedData$MET_Data_Normal
rm(ProcessedData)

METcancer_KIRC=TCGA_GENERIC_CleanUpSampleNames(METcancer_KIRC, 12) 
#n=315
METnormal_KIRC=TCGA_GENERIC_CleanUpSampleNames(METnormal_KIRC, 12) 
#n=160

#
OverlapProbes=intersect(rownames(METnormal_KIRP),rownames(METnormal_KIRC))
METnormal_KIRP=METnormal_KIRP[OverlapProbes,]
METnormal_KIRC=METnormal_KIRC[OverlapProbes,]

METnormal=cbind(METnormal_KIRP, METnormal_KIRC)

commoncpgs=intersect(
     intersect(rownames(METcancer_KIRP), intersect(rownames(METcancer_KIRC), rownames(METcancer_KICH))),
     rownames(METnormal))
commoncpgs=commoncpgs[-c(grep("rs", commoncpgs))]
#n=395143

#restrict meth to just these common CpGs
METcancer_KIRP=METcancer_KIRP[commoncpgs,]
METcancer_KIRC=METcancer_KIRC[commoncpgs,]
METcancer_KICH=METcancer_KICH[commoncpgs,]
METnormal_KIRP=METnormal_KIRP[commoncpgs,]
#METnormal_KIRP=scale(METnormal_KIRP)
METnormal_KIRC=METnormal_KIRC[commoncpgs,]
#METnormal_KIRC=scale(METnormal_KIRC)

#Need to specify which colnames represent cancer and normal
colnames(METcancer_KIRP)=paste(colnames(METcancer_KIRP),"-01", sep="")
colnames(METcancer_KIRC)=paste(colnames(METcancer_KIRC),"-01", sep="")
colnames(METcancer_KICH)=paste(colnames(METcancer_KICH),"-01", sep="")
colnames(METnormal_KIRP)=paste(colnames(METnormal_KIRP),"-11", sep="")
colnames(METnormal_KIRC)=paste(colnames(METnormal_KIRC),"-11", sep="")

Barcodes_METcancer_KIRP=colnames(METcancer_KIRP)
Barcodes_METcancer_KIRC=colnames(METcancer_KIRC)
Barcodes_METcancer_KICH=colnames(METcancer_KICH)
Barcodes_METnormal_KIRP=colnames(METnormal_KIRP)
Barcodes_METnormal_KIRC=colnames(METnormal_KIRC)
Barcodes_METnormal=c(Barcodes_METnormal_KIRC, Barcodes_METnormal_KIRP)

METall=cbind(
     METcancer_KIRP,
     METcancer_KIRC,
     METcancer_KICH,
     METnormal_KIRP,
     METnormal_KIRC)

##################
#Make sample class file to index sample types (KIRC, KIRP, KICH or normal)
####################

SampleClass=rep(NA,ncol(METall))
names(SampleClass)=colnames(METall)

SampleClass[names(SampleClass) %in% Barcodes_METcancer_KIRP]=1
SampleClass[names(SampleClass) %in% Barcodes_METcancer_KIRC]=2
SampleClass[names(SampleClass) %in% Barcodes_METcancer_KICH]=3
SampleClass[names(SampleClass) %in% c(Barcodes_METnormal_KIRP, Barcodes_METnormal_KIRC)]=4

SampleClass=t(as.matrix(SampleClass))
rownames(SampleClass)="Sample_Type"

SampleClass_TCGA=SampleClass

###########################################
#Next run PAM analysis on TCGA samples
##########################################

library(pamr)

dat=load("resultsdir/MethylMix_TCGA_PanCanKidney_AllCpGs.RData")
MethylationStates=MethylMixResults$MethylationStates

MixtureStates=MethylMixResults$MixtureStates

HyperProbes=names(MixtureStates[unlist(lapply(1:length(MixtureStates), function(x) max(MixtureStates[[x]])>0.2))])
HypoProbes=names(MixtureStates[unlist(lapply(1:length(MixtureStates), function(x) min(MixtureStates[[x]])<(-0.2)))])
HyperAndHypoProbes=intersect(HyperProbes, HypoProbes)
#1772
JustHyperProbes=setdiff(HyperProbes, HyperAndHypoProbes)
#38475
JustHypoProbes=setdiff(HypoProbes, HyperAndHypoProbes)
#41785
MethylationStatesHyper=MethylationStates[HyperProbes,]
#40247 probes
MethylationStatesHypo=MethylationStates[HypoProbes,]
#43557 probes 

AbnormalProbes=c(JustHyperProbes, JustHypoProbes, HyperAndHypoProbes)
MethylationStatesAbnormal=MethylationStates[AbnormalProbes,]

###########################################
#running the PAM model
###########################################

METall_select=METall[AbnormalProbes,]

NrFolds=10
#Run with class 4 specified as normal
#PAMresults5=PAM_Analysis_CV_KB(METall,SampleClass,NrFolds,NrLambdaFolds=2, HeteroNorm=4)
#Now compare same without specification of normal 
PAMresults=PAM_Analysis_CV_KB(METall_select,SampleClass_TCGA,NrFolds,NrLambdaFolds=2)
save(PAMresults, file="resultsdir/PAM_PanKidney_CrossValidation_MethylMix0.2_probes.RData")
file=load("resultsdir/PAM_PanKidney_CrossValidation_MethylMix0.2_probes.RData")

tab=PAMresults$NrMisclassificationsPerClass
rownames(tab)=c("KIRP","KIRC","KICH","Normal")
#write.table(tab, file="resultsdir/PAM_predictions_kidney_classes_cross_validation_accuracy.txt", sep="\t")

misclassifications=table(SampleClass_TCGA, PAMresults$Predictions)

rownames(misclassifications)=c("KIRP","KIRC","KICH","Normal")
colnames(misclassifications)=c("KIRP","KIRC","KICH","Normal")
prop.table(misclassifications,1)

write.table(misclassifications, file="resultsdir/PAM_misclassifications_Kidney_MethylMix.txt", sep="\t")
#PAMresults5$OverallPerformance
#63 misclassifications overall
#PAMresults5$NrMisclassificationsPerClass
#tab3=read.table("resultsdir/PAM_predictions_kidney_classes_cross_validation_accuracy.txt", sep="\t", header=T)
#classification improved by just using MethylMix CpG sites 
misclassifications=read.table("resultsdir/PAM_misclassifications_Kidney_MethylMix.txt", sep="\t", header=T)

#posterior here is the maximum posterior probability
misclass=as.data.frame(cbind(SampleClass[1,],PAMresults$Predictions, apply(PAMresults$PosteriorProbabilities,1,max)))
colnames(misclass)=c("actual","predicted","Posterior")
misclassLowProb=rownames(misclass[misclass$Posterior<0.6,])

classified=misclass[misclass$Posterior>0.6,]
classified$actual=revalue(as.factor(classified$actual), c("1"="pRCC","2"="ccRCC","3"="chRCC","4"="NKP"))
classified$predicted=revalue(as.factor(classified$predicted), c("1"="pRCC","2"="ccRCC","3"="chRCC","4"="NKP"))

write.table(classified, file="resultsdir/PAM_Allclassified_PanRCC_HighPosteriorProbabilities.txt", sep="\t")
#classified 855/859 samples with posterior probabilities >0.6
#Of these 

#four samples were excluded due to low posterior probability <0.6 
#three of these samples had low classification posterior probabilities, including one KIRP misclassified as KIRC (TCGA-WN-A9G9-01), and 2 KIRCs misclassified as KIRP (TCGA-B0-5098-01, TCGA-B4-5832-01), 
misclassified=misclass[misclass$actual!=misclass$predicted & misclass$Posterior>0.6,]
misclassified$actual=revalue(as.factor(misclassified$actual), c("1"="pRCC","2"="ccRCC","3"="chRCC","4"="NKP"))
misclassified$predicted=revalue(as.factor(misclassified$predicted), c("1"="pRCC","2"="ccRCC","3"="chRCC","4"="NKP"))

write.table(misclassified, file="resultsdir/PAM_misclassified_PanRCC_HighPosteriorProbabilities.txt", sep="\t")
misclassified=read.table("resultsdir/PAM_misclassified_PanRCC_HighPosteriorProbabilities.txt", sep="\t", header=T)

TestPredictions=as.data.frame(PAMresults$Predictions)

misclassifications=table(SampleClass_TCGA[1,setdiff(colnames(SampleClass_TCGA),misclassLowProb)], TestPredictions[setdiff(rownames(TestPredictions),misclassLowProb),])
rownames(misclassifications)=c("pRCC/KIRP","ccRCC/KIRC","chRCC/KICH","NKP")
colnames(misclassifications)=c("pRCC/KIRP","ccRCC/KIRC","chRCC/KICH","NKP")
prop.table(misclassifications,1)

write.table(misclassifications, file="resultsdir/PAM_misclassifications_Kidney_MethylMix_removed_LowPosteriorProbabilities.txt", sep="\t")
#
misclassifications=read.table("resultsdir/PAM_misclassifications_Kidney_MethylMix_removed_LowPosteriorProbabilities.txt", sep="\t")

counts=rowSums(misclassifications[,1:ncol(misclassifications)])
props=apply(sweep(misclassifications, MARGIN = 2, counts,"/"), 2, max)
misclassifications=cbind(misclassifications, props)
colnames(misclassifications)[5]="Accuracy"

write.table(misclassifications, file="resultsdir/PAM_misclassifications_Kidney_MethylMix_removed_LowPosteriorProbabilities.txt", sep="\t")
misclassifications=read.table("resultsdir/PAM_misclassifications_Kidney_MethylMix_removed_LowPosteriorProbabilities.txt", sep="\t", header=T)
colnames(misclassifications)=c("Sample type",colnames(misclassifications)[1:5])
misclassifications$Accuracy=round(misclassifications$Accuracy,2)

############################################
#PAM model I used to classify the Stanford study samples
###########################################

OverlapProbes=intersect(rownames(beta), rownames(METall_select))
beta_METall=beta[OverlapProbes,]
METall_beta=METall_select[OverlapProbes,]

NrFolds=10
PAMresults=PAM_Analysis_KB(METall_beta, SampleClass, beta_METall, NrLambdaFolds=10)
save(PAMresults, file="resultsdir/PAM_results_all_classasignments_NewData_MethylMix_probes.RData")
dat=load("resultsdir/PAM_results_all_classasignments_NewData_MethylMix_probes.RData")

TestPredictions=PAMresults$TestPredictions
#rownames(TestPredictions)=TestPredictions[,1]
#TestPredictions=TestPredictions[,2:ncol(TestPredictions)]
TestPredictions=as.data.frame(TestPredictions, drop=FALSE)
TestPredictions$Sample_Type=rep(NA, nrow(TestPredictions))

for(i in 1:nrow(TestPredictions)){
  TestPredictions$Sample_Type[i]=as.character(Sample_file[rownames(Sample_file) %in% rownames(TestPredictions)[i],"Sample_Type"])
}

tab=table(TestPredictions$PAMprediction, TestPredictions$Sample_Type)
rownames(tab)=c("ccRCC/KIRC","chRCC/KICH","NKP")

#
write.table(tab, file="resultsdir/PAM_predictions_new_samples_NoMaxNrFeatures_MethylMix_probes.txt", sep="\t")
#Next do PAM to classify differences between oncocytomas and chromophobes and apply to KICH as validation set

#
misclassifications=read.table("resultsdir/PAM_predictions_new_samples_NoMaxNrFeatures_MethylMix_probes.txt", sep="\t", header=T)

###########################################################
#Visualize the posterior probabilities for assignments, both for TCGA classifications and new samples
###########################################################

library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

#file=load("resultsdir/PAM_CrossValidation.RData")
file=load("resultsdir/PAM_PanKidney_CrossValidation_MethylMix0.2_probes.RData")

#
PosteriorProbabilities=PAMresults$PosteriorProbabilities
colnames(PosteriorProbabilities)=c("pRCC/KIRP","ccRCC/KIRC","ChRCC/KICH","NKP")

#Colorblind palette
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#FD61D1", "#A3A500", "#D55E00", "#CC79A7")

#make colorblind version 
rlab=as.matrix(rownames(PosteriorProbabilities))
rownames(rlab)=rlab[,1]
colnames(rlab)=""
rlab[rownames(rlab) %in% Barcodes_METcancer_KIRP,1]=cbp2[1]
rlab[rownames(rlab) %in% Barcodes_METcancer_KIRC,1]=cbp2[2]
rlab[rownames(rlab) %in% Barcodes_METcancer_KICH,1]=cbp2[3]
rlab[rownames(rlab) %in% Barcodes_METnormal,1]=cbp2[4]

#Add misclassLowProb to this
rlab=as.data.frame(rlab)
rlab$misclassLowProb=rep("white", nrow(rlab))
rlab[intersect(misclassLowProb, Barcodes_METcancer_KIRP),]=cbp2[1]
rlab[intersect(misclassLowProb, Barcodes_METcancer_KIRC),]=cbp2[2]
rlab[intersect(misclassLowProb, Barcodes_METcancer_KICH),]=cbp2[3]
colnames(rlab)=c("","")
rlab=as.matrix(rlab)
rlab=rlab[,c(2,1)]

file="resultsdir/heatmap_posterior_possibility_PAM_PanKidney_colorblind_MethylMix"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
par(oma=c(5,0,0,0))
heatmap.3(PosteriorProbabilities, Rowv=FALSE, Colv=FALSE, dendrogram = "none", labRow = "", col=brewer.pal(9,"Reds"), RowSideColors=t(rlab), KeyValueName="Classification probability")
dev.off()

rlab=as.data.frame(rlab)
rlab[rlab[,1]!="white",]

file="resultsdir/heatmap_posterior_possibility_PAM_PanKidney_colorblind_legend"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8, res=600)
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("center", c(colnames(PosteriorProbabilities)), col=c(cbp2[1:4]), pch=15, title="TCGA sample type", cex=2, bty="n")
dev.off()

file="resultsdir/heatmap_posterior_possibility_PAM_PanKidney_colorblind_legend"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8, res=600)
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("center", c(colnames(PosteriorProbabilities)[1:2],c("chRCC/KICH","NKP","oncocytoma","hybrid_oncocytoma_chRCC","hybrid_oncocytoma_RN")), col=c(cbp2[1:7]), pch=15, title="Sample type", cex=2, bty="n")
dev.off()

#Make table showing classification of Stanford study samples
BarcodesHighConfidence=names(which(apply(PosteriorProbabilities, 1, max)>=0.6))
Predictions=as.data.frame(PAMresults6$Predictions)
colnames(Predictions)[1]="Predicted"
Predictions$Sample_type=rep(NA, nrow(Predictions))
Predictions[Barcodes_METcancer_KIRP,"Sample_type"]="KIRP"
Predictions[Barcodes_METcancer_KIRC,"Sample_type"]="KIRC"
Predictions[Barcodes_METcancer_KICH,"Sample_type"]="KICH"
Predictions[Barcodes_METnormal,"Sample_type"]="NKP"

Predictions$Predicted=revalue(as.factor(Predictions$Predicted), c("1"="KIRP","2"="KIRC","3"="KICH","4"="NKP"))
Predictions=Predictions[rownames(Predictions)!= "TCGA-CZ-5986-01",]

write.table(Predictions, file="resultsdir/TCGA_PanKidney_predictions_misclassifications.txt", sep="\t")

#Look at misclassified TCGA KICHs, how are they misclassified 
summary(as.factor(Predictions[Predictions$Sample_type=="KICH","Predicted"]))

######################################################
#Predictions heatmap for newly assigned samples
####################################################

#using more refined color vecotr
n=12
install.packages("viridis")
library(viridis)
viridis_pal(option = "D")(n)  # n = number of colors seeked

cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#FD61D1", "#A3A500", "#D55E00", "#CC79A7")

dat=load("resultsdir/PAM_results_all_classasignments_NewData_MethylMix_probes.RData")

PosteriorProbabilities=PAMresults$PosteriorProbabilites
colnames(PosteriorProbabilities)=c("pRCC/KIRP","ccRCC/KIRC","ChRCC/KICH","NKP")

PosteriorProbabilities=as.data.frame(PosteriorProbabilities)

PosteriorProbabilities$Sample_ID=Sample_file[rownames(PosteriorProbabilities),"Sample_ID2"]
PosteriorProbabilities$Sample_Type2=Sample_file[rownames(PosteriorProbabilities),"Sample_Type2"]
PosteriorProbabilities$col=rep(NA, nrow(PosteriorProbabilities))
PosteriorProbabilities$col[PosteriorProbabilities$Sample_Type2=="RO"]=as.character(cbp2[5])
PosteriorProbabilities$col[PosteriorProbabilities$Sample_Type2=="chRCC"]=as.character(cbp2[3])
PosteriorProbabilities$col[PosteriorProbabilities$Sample_Type2=="NKP"]=as.character(cbp2[4])
PosteriorProbabilities$col[PosteriorProbabilities$Sample_Type2=="ccRCC"]=as.character(cbp2[2])
PosteriorProbabilities$col[PosteriorProbabilities$Sample_Type2=="Hybrid_RO_chRCC_type"]=as.character(cbp2[7])
PosteriorProbabilities$col[PosteriorProbabilities$Sample_Type2=="Hybrid_RO_RN"]=as.character(cbp2[6])

PosteriorProbabilities=PosteriorProbabilities[order(PosteriorProbabilities$Sample_Type2),]

y=c("ccRCC","chRCC","NKP","RO","Hybrid_RO_chRCC_type","Hybrid_RO_RN")
PosteriorProbabilities=PosteriorProbabilities[order(match(PosteriorProbabilities$Sample_Type2, y)),]

#make colorblind version 
rlab=as.matrix(PosteriorProbabilities$col)
rownames(rlab)=rownames(PosteriorProbabilities)
colnames(rlab)=""

file="resultsdir/heatmap_posterior_possibility_PAM_PanKidney_new_sample_assignments_colorblind_MethylMix_probes"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
par(oma=c(5,0,0,5))
heatmap.3(PosteriorProbabilities[,1:4], Rowv=FALSE, Colv=FALSE, symbreaks=F, scale="none", dendrogram = "none", col=brewer.pal(9,"Reds"), KeyValueName="Classification probability", labRow=PosteriorProbabilities$Sample_ID, RowSideColors=t(rlab))
dev.off()

################################################################################
#Make multidimensional scaling plots showing TCGA and Stanford samples 
#################################################################################

OverlapProbes=intersect(rownames(METall_prescale), rownames(ClassCentroids_beta))
METall_CD=METall_prescale[OverlapProbes,]
ClassCentroids_METall=ClassCentroids_beta[OverlapProbes,]

ClassCentroids_Hand=cbind(rowMeans(METall_CD[,Barcodes_METcancer_KIRP]),
                          rowMeans(METall_CD[,Barcodes_METcancer_KIRC]),
                          rowMeans(METall_CD[,Barcodes_METcancer_KICH]),
                          rowMeans(METall_CD[,Barcodes_METnormal]))

BarcodesClearCell=rownames(Sample_file[grep("Clear cell RCC",Sample_file$Sample_Type),])
#There are no KIRPs in newdata
BarcodesChromophobe=rownames(Sample_file[grep("RCC-Chromophobe",Sample_file$Sample_Type),])
BarcodesNormal=rownames(Sample_file[grep("Normal kidney parenchyma",Sample_file$Sample_Type),])
BarcodesOncocytoma=rownames(Sample_file[grep("Oncocytoma",Sample_file$Sample_Type),])
BarcodesHybridoncocyticrenal=rownames(Sample_file[grep("Hybrid oncocytic renal neoplasm",Sample_file$Sample_Type),])
BarcodesHybridoncocyticchromophobe=rownames(Sample_file[grep("Hybrid oncocytic/Chromophobe type",Sample_file$Sample_Type),])

OverlapProbes=intersect(rownames(beta_prescale), intersect(rownames(METall_prescale), rownames(ClassCentroids_Hand)))
m1=ClassCentroids_Hand[OverlapProbes,]
m2=METall_prescale[OverlapProbes,]
m3=beta_prescale[OverlapProbes,]

d1=t(cbind(m1, m2, m3))
d2=dist(d1)

fit <- cmdscale(d2,eig=TRUE, k=2) # k is the number of dim

x <- fit$points[,1]
y <- fit$points[,2]

#Just TCGA
file="resultsdir/MSD_plot_allTCGA"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8.5, res=600)
plot(x[Barcodes_METcancer_KIRP], y[Barcodes_METcancer_KIRP], xlab="Coordinate 1", ylab="Coordinate 2", col=1, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="Multidomensional scaling plot TCGA")
par(new=TRUE)
plot(x[Barcodes_METcancer_KIRC], y[Barcodes_METcancer_KIRC], xlab="Coordinate 1", ylab="Coordinate 2", col=2, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[Barcodes_METcancer_KICH], y[Barcodes_METcancer_KICH], xlab="Coordinate 1", ylab="Coordinate 2", col=3, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[Barcodes_METnormal], y[Barcodes_METnormal], xlab="Coordinate 1", ylab="Coordinate 2", col=4, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[1:4],y[1:4], xlab="Coordinate 1", ylab="Coordinate 2", col=c(1:4), pch=15, cex=1.5, , xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), main="")
#legend(-0.5,-1.5, c("Clear cell RCC", "Papillary RCC", "RCC-Chromophobe", "Normal","Oncocytoma","Hybrid oncocytic renal neoplasm","Hybrid oncocytic/Chromophobe type"), col=c(1:6), title="Sample type", pch=19, cex=1)
legend("topright", c("Clear cell RCC", "Papillary RCC", "RCC-Chromophobe", "Normal","Oncocytoma","Hybrid oncocytic renal neoplasm","Hybrid oncocytic/Chromophobe type"), col=c(1:6), title="Sample type", pch=19, cex=1)
dev.off()

#Just new samples
file="resultsdir/MSD_plot_allNewsamples"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8.5, res=600)
plot(x[BarcodesClearCell], y[BarcodesClearCell], xlab="Coordinate 1", ylab="Coordinate 2", col=1, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="Multidomensional scaling plot new samples")
par(new=TRUE)
plot(x[BarcodesChromophobe], y[BarcodesChromophobe], xlab="Coordinate 1", ylab="Coordinate 2", col=3, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[BarcodesNormal], y[BarcodesNormal], xlab="Coordinate 1", ylab="Coordinate 2", col=4, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[BarcodesOncocytoma], y[BarcodesOncocytoma], xlab="Coordinate 1", ylab="Coordinate 2", col=5, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[BarcodesHybridoncocyticrenal], y[BarcodesHybridoncocyticrenal], xlab="Coordinate 1", ylab="Coordinate 2", col=6, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[BarcodesHybridoncocyticchromophobe], y[BarcodesHybridoncocyticchromophobe], xlab="Coordinate 1", ylab="Coordinate 2", col=7, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[1:4],y[1:4], xlab="Coordinate 1", ylab="Coordinate 2", col=c(1:4), pch=15, cex=1.5, , xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), main="")
#legend(-0.5,-1.5, c("Clear cell RCC", "Papillary RCC", "RCC-Chromophobe", "Normal","Oncocytoma","Hybrid oncocytic renal neoplasm","Hybrid oncocytic/Chromophobe type"), col=c(1:6), title="Sample type", pch=19, cex=1)
legend("topright", c("Clear cell RCC", "Papillary RCC", "RCC-Chromophobe", "Normal","Oncocytoma","Hybrid oncocytic renal neoplasm","Hybrid oncocytic/Chromophobe type"), col=c(1:6), title="Sample type", pch=19, cex=1)
dev.off()

#This version has reduced color saturation for TCGA points
file="resultsdir/MSD_plot_allNewsamples_no_normalization"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8.5, res=600)
plot(x[Barcodes_METcancer_KIRP], y[Barcodes_METcancer_KIRP], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(1,alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="Multidomensional scaling plot TCGA")
par(new=TRUE)
plot(x[Barcodes_METcancer_KIRC], y[Barcodes_METcancer_KIRC], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(2,alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[Barcodes_METcancer_KICH], y[Barcodes_METcancer_KICH], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(3,alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[Barcodes_METnormal], y[Barcodes_METnormal], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(4,alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[1:4],y[1:4], xlab="Coordinate 1", ylab="Coordinate 2", col=c(1:4), pch=15, cex=1.5, , xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), main="")
par(new=TRUE)
plot(x[BarcodesClearCell], y[BarcodesClearCell], xlab="Coordinate 1", ylab="Coordinate 2", col=2, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19)
par(new=TRUE)
plot(x[BarcodesChromophobe], y[BarcodesChromophobe], xlab="Coordinate 1", ylab="Coordinate 2", col=3, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[BarcodesNormal], y[BarcodesNormal], xlab="Coordinate 1", ylab="Coordinate 2", col=4, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[BarcodesOncocytoma], y[BarcodesOncocytoma], xlab="Coordinate 1", ylab="Coordinate 2", col=5, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[BarcodesHybridoncocyticrenal], y[BarcodesHybridoncocyticrenal], xlab="Coordinate 1", ylab="Coordinate 2", col=6, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[BarcodesHybridoncocyticchromophobe], y[BarcodesHybridoncocyticchromophobe], xlab="Coordinate 1", ylab="Coordinate 2", col=7, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[1:4],y[1:4], xlab="Coordinate 1", ylab="Coordinate 2", col=c(1:4), pch=15, cex=1.5, , xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), main="")
dev.off()

legend("topright", c("Clear cell RCC", "Papillary RCC", "RCC-Chromophobe", "Normal","Oncocytoma","Hybrid oncocytic renal neoplasm","Hybrid oncocytic/Chromophobe type"), col=c(1:7), title="Sample type", pch=19, cex=1)
dev.off()

file="resultsdir/MDS_sample_type_legend"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8, res=600)
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("center", c("KIRP","KIRC","KICH","Normal","Oncocytoma", "Hybrid oncocytic/renal","Hybrid oncocytic/chromophobe"), col=c(1:7), pch=18, title="Sample type", cex=2, bty="n")
dev.off()

file="resultsdir/MDS_point_type_legend"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8, res=600)
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("center", c("TCGA study","TCGA study centroid","Stanford study"), col=c(adjustcolor(1,alpha.f=0.5),1,1), pch=c(1,15,19), title="Sample type", cex=2, bty="n")
dev.off()

##################################################################
#Do an MDS plot with just KICH and oncocytoma, but splitting TCGA KICH into clssic and eosinophillic
#Supplementary figure 6
################################################################

dat=load("datadir/MET_KICH_Processed.Rdata")
METcancer_KICH=ProcessedData$MET_Data_Cancer
#No normal data for KICH
rm(ProcessedData)
METcancer_KICH=TCGA_GENERIC_CleanUpSampleNames(METcancer_KICH, 12) 

METcancer_KICH_Var=TCGA_GENERIC_GeneFiltering('MAD',METcancer_KICH,25)
#METcancer_KICH_Var=scale(METcancer_KICH_Var)

ann=read.table("datadir/KICH_TCGA_annotation.txt", sep="\t", header=T)
Barcodes_METcancer_KICH_Eosinophillic=as.character(ann[ann$Eosinophilic.vs.Classic=="Eosinophilic","TCGA.patient.code"])
Barcodes_METcancer_KICH_Eosinophillic=Barcodes_METcancer_KICH_Eosinophillic[Barcodes_METcancer_KICH_Eosinophillic %in% colnames(METcancer_KICH)]

Barcodes_METcancer_KICH_Classic=as.character(ann[ann$Eosinophilic.vs.Classic=="Classic","TCGA.patient.code"])
Barcodes_METcancer_KICH_OncFeatures=as.character(ann[grep("oncocytic features noted", ann$Path.review.notes),"TCGA.patient.code"])

Barcodes_METcancer_KICH_Eosinophillic_oncocytic=intersect(Barcodes_METcancer_KICH_OncFeatures, Barcodes_METcancer_KICH_Eosinophillic)
Barcodes_METcancer_KICH_Eosinophillic_non_oncocytic=setdiff(Barcodes_METcancer_KICH_Eosinophillic, Barcodes_METcancer_KICH_OncFeatures)

#intersect(Barcodes_METcancer_KICH_OncFeatures, Barcodes_METcancer_KICH_Eosinophillic)
#All five ChRCCs with oncocytic features noted are eosinophillic


OverlapProbes=intersect(rownames(METcancer_KICH), rownames(beta))
METcancer_KICH=METcancer_KICH[OverlapProbes,]
beta_prescale=beta[OverlapProbes, c(BarcodesChromophobe, BarcodesOncocytoma)]

METcancer_KICH_Eosinophillic=METcancer_KICH[,colnames(METcancer_KICH) %in% Barcodes_METcancer_KICH_Eosinophillic]
METcancer_KICH_Classic=METcancer_KICH[,colnames(METcancer_KICH) %in% Barcodes_METcancer_KICH_Classic]

ClassCentroids_Hand=cbind(rowMeans(METcancer_KICH_Eosinophillic),
                          rowMeans(METcancer_KICH_Classic))
     
#OverlapProbes=intersect(rownames(beta_prescale), intersect(rownames(METall_prescale), rownames(ClassCentroids_Hand)))
m1=ClassCentroids_Hand
m2=METcancer_KICH[OverlapProbes,]
m3=beta_prescale[OverlapProbes,]

d1=t(cbind(m1, m2, m3))
d2=dist(d1)

fit <- cmdscale(d2,eig=TRUE, k=2) # k is the number of dim

x <- fit$points[,1]
y <- fit$points[,2]

m1=ClassCentroids_Hand
m2=METcancer_KICH
m3=beta_prescale

d1=t(cbind(m1, m2, m3))
d2=dist(d1)

fit <- cmdscale(d2,eig=TRUE, k=2) # k is the number of dim

x <- fit$points[,1]
y <- fit$points[,2]

file="resultsdir/MSD_plot_ChRCC_eosinophillicVnon_eosinophillic_colorblind"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8.5, res=600)
plot(x[Barcodes_METcancer_KICH_Eosinophillic_non_oncocytic], y[Barcodes_METcancer_KICH_Eosinophillic_non_oncocytic], xlab="Coordinate 1", ylab="Coordinate 2", col=cbp2[3], xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=0, main="Multidomensional scaling plot ChRCC and Oncocytoma")
par(new=TRUE)
plot(x[Barcodes_METcancer_KICH_Eosinophillic_oncocytic], y[Barcodes_METcancer_KICH_Eosinophillic_oncocytic], xlab="Coordinate 1", ylab="Coordinate 2", col=cbp2[3], xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=7, main="Multidomensional scaling plot ChRCC and Oncocytoma")
par(new=TRUE)
plot(x[Barcodes_METcancer_KICH_Classic], y[Barcodes_METcancer_KICH_Classic], xlab="Coordinate 1", ylab="Coordinate 2", col=cbp2[3], xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[BarcodesChromophobe], y[BarcodesChromophobe], xlab="Coordinate 1", ylab="Coordinate 2", col=cbp2[3], xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[BarcodesOncocytoma], y[BarcodesOncocytoma], xlab="Coordinate 1", ylab="Coordinate 2", col=cbp2[5], xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[1:2],y[1:2], xlab="Coordinate 1", ylab="Coordinate 2", col=cbp2[3], pch=c(3,4), cex=1.5, , xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), main="")
par(new=TRUE)
dev.off()

###
file="resultsdir/MDS_sample_type_legend_ChRCC_colorblind"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8, res=600)
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("center", c("TCGA chRCC eosinophillic","TCGA chRCC eosinophillic (oncocytic)","TCGA chRCC classic","Stanford chRCC","Stanford oncocytoma", "TCGA chRCC eosinophillic centroid","TCGA chRCC classic centroid"), col=c(cbp2[3],cbp2[3],cbp2[3],cbp2[3],cbp2[5],cbp2[3],cbp2[3]), pch=c(0,7,1,19,19,3,4), title="Sample type", cex=2, bty="n")
dev.off()

#################################################################
########################################################
#Analyzing Chopra et al data 
###########################################################
#################################################################

#Get Chopra biopsy data 
Chop=load("datadir/Archive_of_OSF_Storage/data/pdata.testing.rda")
pdata.testing=as.data.frame(pdata.testing)
pdata.testing$TissueType
#26 oncocytomas 6 KICH

#26 oncocytomas 6 KICH
#There are some patients with multiple tumor biopsies etc. 

#Oncocytoma
pdata.onc.chopra=pdata.testing[which(pdata.testing$Diagnosis=="oncocytoma" & pdata.testing$TissueType=="oncocytoma"),]
Barcodes.pdata.onc.chopra=rownames(pdata.onc.chopra)

pdata.onc.chopra.norm.adj=pdata.testing[which(pdata.testing$Diagnosis=="oncocytoma" & pdata.testing$TissueType=="Normal"),]
Barcodes.pdata.onc.chopra.norm.adj=rownames(pdata.onc.chopra.norm.adj)

#KICH
pdata.KICH.chopra=pdata.KICH.chopra=pdata.testing[which(pdata.testing$Diagnosis=="KICH" & pdata.testing$TissueType=="KICH"),]
Barcodes.pdata.KICH.chopra=rownames(pdata.KICH.chopra)

pdata.KICH.chopra.norm.adj=pdata.testing[which(pdata.testing$Diagnosis=="KICH" & pdata.testing$TissueType=="Normal"),]
Barcodes.pdata.KICH.chopra.norm.adj=rownames(pdata.KICH.chopra.norm.adj)

pdata.all.norm=pdata.testing[which(pdata.testing$TissueType=="Normal"),]
Barcodes.pdata.biopsy.norm=rownames(pdata.all.norm)

#now get beta value data 
#check if there's missing data 

Chop.beta=load("datadir/Archive_of_OSF_Storage/data/beta.testing.rda")
#Is 450k
#need to impute missing values 
#Don't process the standard way, just impute, as full processing algorithm removes too many samples (where these is more than one biopsy per patient, I think)
MET_Data=beta.testing
#2371 probes with 10% missing values removed, but all samples fine
Genes=rownames(MET_Data)
MET_Data=apply(MET_Data,2,as.numeric)
rownames(MET_Data)=Genes
     
SampleNames=colnames(MET_Data)
     
# more than 10% = removal
MissingValueThreshold=0.1
     
# removing samples with too many missing values
NrMissingsPerGene=apply(MET_Data,1,function(x) sum(is.na(x))) /ncol(MET_Data)
cat("Removing",sum(NrMissingsPerGene>MissingValueThreshold),"genes with more than 10% missing values.\n")
if (sum(NrMissingsPerGene>MissingValueThreshold)>0) MET_Data=MET_Data[NrMissingsPerGene<MissingValueThreshold,]
#removed 2371 CpG sites  
    
# removing patients with too many missings values     
NrMissingsPerSample=apply(MET_Data,2,function(x) sum(is.na(x))) /nrow(MET_Data)
cat("Removing",sum(NrMissingsPerSample>MissingValueThreshold),"patients with more than 10% missing values.\n")
if (sum(NrMissingsPerSample>MissingValueThreshold)>0) MET_Data=MET_Data[,NrMissingsPerSample<MissingValueThreshold]
     
# knn samples using Tibshirani's method     
k=15
KNNresults=impute.knn(as.matrix(MET_Data),k) # this does not work, gives storage error. ???
MET_Data_KNN=KNNresults$data
beta.testing.impute=MET_Data_KNN

saveRDS(beta.testing.impute, "datadir/Archive_of_OSF_Storage/data/beta.testing.rds")
beta.testing.impute=readRDS("datadir/Archive_of_OSF_Storage/data/beta.testing.rds")

Barcodes.Chopra.oncocytoma.KICH=c(Barcodes.pdata.onc.chopra, Barcodes.pdata.KICH.chopra)
SampleType.Chopra.onc.KICH=c(rep(2, length(Barcodes.pdata.onc.chopra)), rep(1, length(Barcodes.pdata.KICH.chopra)))
names(SampleType.Chopra.onc.KICH)=Barcodes.Chopra.oncocytoma.KICH

#############################################
#PAM analysis trained on Stanford data
################################################

Sample_file=read.table("datadir/Sample_file.txt", sep="\t", header=T)

Barcodes=c(BarcodesNormal, BarcodesOncocytoma, BarcodesChromophobe)

SampleClass=rep(NA, length(Barcodes))
names(SampleClass)=Barcodes
SampleClass[names(SampleClass) %in% BarcodesNormal]=1
SampleClass[names(SampleClass) %in% BarcodesChromophobe]=2
SampleClass[names(SampleClass) %in% BarcodesOncocytoma]=3
SampleClass=t(as.matrix(SampleClass))
rownames(SampleClass)="Sample_type"

beta_genes_Barcodes=beta[genes,colnames(SampleClass)]
beta_genes_Barcodes=scale(beta_genes_Barcodes)

NrFolds=10

PAMresults.tumor.trained.with.norm=PAM_Analysis_KB(beta_genes_Barcodes, SampleClass, beta.training.impute.onc.ChRCC.genes.s, NrLambdaFolds=10)
PAMresults.tumor.trained.with.norm$TestPredictions

#Save the PAM we used
saveRDS(PAMresults.tumor.trained.with.norm, "datadir/Archive_of_OSF_Storage/data/PAM_Chopra_tumors.rds")
PAMresults.tumor.trained.with.norm=readRDS("datadir/Archive_of_OSF_Storage/data/PAM_Chopra_tumors.rds")
save(PAMresults.tumor.trained.with.norm, )

tab=table(PAMresults.tumor.trained.with.norm$TestPredictions, SampleType.Chopra.onc.KICH.training)
tab=as.matrix(tab)
colnames(tab)=c("chRCC_TCGA","oncocytoma_Chopra")
rownames(tab)=c("Classified_as_normal","Classified_as_chRCC","Classified_as_oncocytoma")
write.table(tab, file="resultsdir/Calssifications_Chopra_tumor_model_trained_with_normal", sep="\t")

#calculate auc for classification of tumors using multi-class auc
library(pROC)
auc=multiclass.roc(as.factor(SampleType.Chopra.onc.KICH.training), PAMresults.tumor.trained.with.norm$TestPredictions[,1])
#0.8657

#PAM analsyis trained with normal tissue applied to needle biopsies 
rlab=as.matrix(rownames(PAMresults.tumor.trained.with.norm$TestPredictions))
rownames(rlab)=rlab[,1]
colnames(rlab)=""
rlab[rownames(rlab) %in% Barcodes_chRCC_training_Chopra,1]=cbp2[4]
rlab[rownames(rlab) %in% Barcodes_onc_training_Chopra,1]=cbp2[5]

#Need to also 
file="resultsdir/PAM_probabilities_analysis_onc_chRCC_Chopra_TCGA_tumors_models_trained_with_normal"
png(file=paste(file,'.png',sep=''), units="in", width=20, height=15, res=600)
heatmap.3(PAMresults.tumor.trained.with.norm$PosteriorProbabilites,  Rowv=FALSE, Colv=FALSE, dendrogram = "none",  col=brewer.pal(9,"Reds"), RowSideColors=t(rlab), KeyValueName="Classification probability", labRow="", labCol="")
dev.off()
#need to go back to this and label misclassified samples

#Need to change names of samples
#Need to deal with multiple of each tumor sample
###################################################
#select the top ranking CpG sites 
####################################################

genranks=read.table("resultsdir/PAM_oncVchRCC_predictive_CpGs_iterative.txt",  sep="\t", header=T,  fill=T, quote="", comment.char="")
colnames(genranks)=gsub("X.","",colnames(genranks))
genranks=genranks[order(genranks$freq.20., -genranks$av.rank.in.CV.20., decreasing=T),]
library(stringr)
rownames(genranks)=str_replace_all(rownames(genranks), '\"', "")

genes=rownames(genranks)[1:20]

beta_genes_Barcodes=beta[genes,colnames(SampleClass)]
beta_genes_Barcodes=scale(beta_genes_Barcodes)

beta.training.impute.onc.ChRCC.genes=beta.training.impute.onc.ChRCC[genes,]
beta.training.impute.onc.ChRCC.genes.s=scale(beta.training.impute.onc.ChRCC.genes)

NrFolds=10

PAMresults.tumor.trained.with.norm.20g=PAM_Analysis_KB(beta_genes_Barcodes, SampleClass, beta.training.impute.onc.ChRCC.genes.s, NrLambdaFolds=10)
PAMresults.tumor.trained.with.norm$TestPredictions

tab=table(PAMresults.tumor.trained.with.norm.20g$TestPredictions, SampleType.Chopra.onc.KICH.training)
tab=as.matrix(tab)
colnames(tab)=c("chRCC_TCGA","oncocytoma_Chopra")
rownames(tab)=c("Classified_as_normal","Classified_as_chRCC","Classified_as_oncocytoma")
#write.table(tab, file="resultsdir/Calssifications_Chopra_tumor_model_trained_with_normal", sep="\t")

#This depreciates the predictive value of the model
#same number of KICHs in my data as for Chopra study 

#Need to figure out which CpG sites represent seperate regions, to see if selecting unique regions does better 
#find R code for finding the distance of each CpG site to each other CpG site 
#bedtools window
#make bedfiles for CpG sites and find any that are within 5kb of that CpG site
bed=cbind(as.character(genranks$chr.), as.numeric(as.character(genranks$pos.)), as.numeric(as.character(genranks$pos.)), rownames(genranks))
write.table(bed, file="resultsdir/chRCC_oncocytoma_predictive_loci_ranked.bed", sep="\t", row.names=F, col.names=F, quote=F)

dir="resultsdir/"
inputfile="chRCC_oncocytoma_predictive_loci_ranked.bed"

Chop=load("datadir/pdata.testing.rda")
pdata.testing=as.data.frame(pdata.testing)
pdata.testing$TissueType
#26 oncocytomas 6 KICH
summary(as.factor(pdata.testing$TissueType))

pdata.testing[pdata.testing$TissueType=="oncocytoma",]
#They only have oncoytoma tumors 

#unique oncocytomas
summary(as.factor(pdata.testing[pdata.testing$TissueType=="oncocytoma","Pt..No."]))
#19 unique oncocytoma patients 
#19 tumors, of which seven had two tumor samples 
#12 with one tumor sample, seven with two tumor samples 

#Is the training data just TCGA
#26 oncocytomas 6 KICH

length(unique(pdata.testing$Alt.Tumor.id[pdata.testing$Type=="T"]))
#There are some patients with multiple tumor biopsies etc. Not certain these are needle biopsies 

#for oncoytoma patients, what proportion were tumor and what normal
#Oncocytoma
pdata.onc.chopra=pdata.testing[which(pdata.testing$Diagnosis=="oncocytoma" & pdata.testing$TissueType=="oncocytoma"),]
Barcodes.pdata.onc.chopra=rownames(pdata.onc.chopra)

pdata.onc.chopra.norm.adj=pdata.testing[which(pdata.testing$Diagnosis=="oncocytoma" & pdata.testing$TissueType=="Normal"),]
Barcodes.pdata.onc.chopra.norm.adj=rownames(pdata.onc.chopra.norm.adj)

#KICH
pdata.KICH.chopra=pdata.KICH.chopra=pdata.testing[which(pdata.testing$Diagnosis=="KICH" & pdata.testing$TissueType=="KICH"),]
Barcodes.pdata.KICH.chopra=rownames(pdata.KICH.chopra)

pdata.KICH.chopra.norm.adj=pdata.testing[which(pdata.testing$Diagnosis=="KICH" & pdata.testing$TissueType=="Normal"),]
Barcodes.pdata.KICH.chopra.norm.adj=rownames(pdata.KICH.chopra.norm.adj)

pdata.all.norm=pdata.testing[which(pdata.testing$TissueType=="Normal"),]
Barcodes.pdata.biopsy.norm=rownames(pdata.all.norm)

#now get beta value data 
#check if there's missing data 

Chop.beta=load("datadir/Archive_of_OSF_Storage/data/beta.testing.rda")
#Is 450k
#need to impute missing values 
#Don't process the standard way, just impute, as full processing algorithm removes too many samples (where these is more than one biopsy per patient, I think)
MET_Data=beta.testing
#2371 probes with 10% missing values removed, but all samples fine
Genes=rownames(MET_Data)
MET_Data=apply(MET_Data,2,as.numeric)
rownames(MET_Data)=Genes
     
SampleNames=colnames(MET_Data)
     
# more than 10% = removal
MissingValueThreshold=0.1
     
# removing clones with too many missing values
NrMissingsPerGene=apply(MET_Data,1,function(x) sum(is.na(x))) /ncol(MET_Data)
cat("Removing",sum(NrMissingsPerGene>MissingValueThreshold),"genes with more than 10% missing values.\n")
if (sum(NrMissingsPerGene>MissingValueThreshold)>0) MET_Data=MET_Data[NrMissingsPerGene<MissingValueThreshold,]
     
# removing patients with too many missings values     
NrMissingsPerSample=apply(MET_Data,2,function(x) sum(is.na(x))) /nrow(MET_Data)
cat("Removing",sum(NrMissingsPerSample>MissingValueThreshold),"patients with more than 10% missing values.\n")
if (sum(NrMissingsPerSample>MissingValueThreshold)>0) MET_Data=MET_Data[,NrMissingsPerSample<MissingValueThreshold]
     
# knn impute using Tibshirani's method     
k=15
KNNresults=impute.knn(as.matrix(MET_Data),k) # this does not work, gives storage error. ???
MET_Data_KNN=KNNresults$data
beta.testing.impute=MET_Data_KNN

saveRDS(beta.testing.impute, "datadir/Archive_of_OSF_Storage/data/beta.testing.rds")
beta.testing.impute=readRDS("datadir/Archive_of_OSF_Storage/data/beta.testing.rds")

Barcodes.Chopra.oncocytoma.KICH=c(Barcodes.pdata.onc.chopra, Barcodes.pdata.KICH.chopra)
SampleType.Chopra.onc.KICH=c(rep(2, length(Barcodes.pdata.onc.chopra)), rep(1, length(Barcodes.pdata.KICH.chopra)))
names(SampleType.Chopra.onc.KICH)=Barcodes.Chopra.oncocytoma.KICH

#############################################
#PAM analysis trained on Stanford data
################################################

Barcodes=c(BarcodesNormal, BarcodesOncocytoma, BarcodesChromophobe)

SampleClass=rep(NA, length(Barcodes))
names(SampleClass)=Barcodes
SampleClass[names(SampleClass) %in% BarcodesNormal]=1
SampleClass[names(SampleClass) %in% BarcodesChromophobe]=2
SampleClass[names(SampleClass) %in% BarcodesOncocytoma]=3
SampleClass=t(as.matrix(SampleClass))
rownames(SampleClass)="Sample_type"

beta_genes_Barcodes=beta[genes,colnames(SampleClass)]
beta_genes_Barcodes=scale(beta_genes_Barcodes)

beta.training.impute.onc.ChRCC.genes=beta.training.impute.onc.ChRCC[genes,]
beta.training.impute.onc.ChRCC.genes.s=scale(beta.training.impute.onc.ChRCC.genes)

NrFolds=10

PAMresults.tumor.trained.with.norm=PAM_Analysis_KB(beta_genes_Barcodes, SampleClass, beta.training.impute.onc.ChRCC.genes.s, NrLambdaFolds=10)
PAMresults.tumor.trained.with.norm$TestPredictions

tab=table(PAMresults.tumor.trained.with.norm$TestPredictions, SampleType.Chopra.onc.KICH.training)
tab=as.matrix(tab)
colnames(tab)=c("chRCC_TCGA","oncocytoma_Chopra")
rownames(tab)=c("Classified_as_normal","Classified_as_chRCC","Classified_as_oncocytoma")
write.table(tab, file="resultsdir/Calssifications_Chopra_tumor_model_trained_with_normal", sep="\t")
read.table("resultsdir/Calssifications_Chopra_tumor_model_trained_with_normal", sep="\t")


#PAM analsyis trained with normal tissue applied to needle biopsies 
rlab=as.matrix(rownames(PAMresults.tumor.trained.with.norm$TestPredictions))
rownames(rlab)=rlab[,1]
colnames(rlab)=""
rlab[rownames(rlab) %in% Barcodes_chRCC_training_Chopra,1]=cbp2[4]
rlab[rownames(rlab) %in% Barcodes_onc_training_Chopra,1]=cbp2[5]

#Need to also 
file="resultsdir/PAM_probabilities_analysis_onc_chRCC_Chopra_TCGA_tumors_models_trained_with_normal"
png(file=paste(file,'.png',sep=''), units="in", width=20, height=15, res=600)
heatmap.3(PAMresults.tumor.trained.with.norm$PosteriorProbabilites,  Rowv=FALSE, Colv=FALSE, dendrogram = "none",  col=brewer.pal(9,"Reds"), RowSideColors=t(rlab), KeyValueName="Classification probability", labRow="", labCol="")
dev.off()
#need to go back to this and label misclassified samples

##############################################
#Model trained in with normal tissue applied to biopsies. 
##################################################

Barcodes.Chopra.oncocytoma.KICH=c(Barcodes.pdata.biopsy.norm, Barcodes.pdata.KICH.chopra, Barcodes.pdata.onc.chopra)
SampleType.Chopra.onc.KICH=c(rep(1, length(Barcodes.pdata.biopsy.norm)), rep(2, length(Barcodes.pdata.KICH.chopra)), rep(3, length(Barcodes.pdata.onc.chopra)))
names(SampleType.Chopra.onc.KICH)=Barcodes.Chopra.oncocytoma.KICH

beta.testing.impute=readRDS("datadir/Archive_of_OSF_Storage/data/beta.testing.rds")
beta.Chopra.test=beta.testing.impute[genes,Barcodes.Chopra.oncocytoma.KICH]

beta_genes_Barcodes.s=scale(beta_genes_Barcodes)
beta.Chopra.test.s=scale(beta.Chopra.test)

PAMresults.biopsies.trained.with.norm=PAM_Analysis_KB(beta_genes_Barcodes, SampleClass, beta.Chopra.test.s, NrLambdaFolds=10)

#save(PAMresults.biopsies.trained.with.norm, file="resultsdir/PAM_Chopra_needle_core_biopsies_with_normal.RData")
dat=load("resultsdir/PAM_Chopra_needle_core_biopsies_with_normal.RData")

#Need to calculate multiclass AUC
library(pROC)
auc=multiclass.roc(as.factor(SampleType.Chopra.onc.KICH), PAMresults.biopsies.trained.with.norm$TestPredictions[,1])
#Multi-class area under the curve: 0.7337

#Multi class auc if restricted to classification of oncoytomas and chRCC biopsies 
SampleType.Chopra.onc.KICH_tumor_biopsy=SampleType.Chopra.onc.KICH[SampleType.Chopra.onc.KICH=="2"|SampleType.Chopra.onc.KICH=="3"]
TestPredictions_onc.KICH_tumor_biopsy=PAMresults.biopsies.trained.with.norm$TestPredictions[,1][names(SampleType.Chopra.onc.KICH_tumor_biopsy)]
auc=multiclass.roc(as.factor(SampleType.Chopra.onc.KICH_tumor_biopsy), TestPredictions_onc.KICH_tumor_biopsy)

table(SampleType.Chopra.onc.KICH_tumor_biopsy, TestPredictions_onc.KICH_tumor_biopsy)

#Can't do Multi class auc if restricted to classification of NKP biopsies, as there is only one class of sample, and only two classes in predictions

#

misclassifications=as.data.frame(cbind(PAMresults.biopsies.trained.with.norm$TestPredictions, SampleType.Chopra.onc.KICH))
BarcodesMiclassificationsChopraBiopsyNorm=rownames(misclassifications[misclassifications$SampleType.Chopra.onc.KICH==1 & misclassifications$PAMprediction!=1,])
BarcodesMiclassificationsChopraBiopsychRCC=rownames(misclassifications[misclassifications$SampleType.Chopra.onc.KICH==2 & misclassifications$PAMprediction!=2,])
BarcodesMiclassificationsChopraBiopsyonc=rownames(misclassifications[misclassifications$SampleType.Chopra.onc.KICH==3 & misclassifications$PAMprediction!=3,])

tab=table(PAMresults.biopsies.trained.with.norm$TestPredictions, SampleType.Chopra.onc.KICH)
tab=as.matrix(tab)
colnames(tab)=c("biopsy NKP","biopsy chRCC","biopsy oncocytoma")
rownames(tab)=c("Classified_as_normal","Classified_as_chRCC","Classified_as_oncocytoma")
#write.table(tab, file="resultsdir/Calssifications_Chopra_biopsies_model_trained_with_normal", sep="\t")

rlab=as.matrix(rownames(PAMresults.biopsies.trained.with.norm$TestPredictions))
rownames(rlab)=rlab[,1]
colnames(rlab)=""
rlab[rownames(rlab) %in% Barcodes.pdata.KICH.chopra,1]=cbp2[4]
rlab[rownames(rlab) %in% Barcodes.pdata.onc.chopra,1]=cbp2[5]
rlab[rownames(rlab) %in% Barcodes.pdata.biopsy.norm,1]=cbp2[3]

#Adding patient number so that we can determine if at least one sample was correctly classified 
#Also def need to add normal adjacent tissue
PosteriorProbabilites.biopsy=as.data.frame(PAMresults.biopsies.trained.with.norm$PosteriorProbabilites)
PosteriorProbabilites.biopsy$patient.no=rep(NA, nrow(PosteriorProbabilites.biopsy))
PosteriorProbabilites.biopsy$patient.no=unlist(lapply(rownames(PosteriorProbabilites.biopsy), function(x) pdata.testing[x,"Pt..No."]))
PosteriorProbabilites.biopsy$TissueType=unlist(lapply(rownames(PosteriorProbabilites.biopsy), function(x) pdata.testing[x,"TissueType"]))
PosteriorProbabilites.biopsy$TissueType=as.factor(PosteriorProbabilites.biopsy$TissueType)
levels(PosteriorProbabilites.biopsy$TissueType)=c(2,1,3)
PosteriorProbabilites.biopsy$TissueType=as.character(PosteriorProbabilites.biopsy$TissueType)
PosteriorProbabilites.biopsy$TestPredictions=PAMresults.biopsies.trained.with.norm$TestPredictions[,"PAMprediction"]

PosteriorProbabilites.biopsy=PosteriorProbabilites.biopsy[order(PosteriorProbabilites.biopsy$TissueType, PosteriorProbabilites.biopsy$patient.no),]
rlab=rlab[match(rownames(PosteriorProbabilites.biopsy), rownames(rlab))]

#Need to also 
#file="resultsdir/PAM_probabilities_analysis_onc_chRCC_Chopra_TCGA_biopsies_models_trained_with_normal"
#png(file=paste(file,'.png',sep=''), units="in", width=3, height=9, res=600)
#heatmap.3(PosteriorProbabilites.biopsy[,1:3],  Rowv=FALSE, Colv=FALSE, dendrogram = "none",  col=brewer.pal(9,"Reds"), RowSideColors=t(rlab), KeyValueName="Classification probability", RowSideColorsSize=4, labRow=PosteriorProbabilites.biopsy$patient.no, labCol="")
#dev.off()

file="resultsdir/PAM_probabilities_analysis_onc_chRCC_Chopra_TCGA_biopsies_models_trained_with_normal"
pdf(file=paste(file,'.pdf',sep=''), height = 9.5,  width = 3, family = "Helvetica")
heatmap.3(PosteriorProbabilites.biopsy[,1:3],  Rowv=FALSE, Colv=FALSE, dendrogram = "none",  col=brewer.pal(9,"Reds"), RowSideColors=t(rlab), KeyValueName="Classification probability", RowSideColorsSize=4, labRow=PosteriorProbabilites.biopsy$patient.no, labCol="")
dev.off()

###################################
#figure out how many samples with multiple biopsies, and how discordant
######################################################

#oncocytoma
PosteriorProbabilites.biopsy.onc=PosteriorProbabilites.biopsy[PosteriorProbabilites.biopsy$TissueType==3,]

length(unique(PosteriorProbabilites.biopsy.onc$patient.no))
#19 unique biopsies 
length(PosteriorProbabilites.biopsy.onc$patient.no)
#26 biopsies

PosteriorProbabilites.biopsy.onc[duplicated(PosteriorProbabilites.biopsy.onc$patient.no),"patient.no"]
#23, 34, 53, 56, 18, 47, 63

#only one sample
setdiff(unique(PosteriorProbabilites.biopsy.onc$patient.no), PosteriorProbabilites.biopsy.onc[duplicated(PosteriorProbabilites.biopsy.onc$patient.no),"patient.no"])
#"10"      "105-339" "105-340" "106-341" "111-358" "116-372" "116-373" "119-381" "17"      "81-258"  "93-303"  "93-304"

#onc
PosteriorProbabilites.biopsy.onc.correct=PosteriorProbabilites.biopsy.onc[PosteriorProbabilites.biopsy.onc$TissueType==3 & PosteriorProbabilites.biopsy.onc$TestPredictions==3,"patient.no"]
PosteriorProbabilites.biopsy.onc.incorrect=PosteriorProbabilites.biopsy.onc[PosteriorProbabilites.biopsy.onc$TissueType==3 & PosteriorProbabilites.biopsy.onc$TestPredictions!=3,"patient.no"]
#14/26 correct

discordant=intersect(PosteriorProbabilites.biopsy.onc.correct, PosteriorProbabilites.biopsy.onc.incorrect)
PosteriorProbabilites.biopsy.onc[PosteriorProbabilites.biopsy.onc$patient.no %in% discordant,]
#Patient no 18, 47, 63

PosteriorProbabilites.biopsy.onc.correct[duplicated(PosteriorProbabilites.biopsy.onc.correct)]
#Patient no 23, 34, 53, 56
PosteriorProbabilites.biopsy.onc.incorrect[duplicated(PosteriorProbabilites.biopsy.onc.incorrect)]

#of 12 incorrect biopsies, three were belonged to a patient with two biopsies; in all three cases the other biopsy was correcttly classified
#Remove the 47 label in MDS plot that is near tumors. This one is the correct one

#chRCC
PosteriorProbabilites.biopsy.chRCC=PosteriorProbabilites.biopsy[PosteriorProbabilites.biopsy$TissueType==2,]
length(PosteriorProbabilites.biopsy.chRCC$patient.no)
#6 altogether
length(unique(PosteriorProbabilites.biopsy.chRCC$patient.no))
#4 unique biopsies 

#4/6 correct

PosteriorProbabilites.biopsy.chRCC.correct=PosteriorProbabilites.biopsy.chRCC[PosteriorProbabilites.biopsy.chRCC$TissueType==2 & PosteriorProbabilites.biopsy.chRCC$TestPredictions==2,"patient.no"]
PosteriorProbabilites.biopsy.chRCC.incorrect=PosteriorProbabilites.biopsy.chRCC[PosteriorProbabilites.biopsy.chRCC$TissueType==2 & PosteriorProbabilites.biopsy.chRCC$TestPredictions!=2,"patient.no"]

discordant=intersect(PosteriorProbabilites.biopsy.chRCC.correct, PosteriorProbabilites.biopsy.chRCC.incorrect)
PosteriorProbabilites.biopsy.chRCC[PosteriorProbabilites.biopsy.chRCC$patient.no %in% discordant,]
#49

PosteriorProbabilites.biopsy.chRCC.correct[duplicated(PosteriorProbabilites.biopsy.chRCC.correct)]
#20

#Normal tissue 
PosteriorProbabilites.biopsy.normal=PosteriorProbabilites.biopsy[PosteriorProbabilites.biopsy$TissueType==1,]
length(PosteriorProbabilites.biopsy.normal$patient.no)
#101 normals

length(unique(PosteriorProbabilites.biopsy.normal$patient.no))
#all are unique
PosteriorProbabilites.biopsy.normal$patient.no[duplicated(PosteriorProbabilites.biopsy.normal$patient.no),]

#Get the size of the FFPE oncocytoma tumors 
onc.size=read.table("datadir/Archive_of_OSF_Storage/data/Chopra_FFPE_sample_info_tumor_size.txt", sep="\t", header=T)
#This is for the FFPE samples
#
PosteriorProbabilites.biopsy.chRCC=PosteriorProbabilites.biopsy[PosteriorProbabilites.biopsy$TissueType==2,]
patient.no.dup=PosteriorProbabilites.biopsy.chRCC$patient.no[duplicated(PosteriorProbabilites.biopsy.chRCC$patient.no)]
PosteriorProbabilites.biopsy.chRCC.dup=PosteriorProbabilites.biopsy.chRCC[PosteriorProbabilites.biopsy.chRCC$patient.no %in% patient.no.dup,]

chRCC.dup.discordant=PosteriorProbabilites.biopsy.chRCC.dup[PosteriorProbabilites.biopsy.chRCC.dup$TestPredictions!="2",]
#Patient no 49
chRCC.dup.concordant=setdiff(patient.no.dup, chRCC.dup.discordant$patient.no)
#patient 20

PosteriorProbabilites.biopsy.normal=PosteriorProbabilites.biopsy[PosteriorProbabilites.biopsy$TissueType==1,]
patient.no.dup=PosteriorProbabilites.biopsy.normal$patient.no[duplicated(PosteriorProbabilites.biopsy.normal$patient.no)]
PosteriorProbabilites.biopsy.normal.dup=PosteriorProbabilites.biopsy.normal[PosteriorProbabilites.biopsy.normal$patient.no %in% patient.no.dup,]

normal.dup.discordant=PosteriorProbabilites.biopsy.normal.dup[PosteriorProbabilites.biopsy.normal.dup$TestPredictions!="1",]
#None
normal.dup.concordant=setdiff(patient.no.dup, normal.dup.discordant$patient.no)
#None

#############################################################################
#MDS plot that includes Chopra samples
#########################################################################

#Load methylation data for TCGA RCCs KIRC, KIRP, KICH

#Order KIRC, KIRP, KICH
dat=load("datadir/MET_KIRP_Processed.Rdata")
METcancer_KIRP=ProcessedData$MET_Data_Cancer
METnormal_KIRP=ProcessedData$MET_Data_Normal
rm(ProcessedData)
METcancer_KIRP=TCGA_GENERIC_CleanUpSampleNames(METcancer_KIRP, 12) 
#n=275
METnormal_KIRP=TCGA_GENERIC_CleanUpSampleNames(METnormal_KIRP, 12) 
#n=44

#KICH
dat=load("datadir/MET_KICH_Processed.Rdata")
METcancer_KICH=ProcessedData$MET_Data_Cancer
#No normal data for KICH
rm(ProcessedData)
METcancer_KICH=TCGA_GENERIC_CleanUpSampleNames(METcancer_KICH, 12) 
#n=65 

#KIRC
dat=load("datadir/ProcessedData_KIRC.RData")
METcancer_KIRC=ProcessedData$MET_Data_Cancer
METnormal_KIRC=ProcessedData$MET_Data_Normal
rm(ProcessedData)

METcancer_KIRC=TCGA_GENERIC_CleanUpSampleNames(METcancer_KIRC, 12) 
#n=315
METnormal_KIRC=TCGA_GENERIC_CleanUpSampleNames(METnormal_KIRC, 12) 
#n=160

#Combine normal tissues from the two cancers 
OverlapProbes=intersect(rownames(METnormal_KIRP),rownames(METnormal_KIRC))
METnormal_KIRP=METnormal_KIRP[OverlapProbes,]
METnormal_KIRC=METnormal_KIRC[OverlapProbes,]

METnormal=cbind(METnormal_KIRP, METnormal_KIRC)

#Get probes that are shared across all of the TCGA cancers
commonvarcpgs=intersect(
     intersect(rownames(METcancer_KIRP), intersect(rownames(METcancer_KIRC), rownames(METcancer_KICH))),
     rownames(METnormal))
commonvarcpgs=commonvarcpgs[-c(grep("rs", commonvarcpgs))]

#restrict meth to just these common variable CpGs
METcancer_KIRP=METcancer_KIRP[commonvarcpgs,]
METcancer_KIRC=METcancer_KIRC[commonvarcpgs,]
METcancer_KICH=METcancer_KICH[commonvarcpgs,]
METnormal_KIRP=METnormal_KIRP[commonvarcpgs,]
METnormal_KIRC=METnormal_KIRC[commonvarcpgs,]

#MET
colnames(METcancer_KIRP)=paste(colnames(METcancer_KIRP),"-01", sep="")
colnames(METcancer_KIRC)=paste(colnames(METcancer_KIRC),"-01", sep="")
colnames(METcancer_KICH)=paste(colnames(METcancer_KICH),"-01", sep="")
colnames(METnormal_KIRP)=paste(colnames(METnormal_KIRP),"-11", sep="")
colnames(METnormal_KIRC)=paste(colnames(METnormal_KIRC),"-11", sep="")

Barcodes_METcancer_KIRP=colnames(METcancer_KIRP)
Barcodes_METcancer_KIRC=colnames(METcancer_KIRC)
Barcodes_METcancer_KICH=colnames(METcancer_KICH)
Barcodes_METnormal_KIRP=colnames(METnormal_KIRP)
Barcodes_METnormal_KIRC=colnames(METnormal_KIRC)
Barcodes_METnormal=c(Barcodes_METnormal_KIRC, Barcodes_METnormal_KIRP)

METall=cbind(
     METcancer_KIRP,
     METcancer_KIRC,
     METcancer_KICH,
     METnormal_KIRP,
     METnormal_KIRC)

#
#Get imputed Chopra needle biopsy test set data
beta.testing.impute=readRDS("datadir/Archive_of_OSF_Storage/data/beta.testing.rds")

#Using 46530# 247624 CpG probes
Barcodes.Chopra.biopsy.onc=rownames(pdata.testing[pdata.testing$TissueType=="oncocytoma",])
Barcodes.Chopra.biopsy.KICH=rownames(pdata.testing[pdata.testing$TissueType=="KICH",])
Barcodes.Chopra.biopsy.KIRC=rownames(pdata.testing[pdata.testing$TissueType=="KIRC",])
Barcodes.Chopra.biopsy.KIRP=rownames(pdata.testing[pdata.testing$TissueType=="KIRP",])
Barcodes.Chopra.biopsy.Normal=rownames(pdata.testing[pdata.testing$TissueType=="Normal",])

#tumor data for onc from Chopra
Barcodes.Chopra.tumor.onc=colnames(beta.training.impute.onc)
#Here's where 

#########################################################################
#Version of MDS plot from needle biopsies (testing)
##########################################################################

cbp2=c("#000000", "#E69F00", "#56B4E9", "#009E73","#FD61D1", "#A3A500", "#D55E00", "#CC79A7")

BarcodesMiclassificationsChopraBiopsyNorm=rownames(misclassifications[misclassifications$SampleType.Chopra.onc.KICH==1 & misclassifications$PAMprediction!=1,])
BarcodesMiclassificationsChopraBiopsychRCC=rownames(misclassifications[misclassifications$SampleType.Chopra.onc.KICH==2 & misclassifications$PAMprediction!=2,])
BarcodesMiclassificationsChopraBiopsyonc=rownames(misclassifications[misclassifications$SampleType.Chopra.onc.KICH==3 & misclassifications$PAMprediction!=3,])
#

patients.dup.onc=pdata.testing$Pt..No.[pdata.testing$TissueType=="oncocytoma"][duplicated(pdata.testing$Pt..No.[pdata.testing$TissueType=="oncocytoma"])]
patients.dup.chRCC=pdata.testing$Pt..No.[pdata.testing$TissueType=="KICH"][duplicated(pdata.testing$Pt..No.[pdata.testing$TissueType=="KICH"])]
#patients.dup.NKP=pdata.testing$Pt..No.[pdata.testing$TissueType=="Normal"][duplicated(pdata.testing$Pt..No.[pdata.testing$TissueType=="Normal"])]
#No duplicated normal tissue biopsies 
patients.dup.KIRP=pdata.testing$Pt..No.[pdata.testing$TissueType=="KIRP"][duplicated(pdata.testing$Pt..No.[pdata.testing$TissueType=="KIRP"])]
patients.dup.KIRC=pdata.testing$Pt..No.[pdata.testing$TissueType=="KIRC"][duplicated(pdata.testing$Pt..No.[pdata.testing$TissueType=="KIRC"])]

summary(as.factor(pdata.testing[pdata.testing$TissueType=="oncocytoma" & pdata.testing$Pt..No. %in% patients.dup.onc,"Pt..No."]))
#Seven oncocytoma patients with duplicate biopsies 
summary(as.factor(pdata.testing[pdata.testing$TissueType=="KICH" & pdata.testing$Pt..No. %in% patients.dup.chRCC,"Pt..No."]))
#Two chromophobe biopsies with duplicate biopsies
summary(as.factor(pdata.testing[pdata.testing$TissueType=="KIRP" & pdata.testing$Pt..No. %in% patients.dup.KIRP,"Pt..No."]))
#Four KIRP patients with duplicate biopsies
summary(as.factor(pdata.testing[pdata.testing$TissueType=="KIRC" & pdata.testing$Pt..No. %in% patients.dup.KIRC,"Pt..No."]))
#20 KIRP patients with duplicate biopsies

#
pdata.dup.onc=pdata.testing[pdata.testing$TissueType=="oncocytoma" & pdata.testing$Pt..No. %in% patients.dup.onc,c("Pt..No.","Diagnosis")]
pdata.dup.onc$x=x[rownames(pdata.dup.onc)]
pdata.dup.onc$y=y[rownames(pdata.dup.onc)]
pdata.dup.onc$FirstSecond=rep(c(1:2), length(unique(pdata.dup.onc$Pt..No.)))

#make matrix of coordinates to draw lines between duplicate biopsies for each patient
library(abind)
library(graphics)

patients=unique(pdata.dup.onc$Pt..No.)
coords=lapply(patients, function(pat)
#First biopsy for patient x,y
c(pdata.dup.onc[pdata.dup.onc$Pt..No. %in% pat & pdata.dup.onc$FirstSecond==1,"x"],
pdata.dup.onc[pdata.dup.onc$Pt..No. %in% pat & pdata.dup.onc$FirstSecond==1,"y"],
#Second biopsy for patient x,y
pdata.dup.onc[pdata.dup.onc$Pt..No. %in% pat & pdata.dup.onc$FirstSecond==2,"x"],
pdata.dup.onc[pdata.dup.onc$Pt..No. %in% pat & pdata.dup.onc$FirstSecond==2,"y"]))
names(coords)=patients

mat1=as.data.frame(t(abind(coords,along=2)))
colnames(mat1)=c("X0","Y0","X1","Y1")
mat1$col=rep(cbp2[5], nrow(mat1))

#
pdata.dup.chRCC=pdata.testing[pdata.testing$TissueType=="KICH" & pdata.testing$Pt..No. %in% patients.dup.chRCC,c("Pt..No.","Diagnosis")]
pdata.dup.chRCC$x=x[rownames(pdata.dup.chRCC)]
pdata.dup.chRCC$y=y[rownames(pdata.dup.chRCC)]
pdata.dup.chRCC$FirstSecond=rep(c(1:2), length(unique(pdata.dup.chRCC$Pt..No.)))

patients=unique(pdata.dup.chRCC$Pt..No.)
coords=lapply(patients, function(pat)
#First biopsy for patient x,y
c(pdata.dup.chRCC[pdata.dup.chRCC$Pt..No. %in% pat & pdata.dup.chRCC$FirstSecond==1,"x"],
pdata.dup.chRCC[pdata.dup.chRCC$Pt..No. %in% pat & pdata.dup.chRCC$FirstSecond==1,"y"],
#Second biopsy for patient x,y
pdata.dup.chRCC[pdata.dup.chRCC$Pt..No. %in% pat & pdata.dup.chRCC$FirstSecond==2,"x"],
pdata.dup.chRCC[pdata.dup.chRCC$Pt..No. %in% pat & pdata.dup.chRCC$FirstSecond==2,"y"]))
names(coords)=patients

mat2=as.data.frame(t(abind(coords,along=2)))
colnames(mat2)=c("X0","Y0","X1","Y1")
mat2$col=rep(cbp2[3], nrow(mat2))

#
pdata.dup.KIRC=pdata.testing[pdata.testing$TissueType=="KIRC" & pdata.testing$Pt..No. %in% patients.dup.KIRC,c("Pt..No.","Diagnosis")]
pdata.dup.KIRC$x=x[rownames(pdata.dup.KIRC)]
pdata.dup.KIRC$y=y[rownames(pdata.dup.KIRC)]
pdata.dup.KIRC$FirstSecond=rep(c(1:2), length(unique(pdata.dup.KIRC$Pt..No.)))

patients=unique(pdata.dup.KIRC$Pt..No.)
coords=lapply(patients, function(pat)
#First biopsy for patient x,y
c(pdata.dup.KIRC[pdata.dup.KIRC$Pt..No. %in% pat & pdata.dup.KIRC$FirstSecond==1,"x"],
pdata.dup.KIRC[pdata.dup.KIRC$Pt..No. %in% pat & pdata.dup.KIRC$FirstSecond==1,"y"],
#Second biopsy for patient x,y
pdata.dup.KIRC[pdata.dup.KIRC$Pt..No. %in% pat & pdata.dup.KIRC$FirstSecond==2,"x"],
pdata.dup.KIRC[pdata.dup.KIRC$Pt..No. %in% pat & pdata.dup.KIRC$FirstSecond==2,"y"]))
names(coords)=patients

mat3=as.data.frame(t(abind(coords,along=2)))
colnames(mat3)=c("X0","Y0","X1","Y1")
mat3$col=rep(cbp2[2], nrow(mat3))

#
pdata.dup.KIRP=pdata.testing[pdata.testing$TissueType=="KIRP" & pdata.testing$Pt..No. %in% patients.dup.KIRP,c("Pt..No.","Diagnosis")]
pdata.dup.KIRP$x=x[rownames(pdata.dup.KIRP)]
pdata.dup.KIRP$y=y[rownames(pdata.dup.KIRP)]
pdata.dup.KIRP$FirstSecond=rep(c(1:2), length(unique(pdata.dup.KIRP$Pt..No.)))

patients=unique(pdata.dup.KIRP$Pt..No.)
coords=lapply(patients, function(pat)
#First biopsy for patient x,y
c(pdata.dup.KIRP[pdata.dup.KIRP$Pt..No. %in% pat & pdata.dup.KIRP$FirstSecond==1,"x"],
pdata.dup.KIRP[pdata.dup.KIRP$Pt..No. %in% pat & pdata.dup.KIRP$FirstSecond==1,"y"],
#Second biopsy for patient x,y
pdata.dup.KIRP[pdata.dup.KIRP$Pt..No. %in% pat & pdata.dup.KIRP$FirstSecond==2,"x"],
pdata.dup.KIRP[pdata.dup.KIRP$Pt..No. %in% pat & pdata.dup.KIRP$FirstSecond==2,"y"]))
names(coords)=patients

mat4=as.data.frame(t(abind(coords,along=2)))
colnames(mat4)=c("X0","Y0","X1","Y1")
mat4$col=rep(cbp2[1], nrow(mat4))

matall=rbind(mat1, mat2, mat3, mat4)

#############
#run MDS
#############

#restrict to probes that were used to classify TCGA cancer type in TCGA
file=load("resultsdir/PAM_PanKidney_CrossValidation_MethylMix0.2_probes.RData")
NonZeroGenes=unique(unlist(PAMresults$Genes))
#6725 CpG sites

#try it just with the probes that are also MethylMix genes
METall_NonZeroGenes=METall[NonZeroGenes,]

ClassCentroids_Hand=cbind(rowMeans(METall_NonZeroGenes[,Barcodes_METcancer_KIRP]),
                          rowMeans(METall_NonZeroGenes[,Barcodes_METcancer_KIRC]),
                          rowMeans(METall_NonZeroGenes[,Barcodes_METcancer_KICH]),
                          rowMeans(METall_NonZeroGenes[,Barcodes_METnormal]))

OverlapProbes=intersect(rownames(METall), rownames(ClassCentroids_Hand))
METall_CD=METall[OverlapProbes,]
ClassCentroids_METall=ClassCentroids_Hand[OverlapProbes,]

OverlapProbes=intersect(rownames(beta.testing.impute),intersect(rownames(beta.training.impute.onc),intersect(rownames(beta), intersect(rownames(METall), rownames(ClassCentroids_Hand)))))
m1=ClassCentroids_Hand[OverlapProbes,]
m2=METall[OverlapProbes,]
m3=beta[OverlapProbes,]
m4=beta.training.impute.onc[OverlapProbes,]
m5=beta.testing.impute[OverlapProbes,]

d1=t(cbind(m1, m2, m3, m4, m5))
d2=dist(d1)

fit <- cmdscale(d2,eig=TRUE, k=2) # k is the number of dim

x <- fit$points[,1]
y <- fit$points[,2]

file="resultsdir/MDS_Pan_kidney_with_Chopra_biopsy_label_misclassifications"
png(file=paste(file,'.png',sep=''), units="in", width=10, height=10, res=600)
plot(x[Barcodes_METcancer_KIRP], y[Barcodes_METcancer_KIRP], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[1],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1)
par(new=TRUE)
plot(x[Barcodes_METcancer_KIRC], y[Barcodes_METcancer_KIRC], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[2],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[Barcodes_METcancer_KICH], y[Barcodes_METcancer_KICH], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[3],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[Barcodes_METnormal], y[Barcodes_METnormal], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[4],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[1:4],y[1:4], xlab="Coordinate 1", ylab="Coordinate 2", col=c(cbp2[1:4]), pch=15, cex=1.5, , xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), main="")
par(new=TRUE)
plot(x[BarcodesClearCell], y[BarcodesClearCell], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[2],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19)
par(new=TRUE)
plot(x[BarcodesChromophobe], y[BarcodesChromophobe], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[3],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[BarcodesNormal], y[BarcodesNormal], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[4],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[BarcodesOncocytoma], y[BarcodesOncocytoma], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[5],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[Barcodes.Chopra.biopsy.Normal], y[Barcodes.Chopra.biopsy.Normal], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[4],alpha.f=1), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=8, main="")
par(new=TRUE)
plot(x[Barcodes.Chopra.biopsy.onc], y[Barcodes.Chopra.biopsy.onc], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[5],alpha.f=1), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=8, main="")
par(new=TRUE)
plot(x[Barcodes.Chopra.biopsy.KICH], y[Barcodes.Chopra.biopsy.KICH], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[3],alpha.f=1), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=8, main="")
par(new=TRUE)
plot(x[Barcodes.Chopra.biopsy.KIRP], y[Barcodes.Chopra.biopsy.KIRP], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[1],alpha.f=1), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=8, main="")
par(new=TRUE)
plot(x[Barcodes.Chopra.biopsy.KIRC], y[Barcodes.Chopra.biopsy.KIRC], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[2],alpha.f=1), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=8, main="")
par(new=TRUE)
text(x[intersect(BarcodesMiclassificationsChopraBiopsyonc, Barcodes.Chopra.biopsy.onc)],y[intersect(BarcodesMiclassificationsChopraBiopsyonc, Barcodes.Chopra.biopsy.onc)], labels=PosteriorProbabilites.biopsy[intersect(BarcodesMiclassificationsChopraBiopsyonc, Barcodes.Chopra.biopsy.onc),"patient.no"], cex=0.9, font=2, col=cbp2[5], adj=c(1,1))
text(x[intersect(BarcodesMiclassificationsChopraBiopsyNorm, Barcodes.Chopra.biopsy.Normal)],y[intersect(BarcodesMiclassificationsChopraBiopsyNorm, Barcodes.Chopra.biopsy.Normal)], labels=PosteriorProbabilites.biopsy[intersect(BarcodesMiclassificationsChopraBiopsyNorm, Barcodes.Chopra.biopsy.Normal),"patient.no"], cex=0.9, font=2, col=cbp2[4], adj=c(1,1))
text(x[intersect(BarcodesMiclassificationsChopraBiopsychRCC, Barcodes.Chopra.biopsy.KICH)],y[intersect(BarcodesMiclassificationsChopraBiopsychRCC, Barcodes.Chopra.biopsy.KICH)], labels=PosteriorProbabilites.biopsy[intersect(BarcodesMiclassificationsChopraBiopsychRCC, Barcodes.Chopra.biopsy.KICH),"patient.no"], cex=0.9, font=2, col=cbp2[3], adj=c(1,1))
for(i in 1:nrow(matall)){
segments(matall[i,"X0"], matall[i,"Y0"], matall[i,"X1"], matall[i,"Y1"], col=as.character(matall$col[i]), lty=2)
}
dev.off()

##################
#Make table showing classification of all Chopra samples, with patient identifier for each patient 
#################

#get PAM model results for biopsies 
dat=load("resultsdir/PAM_Chopra_needle_core_biopsies_with_normal.RData")
TestPredictions=as.data.frame(PAMresults.biopsies.trained.with.norm$TestPredictions, drop=F)
TestPredictions$Patient.No=unlist(lapply(rownames(TestPredictions), function(x) pdata.testing[x,"Pt..No."]))
TestPredictions$Diagnosis=unlist(lapply(rownames(TestPredictions), function(x) pdata.testing[x,"Diagnosis"]))
TestPredictions$TissueType=unlist(lapply(rownames(TestPredictions), function(x) pdata.testing[x,"TissueType"]))
TestPredictions$PAMprediction=revalue(as.factor(TestPredictions$PAMprediction), c("1"="NKP","2"="chRCC/KICH","3"="oncocytoma"))
#TestPredictions$Diagnosis=revalue(as.factor(TestPredictions$Diagnosis), c("AML"="Angiomylolipoma", "KICH"="chRCC/KICH", "KIRC"="ccRCC/KIRC", "KIRP"="pRCC/KIRP","Normal"="NKP","oncocytoma"="oncocytoma"))
TestPredictions$TissueType=revalue(as.factor(TestPredictions$TissueType), c("AML"="Angiomylolipoma", "KICH"="chRCC/KICH", "KIRC"="ccRCC/KIRC", "KIRP"="pRCC/KIRP","Normal"="NKP","oncocytoma"="oncocytoma"))
TestPredictions$Diagnosis=revalue(as.factor(TestPredictions$Diagnosis), c("NKP"="Normal/benign lesion"))

#Indicate which patients are duplicated for the same tissue type 
TestPredictions$Duplicate=rep("No", nrow(TestPredictions))

which(duplicated(TestPredictions[TestPredictions$TissueType=="NKP","Patient.No"]))
TestPredictions[which(duplicated(TestPredictions[TestPredictions$TissueType=="oncocytoma","Patient.No"])),]

patients.dup.onc=TestPredictions[TestPredictions$TissueType=="oncocytoma",][duplicated(TestPredictions[TestPredictions$TissueType=="oncocytoma","Patient.No"]),"Patient.No"]
TestPredictions[TestPredictions$Patient.No %in% patients.dup.onc & TestPredictions$TissueType=="oncocytoma", "Duplicate"]="Yes"

patients.dup.chRCC=TestPredictions[TestPredictions$TissueType=="chRCC/KICH",][duplicated(TestPredictions[TestPredictions$TissueType=="chRCC/KICH","Patient.No"]),"Patient.No"]
TestPredictions[TestPredictions$Patient.No %in% patients.dup.chRCC & TestPredictions$TissueType=="chRCC/KICH", "Duplicate"]="Yes"

#Is the classification correct 
TestPredictions$Misclassified=rep("No", nrow(TestPredictions))
TestPredictions[TestPredictions$TissueType=="oncocytoma" & TestPredictions$PAMprediction!="oncocytoma","Misclassified"]="Yes"
TestPredictions[TestPredictions$TissueType=="oncocytoma" & TestPredictions$PAMprediction!="oncocytoma",]

TestPredictions[TestPredictions$TissueType=="chRCC/KICH" & TestPredictions$PAMprediction!="chRCC/KICH","Misclassified"]="Yes"
TestPredictions[TestPredictions$TissueType=="chRCC/KICH" & TestPredictions$PAMprediction!="chRCC/KICH",]

TestPredictions[TestPredictions$TissueType=="NKP" & TestPredictions$PAMprediction!="NKP","Misclassified"]="Yes"
TestPredictions[TestPredictions$TissueType=="NKP" & TestPredictions$PAMprediction!="NKP",]

#Add in the posterior probabilities for classification of each sample
PosteriorProbabilites=PAMresults.biopsies.trained.with.norm$PosteriorProbabilites
TestPredictions=cbind(TestPredictions, PosteriorProbabilites)
#Add rownames as sample ID, since this is the only sample-specific ID
TestPredictions=cbind(TestPredictions, rownames(TestPredictions))
#Where is the sample number that distinguished the two biopsies? 
TestPredictions=TestPredictions[,c(10,2,3,4,5,1,6,7,8,9)]
colnames(TestPredictions)=c("Sample ID","Patient ID","Diagnosis","Tissue type", "Dupilicate biopsy","Diagnostic model prediction", "Misclassified by diagnostic model?","Posterior probability NKP", "Posterior probability chRCC","Posterior probability oncocytoma")

write.table(TestPredictions, file="resultsdir/PAM_Chopra_needle_core_biopsies.txt", sep="\t", row.names = FALSE)

########################################
#Analyzing Chopra tumor biopsy classifications
##########################################

dim(TestPredictions[TestPredictions$TissueType=="NKP",])
dim(TestPredictions[TestPredictions$TissueType=="chRCC/KICH",])

#101 NKP were available, which were collected from patients with RCC, oncoytoma, angiomylolipoma, or from a cancer-free individual with an unspecified benign lesion. 
#6 chRCC biosies from four patients (duplicate biopsies from two patients)
#26 oncoytoma biopsies from 19 unique patient tumors (duplicate biopsies for seven patients)

TestPredictions[TestPredictions$TissueType=="NKP",]

dim(TestPredictions[TestPredictions$TissueType=="oncocytoma",])

dim(TestPredictions[TestPredictions$TissueType=="chRCC/KICH" & TestPredictions$Duplicate=="Yes",])

length(unique(TestPredictions[TestPredictions$TissueType=="chRCC/KICH","Patient.No"]))
length(unique(TestPredictions[TestPredictions$TissueType=="chRCC/KICH" & TestPredictions$Duplicate=="Yes","Patient.No"]))

length(unique(TestPredictions[TestPredictions$TissueType=="oncocytoma","Patient.No"]))
length(unique(TestPredictions[TestPredictions$TissueType=="oncocytoma" & TestPredictions$Duplicate=="Yes","Patient.No"]))

#misclassification of oncoytoma
dim(TestPredictions[TestPredictions$TissueType=="oncocytoma" & TestPredictions$Misclassified=="Yes",])
length(unique(TestPredictions[TestPredictions$TissueType=="oncocytoma" & TestPredictions$Misclassified=="Yes","Patient.No"]))

TestPredictions[TestPredictions$TissueType=="oncocytoma" & TestPredictions$Misclassified=="Yes"  & TestPredictions$Duplicate=="Yes",]

TestPredictions[TestPredictions$TissueType=="oncocytoma" & TestPredictions$Duplicate=="Yes",]

#Two of six chRCCs, duplicate biopsies from the same patient tumor, were misclassified as oncocytoma. 
#Three of 101 chRCCs were misclassified, all as chRCC. 
#12 of 26 oncocytoma biopsies were misclassified, all as NKP. The 12 misclassified were from 9 unique oncoytoma patients.
#Of seven patients with duplicate oncoytoma biopsies, two had both biopsies correctly classified, three had both biopsies misclassified as NKP, while two more were correctly clasified for one biopsy and misclassified as NKP for the other
#Therefore, it appears that collection of duplicate needle biopsies slightly increases the chance of, but does not guarantee detection of oncoytoma, although the sample size for duplicate sampeles is isufficient to draw strong conclusions.

########################################
#Make supplementary table showing classificiation of tumor samples from the Chopra et al study
##########################################

labs=read.table("datadir/Chopra_oncocytoma_tumor_labels.txt", sep="\t", header=T)
rownames(labs)=labs[,1]
labs$Sample=substr(labs$Sample,1,12)


Chop2=load("datadir/Archive_of_OSF_Storage/data/pdata.training.rda")
pdata.training=as.data.frame(pdata.training)
pdata.training$Type
summary(as.factor(pdata.training$TissueType))
pdata.training[pdata.training$Diagnosis=="KICH",]

PAMresults.tumor.trained.with.norm=readRDS("datadir/Archive_of_OSF_Storage/data/PAM_Chopra_tumors.rds")

TestPredictions=as.data.frame(PAMresults.tumor.trained.with.norm$TestPredictions, drop=F)
TestPredictions$Patient.No=unlist(lapply(rownames(TestPredictions), function(x) pdata.training[x,"PatientID"]))
TestPredictions$Diagnosis=unlist(lapply(rownames(TestPredictions), function(x) pdata.training[x,"Diagnosis"]))
TestPredictions$TissueType=unlist(lapply(rownames(TestPredictions), function(x) pdata.training[x,"TissueType"]))
TestPredictions$PAMprediction=revalue(as.factor(TestPredictions$PAMprediction), c("1"="NKP","2"="chRCC/KICH","3"="oncocytoma"))
#TestPredictions$Diagnosis=revalue(as.factor(TestPredictions$Diagnosis), c("AML"="Angiomylolipoma", "KICH"="chRCC/KICH", "KIRC"="ccRCC/KIRC", "KIRP"="pRCC/KIRP","Normal"="NKP","oncocytoma"="oncocytoma"))
TestPredictions$TissueType=revalue(as.factor(TestPredictions$TissueType), c("AML"="Angiomylolipoma", "KICH"="chRCC/KICH", "KIRC"="ccRCC/KIRC", "KIRP"="pRCC/KIRP","Normal"="NKP","oncocytoma"="oncocytoma"))
TestPredictions$Diagnosis=revalue(as.factor(TestPredictions$Diagnosis), c("NKP"="Normal/benign lesion"))

#####
TestPredictions$Misclassified=rep("No", nrow(TestPredictions))
TestPredictions[TestPredictions$TissueType=="oncocytoma" & TestPredictions$PAMprediction!="oncocytoma","Misclassified"]="Yes"
TestPredictions[TestPredictions$TissueType=="chRCC/KICH" & TestPredictions$PAMprediction!="chRCC/KICH","Misclassified"]="Yes"
TestPredictions[TestPredictions$TissueType=="chRCC/KICH" & TestPredictions$PAMprediction!="chRCC/KICH",]
#Where is the sample number that distinguished the two biopsies? 

PosteriorProbabilites=PAMresults.tumor.trained.with.norm$PosteriorProbabilites
TestPredictions=cbind(TestPredictions, PosteriorProbabilites)
#Add rownames as sample ID, since this is the only sample-specific ID
TestPredictions=cbind(TestPredictions, rownames(TestPredictions))
#Where is the sample number that distinguished the two biopsies? 
TestPredictions=TestPredictions[,c(9,2,3,4,1,5,6,7,8)]

colnames(TestPredictions)=c("Sample ID","Patient ID","Diagnosis","Tissue type","Diagnostic model prediction", "Misclassified by diagnostic model?","Posterior probability NKP", "Posterior probability chRCC","Posterior probability oncocytoma")

write.table(TestPredictions, file="resultsdir/PAM_Chopra_tumors.txt", sep="\t", row.names = FALSE)

#########################
#Make MDS plot labeling Just tumors
##########################

PAMresults.tumor.trained.with.norm=readRDS("datadir/Archive_of_OSF_Storage/data/PAM_Chopra_tumors.rds")

#PAM analsyis trained with normal tissue applied to needle biopsies 
rlab=as.matrix(rownames(PAMresults.tumor.trained.with.norm$TestPredictions))
rownames(rlab)=rlab[,1]
colnames(rlab)=""
rlab[rownames(rlab) %in% Barcodes_chRCC_training_Chopra,1]=cbp2[4]
rlab[rownames(rlab) %in% Barcodes_onc_training_Chopra,1]=cbp2[5]

#Adding patient number so that we can determine if at least one sample was correctly classified 
#Also def need to add normal adjacent tissue
PosteriorProbabilites=as.data.frame(PAMresults.tumor.trained.with.norm$PosteriorProbabilites)
PosteriorProbabilites$patient.no=rep(NA, nrow(PosteriorProbabilites))
PosteriorProbabilites$patient.no=unlist(lapply(rownames(PosteriorProbabilites), function(x) labs[x, "Sample"]))
PosteriorProbabilites$Diagnosis=unlist(lapply(rownames(PosteriorProbabilites), function(x) labs[x,"Diagnosis"]))
PosteriorProbabilites$Diagnosis=as.factor(PosteriorProbabilites$Diagnosis)
PosteriorProbabilites$TissueType=unlist(lapply(rownames(PosteriorProbabilites), function(x) labs[x,"TissueType"]))
PosteriorProbabilites$TissueType=as.factor(PosteriorProbabilites$TissueType)
PosteriorProbabilites$TestPredictions=PAMresults.tumor.trained.with.norm$TestPredictions[,"PAMprediction"]
PosteriorProbabilites$Diagnosis=factor(PosteriorProbabilites$Diagnosis)
PosteriorProbabilites=PosteriorProbabilites[order(PosteriorProbabilites$Diagnosis, PosteriorProbabilites$patient.no),]
rlab=rlab[match(rownames(PosteriorProbabilites), rownames(rlab))]

#get labels for misclassified patients 
misclassifications.oncocytoma=PosteriorProbabilites[PosteriorProbabilites$TissueType=="oncocytoma" & PosteriorProbabilites$TestPredictions!=3,"patient.no"]
names(misclassifications.oncocytoma)=rownames(PosteriorProbabilites[PosteriorProbabilites$TissueType=="oncocytoma" & PosteriorProbabilites$TestPredictions!=3,])

misclassifications.chRCC=PosteriorProbabilites[PosteriorProbabilites$TissueType=="KICH" & PosteriorProbabilites$TestPredictions!=2,"patient.no"]
names(misclassifications.chRCC)=paste(PosteriorProbabilites[PosteriorProbabilites$TissueType=="KICH" & PosteriorProbabilites$TestPredictions!=2,"patient.no"],"-01", sep="")

file="resultsdir/MDS_Pan_kidney_with_Chopra_tumors_label_misclassifications"
png(file=paste(file,'.png',sep=''), units="in", width=10, height=10, res=600)
plot(x[Barcodes_METcancer_KIRP], y[Barcodes_METcancer_KIRP], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[1],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1)
par(new=TRUE)
plot(x[Barcodes_METcancer_KIRC], y[Barcodes_METcancer_KIRC], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[2],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[Barcodes_METcancer_KICH], y[Barcodes_METcancer_KICH], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[3],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[Barcodes_METnormal], y[Barcodes_METnormal], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[4],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[1:4],y[1:4], xlab="Coordinate 1", ylab="Coordinate 2", col=c(cbp2[1:4]), pch=15, cex=1.5, , xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), main="")
par(new=TRUE)
plot(x[BarcodesClearCell], y[BarcodesClearCell], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[2],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19)
par(new=TRUE)
plot(x[BarcodesChromophobe], y[BarcodesChromophobe], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[3],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[BarcodesNormal], y[BarcodesNormal], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[4],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[BarcodesOncocytoma], y[BarcodesOncocytoma], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[5],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[Barcodes.Chopra.tumor.onc], y[Barcodes.Chopra.tumor.onc], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[5],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=23, main="")
text(x[names(misclassifications.oncocytoma)], y[names(misclassifications.oncocytoma)], labels=misclassifications.oncocytoma, cex=0.9, font=2, col=cbp2[5], adj=c(1,1))
text(x[names(misclassifications.chRCC)], y[names(misclassifications.chRCC)], labels=misclassifications.chRCC, cex=0.9, font=2, col=cbp2[3], adj=c(1,1))
dev.off()

#There are no duplicated tumor samples
file="resultsdir/MDS_point_type_legend_with_Chopra"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8, res=600)
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("center", c("TCGA study","TCGA study centroid","Stanford study", "Chopra study FFPE tumor", "Chopra study needle biopsy"), col=c(adjustcolor(1,alpha.f=0.5),1,1), pch=c(1,15,19, 23, 8), title="Study/sample type", cex=2, bty="n")
dev.off()

#Make legend that excludes the samples that we couldn't classify
file="resultsdir/heatmap_posterior_possibility_PAM_PanKidney_colorblind_legend_NoUnclassifiedRCCTypes"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8, res=600)
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("center", c("pRCC/KIRP", "ccRCC/KIRC","chRCC/KICH","NKP","oncocytoma"), col=c(cbp2[1:5]), pch=15, title="Sample type", cex=2, bty="n")
dev.off()


############################
#MDS plot based on the 79 CpG sites that distinguish onc form chRCC, first showing biopsies
###########################

SAM_genes=read.table('resultsdir/SAM_analysis_OncVChRCC.txt', sep="\t", header=T)
SAM_genes=SAM_genes[SAM_genes[,4]=="0" & abs(SAM_genes[,3])>0.2,]
#79 genes
NonZeroGenes=as.character(SAM_genes$Gene.Name)

#try it just with the probes that are also MethylMix genes
METall_NonZeroGenes=METall[rownames(METall) %in% NonZeroGenes,]
#65 available in all samples

ClassCentroids_Hand=cbind(rowMeans(METall_NonZeroGenes[,Barcodes_METcancer_KIRP]),
                          rowMeans(METall_NonZeroGenes[,Barcodes_METcancer_KIRC]),
                          rowMeans(METall_NonZeroGenes[,Barcodes_METcancer_KICH]),
                          rowMeans(METall_NonZeroGenes[,Barcodes_METnormal]))

OverlapProbes=intersect(rownames(METall), rownames(ClassCentroids_Hand))
METall_CD=METall[OverlapProbes,]
ClassCentroids_METall=ClassCentroids_Hand[OverlapProbes,]

OverlapProbes=intersect(rownames(beta.testing.impute),intersect(rownames(beta.training.impute.onc),intersect(rownames(beta), intersect(rownames(METall), rownames(ClassCentroids_Hand)))))
m1=ClassCentroids_Hand[OverlapProbes,]
m2=METall[OverlapProbes,]
m3=beta[OverlapProbes,]
m4=beta.training.impute.onc[OverlapProbes,]
m5=beta.testing.impute[OverlapProbes,]

d1=t(cbind(m1, m2, m3, m4, m5))
d2=dist(d1)

fit <- cmdscale(d2,eig=TRUE, k=2) # k is the number of dim

x <- fit$points[,1]
y <- fit$points[,2]

file="resultsdir/MDS_Pan_kidney_with_Chopra_biopsy_label_misclassifications_chRCC_v_onc_probes"
png(file=paste(file,'.png',sep=''), units="in", width=10, height=10, res=600)
plot(x[Barcodes_METcancer_KIRP], y[Barcodes_METcancer_KIRP], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[1],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1)
par(new=TRUE)
plot(x[Barcodes_METcancer_KIRC], y[Barcodes_METcancer_KIRC], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[2],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[Barcodes_METcancer_KICH], y[Barcodes_METcancer_KICH], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[3],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[Barcodes_METnormal], y[Barcodes_METnormal], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[4],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[1:4],y[1:4], xlab="Coordinate 1", ylab="Coordinate 2", col=c(cbp2[1:4]), pch=15, cex=1.5, , xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), main="")
par(new=TRUE)
plot(x[BarcodesClearCell], y[BarcodesClearCell], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[2],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19)
par(new=TRUE)
plot(x[BarcodesChromophobe], y[BarcodesChromophobe], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[3],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[BarcodesNormal], y[BarcodesNormal], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[4],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[BarcodesOncocytoma], y[BarcodesOncocytoma], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[5],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[Barcodes.Chopra.biopsy.Normal], y[Barcodes.Chopra.biopsy.Normal], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[4],alpha.f=1), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=8, main="")
par(new=TRUE)
plot(x[Barcodes.Chopra.biopsy.onc], y[Barcodes.Chopra.biopsy.onc], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[5],alpha.f=1), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=8, main="")
par(new=TRUE)
plot(x[Barcodes.Chopra.biopsy.KICH], y[Barcodes.Chopra.biopsy.KICH], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[3],alpha.f=1), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=8, main="")
par(new=TRUE)
plot(x[Barcodes.Chopra.biopsy.KIRP], y[Barcodes.Chopra.biopsy.KIRP], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[1],alpha.f=1), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=8, main="")
par(new=TRUE)
plot(x[Barcodes.Chopra.biopsy.KIRC], y[Barcodes.Chopra.biopsy.KIRC], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[2],alpha.f=1), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=8, main="")
par(new=TRUE)
text(x[intersect(BarcodesMiclassificationsChopraBiopsyonc, Barcodes.Chopra.biopsy.onc)],y[intersect(BarcodesMiclassificationsChopraBiopsyonc, Barcodes.Chopra.biopsy.onc)], labels=PosteriorProbabilites.biopsy[intersect(BarcodesMiclassificationsChopraBiopsyonc, Barcodes.Chopra.biopsy.onc),"patient.no"], cex=0.9, font=2, col=cbp2[5], adj=c(1,1))
text(x[intersect(BarcodesMiclassificationsChopraBiopsyNorm, Barcodes.Chopra.biopsy.Normal)],y[intersect(BarcodesMiclassificationsChopraBiopsyNorm, Barcodes.Chopra.biopsy.Normal)], labels=PosteriorProbabilites.biopsy[intersect(BarcodesMiclassificationsChopraBiopsyNorm, Barcodes.Chopra.biopsy.Normal),"patient.no"], cex=0.9, font=2, col=cbp2[4], adj=c(1,1))
text(x[intersect(BarcodesMiclassificationsChopraBiopsychRCC, Barcodes.Chopra.biopsy.KICH)],y[intersect(BarcodesMiclassificationsChopraBiopsychRCC, Barcodes.Chopra.biopsy.KICH)], labels=PosteriorProbabilites.biopsy[intersect(BarcodesMiclassificationsChopraBiopsychRCC, Barcodes.Chopra.biopsy.KICH),"patient.no"], cex=0.9, font=2, col=cbp2[3], adj=c(1,1))
for(i in 1:nrow(matall)){
segments(matall[i,"X0"], matall[i,"Y0"], matall[i,"X1"], matall[i,"Y1"], col=as.character(matall$col[i]), lty=2)
}
dev.off()


#############################################
#MDS plot for just 79 CpG sites with just 79 CpGs that overlap with all genes 
########################################

file="resultsdir/MDS_Pan_kidney_with_Chopra_tumors_label_misclassifications_chRCC_v_onc_probes"
png(file=paste(file,'.png',sep=''), units="in", width=10, height=10, res=600)
plot(x[Barcodes_METcancer_KIRP], y[Barcodes_METcancer_KIRP], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[1],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1)
par(new=TRUE)
plot(x[Barcodes_METcancer_KIRC], y[Barcodes_METcancer_KIRC], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[2],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[Barcodes_METcancer_KICH], y[Barcodes_METcancer_KICH], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[3],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[Barcodes_METnormal], y[Barcodes_METnormal], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[4],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=1, main="")
par(new=TRUE)
plot(x[1:4],y[1:4], xlab="Coordinate 1", ylab="Coordinate 2", col=c(cbp2[1:4]), pch=15, cex=1.5, , xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), main="")
par(new=TRUE)
plot(x[BarcodesClearCell], y[BarcodesClearCell], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[2],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19)
par(new=TRUE)
plot(x[BarcodesChromophobe], y[BarcodesChromophobe], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[3],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[BarcodesNormal], y[BarcodesNormal], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[4],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[BarcodesOncocytoma], y[BarcodesOncocytoma], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[5],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=19, main="")
par(new=TRUE)
plot(x[Barcodes.Chopra.tumor.onc], y[Barcodes.Chopra.tumor.onc], xlab="Coordinate 1", ylab="Coordinate 2", col=adjustcolor(cbp2[5],alpha.f=0.5), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), pch=23, main="")
text(x[names(misclassifications.oncocytoma)], y[names(misclassifications.oncocytoma)], labels=misclassifications.oncocytoma, cex=0.9, font=2, col=cbp2[5], adj=c(1,1))
text(x[names(misclassifications.chRCC)], y[names(misclassifications.chRCC)], labels=misclassifications.chRCC, cex=0.9, font=2, col=cbp2[3], adj=c(1,1))
dev.off()

#
write.table(OverlapProbes, file="resultsdir/Chrcc_Vs_onc_probes_overlapping_with_all_Chopra_MDS.txt", sep="\t", row.names=FALSE)

####################################################################
#Applying PAM analysis to validation cohort Chopra and TCGA chRCC, with incrementally smaller numbers of CpG sites
######################################################################

genranks=read.table("resultsdir/PAM_oncVchRCC_predictive_CpGs_iterative.txt", sep="\t"", header=T)
AllCpGs=rownames(genranks[!duplicated(genranks$region),])
#73 CpGs 

Barcodes=c(BarcodesNormal, BarcodesOncocytoma, BarcodesChromophobe)

SampleClass=rep(NA, length(Barcodes))
names(SampleClass)=Barcodes
SampleClass[names(SampleClass) %in% BarcodesNormal]=1
SampleClass[names(SampleClass) %in% BarcodesChromophobe]=2
SampleClass[names(SampleClass) %in% BarcodesOncocytoma]=3
SampleClass=t(as.matrix(SampleClass))
rownames(SampleClass)="Sample_type"

Ns=seq(10,70,5)

NrFolds=10

PAMs=list()

Aucs=array(NA, c(length(Ns),3))

for(i in 1:length(Ns)){
  N=Ns[i]
  genes=AllCpGs[1:N]
  
  beta_genes_Barcodes=beta[genes,colnames(SampleClass)]
  #scale training set 
  beta_genes_Barcodes=scale(beta_genes_Barcodes)

  beta.training.impute.onc.ChRCC.genes=beta.training.impute.onc.ChRCC[genes,]
  #scale test set
  beta.training.impute.onc.ChRCC.genes.s=scale(beta.training.impute.onc.ChRCC.genes)

  PAMresults=PAM_Analysis_KB(beta_genes_Barcodes, SampleClass, beta.training.impute.onc.ChRCC.genes.s, NrLambdaFolds=10)
  PAMs[[i]]=PAMresults
  ROCcurve=multiclass.roc(as.factor(SampleType.Chopra.onc.KICH.training-1), PAMresults$PosteriorProbabilites[,2])
  Aucs[i,1]=round(as.numeric(as.character(gsub(" .*","", auc(ROCcurve)))),2)
}
names(PAMs)=Ns
save(PAMs, file="resultsdir/PAM_oncVchRCC_Chopra_tumors.RData")

NrMisclassificationsPerClassList=list()
for(i in 1:length(PAMs)){
     tab=table(PAMs[[i]]$TestPredictions, SampleType.Chopra.onc.KICH.training)
     tab=as.data.frame.matrix(tab) 
     colnames(tab)=c("chRCC_TCGA","oncocytoma_Chopra")
     rownames(tab)=c("Classified_as_normal","Classified_as_chRCC","Classified_as_oncocytoma")
     tab["accuracy",]=c(round(tab["Classified_as_chRCC","chRCC_TCGA"]/sum(tab[,"chRCC_TCGA"]),2),round(tab["Classified_as_oncocytoma","oncocytoma_Chopra"]/sum(tab[,"oncocytoma_Chopra"]),2))
     tab["NrMisclassifications",]=c(sum(tab[c("Classified_as_normal","Classified_as_oncocytoma"),"chRCC_TCGA"]), sum(tab[c("Classified_as_normal","Classified_as_chRCC"),"oncocytoma_Chopra"]))
     tab$classifications=rownames(tab)
     tab$Ngenes=rep(names(PAMs)[i], nrow(tab))
     NrMisclassificationsPerClassList[[i]]=tab
}
names(NrMisclassificationsPerClassList)=names(PAMs)

library(abind)
NrMisclassificationsPerClassList_iterative=abind(NrMisclassificationsPerClassList, along=1)
NrMisclassificationsPerClassList_iterative=as.data.frame(NrMisclassificationsPerClassList_iterative)
NrMisclassificationsPerClassList_iterative$Ngenes=as.numeric(as.character(NrMisclassificationsPerClassList_iterative$Ngenes))

file="resultsdir/PAM_oncVchRCC_interative_misclassificationsperclass_Chopra_tumors"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
plot(NrMisclassificationsPerClassList_iterative[NrMisclassificationsPerClassList_iterative$classifications=="accuracy","Ngenes"], as.numeric(as.character(NrMisclassificationsPerClassList_iterative[NrMisclassificationsPerClassList_iterative$classifications=="accuracy","chRCC_TCGA"])), col=cbp2[3], pch=19, cex.axis=1.3, ylab="% accuracy", xlab="N genes", cex.lab=1.3)
points(NrMisclassificationsPerClassList_iterative[NrMisclassificationsPerClassList_iterative$classifications=="accuracy","Ngenes"], as.numeric(as.character(NrMisclassificationsPerClassList_iterative[NrMisclassificationsPerClassList_iterative$classifications=="accuracy","oncocytoma_Chopra"])), col=cbp2[5], pch=19)
dev.off()

file="resultsdir/PAM_oncVchRCC_interative_NrMisclassifications_Chopra_tumors"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
plot(NrMisclassificationsPerClassList_iterative[NrMisclassificationsPerClassList_iterative$classifications=="NrMisclassifications","Ngenes"], as.numeric(as.character(NrMisclassificationsPerClassList_iterative[NrMisclassificationsPerClassList_iterative$classifications=="NrMisclassifications","chRCC_TCGA"])), col=cbp2[3], pch=19, cex.axis=1.3, ylab="n misclassifications", xlab="N genes", cex.lab=1.3, ylim=c(0,40))
points(NrMisclassificationsPerClassList_iterative[NrMisclassificationsPerClassList_iterative$classifications=="NrMisclassifications","Ngenes"], as.numeric(as.character(NrMisclassificationsPerClassList_iterative[NrMisclassificationsPerClassList_iterative$classifications=="NrMisclassifications","oncocytoma_Chopra"])), col=cbp2[5], pch=19, ylim=c(0,40))
dev.off()

tab=NrMisclassificationsPerClassList_iterative[NrMisclassificationsPerClassList_iterative$classifications %in% c("accuracy","NrMisclassifications"),]

accuracy_chRCC=tab[tab$classifications=="accuracy",c("chRCC_TCGA","Ngenes")]
accuracy_chRCC$sample_type=rep("chRCC",nrow(accuracy_chRCC))
colnames(accuracy_chRCC)=c("accuracy","Ngenes","sample_type")

accuracy_onc=tab[tab$classifications=="accuracy",c("oncocytoma_Chopra","Ngenes")]
accuracy_onc$sample_type=rep("oncocytoma",nrow(accuracy_onc))
colnames(accuracy_onc)=c("accuracy","Ngenes","sample_type")

accuracy=rbind(accuracy_chRCC, accuracy_onc)

NrMisclassifications_chRCC=tab[tab$classifications=="NrMisclassifications",c("chRCC_TCGA","Ngenes")]
NrMisclassifications_chRCC$sample_type=rep("chRCC",nrow(NrMisclassifications_chRCC))
colnames(NrMisclassifications_chRCC)=c("NrMisclassifications","Ngenes","sample_type")

NrMisclassifications_onc=tab[tab$classifications=="NrMisclassifications",c("oncocytoma_Chopra","Ngenes")]
NrMisclassifications_onc$sample_type=rep("oncocytoma",nrow(NrMisclassifications_onc))
colnames(NrMisclassifications_onc)=c("NrMisclassifications","Ngenes","sample_type")

NrMisclassifications=rbind(NrMisclassifications_chRCC, NrMisclassifications_onc)
both=cbind(accuracy, NrMisclassifications)
both=both[,c(1:4)]
both=both[order(both$Ngenes),]
both$overall_n=rep(NA, nrow(both))
both[both$sample_type=="chRCC","overall_n"]=length(Barcodes_chRCC_training_Chopra)
both[both$sample_type=="oncocytoma","overall_n"]=length(Barcodes_onc_training_Chopra)
both=both[,c("Ngenes","sample_type","NrMisclassifications","overall_n","accuracy")]
both$NrMisclassifications=as.numeric(as.character(both$NrMisclassifications))
both$NrMisclassificationsPerc=round(both$NrMisclassifications/both$overall_n,2)
colnames(both)[3]="NrMistakes"
both=both[,c("NrMistakes","overall_n","NrMisclassificationsPerc","accuracy","Ngenes","sample_type")]

write.table(both, file="resultsdir/PAM_oncVchRCC_interative_misclassificationsperclass_Chopra_tumors.txt", sep="\t", row.names = FALSE)

##############################################
#Model trained in with normal tissue applied to biopsies. Making figure 
##################################################

Barcodes.Chopra.oncocytoma.KICH=c(Barcodes.pdata.biopsy.norm, Barcodes.pdata.KICH.chopra, Barcodes.pdata.onc.chopra)
SampleType.Chopra.onc.KICH=c(rep(1, length(Barcodes.pdata.biopsy.norm)), rep(2, length(Barcodes.pdata.KICH.chopra)), rep(3, length(Barcodes.pdata.onc.chopra)))
names(SampleType.Chopra.onc.KICH)=Barcodes.Chopra.oncocytoma.KICH

beta.testing.impute=readRDS("datadir/Archive_of_OSF_Storage/data/beta.testing.rds")
beta.Chopra.test=beta.testing.impute[genes,Barcodes.Chopra.oncocytoma.KICH]

beta_genes_Barcodes.s=scale(beta_genes_Barcodes)
beta.Chopra.test.s=scale(beta.Chopra.test)

PAMresults.biopsies.trained.with.norm=PAM_Analysis_KB(beta_genes_Barcodes, SampleClass, beta.Chopra.test.s, NrLambdaFolds=10)

save(PAMresults.biopsies.trained.with.norm, file="resultsdir/PAM_Chopra_needle_core_biopsies_with_normal.RData")
dat=load("resultsdir/PAM_Chopra_needle_core_biopsies_with_normal.RData")

misclassifications=as.data.frame(cbind(PAMresults.biopsies.trained.with.norm$TestPredictions, SampleType.Chopra.onc.KICH))
BarcodesMiclassificationsChopraBiopsyNorm=rownames(misclassifications[misclassifications$SampleType.Chopra.onc.KICH==1 & misclassifications$PAMprediction!=1,])
BarcodesMiclassificationsChopraBiopsychRCC=rownames(misclassifications[misclassifications$SampleType.Chopra.onc.KICH==2 & misclassifications$PAMprediction!=2,])
BarcodesMiclassificationsChopraBiopsyonc=rownames(misclassifications[misclassifications$SampleType.Chopra.onc.KICH==3 & misclassifications$PAMprediction!=3,])

tab=table(PAMresults.biopsies.trained.with.norm$TestPredictions, SampleType.Chopra.onc.KICH)
tab=as.matrix(tab)
colnames(tab)=c("biopsy NKP","biopsy chRCC","biopsy oncocytoma")
rownames(tab)=c("Classified_as_normal","Classified_as_chRCC","Classified_as_oncocytoma")
write.table(tab, file="resultsdir/Calssifications_Chopra_biopsies_model_trained_with_normal", sep="\t")

rlab=as.matrix(rownames(PAMresults.biopsies.trained.with.norm$TestPredictions))
rownames(rlab)=rlab[,1]
colnames(rlab)=""
rlab[rownames(rlab) %in% Barcodes.pdata.KICH.chopra,1]=cbp2[4]
rlab[rownames(rlab) %in% Barcodes.pdata.onc.chopra,1]=cbp2[5]
rlab[rownames(rlab) %in% Barcodes.pdata.biopsy.norm,1]=cbp2[3]

#Adding patient number so that we can determine if at least one sample was correctly classified 
#Also def need to add normal adjacent tissue
PosteriorProbabilites.biopsy=as.data.frame(PAMresults.biopsies.trained.with.norm$PosteriorProbabilites)
PosteriorProbabilites.biopsy$patient.no=rep(NA, nrow(PosteriorProbabilites.biopsy))
PosteriorProbabilites.biopsy$patient.no=unlist(lapply(rownames(PosteriorProbabilites.biopsy), function(x) pdata.testing[x,"Pt..No."]))
PosteriorProbabilites.biopsy$TissueType=unlist(lapply(rownames(PosteriorProbabilites.biopsy), function(x) pdata.testing[x,"TissueType"]))
PosteriorProbabilites.biopsy$TissueType=as.factor(PosteriorProbabilites.biopsy$TissueType)
levels(PosteriorProbabilites.biopsy$TissueType)=c(2,1,3)
PosteriorProbabilites.biopsy$TissueType=as.character(PosteriorProbabilites.biopsy$TissueType)
PosteriorProbabilites.biopsy$TestPredictions=PAMresults.biopsies.trained.with.norm$TestPredictions[,"PAMprediction"]

PosteriorProbabilites.biopsy=PosteriorProbabilites.biopsy[order(PosteriorProbabilites.biopsy$TissueType, PosteriorProbabilites.biopsy$patient.no),]
rlab=rlab[match(rownames(PosteriorProbabilites.biopsy), rownames(rlab))]

#Need to also 
#file="resultsdir/PAM_probabilities_analysis_onc_chRCC_Chopra_TCGA_biopsies_models_trained_with_normal"
#png(file=paste(file,'.png',sep=''), units="in", width=3, height=9, res=600)
#heatmap.3(PosteriorProbabilites.biopsy[,1:3],  Rowv=FALSE, Colv=FALSE, dendrogram = "none",  col=brewer.pal(9,"Reds"), RowSideColors=t(rlab), KeyValueName="Classification probability", RowSideColorsSize=4, labRow=PosteriorProbabilites.biopsy$patient.no, labCol="")
#dev.off()

file="resultsdir/PAM_probabilities_analysis_onc_chRCC_Chopra_TCGA_biopsies_models_trained_with_normal"
pdf(file=paste(file,'.pdf',sep=''), height = 9.5,  width = 3, family = "Helvetica")
heatmap.3(PosteriorProbabilites.biopsy[,1:3],  Rowv=FALSE, Colv=FALSE, dendrogram = "none",  col=brewer.pal(9,"Reds"), RowSideColors=t(rlab), KeyValueName="Classification probability", RowSideColorsSize=4, labRow=PosteriorProbabilites.biopsy$patient.no, labCol="")
dev.off()

#########
#Make tumor classification posterior probability figures
###############

#Get tumor labels from Chopra paper supplement 
labs=read.table("datadir/Chopra_oncocytoma_tumor_labels.txt", sep="\t", header=T)
rownames(labs)=labs[,1]
labs$Sample=substr(labs$Sample,1,12)
rownames(labs)=labs$Sample
labs$patient.no=rownames(labs)

PAMresults.tumor.trained.with.norm=readRDS("datadir/Archive_of_OSF_Storage/data/PAM_Chopra_tumors.rds")

tab=table(PAMresults.tumor.trained.with.norm$TestPredictions, SampleType.Chopra.onc.KICH.training)
tab=as.matrix(tab)
colnames(tab)=c("chRCC_TCGA","oncocytoma_Chopra")
rownames(tab)=c("Classified_as_normal","Classified_as_chRCC","Classified_as_oncocytoma")
#write.table(tab, file="resultsdir/Calssifications_Chopra_tumor_model_trained_with_normal", sep="\t")

#PAM analsyis trained with normal tissue applied to needle biopsies 
rlab=as.matrix(rownames(PAMresults.tumor.trained.with.norm$TestPredictions))
rownames(rlab)=rlab[,1]
colnames(rlab)=""
rlab[rownames(rlab) %in% Barcodes_chRCC_training_Chopra,1]=cbp2[4]
rlab[rownames(rlab) %in% Barcodes_onc_training_Chopra,1]=cbp2[5]

#Adding patient number so that we can determine if at least one sample was correctly classified 
#Also def need to add normal adjacent tissue
PosteriorProbabilites=as.data.frame(PAMresults.tumor.trained.with.norm$PosteriorProbabilites)
PosteriorProbabilites$patient.no=rep(NA, nrow(PosteriorProbabilites))
PosteriorProbabilites$patient.no=unlist(lapply(rownames(PosteriorProbabilites), function(x) labs[x, "Sample"]))
PosteriorProbabilites$Diagnosis=unlist(lapply(rownames(PosteriorProbabilites), function(x) labs[x,"Diagnosis"]))
PosteriorProbabilites$Diagnosis=as.factor(PosteriorProbabilites$Diagnosis)
PosteriorProbabilites$TissueType=unlist(lapply(rownames(PosteriorProbabilites), function(x) labs[x,"TissueType"]))
PosteriorProbabilites$TissueType=as.factor(PosteriorProbabilites$TissueType)
PosteriorProbabilites$TestPredictions=PAMresults.tumor.trained.with.norm$TestPredictions[,"PAMprediction"]
PosteriorProbabilites$Diagnosis=factor(PosteriorProbabilites$Diagnosis)
PosteriorProbabilites=PosteriorProbabilites[order(PosteriorProbabilites$Diagnosis, PosteriorProbabilites$patient.no),]
rlab=rlab[match(rownames(PosteriorProbabilites), rownames(rlab))]

#Need to also 
file="resultsdir/PAM_probabilities_analysis_onc_chRCC_Chopra_TCGA_tumors_models_trained_with_normal"
pdf(file=paste(file,'.pdf',sep=''), height = 9.5,  width = 3, family = "Helvetica")
par(oma = c(2,0,0,2))
heatmap.3(PosteriorProbabilites[,1:3],  Rowv=FALSE, Colv=FALSE, dendrogram = "none",  col=brewer.pal(9,"Reds"), RowSideColors=t(rlab), KeyValueName="Classification probability", RowSideColorsSize=4, labRow=PosteriorProbabilites$patient.no, labCol="")
dev.off()

#####################################
#Make volcano plots indicating differential methylation between chRCC Vs. NKP, and oncocytoma Vs. NKP
Figure 2a
######################################

#Do Wilcox test for oncocytoma vs normal and chromophobe versus normal

#Normal v chromophobe
var=c(rep(1, length(BarcodesNormal)),rep(2, length(BarcodesChromophobe)))
names(var)=c(BarcodesNormal, BarcodesChromophobe)
var=as.factor(var)
beta2=cbind(beta[,colnames(beta) %in% BarcodesNormal],beta[,colnames(beta) %in% BarcodesChromophobe])

wilc_chromophobe=Rowwise_stat(dataframe=beta2, var=var)

#Make volcano plot for chromophobe Vs. normal
file="resultsdir/volcano_chromophobe_normal"
png(file=paste(file,'.png',sep=''), units="in", width=8.5, height=8.5, res=600)
with(wilc_chromophobe, plot(Diff, -log10(P.value), xlim=c(-0.7, 0.7), pch=20, xlab="DNA methylation difference (chromophobe-normal)", cex.lab=1.5, cex.axis=1.5))
with(subset(wilc_chromophobe, Q.value<0.05), points(Diff, -log10(P.value), pch=20, col="red"))
with(subset(wilc_chromophobe, Q.value<0.05 & abs(Diff)>0.2), points(Diff, -log10(P.value), pch=20, col="green"))
abline(v=-0.5, col="blue",lty=2)
abline(v=0.5, col="blue",lty=2)
dev.off()

#############################################
# Volcano Normal v oncocytoma
################################################

var=c(rep(1, length(BarcodesNormal)),rep(2, length(BarcodesOncocytoma)))
names(var)=c(BarcodesNormal, BarcodesOncocytoma)
var=as.factor(var)
beta2=cbind(beta[,colnames(beta) %in% BarcodesNormal],beta[,colnames(beta) %in% BarcodesOncocytoma])

#just find differential methylation for all CpGs 
wilc_oncocytoma=Rowwise_stat(dataframe=beta2, var=var)

#Make volcano plot for oncocytoma Vs. normal
file="resultsdir/volcano_onc_normal"
png(file=paste(file,'.png',sep=''), units="in", width=8.5, height=8.5, res=600)
with(wilc_oncocytoma, plot(Diff, -log10(P.value), xlim=c(-0.7, 0.7), pch=20, xlab="DNA methylation difference (oncocytoma-normal)", cex.lab=1.5, cex.axis=1.5))
with(subset(wilc_oncocytoma, Q.value<0.05), points(Diff, -log10(P.value), pch=20, col="red"))
with(subset(wilc_oncocytoma, Q.value<0.05 & abs(Diff)>0.2), points(Diff, -log10(P.value), pch=20, col="green"))
abline(v=-0.5, col="blue",lty=2)
abline(v=0.5, col="blue",lty=2)
dev.off()

#######################################################
#Plot delta beta for shared abnormally methylated probes in chRCC and oncocytoma
#Figure 2B
##########################################################

AllResultsChRCC=read.table("resultsdir/SAMPAM_analyses/SAM_analysis_ChRCCVsNormal_formatted.txt", sep="\t", header=T, colClasses=c("character"))
rownames(AllResultsChRCC)=AllResultsChRCC$probe

AllResultsOnc=read.table("resultsdir/SAMPAM_analyses/SAM_analysis_OncVsNormal_formatted.txt", sep="\t", header=T, colClasses=c("character"))
rownames(AllResultsOnc)=AllResultsOnc$probe

#of 517 hyper DMRs in ChRCC, 297 are also hyper in oncocytoma
#962 hyper in oncocytoma altogether 

#3377/5058 (67%) of ChRCC hypo
#3377/8550 (39% of oncocytoma hypo probes intersect)
HyperIntersect=intersect(AllResultsChRCC[AllResultsChRCC$direction=="Up","probe"],AllResultsOnc[AllResultsOnc$direction=="Up","probe"])

HypoIntersect=intersect(AllResultsChRCC[AllResultsChRCC$direction=="Down","probe"],AllResultsOnc[AllResultsOnc$direction=="Down","probe"])

#Adding delta betas for other tumor type versus NKP, if the DMR is shared between tumors (in the same direction)
AllResultsChRCC$Delta.beta.Oncocytoma.versus.NKP.DRMs=unlist(lapply(rownames(AllResultsChRCC), function(x) ifelse(x %in% c(HyperIntersect, HypoIntersect), AllResultsOnc[x,"delta.beta"], NA)))

#Adding delta betas for other tumor type versus NKP, if the DMR is shared between tumors (in the same direction)
AllResultsOnc$Delta.beta.ChRCC.versus.NKP.DRMs=unlist(lapply(rownames(AllResultsOnc), function(x) ifelse(x %in% c(HyperIntersect, HypoIntersect), AllResultsChRCC[x,"delta.beta"], NA)))

delta.betas.compare=cbind(AllResultsOnc$delta.beta, AllResultsOnc$Delta.beta.ChRCC.versus.NKP.DRMs)
delta.betas.compare=delta.betas.compare[complete.cases(delta.betas.compare),]

file="resultsdir/figures/RO_ChRCC_shared_DMRs_sqew"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
plot(as.numeric(as.character(delta.betas.compare[,1])), as.numeric(as.character(delta.betas.compare[,2])), ylab="delta beta (ChRCC Vs NKP)", xlab="delta beta (RO Vs NKP)", cex.axis=1.3, cex.lab=1.3)
dev.off()

#Add methylation levels in TCGA studies
write.table(AllResultsOnc, file="resultsdir/SAMPAM_analyses/SAM_analysis_OncVsNormal_formatted.txt", sep="\t")
write.table(AllResultsChRCC, file="resultsdir/SAMPAM_analyses/SAM_analysis_ChRCCVsNormal_formatted.txt", sep="\t")


#######################################
################
#PAM and SAM analysis in of differential methylation Stanford study 
################
##########################################

#Got cbp2 as a colorblind friendly scheme, but changing a couple of them becase they are light or non-distinct colors
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#FD61D1", "#A3A500", "#D55E00", "#CC79A7")

########
#Set up packages
###########

library(plyr)
library(RColorBrewer)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library("gplots")
library("devtools")
library(ggplot2)
library(gridExtra)
library(abind)
#

BarcodesClearCell=rownames(Sample_file[grep("Clear cell RCC",Sample_file$Sample_Type),])
#There are no KIRPs in newdata
BarcodesChromophobe=rownames(Sample_file[grep("RCC-Chromophobe",Sample_file$Sample_Type),])
BarcodesNormal=rownames(Sample_file[grep("Normal kidney parenchyma",Sample_file$Sample_Type),])
BarcodesOncocytoma=rownames(Sample_file[grep("Oncocytoma",Sample_file$Sample_Type),])
BarcodesHybridoncocyticrenal=rownames(Sample_file[grep("Hybrid oncocytic renal neoplasm",Sample_file$Sample_Type),])
BarcodesHybridoncocyticchromophobe=rownames(Sample_file[grep("Hybrid oncocytic/Chromophobe type",Sample_file$Sample_Type),])

####################################
#Load relevant TCGA methylation data
###############################################

dat=load("datadir/MET_KIRP_Processed.Rdata")
METcancer_KIRP=ProcessedData$MET_Data_Cancer
METnormal_KIRP=ProcessedData$MET_Data_Normal
rm(ProcessedData)
METcancer_KIRP=TCGA_GENERIC_CleanUpSampleNames(METcancer_KIRP, 12) 
#n=275
METnormal_KIRP=TCGA_GENERIC_CleanUpSampleNames(METnormal_KIRP, 12) 
#n=44

#KICH
dat=load("~/Documents/Projects/TCGA_KICH/data/2015-05-11/MET_KICH_Processed.Rdata")
METcancer_KICH=ProcessedData$MET_Data_Cancer
#No normal data for KICH
rm(ProcessedData)
METcancer_KICH=TCGA_GENERIC_CleanUpSampleNames(METcancer_KICH, 12) 
#n=65 

#scp -r kbren@crosswood.stanford.edu:/srv/gevaertlab/data/TCGA/Kevin_Data/MethylationData450k/KIRC/gdac_20160715/ProcessedData_KIRC.RData ~/Documents/Projects/TCGA_KIRC/data/2015-04-15/
dat=load("~/Documents/Projects/TCGA_KIRC/data/2015-04-15/ProcessedData_KIRC.RData")
METcancer_KIRC=ProcessedData$MET_Data_Cancer
METnormal_KIRC=ProcessedData$MET_Data_Normal
rm(ProcessedData)

METcancer_KIRC=TCGA_GENERIC_CleanUpSampleNames(METcancer_KIRC, 12) 
#n=315
METnormal_KIRC=TCGA_GENERIC_CleanUpSampleNames(METnormal_KIRC, 12) 
#n=160

#
OverlapProbes=intersect(rownames(METnormal_KIRP),rownames(METnormal_KIRC))
METnormal_KIRP=METnormal_KIRP[OverlapProbes,]
METnormal_KIRC=METnormal_KIRC[OverlapProbes,]

METnormal=cbind(METnormal_KIRP, METnormal_KIRC)

commoncpgs=intersect(
     intersect(rownames(METcancer_KIRP), intersect(rownames(METcancer_KIRC), rownames(METcancer_KICH))),
     rownames(METnormal))
commoncpgs=commoncpgs[-c(grep("rs", commoncpgs))]
#n=395143

#restrict meth to just these common CpGs
METcancer_KIRP=METcancer_KIRP[commoncpgs,]
METcancer_KIRC=METcancer_KIRC[commoncpgs,]
METcancer_KICH=METcancer_KICH[commoncpgs,]
METnormal_KIRP=METnormal_KIRP[commoncpgs,]
#METnormal_KIRP=scale(METnormal_KIRP)
METnormal_KIRC=METnormal_KIRC[commoncpgs,]
#METnormal_KIRC=scale(METnormal_KIRC)

#Need to specify which colnames represent cancer and normal
colnames(METcancer_KIRP)=paste(colnames(METcancer_KIRP),"-01", sep="")
colnames(METcancer_KIRC)=paste(colnames(METcancer_KIRC),"-01", sep="")
colnames(METcancer_KICH)=paste(colnames(METcancer_KICH),"-01", sep="")
colnames(METnormal_KIRP)=paste(colnames(METnormal_KIRP),"-11", sep="")
colnames(METnormal_KIRC)=paste(colnames(METnormal_KIRC),"-11", sep="")

Barcodes_METcancer_KIRP=colnames(METcancer_KIRP)
Barcodes_METcancer_KIRC=colnames(METcancer_KIRC)
Barcodes_METcancer_KICH=colnames(METcancer_KICH)
Barcodes_METnormal_KIRP=colnames(METnormal_KIRP)
Barcodes_METnormal_KIRC=colnames(METnormal_KIRC)
Barcodes_METnormal=c(Barcodes_METnormal_KIRC, Barcodes_METnormal_KIRP)

METall=cbind(
     METcancer_KIRP,
     METcancer_KIRC,
     METcancer_KICH,
     METnormal_KIRP,
     METnormal_KIRC)

#################################
#Make sets of boxplots for all groups of genes we are interested in
#####################################

#Function to make a list ggplot2 boxplots
MakePlotlist=function(cgs){
 black.bold.16.text <- element_text(face = "bold", color = "black", size = 14)
 
     plotlist=list()    
     for(i in 1:length(cgs)){
          cg=cgs[i]  
          
          v1=c(METcancer_KIRP[cg,], METcancer_KIRC[cg,], METcancer_KICH[cg,], METnormal[cg,],
               beta[cg,colnames(beta) %in% BarcodesChromophobe],
               beta[cg,colnames(beta) %in% BarcodesOncocytoma],
               beta[cg,colnames(beta) %in% BarcodesNormal])
          
          v2=c(rep(1, length(METcancer_KIRP[cg,])), 
               rep(2, length(METcancer_KIRC[cg,])), 
               rep(3, length(METcancer_KICH[cg,])),
               rep(4, length(METnormal[cg,])),
               rep(3, length(beta[cg,colnames(beta) %in% BarcodesChromophobe])),
               rep(5, length(beta[cg,colnames(beta) %in% BarcodesOncocytoma])),
               rep(4, length(beta[cg,colnames(beta) %in% BarcodesNormal])))
          
          v3=c(rep(1, length(METcancer_KIRP[cg,])), 
               rep(1, length(METcancer_KIRC[cg,])), 
               rep(1, length(METcancer_KICH[cg,])),
               rep(1, length(METnormal[cg,])),
               rep(2, length(beta[cg,colnames(beta) %in% BarcodesChromophobe])),
               rep(2, length(beta[cg,colnames(beta) %in% BarcodesOncocytoma])),
               rep(2, length(beta[cg,colnames(beta) %in% BarcodesNormal])))
          
          cols=c(rep(1, length(METcancer_KIRP[cg,])), 
                 rep(2, length(METcancer_KIRC[cg,])), 
                 rep(3, length(METcancer_KICH[cg,])),
                 rep(4, length(METnormal[cg,])),
                 rep(3, length(beta[cg,colnames(beta) %in% BarcodesChromophobe])),
                 rep(5, length(beta[cg,colnames(beta) %in% BarcodesOncocytoma])),
                 rep(4, length(beta[cg,colnames(beta) %in% BarcodesNormal])))
          
          tab=as.data.frame(cbind(v1,v2,v3,cols))
          colnames(tab)=c("gene","class","study","col")
          
          tab$study=revalue(as.factor(tab$study),c("1"="TCGA","2"="Stanford"))
          
          
          p1=ggplot(tab, aes(factor(class), as.numeric(as.character(gene))))
          #p2=p1 + geom_boxplot( outlier.shape=NA) + geom_jitter(aes(colour = factor(class)))+ facet_grid(. ~ study, labeller="label_both")  + scale_fill_manual(name = "Sample type", values = c(1,2,3,4,5), labels = c("1" = "1", "2" = "2", "3" = "3", "4"="4", "5"="5"))+ theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),legend.position="none", legend.justification=c(.5,.5), legend.text=element_text(size=16, face="bold"),legend.title=element_text(size=16, face="bold"), plot.title = element_text(lineheight=.8, face="bold", size=20), axis.title = element_text(lineheight=.8, face="bold", size=16)) + theme(axis.text.y = black.bold.16.text)+ theme(strip.text.x = element_text(size = 14, colour = "black", face="bold")) + ggtitle(paste(cg,ann450k2[cg,"Gene"], sep="\n")) + geom_point() + scale_colour_manual(breaks = class, values = unique(as.character(ggplotColours(n = 6)))) + theme(plot.title = element_text(hjust = 0.5))
          p2=p1 + geom_boxplot(aes(fill = factor(class)), alpha=0.3, outlier.shape=NA) + geom_jitter(aes(colour = factor(class))) + facet_grid(. ~ study, labeller="label_both") + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),legend.position="none", legend.justification=c(.5,.5), 
                                                                                                                                                                                         legend.text=element_text(size=16, face="bold"),legend.title=element_text(size=16, face="bold"), plot.title = element_text(lineheight=.8, face="bold", size=20), axis.title = element_text(lineheight=.8, face="bold", size=16)) + theme(axis.text.y = black.bold.16.text)+ theme(strip.text.x = element_text(size = 14, colour = "black", face="bold")) + ggtitle(paste(cg,ann450k2[cg,"Gene"], sep="\n")) + theme(plot.title = element_text(hjust = 0.5))
                    
          plotlist[[i]]=p2
     }
return(plotlist)     
}

##############################################
############################################
#SAM analyses
###########################################
#############################################

#Make new working directory for all new SAM and PAM analyses
dir.create("resultsdir/SAMPAM_analyses")

#parameters for SAM analyses
OutcomeType='TwoClassUnpaired'
DeltaFDRthreshold=0.05
nrPermutations=100
StatisticalTest='wilcoxon'

library(samr)

#########################################################################
#Onc Vs normal
####################################################################

Barcodes=c(BarcodesOncocytoma, BarcodesNormal)

SampleClass=rep(NA, length(Barcodes))
names(SampleClass)=Barcodes
SampleClass[names(SampleClass) %in% BarcodesOncocytoma]=2
SampleClass[names(SampleClass) %in% BarcodesNormal]=1
SampleClass=t(as.matrix(SampleClass))
rownames(SampleClass)="Sample_type"

beta_Barcodes=beta[,colnames(SampleClass)]
beta_Barcodes=TCGA_GENERIC_GeneFiltering('MAD', beta_Barcodes, 25)

ResultsFile='resultsdir/SAM_analysis_OncVsNormal.txt'
SAMresults=SAM_Analysis_meth(beta_Barcodes, SampleClass, OutcomeType, DeltaFDRthreshold, nrPermutations,StatisticalTest, ResultsFile,'')

##################
#Format onc Vs normal 
#######################

AllResults=read.table("resultsdir/SAM_analysis_OncVsNormal.txt", header = T, sep="\t")
AllResults=AllResults[abs(AllResults$Fold.Change)>0.2 & AllResults[,4]=="0",]
AllResults=AllResults[,2:ncol(AllResults)]
colnames(AllResults)=c("probe","delta beta","q.value")
AllResults$direction=rep(NA, nrow(AllResults))
AllResults[AllResults[,2]>0,"direction"]="Up"
AllResults[AllResults[,2]<0,"direction"]="Down"
rownames(AllResults)=AllResults$probe
AllResults=AllResults[order(AllResults$`delta beta`, decreasing = T),]
AllResults=cbind(AllResults, ann450k2[rownames(AllResults),])

colnames=names(genelists)
mat=AllResults
mat$Gene=unlist(mat$Gene)
for(i in 1:length(genelists)){
     genelist=genelists[[colnames[i]]]
     genelist_Meth=lapply(mat$Gene,  function(x) intersect(as.character(unlist(strsplit(as.character(x),", "))),genelist))
     genelist_Meth[sapply(genelist_Meth, length)==0]=NA
     genelist_Meth=lapply(genelist_Meth,function(x) ifelse(length(x)>1,paste(x, collapse=", "),x))
     mat[,colnames[i]]=unlist(genelist_Meth)
}
AllResults=mat

AllResults=as.data.frame(apply(AllResults, 2, as.character))
write.table(AllResults, file="resultsdir/SAM_analysis_OncVsNormal_formatted.txt", sep="\t")

#Adding columns indicating probes that are below or above thresholds of mean methylation in normal tissues
AllResults=read.table("resultsdir/SAM_analysis_OncVsNormal_formatted.txt", sep="\t", header=T)
colnames=names(AllProbesNormBloodUrinary)
AllResults$probe=as.character(AllResults$probe)

for(i in 1:length(colnames)){
  colname=colnames[i]  
  AllResults[,colname]=ifelse(AllResults$probe %in% AllProbesNormBloodUrinary[[colnames[i]]],1,0)
}

write.table(AllResults, file="resultsdir/SAM_analysis_OncVsNormal_formatted.txt", sep="\t")

getplots=MakePlotlist(cgs=ROHyperLowlyMethProbes)

file="resultsdir/Boxplots_SAM_analysis_ROVsNormal_HyperLowlyMethylatedDMRs"
png(file=paste(file,'.png',sep=''), units="in", width=20, height=15, res=600)
do.call(grid.arrange,c(getplots[1:3], ncol=3))
dev.off()                                                                                                                                                                            

getplots=MakePlotlist(cgs=ROHypoHighlyMethProbes[1:10])

file="resultsdir/Boxplots_SAM_analysis_ROVsNormal_HypoHighlyMethylatedDMRs"
png(file=paste(file,'.png',sep=''), units="in", width=20, height=15, res=600)
do.call(grid.arrange,c(getplots[1:10], ncol=5))
dev.off()                                                                                                                                                                            


#write tables with genes that are shared between RO and ChRCC
GSEgenesAll=as.character(AllResults[which(nchar(as.character(AllResults$Delta.beta.ChRCC.versus.NKP.DRMs))>3),"Gene"])
GSEgenesAll=GSEgenesAll[nchar(GSEgenesAll)!=0]
GSEgenesAll=unique(unlist(strsplit(GSEgenesAll, ", ")))

write.table(GSEgenesAll, "resultsdir/GSE_genes_shared_DMRs_ChRCC_RO_all.txt", sep="\t", row.names = F, quote=F, col.names=F)

AllResults$delta.beta=as.numeric((as.character(AllResults$delta.beta)))
AllResults$Delta.beta.ChRCC.versus.NKP.DRMs=as.numeric((as.character(AllResults$Delta.beta.ChRCC.versus.NKP.DRMs)))

#hyper
GSEgenesAll=as.character(AllResults[which(AllResults$delta.beta>0 & AllResults$Delta.beta.ChRCC.versus.NKP.DRMs>0),"Gene"])
GSEgenesAll=GSEgenesAll[nchar(GSEgenesAll)!=0]
GSEgenesAll=unique(unlist(strsplit(GSEgenesAll, ", ")))

write.table(GSEgenesAll, "resultsdir/GSE_genes_shared_DMRs_ChRCC_RO_hyper.txt", sep="\t", quote=F, col.names=F)

#hypo
GSEgenesAll=as.character(AllResults[which(AllResults$delta.beta<0 & AllResults$Delta.beta.ChRCC.versus.NKP.DRMs<0),"Gene"])
GSEgenesAll=GSEgenesAll[nchar(GSEgenesAll)!=0]
GSEgenesAll=unique(unlist(strsplit(GSEgenesAll, ", ")))

write.table(GSEgenesAll, "resultsdir/GSE_genes_shared_DMRs_ChRCC_RO_hypo.txt", sep="\t", quote=F, col.names=F)

#######################################
#Make boxplots for most differentially methylated onc Vs norm
##################################

AllResults=read.table("resultsdir/SAM_analysis_OncVsNormal.txt", header = T, sep="\t")

cgs=as.character(rownames(AllResults[order(abs(AllResults$delta.beta), decreasing = T)[1:20],1:5]))
cgs=cgs[cgs %in% rownames(METall) & cgs %in% rownames(beta)]

getplots=MakePlotlist(cgs)

file="resultsdir/Boxplots_SAM_analysis_OncVsNormal"
png(file=paste(file,'.png',sep=''), units="in", width=20, height=15, res=600)
do.call(grid.arrange,c(getplots[1:10], ncol=5))
dev.off()

#########################################################################
#ChRCC Vs normal
####################################################################

Barcodes=c(BarcodesChromophobe, BarcodesNormal)

SampleClass=rep(NA, length(Barcodes))
names(SampleClass)=Barcodes
SampleClass[names(SampleClass) %in% BarcodesChromophobe]=2
SampleClass[names(SampleClass) %in% BarcodesNormal]=1
SampleClass=t(as.matrix(SampleClass))
rownames(SampleClass)="Sample_type"

beta_Barcodes=beta[,colnames(SampleClass)]
beta_Barcodes=TCGA_GENERIC_GeneFiltering('MAD', beta_Barcodes, 25)

ResultsFile='resultsdir/SAM_analysis_ChRCCVsNormal.txt'
SAMresults=SAM_Analysis_meth(beta_Barcodes, SampleClass, OutcomeType, DeltaFDRthreshold, nrPermutations,StatisticalTest, ResultsFile,'')

######################
#format ChRCC vs normal file
######################

AllResults=read.table("resultsdir/SAM_analysis_ChRCCVsNormal.txt", header = T, sep="\t")
AllResults=AllResults[abs(AllResults$Fold.Change)>0.2 & AllResults[,4]=="0",]
AllResults=AllResults[,2:ncol(AllResults)]
colnames(AllResults)=c("probe","delta beta","q.value")
AllResults$direction=rep(NA, nrow(AllResults))
AllResults[AllResults[,2]>0,"direction"]="Up"
AllResults[AllResults[,2]<0,"direction"]="Down"
rownames(AllResults)=AllResults$probe
AllResults=AllResults[order(AllResults$`delta beta`, decreasing = T),]
AllResults=cbind(AllResults, ann450k2[rownames(AllResults),])

colnames=names(genelists)
mat=AllResults
mat$Gene=unlist(mat$Gene)
for(i in 1:length(genelists)){
     genelist=genelists[[colnames[i]]]
     genelist_Meth=lapply(mat$Gene,  function(x) intersect(as.character(unlist(strsplit(as.character(x),", "))),genelist))
     genelist_Meth[sapply(genelist_Meth, length)==0]=NA
     genelist_Meth=lapply(genelist_Meth,function(x) ifelse(length(x)>1,paste(x, collapse=", "),x))
     mat[,colnames[i]]=unlist(genelist_Meth)
}
AllResults=mat

AllResults=as.data.frame(apply(AllResults, 2, as.character))
write.table(AllResults, file="resultsdir/SAM_analysis_ChRCCVsNormal_formatted.txt", sep="\t")

#Adding columns indicating probes that are below or above thresholds of mean methylation in normal tissues
AllResults=read.table("resultsdir/SAM_analysis_ChRCCVsNormal_formatted.txt", sep="\t", header=T)
colnames=names(AllProbesNormBloodUrinary)
AllResults$probe=as.character(AllResults$probe)

for(i in 1:length(colnames)){
  colname=colnames[i]  
  AllResults[,colname]=ifelse(AllResults$probe %in% AllProbesNormBloodUrinary[[colnames[i]]],1,0)
}

write.table(AllResults, file="resultsdir/SAM_analysis_ChRCCVsNormal_formatted.txt", sep="\t")

apply(AllResults[,colnames],2,function(x) summary(as.factor(x)))

#hyper DMRs
apply(AllResults[AllResults$direction=="Up",colnames],2,function(x) summary(as.factor(x)))

apply(AllResults[AllResults$direction=="Down",colnames],2,function(x) summary(as.factor(x)))

summary(as.factor(AllResults$direction))

length(which(AllResults$direction=="Up" & AllResults$LessThan0.1_normal==1))
#0

length(which(AllResults$direction=="Down" & AllResults$GreaterThan0.9_normal==1))
#25 hypo probes highly methylated across normal tissues
AllResults[which(AllResults$direction=="Down" & AllResults$GreaterThan0.9_normal==1),"probe"]

ChRCCHypoHighlyMethProbes=AllResults[which(AllResults$direction=="Down" & AllResults$GreaterThan0.9_normal==1),]
ChRCCHypoHighlyMethProbes=ChRCCHypoHighlyMethProbes[order(as.numeric(as.character(ChRCCHypoHighlyMethProbes$delta.beta)), decreasing=F),]
ChRCCHypoHighlyMethProbes=as.character(ChRCCHypoHighlyMethProbes$probe)

getplots=MakePlotlist(cgs=ChRCCHypoHighlyMethProbes[1:10])

file="resultsdir/Boxplots_SAM_analysis_ChRCCVsNormal_HypoHighlyMethylatedDMRs"
png(file=paste(file,'.png',sep=''), units="in", width=20, height=15, res=600)
do.call(grid.arrange,c(getplots[1:10], ncol=5))
dev.off()                                                                                                                                                                            

#######################################
#Make boxplots for most differentially methylated onc Vs norm
##################################

AllResults=read.table("resultsdir/SAM_analysis_ChRCCVsNormal_formatted.txt", sep="\t", header=T)

cgs=as.character(rownames(AllResults[order(abs(AllResults$delta.beta), decreasing = T)[1:20],1:5]))
cgs=cgs[cgs %in% rownames(METall) & cgs %in% rownames(beta)]

library(ggplot2)
library(gridExtra)
getplots=MakePlotlist(cgs)

file="resultsdir/Boxplots_SAM_analysis_ChRCCVsNormal"
png(file=paste(file,'.png',sep=''), units="in", width=20, height=15, res=600)
do.call(grid.arrange,c(getplots[1:10], ncol=5))
dev.off()

file="resultsdir/Boxplots_SAM_analysis_ChRCCVsNormal_legend"
png(file=paste(file,'.png',sep=''), units="in", width=20, height=15, res=600)
p1 + geom_boxplot(aes(fill = factor(class)), alpha=0.3, outlier.shape=NA) + geom_jitter(aes(colour = factor(class))) + facet_grid(. ~ study, labeller="label_both") + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), legend.text=element_blank(),legend.title=element_blank()) + theme(legend.key.size = unit(2, "cm"))
dev.off()                                                                                                                                                                            
 
#########################################################################
#Onc V ChRCC
####################################################################

#First do analysis fo difference between onc and chr
Barcodes=c(BarcodesOncocytoma, BarcodesChromophobe)

SampleClass=rep(NA, length(Barcodes))
names(SampleClass)=Barcodes
SampleClass[names(SampleClass) %in% BarcodesChromophobe]=1
SampleClass[names(SampleClass) %in% BarcodesOncocytoma]=2
SampleClass=t(as.matrix(SampleClass))
rownames(SampleClass)="Sample_type"

beta_Barcodes=beta[,colnames(SampleClass)]
beta_Barcodes=TCGA_GENERIC_GeneFiltering('MAD', beta_Barcodes, 25)

#ResultsFile='resultsdir/SAM_analysis_OncVChRCC.txt'
#SAMresults=SAM_Analysis_meth(beta_Barcodes, SampleClass, OutcomeType, DeltaFDRthreshold, nrPermutations,StatisticalTest, ResultsFile,'')

ResultsFile='resultsdir/SAM_analysis_OncVChRCC_retest.txt'
SAMresults=SAM_Analysis_meth(beta_Barcodes, SampleClass, OutcomeType, DeltaFDRthreshold, nrPermutations,StatisticalTest, ResultsFile,'')

############
#Format Onv Vs ChRCC file
##################

AllResults=read.table("resultsdir/SAM_analysis_OncVChRCC.txt", header = T, sep="\t")
AllResults=AllResults[abs(AllResults$Fold.Change)>0.2 & AllResults[,4]=="0",]
AllResults=AllResults[,2:ncol(AllResults)]
colnames(AllResults)=c("probe","delta beta","q.value")
AllResults$direction=rep(NA, nrow(AllResults))
AllResults[AllResults[,2]>0,"direction"]="Up"
AllResults[AllResults[,2]<0,"direction"]="Down"
rownames(AllResults)=AllResults$probe
AllResults=AllResults[order(AllResults$`delta beta`, decreasing = T),]
AllResults=cbind(AllResults, ann450k2[rownames(AllResults),])

colnames=names(genelists)
mat=AllResults
mat$Gene=unlist(mat$Gene)
for(i in 1:length(genelists)){
     genelist=genelists[[colnames[i]]]
     genelist_Meth=lapply(mat$Gene,  function(x) intersect(as.character(unlist(strsplit(as.character(x),", "))),genelist))
     genelist_Meth[sapply(genelist_Meth, length)==0]=NA
     genelist_Meth=lapply(genelist_Meth,function(x) ifelse(length(x)>1,paste(x, collapse=", "),x))
     mat[,colnames[i]]=unlist(genelist_Meth)
}
AllResults=mat

AllResults=as.data.frame(apply(AllResults, 2, as.character))
write.table(AllResults, file="resultsdir/SAM_analysis_OncVChRCC_formatted.txt", sep="\t")

#Adding columns indicating probes that are below or above thresholds of mean methylation in normal tissues
AllResults=read.table("resultsdir/SAM_analysis_OncVChRCC_formatted.txt", sep="\t", header=T)
colnames=names(AllProbesNormBloodUrinary)
AllResults$probe=as.character(AllResults$probe)

for(i in 1:length(colnames)){
  colname=colnames[i]  
  AllResults[,colname]=ifelse(AllResults$probe %in% AllProbesNormBloodUrinary[[colnames[i]]],1,0)
}

write.table(AllResults, file="resultsdir/SAM_analysis_OncVChRCC_formatted.txt", sep="\t")

#check what number of DMRs are within these methylation categories in normal tissue
apply(AllResults[,colnames],2,function(x) summary(as.factor(x)))

#get probes that differ between onc and norm, that also differ between eosinophillic and classic ChRCC
AllResultsEoVCl=read.table('resultsdir/SAM_analysis_EosinophillicVClassic.txt', sep="\t", header=T)

AllResults$delta.beta.ChRCC.eosinophillicVclassic=ifelse(AllResults$probe %in% AllResultsEoVCl$Gene.Name, AllResultsEoVCl$Fold.Change, NA)


AllResults$delta.beta.ChRCC.eosinophillicVclassic[AllResults$direction=="Up"]
#Six hypermethylted in RO are differentially methylated in Eosinophillic ChRCC, four hyper in eosinophillic
AllResults$delta.beta.ChRCC.eosinophillicVclassic[AllResults$direction=="Down"]
#two hypomethylted in RO are differentially methylated in Eosinophillic ChRCC, one hyper in eosinophillic

write.table(AllResults, file="resultsdir/SAM_analysis_OncVChRCC_formatted.txt", sep="\t")


AllResults=read.table("resultsdir/SAM_analysis_OncVChRCC_formatted.txt", sep="\t", header=T)

summary(as.factor(unlist(AllResults[AllResults$direction=="Up","Gene"])))
AllResults$Gene

unique(as.character(AllResults[AllResults$direction=="Up","Gene"]))

########################
#Gene set enrichments for onc V Chr
#############################

AllResults=read.table("resultsdir/SAM_analysis_OncVChRCC_formatted.txt", sep="\t", header=T)

AllResults$delta.beta=as.numeric((as.character(AllResults$delta.beta)))

#hyper
GSEgenesAll=as.character(AllResults[which(AllResults$delta.beta>0),"Gene"])
GSEgenesAll=GSEgenesAll[nchar(GSEgenesAll)!=0]
GSEgenesAll=unique(unlist(strsplit(GSEgenesAll, ", ")))
GSEgenesAllHyper=GSEgenesAll

write.table(GSEgenesAll, "resultsdir/GSE_genes_RO_Vs_ChRCC_hyper.txt", sep="\t", quote=F, col.names=F)

#hypo
GSEgenesAll=as.character(AllResults[which(AllResults$delta.beta<0),"Gene"])
GSEgenesAll=GSEgenesAll[nchar(GSEgenesAll)!=0]
GSEgenesAll=unique(unlist(strsplit(GSEgenesAll, ", ")))
GSEgenesAllHypo=GSEgenesAll

write.table(GSEgenesAll, "resultsdir/GSE_genes_RO_Vs_ChRCC_hypo.txt", sep="\t", quote=F, col.names=F)

#############################################
#Hypergeometric test onch chrcc
#####################################

GSEgenesAllHyper=checkGeneSymbols.KB(GSEgenesAllHyper)
#six unique genes hyper in onc
GSEgenesAllHypo=checkGeneSymbols.KB(GSEgenesAllHypo)

#Rohan differentially expressed between onc and chr
Rohan=read.table("datadir/Rohan_oncocytomaVchromophobe_gene_expr.txt", sep="\t", header=T)
Rohan$Gene[Rohan$Gene=="(N.A.)"]=NA
Rohan$Gene=as.character(Rohan$Gene)
Rohan=Rohan[nchar(Rohan$Gene)!=0,]
Rohan=Rohan[!is.na(Rohan$Gene),]
Rohan$Gene.HGNC=checkGeneSymbols.KB(Rohan$Gene)

Rohan_down=Rohan[!is.na(Rohan$Gene.HGNC) & Rohan$Direction.in.onc=="down","Gene.HGNC"]
Rohan_up=Rohan[!is.na(Rohan$Gene.HGNC) & Rohan$Direction.in.onc=="up","Gene.HGNC"]

intersect(GSEgenesAllHyper, Rohan_down)
intersect(GSEgenesAllHyper, Rohan_up)
intersect(GSEgenesAllHypo, Rohan_up)
intersect(GSEgenesAllHyper, Rohan_up)
#No intesections for any of the above

#######################################
#Make boxplots for most differentially methylated onc Vs ChRCC
##################################

AllResults=read.table("resultsdir/SAM_analysis_OncVChRCC_formatted.txt", sep="\t", header=T)
summary(as.factor(AllResults$direction))

cgs=as.character(AllResults[order(abs(AllResults$delta.beta), decreasing = T)[1:20],"probe"])
cgs=cgs[cgs %in% rownames(METall) & cgs %in% rownames(beta)]

getplots=MakePlotlist(cgs)

file="resultsdir/Boxplots_SAM_analysis_OncVChRCC"
png(file=paste(file,'.png',sep=''), units="in", width=20, height=15, res=600)
do.call(grid.arrange,c(getplots[1:10], ncol=5))
dev.off()

###############################
#PAM oncVs chrcc Max nr features 5,10,20,30,50,100
##################################

Barcodes=c(BarcodesOncocytoma, BarcodesChromophobe)

SampleClass=rep(NA, length(Barcodes))
names(SampleClass)=Barcodes
SampleClass[names(SampleClass) %in% BarcodesChromophobe]=1
SampleClass[names(SampleClass) %in% BarcodesOncocytoma]=2
SampleClass=t(as.matrix(SampleClass))
rownames(SampleClass)="Sample_type"

beta_Barcodes=beta[,colnames(SampleClass)]
beta_Barcodes=TCGA_GENERIC_GeneFiltering('MAD', beta_Barcodes, 25)

Ns=c(100,50,30,20,10,9,8,7,6,5)
#Ns=c(10,5,4)
names(PAMs)=Ns

PAMs=list()

#Running into issue where the number of genes is can be as few as one
#changes some bugs in M
Aucs=array(NA, c(length(Ns),3))
for(i in 1:length(Ns)){
  N=Ns[i]
  print(paste(N, " genes", sep=""))
  PAMresults=PAM_Analysis_MaxNrFeatures_CV_KB2(beta_genes_Barcodes,SampleClass,NrLambdaFolds=10,MaxNrFeatures=N)
  PAMs[[i]]=PAMresults
  ROCcurve=roc(as.factor(SampleClass-1),PAMresults$PosteriorProbabilities[,2])
  Aucs[i,1]=round(as.numeric(as.character(ci(ROCcurve))[[1]]),2)
  Aucs[i,2]=round(as.numeric(as.character(ci(ROCcurve))[[2]]),2)
  Aucs[i,3]=round(as.numeric(as.character(ci(ROCcurve))[[3]]),2)
}

rownames(Aucs)=Ns
colnames(Aucs)=c("95% CI (lower)","auc", "95% CI (upper)")

write.table(Aucs, file="resultsdir/PAM_MaxNrFeatures_Multiple_ChRCCVsnormal_Aucs.txt", sep="\t", row.names = F)

Missclassifications=lapply(PAMs, function(x) x$NrMisclassificationsPerClass)
Missclassifications=abind(Missclassifications, along=1)
rownames(Missclassifications)=c(1:nrow(Missclassifications))
Missclassifications=as.data.frame(Missclassifications)
Missclassifications$n.genes=rep(Ns,each=2)
Missclassifications$sample=rep(c("NKP","ChRCC"), length(Ns))

write.table(Missclassifications, file="resultsdir/PAM_MaxNrFeatures_Multiple_ChRCCVsnormal_Missclassifications.txt", sep="\t")

save(PAMs, file="resultsdir/PAM_MaxNrFeatures_Multiple_ChRCCVsnormal.RData")
file=load("resultsdir/PAM_MaxNrFeatures_Multiple_ChRCCVsnormal.RData")

##################################
#ChRCC Eosinophillic Vs. ChRCC classic (OncVChRCC genes)
#####################################

ann=read.table("datadir/KICH_TCGA_annotation.txt", sep="\t", header=T)
Barcodes_METcancer_KICH_Eosinophillic=as.character(ann[ann$Eosinophilic.vs.Classic=="Eosinophilic","TCGA.patient.code"])
Barcodes_METcancer_KICH_Eosinophillic=Barcodes_METcancer_KICH_Eosinophillic[Barcodes_METcancer_KICH_Eosinophillic %in% colnames(METcancer_KICH)]

Barcodes_METcancer_KICH_Classic=as.character(ann[ann$Eosinophilic.vs.Classic=="Classic","TCGA.patient.code"])
Barcodes_METcancer_KICH_Classic=Barcodes_METcancer_KICH_Classic[Barcodes_METcancer_KICH_Classic %in% colnames(METcancer_KICH)]

AllResults=read.table("resultsdir/SAM_analysis_OncVChRCC_formatted.txt", sep="\t", header=T)

genes=as.character(AllResults$probe)

Barcodes=c(Barcodes_METcancer_KICH_Eosinophillic, Barcodes_METcancer_KICH_Classic)

dat=load("~/Documents/Projects/TCGA_KICH/data/2015-05-11/MET_KICH_Processed.Rdata")
METcancer_KICH=ProcessedData$MET_Data_Cancer
#No normal data for KICH
rm(ProcessedData)
METcancer_KICH=TCGA_GENERIC_CleanUpSampleNames(METcancer_KICH, 12) 
METcancer_Barcodes=METcancer_KICH[rownames(METcancer_KICH) %in% genes,Barcodes]

SampleClass=rep(NA, length(Barcodes))
names(SampleClass)=Barcodes
SampleClass[names(SampleClass) %in% Barcodes_METcancer_KICH_Classic]=1
SampleClass[names(SampleClass) %in% Barcodes_METcancer_KICH_Eosinophillic]=2
SampleClass=t(as.matrix(SampleClass))
rownames(SampleClass)="Sample_type"

ResultsFile='resultsdir/SAM_analysis_EosinophillicVClassic.txt'
SAMresults=SAM_Analysis_meth(METcancer_Barcodes, SampleClass, OutcomeType, DeltaFDRthreshold, nrPermutations,StatisticalTest, ResultsFile,'')
#of 65 probes, 8 were differentially methylated between eosinophillic and classic 

wilc_METcancer_KICH_eosin_class=Rowwise_stat(dataframe=METcancer_Barcodes, var=as.factor(SampleClass[1,]))
wilc_METcancer_KICH_eosin_class=wilc_METcancer_KICH_eosin_class[order(wilc_METcancer_KICH_eosin_class$Diff),]
wilc_METcancer_KICH_eosin_class$Direction_in_onc_Vs_ChRCC=ifelse(rownames(wilc_METcancer_KICH_eosin_class) %in% as.character(AllResults[AllResults$direction=="Up","probe"]), "Hyper","Hypo")
#One of 65 probes, cg04601957 was significantly differentially methylated

write.table(wilc_METcancer_KICH_eosin_class, file="resultsdir/Wilcox_KICHEosinophillicVsClassic_.txt", sep="\t")


#############################################################
#######################################################
#Comparison of SAM analyses results between Onv Vs norm abd ChRCC and norm
##########################################################
##############################################################

AllResultsChRCC=read.table("resultsdir/SAM_analysis_ChRCCVsNormal_formatted.txt", sep="\t", header=T, colClasses=c("character"))
rownames(AllResultsChRCC)=AllResultsChRCC$probe

AllResultsOnc=read.table("resultsdir/SAM_analysis_OncVsNormal_formatted.txt", sep="\t", header=T, colClasses=c("character"))
rownames(AllResultsOnc)=AllResultsOnc$probe

#of 517 hyper DMRs in ChRCC, 297 are also hyper in oncocytoma
#962 hyper in oncocytoma altogether 

#3377/5058 (67%) of ChRCC hypo
#3377/8550 (39% of oncocytoma hypo probes intersect)
HyperIntersect=intersect(AllResultsChRCC[AllResultsChRCC$direction=="Up","probe"],AllResultsOnc[AllResultsOnc$direction=="Up","probe"])

HypoIntersect=intersect(AllResultsChRCC[AllResultsChRCC$direction=="Down","probe"],AllResultsOnc[AllResultsOnc$direction=="Down","probe"])

#Adding delta betas for other tumor type versus NKP, if the DMR is shared between tumors (in the same direction)
AllResultsChRCC$Delta.beta.Oncocytoma.versus.NKP.DRMs=unlist(lapply(rownames(AllResultsChRCC), function(x) ifelse(x %in% c(HyperIntersect, HypoIntersect), AllResultsOnc[x,"delta.beta"], NA)))

#Adding delta betas for other tumor type versus NKP, if the DMR is shared between tumors (in the same direction)
AllResultsOnc$Delta.beta.ChRCC.versus.NKP.DRMs=unlist(lapply(rownames(AllResultsOnc), function(x) ifelse(x %in% c(HyperIntersect, HypoIntersect), AllResultsChRCC[x,"delta.beta"], NA)))


delta.betas.compare=cbind(AllResultsOnc$delta.beta, AllResultsOnc$Delta.beta.ChRCC.versus.NKP.DRMs)
delta.betas.compare=delta.betas.compare[complete.cases(delta.betas.compare),]

file="resultsdir/RO_ChRCC_shared_DMRs_sqew"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
plot(as.numeric(as.character(delta.betas.compare[,1])), as.numeric(as.character(delta.betas.compare[,2])), ylab="delta beta (ChRCC Vs NKP)", xlab="delta beta (RO Vs NKP)", cex.axis=1.3, cex.lab=1.3)
dev.off()

#Add methylation levels in TCGA studies

write.table(AllResultsOnc, file="resultsdir/SAM_analysis_OncVsNormal_formatted.txt", sep="\t")
write.table(AllResultsChRCC, file="resultsdir/SAM_analysis_ChRCCVsNormal_formatted.txt", sep="\t")

############################################                          
#################################
#PAM analyses
#####################################              
################################################              

#PAM settings
NrFolds=10

#################
#PAM ChRCC Vs norm
#########################
              
Barcodes=c(BarcodesChromophobe, BarcodesNormal)

SampleClass=rep(NA, length(Barcodes))
names(SampleClass)=Barcodes
SampleClass[names(SampleClass) %in% BarcodesChromophobe]=2
SampleClass[names(SampleClass) %in% BarcodesNormal]=1
SampleClass=t(as.matrix(SampleClass))
rownames(SampleClass)="Sample_type"

SAM_genes=read.table('resultsdir/SAM_analysis_ChRCCVsNormal.txt', sep="\t", header=T)
SAM_genes=SAM_genes[SAM_genes[,4]=="0" & abs(SAM_genes[,3])>0.2,]
#5575 CpG probes

#First apply to all genes
#Next maxnrfeatures 5,10,20,30 genes
#Then try support vector machines with 1,2,3 genes
#Should PAM analysis be scaled? #Not scaling for moment

genes=as.character(SAM_genes[,2])
beta_genes_Barcodes=beta[genes,colnames(SampleClass)]
beta_genes_Barcodes=scale(beta_genes_Barcodes)

PAMresults=PAM_Analysis_CV_KB(beta_genes_Barcodes,SampleClass,NrFolds,NrLambdaFolds=2)
save(PAMresults, file="resultsdir/PAM_analysis_CV_ChRCCVsNormal_all_genes.RData")
 
file=load("resultsdir/PAM_analysis_CV_ChRCCVsNormal_all_genes.RData")

PAMresults$OverallPerformance
     NrMisclassifications NrMisclassificationsPerc  Accuracy
[1,]                    2               0.08695652 0.9130435
     
PAMresults$NrMisclassificationsPerClass
       NrMistakes TotalNrCasesInClass NrMisclassificationsPerc Accuracy
Class1          0                  15                     0.00     1.00
Class2          2                   8                     0.25     0.75
 
ROCcurve=roc(as.factor(SampleClass-1),PAMresults$Predictions[,1])
           
file="resultsdir/roc_PAM_ChRCCVsnormal"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
plot(ROCcurve, cex.lab=1.3, cex.axis=1.3, main=paste("AUC = ",gsub(".* ","",auc(ROCcurve))," (", gsub("[(DeLong)]","", capture.output(suppressWarnings(ci(ROCcurve)))),")", sep=""), cex.main=1.3)
dev.off()

###############################
##PAM ChRCC Vs norm heatmap classifications
################################

library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
library(RColorBrewer)

file=load("resultsdir/PAM_analysis_CV_ChRCCVsNormal_all_genes.RData")

PosteriorProbabilities=PAMresults$PosteriorProbabilities
colnames(PosteriorProbabilities)=c("NKP","ChRCC")

PosteriorProbabilities=as.data.frame(PosteriorProbabilities)

PosteriorProbabilities$Sample_ID=Sample_file[rownames(PosteriorProbabilities),"Sample_ID2"]
PosteriorProbabilities$Sample_Type2=Sample_file[rownames(PosteriorProbabilities),"Sample_Type2"]

#Colorblind palette
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#FD61D1", "#A3A500", "#D55E00", "#CC79A7")

#make colorblind version 
rlab=as.matrix(rownames(PosteriorProbabilities))
rownames(rlab)=rlab[,1]
colnames(rlab)=""
rlab[rownames(rlab) %in% BarcodesNormal,1]=cbp2[4]
rlab[rownames(rlab) %in% BarcodesChromophobe,1]=cbp2[3]

file="resultsdir/heatmap_posterior_PAM_ChRCCVsnormal"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
par(oma=c(6,0,0,0))
heatmap.3(PosteriorProbabilities[,c(1:2)], Rowv=FALSE, Colv=FALSE, dendrogram = "none", col=brewer.pal(9,"Reds"), RowSideColors=t(rlab), KeyValueName="Classification probability", labRow=PosteriorProbabilities$Sample_ID)
dev.off()

#Make legend 
file="resultsdir/heatmap_posterior_PAM_ChRCCVsnormal_legend"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("center", c("chRCC/KICH","NKP"), col=c(cbp2[3:4]), pch=15, title="Sample type", cex=2, bty="n")
dev.off()

###############################
#PAM ChRCC Vs norm Max nr features 5,10,20,30,50,100
##################################

Ns=c(100,50,30,20,10,9,8,7,6,5,4,3)
#Ns=c(10,5,4)
names(PAMs)=Ns

PAMs=list()

Aucs=array(NA, c(length(Ns),3))
for(i in 1:length(Ns)){
  N=Ns[i]
  PAMresults=PAM_Analysis_MaxNrFeatures_CV_KB(beta_genes_Barcodes,SampleClass,NrLambdaFolds=10,MaxNrFeatures=N)
  PAMs[[i]]=PAMresults
  ROCcurve=roc(as.factor(SampleClass-1),PAMresults$PosteriorProbabilities[,2])
  Aucs[i,1]=round(as.numeric(as.character(ci(ROCcurve))[[1]]),2)
  Aucs[i,2]=round(as.numeric(as.character(ci(ROCcurve))[[2]]),2)
  Aucs[i,3]=round(as.numeric(as.character(ci(ROCcurve))[[3]]),2)
}

rownames(Aucs)=Ns
colnames(Aucs)=c("95% CI (lower)","auc", "95% CI (upper)")

write.table(Aucs, file="resultsdir/PAM_MaxNrFeatures_Multiple_ChRCCVsnormal_Aucs.txt", sep="\t", row.names = F)

Missclassifications=lapply(PAMs, function(x) x$NrMisclassificationsPerClass)
Missclassifications=abind(Missclassifications, along=1)
rownames(Missclassifications)=c(1:nrow(Missclassifications))
Missclassifications=as.data.frame(Missclassifications)
Missclassifications$n.genes=rep(Ns,each=2)
Missclassifications$sample=rep(c("NKP","ChRCC"), length(Ns))

write.table(Missclassifications, file="resultsdir/PAM_MaxNrFeatures_Multiple_ChRCCVsnormal_Missclassifications.txt", sep="\t")

save(PAMs, file="resultsdir/PAM_MaxNrFeatures_Multiple_ChRCCVsnormal.RData")
file=load("resultsdir/PAM_MaxNrFeatures_Multiple_ChRCCVsnormal.RData")

######################
#remake boxplots but with misclassifications indicated
######################

#get CpG sites
AllResults=read.table("resultsdir/SAM_analysis_ChRCCVsNormal_formatted.txt", sep="\t", header=T)

cgs=as.character(rownames(AllResults[order(abs(AllResults$delta.beta), decreasing = T)[1:20],1:5]))
cgs=cgs[cgs %in% rownames(METall) & cgs %in% rownames(beta)]

#define sample class
Barcodes=c(BarcodesChromophobe, BarcodesNormal)

SampleClass=rep(NA, length(Barcodes))
names(SampleClass)=Barcodes
SampleClass[names(SampleClass) %in% BarcodesChromophobe]=2
SampleClass[names(SampleClass) %in% BarcodesNormal]=1
SampleClass=t(as.matrix(SampleClass))
rownames(SampleClass)="Sample_type"

Sample_file$Sample_ID2=gsub("RO_","onc_", Sample_file$Sample_ID2)

#get pam pam classifications
file=load("resultsdir/PAM_analysis_CV_ChRCCVsNormal_all_genes.RData")
Predictions=as.data.frame(PAMresults$Predictions)
BarcodesMisclassified=as.character(names(which(Predictions$V1!=SampleClass[1,])))
Sample_file[BarcodesMisclassified,]

#Make new plotting function but showing misclassifications
MakePlotlist_misclassifications=function(cgs){
black.bold.16.text <- element_text(face = "bold", color = "black", size = 14)
 
     plotlist=list()    
     for(i in 1:length(cgs)){
          cg=cgs[i]  
          
          v1=c(METcancer_KIRP[cg,], METcancer_KIRC[cg,], METcancer_KICH[cg,], METnormal[cg,],
               beta[cg,colnames(beta) %in% BarcodesChromophobe],
               beta[cg,colnames(beta) %in% BarcodesOncocytoma],
               beta[cg,colnames(beta) %in% BarcodesNormal])
          
          v2=c(rep(1, length(METcancer_KIRP[cg,])), 
               rep(2, length(METcancer_KIRC[cg,])), 
               rep(3, length(METcancer_KICH[cg,])),
               rep(4, length(METnormal[cg,])),
               rep(3, length(beta[cg,colnames(beta) %in% BarcodesChromophobe])),
               rep(5, length(beta[cg,colnames(beta) %in% BarcodesOncocytoma])),
               rep(4, length(beta[cg,colnames(beta) %in% BarcodesNormal])))
          
          v3=c(rep(1, length(METcancer_KIRP[cg,])), 
               rep(1, length(METcancer_KIRC[cg,])), 
               rep(1, length(METcancer_KICH[cg,])),
               rep(1, length(METnormal[cg,])),
               rep(2, length(beta[cg,colnames(beta) %in% BarcodesChromophobe])),
               rep(2, length(beta[cg,colnames(beta) %in% BarcodesOncocytoma])),
               rep(2, length(beta[cg,colnames(beta) %in% BarcodesNormal])))
          
          cols=c(rep(1, length(METcancer_KIRP[cg,])), 
                 rep(2, length(METcancer_KIRC[cg,])), 
                 rep(3, length(METcancer_KICH[cg,])),
                 rep(4, length(METnormal[cg,])),
                 rep(3, length(beta[cg,colnames(beta) %in% BarcodesChromophobe])),
                 rep(5, length(beta[cg,colnames(beta) %in% BarcodesOncocytoma])),
                 rep(4, length(beta[cg,colnames(beta) %in% BarcodesNormal])))
          
          #highlight misclassifications
          misclass=c(rep(1, length(METcancer_KIRP[cg,])), 
               rep(1, length(METcancer_KIRC[cg,])), 
               rep(1, length(METcancer_KICH[cg,])),
               rep(1, length(METnormal[cg,])),
               rep(1, length(beta[cg,colnames(beta) %in% BarcodesChromophobe])),
               rep(1, length(beta[cg,colnames(beta) %in% BarcodesOncocytoma])),
               rep(1, length(beta[cg,colnames(beta) %in% BarcodesNormal])))
          names(misclass)=c(names(METcancer_KIRP[cg,]), names(METcancer_KIRC[cg,]), names(METcancer_KICH[cg,]), names(METnormal[cg,]),
               names(beta[cg,colnames(beta) %in% BarcodesChromophobe]),
               names(beta[cg,colnames(beta) %in% BarcodesOncocytoma]),
               names(beta[cg,colnames(beta) %in% BarcodesNormal]))
          misclass[names(misclass) %in% BarcodesMisclassified]=4
          
          tab=as.data.frame(cbind(v1,v2,v3,cols))
          colnames(tab)=c("gene","class","study","col")
          
          tab$study=revalue(as.factor(tab$study),c("1"="TCGA","2"="Stanford"))
          
          
          p1=ggplot(tab, aes(factor(class), as.numeric(as.character(gene))))
          #p2=p1 + geom_boxplot( outlier.shape=NA) + geom_jitter(aes(colour = factor(class)))+ facet_grid(. ~ study, labeller="label_both")  + scale_fill_manual(name = "Sample type", values = c(1,2,3,4,5), labels = c("1" = "1", "2" = "2", "3" = "3", "4"="4", "5"="5"))+ theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),legend.position="none", legend.justification=c(.5,.5), legend.text=element_text(size=16, face="bold"),legend.title=element_text(size=16, face="bold"), plot.title = element_text(lineheight=.8, face="bold", size=20), axis.title = element_text(lineheight=.8, face="bold", size=16)) + theme(axis.text.y = black.bold.16.text)+ theme(strip.text.x = element_text(size = 14, colour = "black", face="bold")) + ggtitle(paste(cg,ann450k2[cg,"Gene"], sep="\n")) + geom_point() + scale_colour_manual(breaks = class, values = unique(as.character(ggplotColours(n = 6)))) + theme(plot.title = element_text(hjust = 0.5))
          #p2=p1 + geom_boxplot(aes(fill = factor(class)), alpha=0.3, outlier.shape=NA) + geom_jitter(aes(colour = factor(class), shape=as.character(as.numeric(misclass)))) + facet_grid(. ~ study, labeller="label_both") + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),legend.position="none", legend.justification=c(.5,.5), legend.text=element_text(size=16, face="bold"),legend.title=element_text(size=16, face="bold"), plot.title = element_text(lineheight=.8, face="bold", size=20), axis.title = element_text(lineheight=.8, face="bold", size=16)) + theme(axis.text.y = black.bold.16.text)+ theme(strip.text.x = element_text(size = 14, colour = "black", face="bold")) + ggtitle(paste(cg,ann450k2[cg,"Gene"], sep="\n")) + theme(plot.title = element_text(hjust = 0.5))
          #this version adds labels to misclassified with transparency of lables controlled by alpha
#p2=p1 + geom_boxplot(aes(fill = factor(class)), alpha=0.3, outlier.shape=NA) + geom_jitter(aes(colour = factor(class), shape=as.character(as.numeric(misclass)))) + facet_grid(. ~ study, labeller="label_both") + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),legend.position="none", legend.justification=c(.5,.5), legend.text=element_text(size=16, face="bold"),legend.title=element_text(size=16, face="bold"), plot.title = element_text(lineheight=.8, face="bold", size=20), axis.title = element_text(lineheight=.8, face="bold", size=16)) + theme(axis.text.y = black.bold.16.text)+ theme(strip.text.x = element_text(size = 14, colour = "black", face="bold")) + ggtitle(paste(cg,ann450k2[cg,"Gene"], sep="\n")) + theme(plot.title = element_text(hjust = 0.5)) + geom_text(aes(label=ifelse(misclass==4, Sample_file[as.character(names(misclass)),"Sample_ID2"],'')),hjust=0, vjust="left", alpha=0.5, size=3)
#p2=p1 + geom_boxplot(aes(fill = tab$class2), alpha=0.3, outlier.shape=NA) + scale_color_manual(values=c(cbp2[1:5])) + scale_fill_manual(values=c(cbp2[1:5])) + geom_jitter(aes(colour = as.character(tab$class2), shape=as.character(as.numeric(misclass)))) + facet_grid(. ~ study, labeller="label_both") + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),legend.position="none", legend.justification=c(.5,.5), legend.text=element_text(size=16, face="bold"),legend.title=element_text(size=16, face="bold"), plot.title = element_text(lineheight=.8, face="bold", size=20), axis.title = element_text(lineheight=.8, face="bold", size=16)) + theme(axis.text.y = black.bold.16.text)+ theme(strip.text.x = element_text(size = 14, colour = "black", face="bold")) + ggtitle(paste(cg,ann450k2[cg,"Gene"], sep="\n")) + theme(plot.title = element_text(hjust = 0.5)) + geom_text(aes(label=ifelse(misclass==4, Sample_file[as.character(names(misclass)),"Sample_ID2"],'')),hjust=0, vjust="left", alpha=0.5, size=3)
p2=p1 + geom_boxplot(aes(fill = factor(class)), alpha=0.3, outlier.shape=NA) + scale_color_manual(values=c(cbp2[1:5])) + scale_fill_manual(values=c(cbp2[1:5])) + geom_jitter(aes(colour = factor(class), shape=as.character(as.numeric(misclass)))) + facet_grid(. ~ study, labeller="label_both") + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),legend.position="none", legend.justification=c(.5,.5), legend.text=element_text(size=16, face="bold"),legend.title=element_text(size=16, face="bold"), plot.title = element_text(lineheight=.8, face="bold", size=20), axis.title = element_text(lineheight=.8, face="bold", size=16)) + theme(axis.text.y = black.bold.16.text)+ theme(strip.text.x = element_text(size = 14, colour = "black", face="bold")) + ggtitle(paste(cg,ann450k2[cg,"Gene"], sep="\n")) + theme(plot.title = element_text(hjust = 0.5)) + geom_text(aes(label=ifelse(misclass==4, as.character(Sample_file[as.character(names(misclass)),"Sample_ID2"]),'')),hjust=0, alpha=0.5, size=4, angle=-135)

          plotlist[[i]]=p2
     }
return(plotlist)     
}

getplots=MakePlotlist_misclassifications(cgs)

file="resultsdir/Boxplots_SAM_analysis_ChRCCVsNormal_TCGA_with_misclassifications"
png(file=paste(file,'.png',sep=''), units="in", width=20, height=15, res=600)
do.call(grid.arrange,c(getplots[1:10], ncol=5))
dev.off()

###################
#PAM ChRCC Vs norm Max nr features: just hyper geens
#########################

genes=as.character(SAM_genes[SAM_genes$Fold.Change>0,"Gene.Name"])
beta_genes_Barcodes=beta[genes,colnames(SampleClass)]
beta_genes_Barcodes=scale(beta_genes_Barcodes)

PAMs=list()
names(PAMs)=Ns

Aucs=array(NA, c(length(Ns),3))
for(i in 1:length(Ns)){
  N=Ns[i]
  PAMresults=PAM_Analysis_MaxNrFeatures_CV_KB(beta_genes_Barcodes,SampleClass,NrLambdaFolds=10,MaxNrFeatures=N)
  PAMs[[i]]=PAMresults
  ROCcurve=roc(as.factor(SampleClass-1),PAMresults$PosteriorProbabilities[,2])
  Aucs[i,1]=round(as.numeric(as.character(ci(ROCcurve))[[1]]),2)
  Aucs[i,2]=round(as.numeric(as.character(ci(ROCcurve))[[2]]),2)
  Aucs[i,3]=round(as.numeric(as.character(ci(ROCcurve))[[3]]),2)
}

rownames(Aucs)=Ns
colnames(Aucs)=c("95% CI (lower)","auc", "95% CI (upper)")

write.table(Aucs, file="resultsdir/PAM_MaxNrFeatures_Multiple_ChRCCVsnormal_JustHyper_Aucs.txt", sep="\t")

save(PAMs, file="resultsdir/PAM_MaxNrFeatures_Multiple_ChRCCVsnormal_JustHyper.RData")NrMisclassificationsPerClass

Missclassifications=lapply(PAMs, function(x) x$NrMisclassificationsPerClass)
Missclassifications=abind(Missclassifications, along=1)
rownames(Missclassifications)=c(1:nrow(Missclassifications))
Missclassifications=as.data.frame(Missclassifications)
Missclassifications$n.genes=rep(Ns,each=2)
Missclassifications$sample=rep(c("NKP","ChRCC"), length(Ns))

write.table(Missclassifications, file="resultsdir/PAM_MaxNrFeatures_Multiple_ChRCCVsnormal_JustHyper_Missclassifications.txt", sep="\t")

beta_genes_Barcodes=beta[genes,colnames(SampleClass)]
probe="cg27630540"
stripchart(beta_genes_Barcodes[probe,]~SampleClass[1,], vertical=T, method="jitter")

#Looks great, but not unmethylated in normal tissue

#####################################
#PAM onc Vs NKP
#####################################

Barcodes=c(BarcodesOncocytoma, BarcodesNormal)

SampleClass=rep(NA, length(Barcodes))
names(SampleClass)=Barcodes
SampleClass[names(SampleClass) %in% BarcodesOncocytoma]=2
SampleClass[names(SampleClass) %in% BarcodesNormal]=1
SampleClass=t(as.matrix(SampleClass))
rownames(SampleClass)="Sample_type"

SAM_genes=read.table('resultsdir/SAM_analysis_OncVsNormal.txt', sep="\t", header=T)
SAM_genes=SAM_genes[SAM_genes[,4]=="0" & abs(SAM_genes[,3])>0.2,]
#9512 CpG probes

#First apply to all genes
#Next maxnrfeatures 5,10,20,30 genes
#Then try support vector machines with 1,2,3 genes
#Should PAM analysis be scaled? #Not scaling for moment

genes=as.character(SAM_genes[,2])
beta_genes_Barcodes=beta[genes,colnames(SampleClass)]
beta_genes_Barcodes=scale(beta_genes_Barcodes)

PAMresults=PAM_Analysis_CV_KB(beta_genes_Barcodes,SampleClass,NrFolds,NrLambdaFolds=2)
save(PAMresults, file="resultsdir/PAM_analysis_CV_OncVsNormal_all_genes.RData")
file=load("resultsdir/PAM_analysis_CV_OncVsNormal_all_genes.RData")

PAMresults$OverallPerformance
     NrMisclassifications NrMisclassificationsPerc Accuracy
[1,]                    1               0.03703704 0.962963
     
PAMresults$NrMisclassificationsPerClass
       NrMistakes TotalNrCasesInClass NrMisclassificationsPerc  Accuracy
Class1          0                  15               0.00000000 1.0000000
Class2          1                  12               0.08333333 0.9166667
     
 
ROCcurve=roc(as.factor(SampleClass-1),PAMresults$Predictions[,1])
           
file="resultsdir/roc_PAM_OncVsnormal"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
plot(ROCcurve, cex.lab=1.3, cex.axis=1.3, main=paste("AUC = ",gsub(".* ","",auc(ROCcurve))," (", gsub("[(DeLong)]","", capture.output(suppressWarnings(ci(ROCcurve)))),")", sep=""), cex.main=1.3)
dev.off()

##############################
#posterior probabilities heatmap RO versus NKP
##################################

file=load("resultsdir/PAM_analysis_CV_OncVsNormal_all_genes.RData")

PosteriorProbabilities=PAMresults$PosteriorProbabilities
colnames(PosteriorProbabilities)=c("NKP","RO")
PosteriorProbabilities=as.data.frame(PosteriorProbabilities)

PosteriorProbabilities$Sample_ID=Sample_file[rownames(PosteriorProbabilities),"Sample_ID2"]
PosteriorProbabilities$Sample_Type2=Sample_file[rownames(PosteriorProbabilities),"Sample_Type2"]
PosteriorProbabilities$Sample_ID=as.factor(gsub("RO","oncocytoma",PosteriorProbabilities$Sample_ID))

rlab=as.matrix(rownames(PosteriorProbabilities))
rownames(rlab)=rlab[,1]
colnames(rlab)=""
rlab[rownames(rlab) %in% BarcodesNormal,1]=cbp2[4]
rlab[rownames(rlab) %in% BarcodesOncocytoma,1]=cbp2[5]

file="resultsdir/heatmap_posterior_PAM_ROVsnormal"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
par(oma=c(5,0,0,5))
heatmap.3(PosteriorProbabilities[,c(1:2)], Rowv=FALSE, Colv=FALSE, dendrogram = "none", col=brewer.pal(9,"Reds"), RowSideColors=t(rlab), KeyValueName="Classification probability", labRow=PosteriorProbabilities$Sample_ID, labCol=c("NKP","oncocytoma"))
dev.off()

file="resultsdir/heatmap_posterior_PAM_ROVsnormal_legend"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("center", c("oncocytoma","NKP"), col=c(cbp2[c(5,4)]), pch=15, title="Sample type", cex=2, bty="n")
dev.off()

##############################
#Boxplots onv v norm, but including misclassification indices
###############################

AllResults=read.table("resultsdir/SAM_analysis_OncVsNormal_formatted.txt", sep="\t", header=T)

cgs=as.character(rownames(AllResults[order(abs(AllResults$delta.beta), decreasing = T)[1:20],1:5]))
cgs=cgs[cgs %in% rownames(METall) & cgs %in% rownames(beta)]

#define sample class
Barcodes=c(BarcodesOncocytoma, BarcodesNormal)

SampleClass=rep(NA, length(Barcodes))
names(SampleClass)=Barcodes
SampleClass[names(SampleClass) %in% BarcodesOncocytoma]=2
SampleClass[names(SampleClass) %in% BarcodesNormal]=1
SampleClass=t(as.matrix(SampleClass))
rownames(SampleClass)="Sample_type"

#get pam pam classifications
file=load("resultsdir/PAM_analysis_CV_OncVsNormal_all_genes.RData")
Predictions=as.data.frame(PAMresults$Predictions)
BarcodesMisclassified=as.character(names(which(Predictions[,1]!=SampleClass[1,])))
#BarcodesMisclassified=gsub(".*_","",BarcodesMisclassified)
#200397540095_R03C01
Sample_file["200397540095_R03C01",]
#34367755-T
Sample_file[Sample_file$Sample_ID=="JL076-T",]
#oncocytoma_22
Sample_file[Sample_file$Sample_ID=="RCC103-T",]
#oncocytoma_23

getplots=MakePlotlist_misclassifications(cgs)

file="resultsdir/Boxplots_SAM_analysis_OncVsNormal_with_misclassifications"
png(file=paste(file,'.png',sep=''), units="in", width=20, height=15, res=600)
do.call(grid.arrange,c(getplots[1:10], ncol=5))
dev.off()


###############################
#OncVsChRCC (all genes that differ between them)
#######################################
  
SAM_genes=read.table('resultsdir/SAM_analysis_OncVChRCC.txt', sep="\t", header=T)
SAM_genes=SAM_genes[SAM_genes[,4]=="0" & abs(SAM_genes[,3])>0.2,]
roc_OncVsChRCC_PAM20_genes_ranks.png#79 genes

Barcodes=c(BarcodesOncocytoma, BarcodesChromophobe)

SampleClass=rep(NA, length(Barcodes))
names(SampleClass)=Barcodes
SampleClass[names(SampleClass) %in% BarcodesChromophobe]=1
SampleClass[names(SampleClass) %in% BarcodesOncocytoma]=2
SampleClass=t(as.matrix(SampleClass))
rownames(SampleClass)="Sample_type"

genes=as.character(SAM_genes[,2])
beta_genes_Barcodes=beta[genes,colnames(SampleClass)]
beta_genes_Barcodes=scale(beta_genes_Barcodes)

PAMresults=PAM_Analysis_CV_KB(beta_genes_Barcodes,SampleClass,NrFolds,NrLambdaFolds=2)
save(PAMresults, file="resultsdir/PAM_analysis_CV_OncVsChRCC_all_genes.RData")
file=load("resultsdir/PAM_analysis_CV_OncVsChRCC_all_genes.RData")


PAMresults$OverallPerformance
     NrMisclassifications NrMisclassificationsPerc Accuracy
[1,]                    1                     0.05     0.95
     
PAMresults$NrMisclassificationsPerClass
       NrMistakes TotalNrCasesInClass NrMisclassificationsPerc  Accuracy
Class1          0                   8               0.00000000 1.0000000
Class2          1                  12               0.08333333 0.9166667


ROCcurve=roc(as.factor(SampleClass-1),PAMresults$Predictions[,1])
           
file="resultsdir/roc_OncVsChRCC_all_genes"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
plot(ROCcurve, cex.lab=1.3, cex.axis=1.3, main=paste("AUC = ",gsub(".* ","",auc(ROCcurve))," (", gsub("[(DeLong)]","", capture.output(suppressWarnings(ci(ROCcurve)))),")", sep=""), cex.main=1.3)
dev.off()

Predictions=as.data.frame(PAMresults$Predictions)
BarcodesMisclassified=as.character(names(which(Predictions[,1]!=SampleClass[1,])))


##############################
#posterior probabilities heatmap ChRCC versus RO
##################################

file=load("resultsdir/PAM_analysis_CV_OncVsChRCC_all_genes.RData")

PosteriorProbabilities=PAMresults$PosteriorProbabilities
colnames(PosteriorProbabilities)=c("ChRCC","RO")
PosteriorProbabilities=as.data.frame(PosteriorProbabilities)

PosteriorProbabilities$Sample_ID=Sample_file[rownames(PosteriorProbabilities),"Sample_ID2"]
PosteriorProbabilities$Sample_Type2=Sample_file[rownames(PosteriorProbabilities),"Sample_Type2"]

PosteriorProbabilities$Sample_ID=gsub("RO","oncocytoma",PosteriorProbabilities$Sample_ID)

rlab=as.matrix(rownames(PosteriorProbabilities))
rownames(rlab)=rlab[,1]
colnames(rlab)=""
rlab[rownames(rlab) %in% BarcodesChromophobe,1]=cbp2[3]
rlab[rownames(rlab) %in% BarcodesOncocytoma,1]=cbp2[5]

file="resultsdir/heatmap_posterior_PAM_ChRCCVsRO"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
#par(oma=c(8,0,0,0))
#heatmap.3(PosteriorProbabilities[,c(1:2)], Rowv=FALSE, Colv=FALSE, dendrogram = "none", col=brewer.pal(9,"Reds"), RowSideColors=t(rlab), KeyValueName="Classification probability", labRow=PosteriorProbabilities$Sample_ID)
par(oma=c(5,0,0,5))
heatmap.3(PosteriorProbabilities[,c(1:2)], Rowv=FALSE, Colv=FALSE, dendrogram = "none", col=brewer.pal(9,"Reds"), RowSideColors=t(rlab), KeyValueName="Classification probability", labRow=PosteriorProbabilities$Sample_ID, labCol="")
dev.off()

#######################################################################
#Make supplementary table showing classification to chRCC and RO samples from ten fold cross-validation
########################################################################

file=load("resultsdir/PAM_analysis_CV_OncVsChRCC_all_genes.RData")

#
TestPredictions=as.data.frame(PAMresults$Predictions, drop=F)
TestPredictions=cbind(TestPredictions, PosteriorProbabilities)
colnames(TestPredictions)[1]="PAMprediction"
TestPredictions$PAMprediction=revalue(as.factor(TestPredictions$PAMprediction), c("1"="chRCC/KICH","2"="oncocytoma"))
TestPredictions$Sample_Type2=revalue(as.factor(TestPredictions$Sample_Type2), c("AML"="Angiomylolipoma", "chRCC"="chRCC/KICH", "KIRC"="ccRCC/KIRC", "KIRP"="pRCC/KIRP","Normal"="NKP","oncocytoma"="oncocytoma"))
TestPredictions$Misclassified=rep("No", nrow(TestPredictions))
TestPredictions[TestPredictions$Sample_Type2=="oncocytoma" & TestPredictions$PAMprediction!="oncocytoma","Misclassified"]="Yes"
TestPredictions[TestPredictions$Sample_Type2=="chRCC/KICH" & TestPredictions$PAMprediction!="chRCC/KICH","Misclassified"]="Yes"
TestPredictions=TestPredictions[,c(4,5,1,6,2,3)]
colnames(TestPredictions)=c("Sample ID","Tumor type","Diagnostic model prediction", "Misclassified by diagnostic model?","Posterior probability chRCC", "Posterior probability oncocytoma")

write.table(TestPredictions, file="resultsdir/PAM_analysis_CV_OncVsChRCC_all_genes_table.txt", sep="\t", row.names = FALSE)


################################

Ns=c(5,10,15,20,25,30,35,40)
#Ns=c(10,5,4)

PAMs=list()

Aucs=array(NA, c(length(Ns),3))
for(i in 1:length(Ns)){
  N=Ns[i]
  PAMresults=PAM_Analysis_MaxNrFeatures_CV_KB2(beta_genes_Barcodes,SampleClass,NrLambdaFolds=10,MaxNrFeatures=N)
  PAMs[[i]]=PAMresults
  ROCcurve=roc(as.factor(SampleClass-1),PAMresults$PosteriorProbabilities[,2])
  Aucs[i,1]=round(as.numeric(as.character(ci(ROCcurve))[[1]]),2)
  Aucs[i,2]=round(as.numeric(as.character(ci(ROCcurve))[[2]]),2)
  Aucs[i,3]=round(as.numeric(as.character(ci(ROCcurve))[[3]]),2)
}
names(PAMs)=Ns
save(PAMs, file="resultsdir/PAM_oncVchRCC_interative.RData")

dat=load("~/Documents/Projects/Oncocytoma/results/PAM_oncVchRCC_interative.RData")

rownames(Aucs)=Ns
colnames(Aucs)=c("95% CI (lower)","auc", "95% CI (upper)")
Aucs=as.data.frame(Aucs)
plot(rownames(Aucs),Aucs$auc)

#Make list of number of misclassifications per class
NrMisclassificationsPerClassList=list()
for(i in 1:length(PAMs)){
     tab=as.data.frame(PAMs[[i]]$NrMisclassificationsPerClass)
     rownames(tab)=c("chRCC", "oncocytoma")
     tab$NGenes=rep(names(PAMs[i]), nrow(tab))
     NrMisclassificationsPerClassList[[i]]=tab
}
names(NrMisclassificationsPerClassList)=names(PAMs)

library(abind)
NrMisclassificationsPerClassList_iterative=abind(NrMisclassificationsPerClassList, along=1)
NrMisclassificationsPerClassList_iterative=as.data.frame(NrMisclassificationsPerClassList_iterative)
NrMisclassificationsPerClassList_iterative$NrMisclassificationsPerc=round(as.numeric(as.character(NrMisclassificationsPerClassList_iterative$NrMisclassificationsPerc)),2)
NrMisclassificationsPerClassList_iterative$Accuracy=round(as.numeric(as.character(NrMisclassificationsPerClassList_iterative$Accuracy)),2)

NrMisclassificationsPerClassList_iterative$Tumor.type=gsub("[.].*","",rownames(NrMisclassificationsPerClassList_iterative))
write.table(NrMisclassificationsPerClassList_iterative, file="resultsdir/PAM_oncVchRCC_interative_misclassificationsperclass.txt", sep="\t")
#Plotting accuracies of accuracies of different PAM models

file="resultsdir/PAM_oncVchRCC_interative_misclassificationsperclass_plot"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
plot(as.numeric(as.character(NrMisclassificationsPerClassList_iterative[NrMisclassificationsPerClassList_iterative$Tumor.type=="chRCC","NGenes"])),
NrMisclassificationsPerClassList_iterative[NrMisclassificationsPerClassList_iterative$Tumor.type=="chRCC","Accuracy"], col=cbp2[3], pch=19, cex.axis=1.3, ylab="% accuracy", xlab="N genes", cex.lab=1.3)
points(as.numeric(as.character(NrMisclassificationsPerClassList_iterative[NrMisclassificationsPerClassList_iterative$Tumor.type=="oncocytoma","NGenes"])),
NrMisclassificationsPerClassList_iterative[NrMisclassificationsPerClassList_iterative$Tumor.type=="oncocytoma","Accuracy"], col=cbp2[5], pch=19)
dev.off()

file="resultsdir/PAM_oncVchRCC_interative_misclassificationsperclass_plot_legend"
png(file=paste(file,'.png',sep=''), units="in", width=6, height=6, res=600)
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("center", c("chRCC","oncocytoma"), col=c(cbp2[3],cbp2[5]), pch=19, title="Sample type", cex=2, bty="n")
dev.off()

##################################################################
#get most frequently used CpG sites for PAM onc versus chRCC
###################################################################

generankslist=list()

for(i in 1:length(PAMs)){
generanks=as.data.frame(abind(PAMs[[i]]$GenesRanks, along=1))

tab=as.data.frame(table(generanks$id)) 
rownames(tab)=tab$Var1

generanks$freq=unlist(lapply(generanks$id, function(x) tab[x,"Freq"]))
generanks=generanks[order(generanks$freq, decreasing = T),]

#get mean for proprtion selected in cross validation 
colnames(generanks)=gsub("-",".",colnames(generanks))
#tab1=as.data.frame(aggregate(as.numeric(as.character(generanks$prop.selected.in.CV))~generanks$id, FUN=mean))

tab2=as.data.frame(aggregate(as.numeric(as.character(generanks$`1.score`))~generanks$id, FUN=mean))
tab3=as.data.frame(aggregate(as.numeric(as.character(generanks$`2.score`))~generanks$id, FUN=mean))
tab4=as.data.frame(aggregate(as.numeric(as.character(generanks$av.rank.in.CV))~generanks$id, FUN=mean))
tab5=as.data.frame(aggregate(as.numeric(as.character(generanks$freq))~generanks$id, FUN=mean))
prop=unlist(lapply(levels(generanks$id), function(k) sum(as.numeric(as.character(abind(PAMs[[i]]$GenesRanks, along=1)[rownames(abind(PAMs[[i]]$GenesRanks, along=1)) %in% k,"prop-selected-in-CV"])))/10))
names(prop)=levels(generanks$id)
tab1=as.data.frame(cbind(names(prop), prop))
colnames(tab1)=colnames(tab2)
tab1[,2]=as.numeric(as.character(tab1[,2]))

gen=Reduce(function(...) merge(..., all=TRUE, by="generanks$id"), list(tab1, tab2, tab3, tab4, tab5))
gen=gen[,c(1,2,5,6,3,4)]
colnames(gen)=c("ID","prop.selected.in.CV","av.rank.in.CV", "freq","score.chRCC","score.oncocytoma")
gene=gen[order(gen$freq, gen$prop.selected.in.CV, -gen$av.rank.in.CV, decreasing=T),]
colnames(gene)[2:ncol(gene)]=gsub(" ", "", paste(colnames(gene)[2:ncol(gene)],".",names(PAMs)[i], by=""))

generankslist[[i]]=gene
}
names(generankslist)=names(PAMs)

generankslist2=Reduce(function(...) merge(..., all=TRUE, by="ID"), generankslist)

#fine to here 
#Make averages across all selected gene numbers and then order by ranking
generankslist2$av.rank.in.CV.overall.average=rowMeans(generankslist2[,grep("av.rank.in.CV", colnames(generankslist2))])
generankslist2$prop.selected.in.CV.overall.average=rowMeans(generankslist2[,grep("prop.selected.in.CV", colnames(generankslist2))])
generankslist2$freq.overall.average=rowMeans(generankslist2[,grep("freq.", colnames(generankslist2))])
generankslist2=generankslist2[order(generankslist2$av.rank.in.CV.overall.average),]
rownames(generankslist2)=generankslist2$ID
generankslist2=generankslist2[,order(colnames(generankslist2))]
generankslist2=generankslist2[,-grep("ID", colnames(generankslist2))]
generankslist2=round(generankslist2, 2)
generankslist2=cbind(ann450k2[rownames(generankslist2),], generankslist2)
generankslist2=as.matrix(generankslist2)
write.table(generankslist2, file="resultsdir/PAM_oncVchRCC_predictive_CpGs_iterative.txt", sep="\t")

##################
#make plots indicating ranks, centroids, meth difference, etc. Making it for the 20 gene PAM model, since this is the minimal optimal model
####################

genranks=read.table("resultsdir/PAM_oncVchRCC_predictive_CpGs_iterative.txt",  sep="\t", header=T,  fill=T, quote="", comment.char="")
genranks=genranks[order(genranks$freq.20., -genranks$av.rank.in.CV.20., decreasing=T),]

colnames(genranks)=gsub("X.","",colnames(genranks))
library(stringr)
rownames(genranks)=str_replace_all(rownames(genranks), '\"', "")

SAM_genes=read.table('resultsdir/SAM_analysis_OncVChRCC.txt', sep="\t", header=T)
SAM_genes=SAM_genes[SAM_genes[,4]=="0" & abs(SAM_genes[,3])>0.2,]
rownames(SAM_genes)=SAM_genes$Gene.Name

#add differential methylation 
genranks$direction=rep(NA, nrow(genranks))
genranks$direction=SAM_genes[rownames(genranks),"Fold.Change"] 
genranks$dircol=ifelse(genranks$direction>0,"hyper","hypo")
genranks$Gene.=as.character(genranks$Gene.)

genranks$probe2=as.factor(ifelse(nchar(as.character(genranks$Gene.))>0,paste(rownames(genranks),"-", genranks$Gene., sep=""), rownames(genranks)))

#########################
########################
library(forcats)

dat1=as.data.frame(cbind(as.character(genranks$probe2), genranks$av.rank.in.CV.20., genranks$dircol))
colnames(dat1)=c("probe","cent","col")
dat1$cent=as.numeric(as.character(dat1$cent))
dat1$probe <- factor(dat1$probe, levels = dat1$probe)
dat1$factor=rep("mean.rank.in.CV", nrow(dat1))

dat2=as.data.frame(cbind(as.character(genranks$probe2), genranks$freq.20., genranks$dircol))
colnames(dat2)=c("probe","cent","col")
dat2$cent=as.numeric(as.character(dat2$cent))
dat2$probe <- factor(dat2$probe, levels = dat2$probe) 
dat2$factor=rep("shrinkage.survival.index", nrow(dat2))

dat3=as.data.frame(cbind(as.character(genranks$probe2), genranks$score.chRCC.20., genranks$dircol))
colnames(dat3)=c("probe","cent","col")
dat3$cent=as.numeric(as.character(dat3$cent))
dat3$probe <- factor(dat3$probe, levels = dat3$probe) 
dat3$factor=rep("centroid.chRCC", nrow(dat3))

dat4=as.data.frame(cbind(as.character(genranks$probe2), genranks$score.oncocytoma.20., genranks$dircol))
colnames(dat4)=c("probe","cent","col")
dat4$cent=as.numeric(as.character(dat4$cent))
dat4$probe <- factor(dat4$probe, levels = dat4$probe) 
dat4$factor=rep("centroid.oncocytoma", nrow(dat4))

dat5=as.data.frame(cbind(as.character(genranks$probe2), genranks$direction, genranks$dircol))
colnames(dat5)=c("probe","cent","col")
dat5$cent=as.numeric(as.character(dat5$cent))
dat5$probe <- factor(dat5$probe, levels = dat5$probe) 
dat5$factor=rep("direction", nrow(dat5))

all.dat=rbind(dat1, dat2, dat3, dat4, dat5)
dat=all.dat
colnames(dat)[4]="Parameter"
colnames(dat)[3]="Direction"
dat$Parameter=gsub("direction","beta.value.difference", dat$Parameter)

dat$Parameter = factor(dat$Parameter, levels=c("shrinkage.survival.index","mean.rank.in.CV","beta.value.difference","centroid.chRCC","centroid.oncocytoma"))
colnames(dat)[4]="Metric"

file="resultsdir/roc_OncVsChRCC_PAM20_genes_ranks"
png(file=paste(file,'.png',sep=''), units="in", width=15, height=9, res=600)
ggplot(dat) + aes(x = cent, y = fct_rev(probe), color = Direction) + geom_point()+ geom_segment(aes(x = 0, y = probe, xend = cent, yend = probe), color = "grey50") + scale_fill_manual(values = c("hypo"="blue", "hyper"="red")) + facet_grid(. ~ Metric, labeller="label_both", scales="free") + theme(strip.text.x = element_text(size = 10, colour = "black", face="bold")) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(legend.key.size = unit(1.5, "cm"))
dev.off()

#########################
#Make data frame of this
################################

genranks=read.table("resultsdir/PAM_oncVchRCC_predictive_CpGs_iterative.txt",  sep="\t", header=T,  fill=T, quote="", comment.char="")

#genranks=read.table("resultsdir/PAM_oncVchRCC_predictive_CpGs_iterative.txt",  sep="\t", header=T)
colnames(genranks)[2:50]=colnames(genranks)[1:49]
colnames(genranks)[1]="probe"
rownames(genranks)=as.character(genranks$probe)

genranks=genranks[order(genranks$freq.20., -genranks$av.rank.in.CV.20., decreasing=T),]

SAM_genes=read.table('resultsdir/SAM_analysis_OncVChRCC.txt', sep="\t", header=T)
SAM_genes=SAM_genes[SAM_genes[,4]=="0" & abs(SAM_genes[,3])>0.2,]
rownames(SAM_genes)=SAM_genes$Gene.Name

#add differential methylation 
genranks$direction=rep(NA, nrow(genranks))
genranks$direction=SAM_genes[rownames(genranks),"Fold.Change"] 
genranks$dircol=ifelse(genranks$direction>0,"hyper","hypo")
genranks$Gene.=as.character(genranks$Gene.)

dat=cbind(as.character(genranks$probe), as.character(genranks$chr.),as.character(genranks$pos.), as.character(genranks$Gene.),as.character(genranks$region), as.character(genranks$UCSC_RefGene_Group.),as.character(genranks$Relation_to_Island.),genranks$av.rank.in.CV.20., genranks$freq.20., genranks$score.chRCC.20., genranks$score.oncocytoma.20., genranks$direction)
colnames(dat)=c("probe","chr","position_hg19","gene","region","UCSC_RefGene_Group","Relation_to_Island","mean.rank.in.CV","shrinkage.survival.index","centroid.chRCC","centroid.oncocytoma","beta.value.difference")
dat=as.data.frame(dat)


#add sequence 
library(BiocManager)
#install(pkgs = "BSgenome.Hsapiens.UCSC.hg19",update = TRUE, ask = TRUE)
library(BSgenome.Hsapiens.UCSC.hg19)

probes=as.character(dat$probe)
get3=getpyroseqs(probes)

dat2=cbind(dat, get3)

write.table(dat2, file="resultsdir/assay_design_info.txt", sep="\t")

###################################################
#use bedtools merge to find CpG sites among the top ranking set that are within 10kb of each other 
####################################################

genranks=read.table("resultsdir/PAM_oncVchRCC_predictive_CpGs_iterative.txt",  sep="\t", header=T,  fill=T, quote="", comment.char="")
colnames(genranks)=gsub("X.","",colnames(genranks))
genranks=genranks[order(genranks$freq.20., -genranks$av.rank.in.CV.20., decreasing=T),]
library(stringr)
rownames(genranks)=str_replace_all(rownames(genranks), '\"', "")

#This depreciates the predictive value of the model
#same number of KICHs in my data as for Chopra study 

#Need to figure out which CpG sites represent seperate regions, to see if selecting unique regions does better 
#find R code for finding the distance of each CpG site to each other CpG site 
#bedtools window
#make bedfiles for CpG sites and find any that are within 5kb of that CpG site
bed=cbind(as.character(genranks$chr.), as.numeric(as.character(genranks$pos.)), as.numeric(as.character(genranks$pos.)), rownames(genranks))
write.table(bed, file="resultsdir/chRCC_oncocytoma_predictive_loci_ranked.bed", sep="\t", row.names=F, col.names=F, quote=F)

#sort bed file
dir="resultsdir/"
inputfile="chRCC_oncocytoma_predictive_loci_ranked.bed"
outputfile="chRCC_oncocytoma_predictive_loci_ranked_sort.bed"

command=paste("/Applications/bedtools2/bin/bedtools sort -i ", paste(dir, inputfile, sep="")," > ", paste(dir, outputfile, sep=""))
system(command)

#merging CpG probes that are within 5kb of another CpG probe. output is ordered by chromosome and positio
inputfile="chRCC_oncocytoma_predictive_loci_ranked_sort.bed"
command=paste("/Applications/bedtools2/bin/bedtools merge -i ", paste(dir, inputfile, sep=""), " -c 4 -o collapse -d 5000", " > ", paste(dir, "chRCC_oncocytoma_predictive_loci_ranked.merged.bed", sep=""), sep="")

command=paste("/Applications/bedtools2/bin/bedtools merge -i ", paste(dir, inputfile, sep=""), " -c 4 -o collapse -d 10000", " > ", paste(dir, "chRCC_oncocytoma_predictive_loci_ranked.merged.bed", sep=""), sep="")
system(command)

window=read.table(paste(dir, "chRCC_oncocytoma_predictive_loci_ranked.merged.bed", sep=""), sep="\t", header=F)
colnames(window)=c("chr","start","end","probes")
window$probes=as.character(window$probes)

lapply(window$probes, function(x) strsplit(x, ","))
windowlist=strsplit(as.character(window$probes),",")
names(windowlist)=paste(window$chr,":", window$start,"-", window$end, sep="")
library(abind)

genranks$region=rep(NA, nrow(genranks))
for(i in 1:nrow(genranks)){
y=rownames(genranks)[i]
genranks$region[i]=names(windowlist[unlist(lapply(windowlist, function(x) y %in% x))])
}
#74 unique regions
dupregs=genranks$region[duplicated(genranks$region)]
genranks[genranks$region %in% dupregs, c("Gene.","region")]
genranks=genranks[order(genranks$freq.20., -genranks$av.rank.in.CV.20., decreasing=T),]

#writing genranks to file so I have the region info
write.table(genranks, file="resultsdir/PAM_oncVchRCC_predictive_CpGs_iterative.txt",  sep="\t")
ord=read.table("datadir/order_genranks_temp.txt", sep="\t")

tab=cbind(rownames(genranks),as.character(genranks[match(ord$V1, rownames(genranks)),"Gene."]), as.character(genranks[match(ord$V1, rownames(genranks)),"region"]))
write.table(tab, file="datadir/genes_order_genranks_temp.txt", sep="\t", quote=F)



















