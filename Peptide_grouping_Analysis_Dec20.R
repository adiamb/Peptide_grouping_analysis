setwd('~/Dropbox/PEPARRAY-75_selected_peptides/')
list.files()
require(readr)
peplist = read_delim('~/Dropbox/PEPARRAY-75_selected_peptides/PEPARRAY-75_provided_peptides.txt',  delim ="\t")
peplist$CYCLES = NULL
peplist$PROBE_SEQUENCE = NULL
require(reshape2)
peplist=transform(peplist, PROBE_ID = colsplit(PROBE_ID, "_", names = c("Pos", "mer", "flu_prot", "subtype", "variant")))
peplist= cbind(peplist, peplist$PROBE_ID)
peplist
peplist$PROBE_ID = NULL
require(data.table)
setDT(peplist)
#read in another text file containing all the probe IDS
all_probes = read_delim('~/Desktop/ProteinIDS_seq.txt', delim = ";", col_names = F)
colnames(all_probes) = c("PROTID", "FULLNAME", "SEQUENCE")
flu_main=all_probes[grep("Influenza", all_probes$FULLNAME),]
require(stringr)

flu_main$subtype=ifelse(str_detect(flu_main$FULLNAME, "California"), "pH1N1",
       ifelse(str_detect(flu_main$FULLNAME, "Solomon"), "gH1N1", 
              ifelse(str_detect(flu_main$FULLNAME, "Puerto|South"), "pastpH1N1",
                     "H3N2")))



flu_main$Protname = ifelse(str_detect(flu_main$FULLNAME, "Polymerase basic protein 2|polymerase PB2"), "PB2",
                           ifelse(str_detect(flu_main$FULLNAME, "Polymerase acidic protein|polymerase PA"), "PA",
                                  ifelse(str_detect(flu_main$FULLNAME, "RNA-directed RNA polymerase catalytic subunit|polymerase PB1"), "PB1",
                                         ifelse(str_detect(flu_main$FULLNAME, "Hemagglutinin|hemagglutinin|HA|haemagglutinin"), "HA",
                                                ifelse(str_detect(flu_main$FULLNAME, "Neuraminidase|neuraminidase|NA"), "Neu",
                                                       ifelse(str_detect(flu_main$FULLNAME, "Matrix protein 2|matrix protein 2"), "M2",
                                                              ifelse(str_detect(flu_main$FULLNAME, "Matrix protein 1|matrix protein 1"),"M1",
                                                                     ifelse(str_detect(flu_main$FULLNAME, "Nuclear export protein|nuclear export protein"),"NEP",
                                                                            ifelse(str_detect(flu_main$FULLNAME, "Non-structural protein 1|non-structural protein 1|nonstructural protein 1"), "NS1",
                                                                                   ifelse(str_detect(flu_main$FULLNAME, "Nucleoprotein|nucleocapsid protein|nucleoprotein"), "NP",
                                                                                          ifelse(str_detect(flu_main$FULLNAME, "PA-X protein|PB1-F2"), "otherPA_PB",
                                                                                                 ifelse(str_detect(flu_main$FULLNAME, "nonstructural protein 2"), "NS2",
                                                                                                 NA))))))))))))
                                                                                            






##select only peptide sequence and the protein names and the subtype from the flu_main data frames 
flu_main1 = select(flu_main, SEQUENCE, subtype, Protname)
peplist2 = select(peplist, PEPTIDE_SEQUENCE, subtype, flu_prot)
names(peplist2) = names(flu_main1) ##match the colnames from each of the data frame
final_pep_df=rbind.data.frame(peplist2, flu_main1)
##aggregate to get common sequences between strains
require(data.table)
setDT(final_pep_df)
final_grouping = final_pep_df[, list(subtype = paste(unique(subtype), collapse = "/")), by = list(SEQUENCE, Protname)]
grouping_list = split.data.frame(final_grouping, final_grouping$subtype)



##read in the peptide array results - using the qunatile normalized data from roche 
total_IgG = fread('~/Documents/IgG_Analysis/normalized_qunatile/final_aggregated_data_Dec12.csv', showProgress = T)
total_IgG = as_data_frame(total_IgG)# convert to normal dataframe
total_IgG$Diagnosis = ifelse(total_IgG$Dx == "K", 0, 1) #add diagnosis in the last column
##check the dims of each dataframe in the list
lapply(grouping_list, FUN = function(x) (dim(x)))
##for loop to match and separate the IgG responses to each peptide grouping 
pep_group = list()
for (i in names(grouping_list)){
  temp = total_IgG[colnames(total_IgG) %in% grouping_list[[i]]$SEQUENCE]
  temp2 = cbind(total_IgG[c(1:7, 142030)], temp)
  name = paste(i)
  pep_group[[name]] = temp2
}


###peptides present only in pH1N1
pH1N1=pep_group$pH1N1
###peptide present only in gH1N1
gH1N1=pep_group$gH1N1
###peptide present only in H3N2
H3N2=pep_group$H3N2
#################################################################################TOTAL SAMPLES ANALYSIS#########################################################################################
##SAMR analysis pH1N1
require(samr)
require(Biobase)
require(limma)
y = ifelse(total_IgG$Dx == "K", 1, 2)
x = log1p(t(pH1N1[9:ncol(pH1N1)]))

d <- list(x=x, y=y, geneid=row.names(x), logged2 = T)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=1000, assay.type = "array")
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.12
samr.plot(samr.obj,delta)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
require(ggplot2)
ggplot(pH1N1, aes(factor(y), ENVETMNSNTLELRSR, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
grouping_list$pH1N1[grep("REILTKTTVDHMSII", grouping_list$pH1N1$SEQUENCE),]

##limma analysis pH1N1
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimensions are right !
design = model.matrix(~Diagnosis+Gender, data = pH1N1) 
fit = eBayes(lmFit(object, design))
topTable(fit, coef="Diagnosis", adjust="BH", number = 15)

###SAMR analysis with gH1N1

y = ifelse(total_IgG$Dx == "K", 1, 2)
x = log1p(t(gH1N1[9:ncol(gH1N1)]))

d <- list(x=x, y=y, geneid=row.names(x), logged2 = T)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=1000, assay.type = "array")
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.2
samr.plot(samr.obj,delta)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
require(ggplot2)
ggplot(gH1N1, aes(factor(y), EAMEVANQARQMVQAM, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
grouping_list$gH1N1[grep("EAMEVANQARQMVQAM", grouping_list$gH1N1$SEQUENCE),]

##limma analysis gH1N1
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimensions are right !
design = model.matrix(~Diagnosis+Gender, data = gH1N1) 
fit = eBayes(lmFit(object, design))
topTable(fit, coef="Diagnosis", adjust="BH", number = 15)

###SAMR analysis with H3N3

y = ifelse(total_IgG$Dx == "K", 1, 2)
x = log1p(t(H3N2[9:ncol(H3N2)]))

d <- list(x=x, y=y, geneid=row.names(x), logged2 = T)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=1000, assay.type = "array")
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.24
samr.plot(samr.obj,delta)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
require(ggplot2)
ggplot(H3N2, aes(factor(y), IRNGTYDPDVYRDEAL, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
grouping_list$H3N2[grep("IRNGTYDPDVYRDEAL", grouping_list$H3N2$SEQUENCE),]

##limma analysis gH1N1
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimensions are right !
design = model.matrix(~Diagnosis+Gender, data = H3N2) 
fit = eBayes(lmFit(object, design))
topTable(fit, coef="Diagnosis", adjust="BH", number = 15)
#################################################################################################split analysis - recent onset and pandemrix ########################################################
##split samples into recent onset and pandemrix
require(dplyr)
pH1N1_recent = filter(pH1N1, sample_group == "early onset")
pH1N1_PX = filter(pH1N1, sample_group == "pandemrix")
##recent onset SAMR pH1N1
y = ifelse(pH1N1_recent$Dx == "K", 1, 2)
x = log(t(pH1N1_recent[9:ncol(pH1N1_recent)]))

d <- list(x=x, y=y, geneid=row.names(x), logged2 = T)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=1000, assay.type = "array")
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.2
samr.plot(samr.obj,delta)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
require(ggplot2)
ggplot(pH1N1, aes(factor(y), ENVETMNSNTLELRSR, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
grouping_list$pH1N1[grep("REILTKTTVDHMSII", grouping_list$pH1N1$SEQUENCE),]

##limma analysis pH1N1
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimensions are right !
design = model.matrix(~Diagnosis+Gender, data = pH1N1_recent) 
fit = eBayes(lmFit(object, design))
topTable(fit, coef="Diagnosis", adjust="BH", number = 15)

##pandemrix pH1N1
y = ifelse(pH1N1_PX$Dx == "K", 1, 2)
x = log1p(t(pH1N1_PX[9:ncol(pH1N1_PX)]))

d <- list(x=x, y=y, geneid=row.names(x), logged2 = T)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=1000, assay.type = "array")
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.2
samr.plot(samr.obj,delta)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
require(ggplot2)
ggplot(pH1N1, aes(factor(Diagnosis), ESVFTGKNTDLEALME, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
grouping_list$pH1N1[grep("NFAAGQSVVSAKLAGN", grouping_list$pH1N1$SEQUENCE),]

##limma analysis pH1N1
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimensions are right !
design = model.matrix(~Diagnosis+Gender, data = pH1N1_PX) 
fit = eBayes(lmFit(object, design))
topTable(fit, coef="Diagnosis", adjust="BH", number = 15)
#################################################################################################split samples gH1N1#############################################################
require(dplyr)
gH1N1_recent = filter(gH1N1, sample_group == "early onset")
gH1N1_PX = filter(gH1N1, sample_group == "pandemrix")
##recent onset SAMR gH1N1
y = ifelse(gH1N1_recent$Dx == "K", 1, 2)
x = log(t(gH1N1_recent[9:ncol(gH1N1_recent)]))

d <- list(x=x, y=y, geneid=row.names(x), logged2 = T)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=1000, assay.type = "array")
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.2
samr.plot(samr.obj,delta)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
require(ggplot2)
ggplot(gH1N1, aes(factor(Diagnosis), GRWTKNTETGAPQLNP, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
grouping_list$gH1N1[grep("GRWTKNTETGAPQLNP", grouping_list$gH1N1$SEQUENCE),]

##limma analysis gH1N1
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimensions are right !
design = model.matrix(~Diagnosis+Gender, data = gH1N1_recent) 
fit = eBayes(lmFit(object, design))
topTable(fit, coef="Diagnosis", adjust="BH", number = 15)

##pandemrix gH1N1
y = ifelse(gH1N1_PX$Dx == "K", 1, 2)
x = log1p(t(gH1N1_PX[9:ncol(gH1N1_PX)]))

d <- list(x=x, y=y, geneid=row.names(x), logged2 = T)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=1000, assay.type = "array")
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.2
samr.plot(samr.obj,delta)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
require(ggplot2)
ggplot(gH1N1, aes(factor(Diagnosis), ESVFTGKNTDLEALME, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
grouping_list$gH1N1[grep("NFAAGQSVVSAKLAGN", grouping_list$gH1N1$SEQUENCE),]

##limma analysis gH1N1
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimensions are right !
design = model.matrix(~Diagnosis+Gender, data = gH1N1_PX) 
fit = eBayes(lmFit(object, design))
topTable(fit, coef="Diagnosis", adjust="BH", number = 15)
#################################################################################################split samples H3N2#############################################################
require(dplyr)
H3N2_recent = filter(H3N2, sample_group == "early onset")
H3N2_PX = filter(H3N2, sample_group == "pandemrix")
##recent onset SAMR H3N2
y = ifelse(H3N2_recent$Dx == "K", 1, 2)
x = log(t(H3N2_recent[9:ncol(H3N2_recent)]))

d <- list(x=x, y=y, geneid=row.names(x), logged2 = T)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=1000, assay.type = "array")
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.2
samr.plot(samr.obj,delta)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
require(ggplot2)
ggplot(H3N2, aes(factor(Diagnosis), GRWTKNTETGAPQLNP, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
grouping_list$H3N2[grep("GRWTKNTETGAPQLNP", grouping_list$H3N2$SEQUENCE),]

##limma analysis H3N2
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimensions are right !
design = model.matrix(~Diagnosis+Gender, data = H3N2_recent) 
fit = eBayes(lmFit(object, design))
topTable(fit, coef="Diagnosis", adjust="BH", number = 15)

##pandemrix H3N2
y = ifelse(H3N2_PX$Dx == "K", 1, 2)
x = log1p(t(H3N2_PX[9:ncol(H3N2_PX)]))

d <- list(x=x, y=y, geneid=row.names(x), logged2 = T)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=1000, assay.type = "array")
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.17
samr.plot(samr.obj,delta)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
require(ggplot2)
ggplot(H3N2, aes(factor(Diagnosis), ESVFTGKNTDLEALME, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
grouping_list$H3N2[grep("NFAAGQSVVSAKLAGN", grouping_list$H3N2$SEQUENCE),]

##limma analysis H3N2
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimensions are right !
design = model.matrix(~Diagnosis+Gender, data = H3N2_PX) 
fit = eBayes(lmFit(object, design))
topTable(fit, coef="Diagnosis", adjust="BH", number = 15)



#####$shared flu epitopes ##### analysis
shared_flu = pep_group$`gH1N1/pH1N1/H3N2/pastpH1N1`

require(dplyr)
shared_flu_recent = filter(shared_flu, sample_group == "early onset")
shared_flu_PX = filter(shared_flu, sample_group == "pandemrix")
##recent onset SAMR shared_flu
y = ifelse(shared_flu_recent$Dx == "K", 1, 2)
x = log(t(shared_flu_recent[9:ncol(shared_flu_recent)]))

d <- list(x=x, y=y, geneid=row.names(x), logged2 = T)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=1000, assay.type = "array")
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.18
samr.plot(samr.obj,delta)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
require(ggplot2)
ggplot(shared_flu_recent, aes(factor(Diagnosis), DNTKWNENQNPRMFLA, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
grouping_list$`gH1N1/pH1N1/H3N2/pastpH1N1`[grep("DNTKWNENQNPRMFLA", grouping_list$`gH1N1/pH1N1/H3N2/pastpH1N1`$SEQUENCE),]

##limma analysis shared_flu
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimensions are right !
design = model.matrix(~Diagnosis+Gender, data = shared_flu_recent) 
fit = eBayes(lmFit(object, design))
topTable(fit, coef="Diagnosis", adjust="BH", number = 15)

##pandemrix shared_flu
y = ifelse(shared_flu_PX$Dx == "K", 1, 2)
x = log1p(t(shared_flu_PX[9:ncol(shared_flu_PX)]))

d <- list(x=x, y=y, geneid=row.names(x), logged2 = T)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=1000, assay.type = "array")
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.25
samr.plot(samr.obj,delta)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
require(ggplot2)
ggplot(shared_flu_PX, aes(factor(Diagnosis), SADPLASLLEMCHSTQ, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
grouping_list$shared_flu[grep("NFAAGQSVVSAKLAGN", grouping_list$shared_flu$SEQUENCE),]

##limma analysis shared_flu
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimensions are right !
design = model.matrix(~Diagnosis+Gender, data = shared_flu_PX) 
fit = eBayes(lmFit(object, design))
topTable(fit, coef="Diagnosis", adjust="BH", number = 15)
#######################################################################################sharedH1#####################################################################################
#####$shared flu epitopes ##### analysis
sharedpH1_flu = pep_group$`pH1N1/pastpH1N1`
require(samr)
require(dplyr)
sharedpH1_flu_recent = filter(sharedpH1_flu, sample_group == "early onset")
sharedpH1_flu_PX = filter(sharedpH1_flu, sample_group == "pandemrix")
##recent onset SAMR sharedpH1_flu
y = ifelse(sharedpH1_flu_recent$Dx == "K", 1, 2)
x = log(t(sharedpH1_flu_recent[9:ncol(sharedpH1_flu_recent)]))

d <- list(x=x, y=y, geneid=row.names(x), logged2 = T)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=1000, assay.type = "array")
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.2
samr.plot(samr.obj,delta)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
require(ggplot2)
ggplot(sharedpH1_flu_recent, aes(factor(Diagnosis), RGVQIASNENMETMES, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
grouping_list$`pH1N1/pastpH1N1`[grep("RGVQIASNENMETMES", grouping_list$`pH1N1/pastpH1N1`$SEQUENCE),]

##limma analysis sharedpH1_flu
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimensions are right !
design = model.matrix(~Diagnosis+Gender, data = sharedpH1_flu_recent) 
fit = eBayes(lmFit(object, design))
topTable(fit, coef="Diagnosis1", adjust="BH", number = 30)

##pandemrix sharedpH1_flu
y = ifelse(sharedpH1_flu_PX$Dx == "K", 0, 1)
x = log1p(t(sharedpH1_flu_PX[9:ncol(sharedpH1_flu_PX)]))

d <- list(x=x, y=y, geneid=row.names(x), logged2 = T)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=1000, assay.type = "array")
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.25
samr.plot(samr.obj,delta)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
require(ggplot2)
ggplot(sharedpH1_flu_PX, aes(factor(Diagnosis), MATKADYTLDEESRAR, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
grouping_list$sharedpH1_flu[grep("NFAAGQSVVSAKLAGN", grouping_list$sharedpH1_flu$SEQUENCE),]

##limma analysis sharedpH1_flu
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimensions are right !
design = model.matrix(~Diagnosis+Gender, data = sharedpH1_flu_PX) 
fit = eBayes(lmFit(object, design))
topTable(fit, coef="Diagnosis", adjust="BH", number = 15)

##glmnet
require(glmnet)
y=shared_flu_PX$Diagnosis
x=as.matrix(scale(shared_flu_PX[9:ncol(shared_flu_PX)]))
fit1 = cv.glmnet(x=x, y=y, family ="binomial", type.measure = "class", nfolds = 3, alpha = 0.8)
### check lambda and lambda.min
plot(fit1)
fit1$lambda
fit1$lambda.min
plot(fit1)
fit1$glmnet.fit$beta[which(fit1$glmnet.fit$beta[,21]!=0),21] %>% sort() 
