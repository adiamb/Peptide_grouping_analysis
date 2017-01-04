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
###peptides present only in pH1N1
pH1N1=total_IgG[colnames(total_IgG) %in% grouping_list$pH1N1$SEQUENCE]
pH1N1 = cbind(total_IgG[c(1:7, 142030)], pH1N1)
###peptide present only in gH1N1
gH1N1=total_IgG[colnames(total_IgG) %in% grouping_list$gH1N1$SEQUENCE]
gH1N1 = cbind(total_IgG[c(1:7, 142030)], gH1N1)
###peptide present only in H3N2
H3N2=total_IgG[colnames(total_IgG) %in% grouping_list$H3N2$SEQUENCE]
H3N2 = cbind(total_IgG[c(1:7, 142030)], H3N2)
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

















