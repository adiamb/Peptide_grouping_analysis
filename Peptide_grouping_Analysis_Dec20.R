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
total_IgG = as_data_frame(total_IgG) # convert to normal dataframe
###peptide present only in pH1N1
pH1N1=total_IgG[colnames(total_IgG) %in% grouping_list$pH1N1$SEQUENCE]
pH1N1 = cbind(total_IgG[1:7], pH1N1)
###peptide present only in gH1N1
gH1N1=total_IgG[colnames(total_IgG) %in% grouping_list$gH1N1$SEQUENCE]
gH1N1 = cbind(total_IgG[1:7], gH1N1)
###peptide present only in H3N2
H3N2=total_IgG[colnames(total_IgG) %in% grouping_list$H3N2$SEQUENCE]
H3N2 = cbind(total_IgG[1:7], H3N2)













