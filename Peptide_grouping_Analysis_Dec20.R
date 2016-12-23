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
peplist2 = peplist[, list(subtype = paste(unique(subtype), collapse = "/")), by = list(PEPTIDE_SEQUENCE, flu_prot, Pos)]
#read in another text file containing all the probe IDS
all_probes = read_delim('~/Desktop/ProteinIDS_seq.txt', delim = ";", col_names = F)
colnames(all_probes) = c("PROTID", "FULLNAME", "SEQUENCE")
flu_main=all_probes[grep("Influenza", all_probes$FULLNAME),]
require(stringr)
flu_list=split.data.frame(flu_main, flu_main$FULLNAME)
flu1=as_data_frame(str_split_fixed(flu_main$FULLNAME, " ", 6))

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
                                                                                            

flu_main_unique=flu_main[, list(PROTID = paste(unique(PROTID), collapse = "/"), FULLNAME = paste(unique(FULLNAME), collapse = "/")), by=list(SEQUENCE, subtype)]
final_grouping=merge.data.frame(peplist2, flu_main_unique, by.x ="PEPTIDE_SEQUENCE", by.y = "SEQUENCE", all= T)

flu_main1=colsplit(flu_main$FULLNAME, " ", names = 1:6)
flu_main2=unique(flu_main1)
