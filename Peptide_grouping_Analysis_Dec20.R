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
#read in another text file containing all the probe IDS
all_probes = read_delim('~/Desktop/ProteinIDS_seq.txt', delim = ";", col_names = F)
colnames(all_probes) = c("PROTID", "FULLNAME", "SEQUENCE")
flu_main=all_probes[grep("Influenza", all_probes$FULLNAME),]
flu_main1=colsplit(flu_main$FULLNAME, " ", names = 1:6)
flu_main2=unique(flu_main1)
