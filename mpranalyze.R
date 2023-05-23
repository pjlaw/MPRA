
library(tidyverse)
library(MPRAnalyze)
library(glue)


cell_line="HCEC-1CT" #SW403 "HCEC-1CT" #HT29"
setwd(paste0("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/mpraflow/", cell_line))

#readin count matices, replace NA with 0
# dna_dat=read_tsv("dna_counts.tsv")
# dna_dat0 = dna_dat %>% replace(is.na(.),0)
# rna_dat=read_tsv("rna_counts.tsv")
# rna_dat0 = rna_dat %>% replace(is.na(.),0)
# 
# 
# #pivot the matrices so column
# dna_dat0 = dna_dat0 %>% separate(seq_id, into=c("rsid", "allele", "dir"))
# #pivot tables to group by rsid
# dna_mat = dna_dat0 %>% pivot_wider(names_from=c(allele, dir), values_from=c(-rsid, -allele, -dir))
# 
# rna_dat0 = rna_dat0 %>% separate(seq_id, into=c("rsid", "allele", "dir"))
# rna_mat = rna_dat0 %>% pivot_wider(names_from=c(allele, dir), values_from=c(-rsid, -allele, -dir))
# 
# #TODO - filter rsid?
# 
# #check that all rownames are the same for RNA and DNA matrix
# all(dna_mat$rsid ==rna_mat$rsid)
# 
# 
# variant_annot=read_tsv("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/CRC_MPRA_design_library_all_variants_v3_newgwas_with_controls")
# control_snps = variant_annot %>% filter(test_control=="control") %>% pull(rsid)
# #snp_data = dna_dat0[,"seq_id"] %>% separate(seq_id, into=c("rsid", "allele", "dir"))
# snp_data = dna_mat$rsid
# 
# #correct for snps missing rsid - should only be one (rs113101564)
# missing_names=which(snp_data=="")
# snp_data[missing_names]=paste0("missing_", missing_names)
# 
# is_control = snp_data %in% control_snps
# 
# # drop name column, convert to matrix (required by MPRAnalyze), readd rsid as rownames
# dna_mat=as.matrix(dna_mat[,-1])
# rna_mat=as.matrix(rna_mat[,-1])
# rownames(dna_mat)=snp_data
# rownames(rna_mat)=snp_data
# 
# #make annotation dataframes from matrix column names
# dna_annot = tibble(sample=colnames(dna_mat))
# dna_annot = dna_annot %>% separate(sample, into=c("type", "cell_line", "replicate", "barcode", "allele", "dir"), remove=F, sep="_")
# dna_annot = dna_annot %>% mutate(barcode_full=glue("{barcode}_{allele}_{dir}"))
# rna_annot = tibble(sample=colnames(rna_mat))
# rna_annot = rna_annot %>% separate(sample, into=c("type", "cell_line", "replicate", "barcode", "allele", "dir"), remove=F, sep="_")
# rna_annot = rna_annot %>% mutate(barcode_full=glue("{barcode}_{allele}_{dir}"))
# 
# #make columns into factors
# dna_annot = dna_annot %>% mutate(across(everything(), as.factor))
# dna_annot = as.data.frame(dna_annot)
# rownames(dna_annot) = dna_annot$sample
# 
# rna_annot = rna_annot %>% mutate(across(everything(), as.factor))
# rna_annot = as.data.frame(rna_annot)
# rownames(rna_annot) = rna_annot$sample
# 
# #create MpraObject
# obj = MpraObject(dnaCounts=dna_mat, rnaCounts=rna_mat, dnaAnnot=dna_annot, rnaAnnot=rna_annot, controls=is_control)
# obj <- estimateDepthFactors(obj, lib.factor = "replicate", which.lib = "both", depth.estimator = "uq")
# saveRDS(obj, paste0(cell_line,"_mpranalyze_obj.rds"))

######

obj = readRDS(paste0(cell_line,"_mpranalyze_obj.rds"))
obj <- estimateDepthFactors(obj, lib.factor = c("replicate","allele","dir"), which.lib = "both", depth.estimator = "uq")

library(BiocParallel)
#obj = analyzeComparative(obj, dnaDesign = ~ barcode_full+replicate+allele+dir, rnaDesign = ~replicate+allele+dir, reducedDesign = ~replicate, BPPARAM=MulticoreParam())
obj = analyzeComparative(obj, dnaDesign = ~ barcode_full+replicate+allele+dir, rnaDesign = ~replicate+allele+dir, reducedDesign = ~replicate+dir, BPPARAM=MulticoreParam(), mode="scale")
res <- testLrt(obj)
write_tsv(res %>% rownames_to_column(), paste0(cell_line,"_mpranalyze_res_all_scaled2.txt"))

########
library(tidyverse)
library(MPRAnalyze)
library(glue)
library(BiocParallel)

args = commandArgs(trailingOnly=TRUE)
cell_line=args[1]
max_barcodes=as.numeric(args[2])

#cell_line="HCEC-1CT" #"HT29"
#TODO determine this programmatically, eg 80% of data has less than x barcodes
#max_barcodes=30 #50#100

setwd(paste0("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/mpraflow/", cell_line))

obj = readRDS(paste0(cell_line,"_mpranalyze_obj.rds"))

dna_mat = dnaCounts(obj)
dna_mat[is.na(dna_mat)]=0
#rename the rsids that are missing info (how to tell them apart?)
missing_names=which(rownames(dna_mat)=="")
rownames(dna_mat)[missing_names]=paste0("missing_", missing_names)

dna_annot = dnaAnnot(obj)

filtered_annot = dna_annot %>% filter(as.numeric(barcode)<=max_barcodes)
wanted_observations_dna = filtered_annot %>% pull(sample) %>% as.character()

filtered_dna_mat = dna_mat[,wanted_observations_dna]

snp_zeros = rowSums(filtered_dna_mat==0)
n_snps = length(snp_zeros)

#snp has reads in at least half the observations
filtered_dna_mat =filtered_dna_mat[snp_zeros<=(n_snps/2),]

#recreate annotation dataframes from matrix column names
dna_annot = tibble(sample=wanted_observations_dna)
dna_annot = dna_annot %>% separate(sample, into=c("type", "cell_line", "replicate", "barcode", "allele", "dir"), remove=F, sep="_")
dna_annot = dna_annot %>% mutate(barcode_full=glue("{barcode}_{allele}_{dir}"))

#do the same for RNA
rna_mat = rnaCounts(obj)
rna_mat[is.na(rna_mat)]=0
missing_names=which(rownames(rna_mat)=="")
rownames(rna_mat)[missing_names]=paste0("missing_", missing_names)

rna_annot = rnaAnnot(obj)

filtered_annot = rna_annot %>% filter(as.numeric(barcode)<=max_barcodes)
wanted_observations_rna = filtered_annot %>% pull(sample) %>% as.character()

filtered_rna_mat = rna_mat[,wanted_observations_rna]

snp_zeros = rowSums(filtered_rna_mat==0)
n_snps = length(snp_zeros)

#snp has reads in at least half the observations
filtered_rna_mat =filtered_rna_mat[snp_zeros<=(n_snps/2),]

rna_annot = tibble(sample=wanted_observations_rna)
rna_annot = rna_annot %>% separate(sample, into=c("type", "cell_line", "replicate", "barcode", "allele", "dir"), remove=F, sep="_")
rna_annot = rna_annot %>% mutate(barcode_full=glue("{barcode}_{allele}_{dir}"))

#make columns into factors
dna_annot = dna_annot %>% mutate(across(everything(), as.factor))
dna_annot = as.data.frame(dna_annot)
rownames(dna_annot) = dna_annot$sample

rna_annot = rna_annot %>% mutate(across(everything(), as.factor))
rna_annot = as.data.frame(rna_annot)
rownames(rna_annot) = rna_annot$sample

variant_annot=read_tsv("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/CRC_MPRA_design_library_all_variants_v3_newgwas_with_controls")
control_snps = variant_annot %>% filter(test_control=="control") %>% pull(rsid)
#snp_data = dna_dat0[,"seq_id"] %>% separate(seq_id, into=c("rsid", "allele", "dir"))
snp_data = rownames(filtered_dna_mat)

is_control = snp_data %in% control_snps

#create MpraObject
obj = MpraObject(dnaCounts=filtered_dna_mat, rnaCounts=filtered_rna_mat, dnaAnnot=dna_annot, rnaAnnot=rna_annot, controls=is_control)
obj <- estimateDepthFactors(obj, lib.factor = "replicate", which.lib = "both", depth.estimator = "uq")
saveRDS(obj, paste0(cell_line,"_filtered_mpranalyze_obj_",max_barcodes,".rds"))


#####

# obj = readRDS(paste0(cell_line,"_filtered_mpranalyze_obj.rds"))


obj_pruned = analyzeComparative(obj, dnaDesign = ~ barcode_full+replicate+allele+dir, rnaDesign = ~replicate+allele+dir, reducedDesign = ~replicate+dir, BPPARAM=MulticoreParam())
saveRDS(obj_pruned, paste0(cell_line,"_filtered_mpranalyze_obj_",max_barcodes,"_fitted.rds"))
resobj_pruned <- testLrt(obj_pruned)
write_tsv(resobj_pruned %>% rownames_to_column(), paste0(cell_line,"_filtered_mpranalyze_res_",max_barcodes,".txt"))


# obj_scale = analyzeComparative(obj, dnaDesign = ~ barcode_full+replicate+allele+dir, rnaDesign = ~replicate+allele+dir, reducedDesign = ~replicate, BPPARAM=MulticoreParam(), mode="scale")
# saveRDS(obj_scale, paste0(cell_line,"_filtered_mpranalyze_obj_scaled_fitted.rds"))
# res_scale <- testLrt(obj_scale)
# write_tsv(res_scale %>% rownames_to_column(), paste0(cell_line,"_filtered_mpranalyze_res_scaled.txt"))

###########

library(tidyverse)
library(MPRAnalyze)
library(glue)

args = commandArgs(trailingOnly=TRUE)
cell_line=args[1]
max_barcodes=as.numeric(args[2])

print(paste("running", cell_line, "filtering to", max_barcodes, "barcodes"))

#cell_line="HCEC-1CT"
setwd(paste0("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/mpraflow/", cell_line))

library(BiocParallel)

# obj = readRDS(paste0(cell_line,"_mpranalyze_obj_extraQC.rds"))
# 
# dna_mat = dnaCounts(obj)
# 
# dna_annot = tibble(sample=colnames(dna_mat))
# dna_annot = dna_annot %>% separate(sample, into=c("type", "replicate", "allele", "dir", "barcode"), remove=F, sep="_", convert=T)
# 
# filtered_dna_annot = dna_annot %>% filter(barcode<max_barcodes)
# filtered_dna_mat = dna_mat[,filtered_dna_annot %>% pull(sample)]
# 
# #remake the annotation file
# filtered_dna_annot = filtered_dna_annot %>% mutate(barcode_full=glue("{barcode}_{allele}_{dir}"))
# filtered_dna_annot = filtered_dna_annot %>% mutate(across(everything(), as.factor))
# filtered_dna_annot = as.data.frame(filtered_dna_annot)
# rownames(filtered_dna_annot) = filtered_dna_annot$sample
# 
# #do the same for RNA
# rna_mat = rnaCounts(obj)
# 
# rna_annot = tibble(sample=colnames(rna_mat))
# rna_annot = rna_annot %>% separate(sample, into=c("type", "replicate", "allele", "dir", "barcode"), remove=F, sep="_", convert=T)
# 
# filtered_rna_annot = rna_annot %>% filter(barcode<max_barcodes)
# filtered_rna_mat = rna_mat[,filtered_rna_annot %>% pull(sample)]
# 
# #remake the annotation file
# filtered_rna_annot = filtered_rna_annot %>% mutate(barcode_full=glue("{barcode}_{allele}_{dir}"))
# filtered_rna_annot = filtered_rna_annot %>% mutate(across(everything(), as.factor))
# filtered_rna_annot = as.data.frame(filtered_rna_annot)
# rownames(filtered_rna_annot) = filtered_rna_annot$sample
# 
# all(rownames(filtered_dna_mat) == rownames(filtered_rna_mat))
# 
# # 
# # ######
# 
# #get control snps
# variant_annot=read_tsv("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/CRC_MPRA_design_library_all_variants_v3_newgwas_with_controls")
# control_snps = variant_annot %>% filter(test_control=="control") %>% pull(rsid)
# is_control = rownames(filtered_dna_mat) %in% control_snps
# 
# # make mpraobject
# obj = MpraObject(dnaCounts=filtered_dna_mat, rnaCounts=filtered_rna_mat, dnaAnnot=filtered_dna_annot, rnaAnnot=filtered_dna_annot, controls=is_control, BPPARAM=MulticoreParam())
# obj <- estimateDepthFactors(obj, lib.factor = c("replicate","allele","dir"), which.lib = "both", depth.estimator = "uq")
# saveRDS(obj, paste0(cell_line,"_mpranalyze_obj_extraQC_",max_barcodes, ".rds"))

#######

obj = readRDS(paste0(cell_line,"_mpranalyze_obj_extraQC_",max_barcodes, ".rds"))

obj = analyzeComparative(obj, dnaDesign = ~ barcode_full+replicate+allele+dir, rnaDesign = ~replicate+allele+dir, reducedDesign = ~replicate+dir, BPPARAM=MulticoreParam())
res <- testLrt(obj)
write_tsv(res %>% rownames_to_column(), paste0(cell_line,"_mpranalyze_res_extraQC_",max_barcodes, ".txt"))


######

labeled_counts=read_tsv(paste0(cell_line, "_final_labeled_counts.txt"))

#exlcude snps with no name/unassigned
clean_counts=labeled_counts %>% filter(label!="no_BC" & !(startsWith(label, "..")) ) %>% select(-Sequence, -Barcode)
clean_counts = clean_counts  %>% replace(is.na(.),0)
colnames(clean_counts) = c("label", "DNA_1", "DNA_2", "DNA_3", "RNA_1", "RNA_2", "RNA_3")

#correct for snps missing rsid - should only be one (rs113101564)
#missing_names=which(snp_data=="")
#snp_data[missing_names]=paste0("missing_", missing_names)

#group by barcode name
n_bc_indiv = clean_counts %>% group_by(label) %>% summarise(n=n())
#separate name into components
n_bc_indiv = n_bc_indiv %>% separate(label, into=c("rsid", "allele", "dir"), remove=F)
#group again by rsid
n_bc_grouped = n_bc_indiv %>% group_by(rsid) %>% summarise(n=n())
print(paste(nrow(n_bc_grouped), "variants mapped to barcodes"))

#find those that have all 4 barcode types (fwd_ref, fwd_alt, rev_ref, rev_alt)
wanted_variants = n_bc_grouped %>% filter(n==4) %>% pull(rsid)
print(paste("there are", length(wanted_variants), "variants with complete information"))

clean_counts = clean_counts %>% separate(label, into=c("rsid", "allele", "dir"), remove=F)
clean_counts = clean_counts %>% filter(rsid %in% wanted_variants)

#clean data further, only keeping dna if corresponding rna is also present in the same replicate (and vice versa)
for (i in 1:3){
  wanted_dna=paste0("DNA_",i)
  wanted_rna=paste0("RNA_",i)
  this_dna=clean_counts %>% pull({{wanted_dna}})
  this_rna=clean_counts %>% pull({{wanted_rna}})
  counts = this_dna*this_rna
  zero_counts = which(counts==0)
  this_dna[zero_counts]=0
  this_rna[zero_counts]=0

  clean_counts[{{wanted_dna}}]=this_dna
  clean_counts[{{wanted_rna}}]=this_rna

}

#remove any barcodes where count is now zero across all replicates
row_counts = clean_counts%>% select(-label, -rsid, -allele, -dir) %>% rowSums()
clean_counts = clean_counts[which(row_counts>0),]


dna_counts=clean_counts %>% select(rsid, allele, dir, starts_with("DNA"))
dna_mat = dna_counts %>% group_by(rsid, allele, dir) %>% mutate(n = 1:n()) %>% pivot_wider(names_from=c(allele, dir, n), values_from=starts_with("DNA"), values_fill = 0)

rna_counts=clean_counts %>% select(rsid, allele, dir, starts_with("RNA"))
rna_mat = rna_counts %>% group_by(rsid, allele, dir) %>% mutate(n = 1:n()) %>% pivot_wider(names_from=c(allele, dir, n), values_from=starts_with("RNA"), values_fill = 0)

# # drop name column, convert to matrix (required by MPRAnalyze), readd rsid as rownames
snp_data = dna_mat %>% pull(rsid)
dna_mat = as.matrix(dna_mat %>% ungroup() %>% select(-rsid))
rownames(dna_mat)=snp_data
rna_mat = as.matrix(rna_mat %>% ungroup() %>% select(-rsid))
rownames(rna_mat)=snp_data

#get control snps
variant_annot=read_tsv("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/CRC_MPRA_design_library_all_variants_v3_newgwas_with_controls")
control_snps = variant_annot %>% filter(test_control=="control") %>% pull(rsid)
is_control = snp_data %in% control_snps

#
#make annotation dataframes from matrix column names
dna_annot = tibble(sample=colnames(dna_mat))
dna_annot = dna_annot %>% separate(sample, into=c("type", "replicate", "allele", "dir", "barcode"), remove=F, sep="_")
dna_annot = dna_annot %>% mutate(barcode_full=glue("{barcode}_{allele}_{dir}"))
rna_annot = tibble(sample=colnames(rna_mat))
rna_annot = rna_annot %>% separate(sample, into=c("type", "replicate", "allele", "dir", "barcode"), remove=F, sep="_")
rna_annot = rna_annot %>% mutate(barcode_full=glue("{barcode}_{allele}_{dir}"))

#make columns into factors
dna_annot = dna_annot %>% mutate(across(everything(), as.factor))
dna_annot = as.data.frame(dna_annot)
rownames(dna_annot) = dna_annot$sample

rna_annot = rna_annot %>% mutate(across(everything(), as.factor))
rna_annot = as.data.frame(rna_annot)
rownames(rna_annot) = rna_annot$sample
#
