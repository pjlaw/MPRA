library(tidyverse)
library(glue)
library(MPRAnalyze)

cell_line = "RPMI_8226" #"KMS_11" #"L363" 
setwd(paste0("//cancgene/cancgene/shared/Mol and Pop Genetics Team/Lab Members/Molly/GWAS Work/MPRA/MPRAnalyze/", cell_line))

dna_dat=read_tsv("dna_counts.tsv")
dna_dat0 = dna_dat %>% replace(is.na(.),0)
rna_dat=read_tsv("rna_counts.tsv")
rna_dat0 = rna_dat %>% replace(is.na(.),0)

dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "1_1", "1:1")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "10_1", "10:1")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "10_6", "10:6")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "11_1", "11:1")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "12_9", "12:9")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "13_3", "13:3")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "15_7", "15:7")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "15_8", "15:8")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "16_9", "16:9")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "17_4", "17:4")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "17_5", "17:5")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "17_6", "17:6")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "19_4", "19:4")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "19_7", "19:7")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "2_1", "2:1")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "2_2", "2:2")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "20_3", "20:3")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "3_1", "3:1")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "3_3", "3:3")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "5_1", "5:1")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "6_3", "6:3")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "7_1", "7:1")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "8_1", "8:1")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "8_7", "8:7")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "8_8", "8:8")
dna_dat0$seq_id <- str_replace(dna_dat0$seq_id, "9_3", "9:3")



#pivot the matrices so column
dna_dat0 = dna_dat0 %>% separate(seq_id, into=c("rsid", "dir", "allele"), sep="_")
dna_dat0$dir <- str_replace(dna_dat0$dir, "fdw", "fwd")
dna_dat0 %>% select(rsid, dir,allele) -> manifest
dna_dat0$dir <- str_replace(dna_dat0$dir, "ctrl", "fwd")
dna_dat0 = dna_dat0 %>% replace(is.null(.),0)

dna_dat0 %>% filter(rsid!="rs1956905" & rsid!="rs3825588") -> dna_dat0
manifest %>% filter(rsid!="rs1956905" & rsid!="rs3825588") -> manifest


#pivot tables to group by rsid
dna_mat = dna_dat0 %>% distinct() %>% pivot_wider(names_from=c(allele, dir), values_from=c(-rsid,-allele,-dir))
dna_dat0 = dna_dat0 %>% replace(is.na(.),0)
#dna_mat <- dna_mat %>% mutate(across(everything(), ~replace(., lengths(.) == 0, 0)))



rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "1_1", "1:1")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "10_1", "10:1")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "10_6", "10:6")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "11_1", "11:1")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "12_9", "12:9")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "13_3", "13:3")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "15_7", "15:7")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "15_8", "15:8")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "16_9", "16:9")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "17_4", "17:4")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "17_5", "17:5")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "17_6", "17:6")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "19_4", "19:4")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "19_7", "19:7")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "2_1", "2:1")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "2_2", "2:2")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "20_3", "20:3")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "3_1", "3:1")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "3_3", "3:3")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "5_1", "5:1")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "6_3", "6:3")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "7_1", "7:1")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "8_1", "8:1")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "8_7", "8:7")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "8_8", "8:8")
rna_dat0$seq_id <- str_replace(rna_dat0$seq_id, "9_3", "9:3")


rna_dat0 = rna_dat0 %>% separate(seq_id, into=c("rsid", "dir", "allele"), sep="_")
rna_dat0$dir <- str_replace(rna_dat0$dir, "fdw", "fwd")
rna_dat0$dir <- str_replace(rna_dat0$dir, "ctrl", "fwd")
rna_dat0 = rna_dat0 %>% replace(is.null(.),0)

rna_dat0 %>% filter(rsid!="rs1956905" & rsid!="rs3825588") -> rna_dat0

rna_mat = rna_dat0 %>% pivot_wider(names_from=c(allele, dir), values_from=c(-rsid, -allele, -dir))
rna_dat0 = rna_dat0 %>% replace(is.na(.),0)

 
all(dna_mat$rsid ==rna_mat$rsid)
snp_data = dna_mat$rsid

#don't want to include control SNPs
#control_snps = manifest %>% filter(dir=="ctrl") %>% pull(rsid)
#is_control = snp_data %in% control_snps

# drop name column, convert to matrix (required by MPRAnalyze), readd rsid as rownames
dna_mat=as.matrix(dna_mat[,-1])
rna_mat=as.matrix(rna_mat[,-1])
rownames(dna_mat)=snp_data
rownames(rna_mat)=snp_data

#make annotation dataframes from matrix column names
dna_annot = tibble(sample=colnames(dna_mat))
dna_annot = dna_annot %>% separate(sample, into=c("type", "cell_line", "replicate", "barcode", "allele", "dir"), remove=F, sep="_")
dna_annot = dna_annot %>% mutate(barcode_full=glue("{barcode}_{allele}_{dir}"))
rna_annot = tibble(sample=colnames(rna_mat))
rna_annot = rna_annot %>% separate(sample, into=c("type", "cell_line", "replicate", "barcode", "allele", "dir"), remove=F, sep="_")
rna_annot = rna_annot %>% mutate(barcode_full=glue("{barcode}_{allele}_{dir}"))

#make columns into factors
dna_annot = dna_annot %>% mutate(across(everything(), as.factor))
dna_annot = as.data.frame(dna_annot)
rownames(dna_annot) = dna_annot$sample

rna_annot = rna_annot %>% mutate(across(everything(), as.factor))
rna_annot = as.data.frame(rna_annot)
rownames(rna_annot) = rna_annot$sample

#create MpraObject
obj = MpraObject(dnaCounts=dna_mat, rnaCounts=rna_mat, dnaAnnot=dna_annot, rnaAnnot=rna_annot)
obj <- estimateDepthFactors(obj, lib.factor = "replicate", which.lib = "both", depth.estimator = "uq")
saveRDS(obj, "RPMI8226_mpranalyze_obj.rds")

#obj = readRDS(paste0(cell_line,"_mpranalyze_obj_extraQC_",max_barcodes, ".rds"))

obj = analyzeComparative(obj, dnaDesign = ~ barcode_full+replicate+allele+dir,
                         rnaDesign = ~replicate+allele+dir,
                         reducedDesign = ~replicate+dir, mode='scale')
res <- testLrt(obj)
write_tsv(res %>% rownames_to_column(), "RPMI8226_mpranalyze_res_testrun.txt")

