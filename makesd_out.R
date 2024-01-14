library(data.table)
library(dplyr)
library(tidyr)

##generate .esd file for SMR tests
args <- commandArgs(T)
tis <- args[1]
n_chr <- args[2]
Phe1 <- args[3]
Phe2 <- args[4]

chr_file <- paste0("~/2SMR/", Phe1, "2", Phe2, "/", tis, "_", Phe1, "2", Phe2, "_pairs/", tis, ".", Phe1, "2", Phe2, "_chr", n_chr, ".txt")
if(Phe2 == "m6A"){
	Phe_file <- paste0("~/2SMR/", Phe1, "2m6A/", tis, "_m6A-QTL_all/", tis, ".m6A-QTL.all_chr", n_chr, ".txt")
}else{
	Phe_file <- paste0("~/2SMR/m6A2", Phe2, "/", tis, "_", Phe2, "QTL_all/", tis, ".", Phe2, "QTL.all_chr", n_chr, ".txt")
}
maf_file <- paste0("~/Database/GTEx_v8/maf/chr", n_chr, ".frq")

smr_pairs <- fread(chr_file, sep = "\t", header = FALSE, col.names = c(Phe1, Phe2)) %>% distinct(Phe2)
maf <- fread(maf_file, sep = "\t", header = TRUE) %>% select(CHR, SNP, A1, A2, MAF) %>% rename(Chr = CHR, Freq = MAF)
if(Phe1 == "me"){
	Phe_dat <- fread(Phe_file, sep = "\t", header = FALSE, col.names = c("Phenotype", "SNP", "Beta", "effect_allele", "other_allele", "p")) %>%
		inner_join(y = maf, by = "SNP")
}else{
	Phe_dat <- fread(Phe_file, sep = "\t", header = FALSE, col.names = c("Phenotype", "SNP", "Beta", "se", "effect_allele", "other_allele", "Freq", "p"))
}
Phe_dat <- Phe_dat %>% separate(SNP, c(NA, "Bp", NA, NA), ":", convert = TRUE, remove = FALSE) %>%
  mutate(Beta = ifelse(effect_allele == A1, Beta, -Beta),
         z = case_when(Beta > 0 ~ qnorm(1-p/2), Beta < 0 ~ -qnorm(1-p/2)), se = Beta/z) %>%
  select(Phenotype, Chr, SNP, Bp, A1, A2, Freq, Beta, se, p)

for(i in 1:nrow(smr_pairs)){
  Phe_esd <- Phe_dat %>% filter(Phenotype == smr_pairs[i,Phe2]) %>% select(-Phenotype)
  if(Phe2 == "m6A"){
	  Phe_esd_file <- paste0("~/SMR/", Phe1, "2m6A/", tis, "_m6A-QTL_out_esd/", smr_pairs[i,"m6A"], ".esd")
  }else{
	  Phe_esd_file <- paste0("~/SMR/m6A2", Phe2, "/", tis, "_", Phe2, "QTL_out_esd/", smr_pairs[i,Phe2], ".esd")
  }
  write.table(Phe_esd, file = Phe_esd_file, quote = FALSE, sep = "\t", row.names = FALSE)
}
