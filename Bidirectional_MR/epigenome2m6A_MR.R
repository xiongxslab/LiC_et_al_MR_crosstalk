library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(plinkbinr)
library(foreach)
library(doParallel)

perform_mr <- function(i){
	Phe2_out_dat <- Phe2_dat %>% filter(Phenotype == mr_pairs[i,Phe2])
	Phe1_exp_dat <- Phe1_dat %>% filter(Phenotype == mr_pairs[i,Phe1]) %>% filter(SNP %in% Phe2_out_dat$SNP)
	check <- 1
	if(length(Phe1_exp_dat$SNP) > 1){
		check <<- try(
			      Phe1_exp_dat <- ld_clump_local(Phe1_exp_dat %>% rename(rsid = SNP), clump_kb = 100, clump_r2 = 0.01, clump_p = 1,
							     bfile = b_file, plink_bin = get_plink_exe()) %>% rename(SNP = rsid),
			      silent = FALSE
		)
	}
	if("try-error" %in% class(check)){
		next
	}
	Phe2_out_dat <- Phe2_out_dat %>% filter(SNP %in% Phe1_exp_dat$SNP)
	Phe1_exp_dat <- format_data(Phe1_exp_dat, type = "exposure", phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta", se_col = "se",
				    effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "pval", samplesize_col = "samplesize")
  	Phe2_out_dat <- format_data(Phe2_out_dat, type = "outcome", phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta", se_col = "se",
				    effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "pval", samplesize_col = "samplesize")
  	mr_dat <- harmonise_data(exposure_dat = Phe1_exp_dat, outcome_dat = Phe2_out_dat) %>% filter(mr_keep == TRUE)
	if(length(mr_dat$SNP) >= 1){
		if(length(mr_dat$SNP) == 1){
			mr_res <- mr(mr_dat, method_list = "mr_wald_ratio")
		}else if(length(mr_dat$SNP) > 1){
			mr_res <- mr(mr_dat, method_list = "mr_ivw")
		}
		mr_dat_file <- paste0("~/2SMR/", Phe1, "2m6A/", tis, "_", Phe1, "2m6A_dat/", mr_pairs[i,Phe1], "-", mr_pairs[i,Phe2], ".mr.dat")
		write.table(mr_dat, file = mr_dat_file, quote = FALSE, sep = "\t", row.names = FALSE)
		mr_res <- mr_res %>% select(-id.exposure, -id.outcome) %>% select(exposure, everything())
		return(mr_res)
	}
}

##Preform epigenome to m6A MR across four tissues
args <- commandArgs(T)
tis <- args[1]
n_chr <- args[2]
Phe1 <- args[3] ##"ha" or "me"

Phe2 <- "m6A"
##list of sample sizes for m6A-QTLs, haQTLs, and mQTLs in four tissues
n <- list("m6A"=c("Brain"=53, "Lung"=32, "Muscle"=40, "Heart"=40),
	  "ha"=c("Brain"=113, "Lung"=66, "Muscle"=108, "Heart"=100),
          "me"=c("Lung"=190, "Muscle"=42))
##list of degrees of freedom for m6A-QTLs and haQTLs in four tissues, used for calculating se
df <- list("m6A"=c("Brain"=44, "Lung"=28, "Muscle"=33, "Heart"=33),
	   "ha"=c("Brain"=96, "Lung"=62, "Muscle"=96, "Heart"=93))

chr_file <- paste0("~/2SMR/", Phe1, "2m6A/", tis, "_", Phe1, "2m6A_pairs/", tis, ".", Phe1, "2m6A_chr", n_chr, ".txt")
Phe1_file <- paste0("~/Database/", tis, "_", Phe1, "QTL_signif/", tis, ".", Phe1, "QTL.signif_chr", n_chr, ".txt")
Phe2_file <- paste0("~/Database/", tis, "_m6A-QTL_all/", tis, ".m6A-QTL.all_chr", n_chr, ".txt")
maf_file <- paste0("~/Database/GTEx_v8/maf/chr", n_chr, ".frq")
b_file <- paste0("~/Database/GTEx_v8/bfile/chr", n_chr)
out_file <- paste0("~/2SMR/", Phe1, "2m6A/", tis, "_", Phe1, "2m6A_mr/chr", n_chr, ".mr.res")

mr_pairs <- fread(chr_file, sep = "\t", header = FALSE, col.names = c(Phe1, Phe2))
maf <- fread(maf_file, sep = "\t", header = TRUE) %>% select(SNP, MAF) %>% rename(eaf = MAF)
if(Phe1 == "ha"){
	Phe1_dat <- fread(Phe1_file, sep = "\t", header = FALSE, col.names = c("Phenotype", "SNP", "beta", "effect_allele", "other_allele", "pval")) %>%
		mutate(se = abs(beta/qt(pval/2, df[[Phe1]][tis])), samplesize = n[[Phe1]][tis]) %>% inner_join(y = maf, by = "SNP")
}else{
	Phe1_dat <- fread(Phe1_file, sep = "\t", header = FALSE, col.names = c("Phenotype", "SNP", "beta", "se", "effect_allele", "other_allele", "eaf", "pval")) %>%
                mutate(samplesize = n[[Phe1]][tis])
}
Phe2_dat <- fread(Phe2_file, sep = "\t", header = FALSE, col.names = c("Phenotype", "SNP", "beta", "effect_allele", "other_allele", "pval")) %>%
	mutate(se = abs(beta/qt(pval/2, df[[Phe2]][tis])), samplesize = n[[Phe2]][tis]) %>% inner_join(y = maf, by = "SNP")

cl <- makeCluster(12)
registerDoParallel(cl)

res <- foreach(i = 1:nrow(mr_pairs), .combine = rbind, .packages = c("dplyr", "TwoSampleMR", "ieugwasr", "plinkbinr"), .errorhandling = "pass") %dopar% {
	perform_mr(i)
}
write.table(res, file = out_file, quote = FALSE, sep = "\t", row.names = FALSE)

stopImplicitCluster()
stopCluster(cl)
