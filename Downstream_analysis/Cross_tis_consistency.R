library(data.table)
library(dplyr)

##Consistency of m6A-H3K27ac bi-directional MR pairs across four tissues
##m6A-to-H3K27ac
#Merge indexes
m6A_idx_m6A2ha <- data.frame()
for(i in c("Brain", "Lung", "Muscle", "Heart")){
	m6A_idx_m6A2ha <- bind_rows(m6A_idx_m6A2ha,
				    fread(paste0("~/Downstream_analysis/Cross_tis_consistency/", i, ".m6A_idx.m6A2ha.txt"), sep = "\t", header = FALSE,
					  col.names = c("m6A_idx", "m6A_peak")) %>% mutate(Tissue = i))
}

ha_idx_m6A2ha <- data.frame()
for(i in c("Brain", "Lung", "Muscle", "Heart")){
	ha_idx_m6A2ha <- bind_rows(ha_idx_m6A2ha,
				   fread(paste0("~/Downstream_analysis/Cross_tis_consistency/", i, ".ha_idx.m6A2ha.txt"), sep = "\t", header = FALSE,
					 col.names = c("ha_idx", "ha_peak")) %>% mutate(Tissue = i))
}

#Merge MR results
MR_results_m6A2ha <- data.frame()
for(i in c("Brain", "Lung", "Muscle", "Heart")){
        MR_results_m6A2ha <- bind_rows(MR_results_m6A2ha,
				       fread(paste0("~/2SMR/Results/", i, ".m6A2ha.mr.fdr.res"), sep = "\t", header = TRUE) %>%
					       select(exposure, outcome, b, se, pval, fdr) %>% mutate(Tissue = i))
}

MR_results_m6A2ha_Tis1 <- MR_results_m6A2ha %>%
	rename(Exposure_Tis1 = exposure, Outcome_Tis1 = outcome, b_Tis1 = b, se_Tis1 = se, pval_Tis1 = pval, fdr_Tis1 = fdr, Tissue1 = Tissue) %>%
	inner_join(y = m6A_idx_m6A2ha, by = c("Exposure_Tis1" = "m6A_peak", "Tissue1" = "Tissue")) %>%
	inner_join(y = ha_idx_m6A2ha, by = c("Outcome_Tis1" = "ha_peak", "Tissue1" = "Tissue"))

MR_results_m6A2ha_Tis2 <- MR_results_m6A2ha %>%
	rename(Exposure_Tis2 = exposure, Outcome_Tis2 = outcome, b_Tis2 = b, se_Tis2 = se, pval_Tis2 = pval, fdr_Tis2 = fdr, Tissue2 = Tissue) %>%
	inner_join(y = m6A_idx_m6A2ha, by = c("Exposure_Tis2" = "m6A_peak", "Tissue2" = "Tissue")) %>%
        inner_join(y = ha_idx_m6A2ha, by = c("Outcome_Tis2" = "ha_peak", "Tissue2" = "Tissue"))

#Pair MR results from Tissue1 and Tissue2 by index combinations of exposures and outcomes
MR_results_m6A2ha_combine <- MR_results_m6A2ha_Tis1 %>%
	inner_join(y = MR_results_m6A2ha_Tis2, by = c("m6A_idx", "ha_idx")) %>%
	filter(!(Tissue1 == Tissue2 & pval_Tis1 != pval_Tis2))

##H3K27ac-to-m6A
ha_idx_ha2m6A <- data.frame()
for(i in c("Brain", "Lung", "Muscle", "Heart")){
        ha_idx_ha2m6A <- bind_rows(ha_idx_ha2m6A,
				   fread(paste0("~/Downstream_analysis/Cross_tis_consistency/", i, ".ha_idx.ha2m6A.txt"), sep = "\t", header = FALSE,
					 col.names = c("ha_idx", "ha_peak")) %>% mutate(Tissue = i))
}

m6A_idx_ha2m6A <- data.frame()
for(i in c("Brain", "Lung", "Muscle", "Heart")){
        m6A_idx_ha2m6A <- bind_rows(m6A_idx_ha2m6A,
				    fread(paste0("~/Downstream_analysis/Cross_tis_consistency/", i, ".m6A_idx.ha2m6A.txt"), sep = "\t", header = FALSE,
					  col.names = c("m6A_idx", "m6A_peak")) %>% mutate(Tissue = i))
}

MR_results_ha2m6A <- data.frame()
for(i in c("Brain", "Lung", "Muscle", "Heart")){
        MR_results_ha2m6A <- bind_rows(MR_results_ha2m6A,
				       fread(paste0("~/2SMR/Results/", i, ".m6A2ha.mr.fdr.res"), sep = "\t", header = TRUE) %>%
					       select(exposure, outcome, b, se, pval, fdr) %>% mutate(Tissue = i))
}

MR_results_ha2m6A_Tis1 <- MR_results_ha2m6A %>%
        rename(Exposure_Tis1 = exposure, Outcome_Tis1 = outcome, b_Tis1 = b, se_Tis1 = se, pval_Tis1 = pval, fdr_Tis1 = fdr, Tissue1 = Tissue) %>%
        inner_join(y = ha_idx_ha2m6A, by = c("Exposure_Tis1" = "ha_peak", "Tissue1" = "Tissue")) %>%
        inner_join(y = m6A_idx1_ha2m6A, by = c("Outcome_Tis1" = "m6A_peak", "Tissue1" = "Tissue"))

MR_results_ha2m6A_Tis2 <- MR_results_ha2m6A %>%
        rename(Exposure_Tis2 = exposure, Outcome_Tis2 = outcome, b_Tis2 = b, se_Tis2 = se, pval_Tis2 = pval, fdr_Tis2 = fdr, Tissue2 = Tissue) %>%
        inner_join(y = ha_idx_ha2m6A, by = c("Exposure_Tis2" = "ha_peak", "Tissue2" = "Tissue")) %>%
        inner_join(y = m6A_idx_ha2m6A, by = c("Outcome_Tis2" = "m6A_peak", "Tissue2" = "Tissue"))

MR_results_ha2m6A_combine <- MR_results_ha2m6A_Tis1 %>%
        inner_join(y = MR_results_ha2m6A_Tis2, by = c("ha_idx", "m6A_idx")) %>%
        filter(!(Tissue1 == Tissue2 & pval_Tis1 != pval_Tis2))

##Consistency of m6A-DNAme bi-directional MR pairs across lung and muscle
##m6A-to-DNAme
#Merge indexes, only m6A peaks
m6A_idx_m6A2me <- data.frame()
for(i in c("Lung", "Muscle")){
        m6A_idx_m6A2me <- bind_rows(m6A_idx_m6A2me,
                                    fread(paste0("~/Downstream_analysis/Cross_tis_consistency/", i, ".m6A_idx.m6A2me.txt"), sep = "\t", header = FALSE,
                                          col.names = c("m6A_idx", "m6A_peak")) %>% mutate(Tissue = i))
}

#Merge MR results
MR_results_m6A2me <- data.frame()
for(i in c("Lung", "Muscle")){
        MR_results_m6A2me <- bind_rows(MR_results_m6A2me,
                                       fread(paste0("~/2SMR/Results/", i, ".m6A2me.mr.fdr.res"), sep = "\t", header = TRUE) %>%
					       select(exposure, outcome, b, se, pval, fdr) %>% mutate(Tissue = i))
}

MR_results_m6A2me_Tis1 <- MR_results_m6A2me %>%
        rename(Exposure_Tis1 = exposure, Outcome_Tis1 = outcome, b_Tis1 = b, se_Tis1 = se, pval_Tis1 = pval, fdr_Tis1 = fdr, Tissue1 = Tissue) %>%
        inner_join(y = m6A_idx_m6A2me, by = c("Exposure_Tis1" = "m6A_peak", "Tissue1" = "Tissue"))

MR_results_m6A2me_Tis2 <- MR_results_m6A2me %>%
        rename(Exposure_Tis2 = exposure, Outcome_Tis2 = outcome, b_Tis2 = b, se_Tis2 = se, pval_Tis2 = pval, fdr_Tis2 = fdr, Tissue2 = Tissue) %>%
        inner_join(y = m6A_idx_m6A2me, by = c("Exposure_Tis2" = "m6A_peak", "Tissue2" = "Tissue"))

#Pair MR results from Tissue1 and Tissue2 by index combinations of exposures and outcomes
MR_results_m6A2me_combine <- MR_results_m6A2me_Tis1 %>%
        inner_join(y = MR_results_m6A2me_Tis2, by = c("m6A_idx", "Outcome_Tis1" = "Outcome_Tis2")) %>%
        filter(!(Tissue1 == Tissue2 & pval_Tis1 != pval_Tis2))

##DNAme-to-m6A
m6A_idx_me2m6A <- data.frame()
for(i in c("Lung", "Muscle")){
        m6A_idx_me2m6A <- bind_rows(m6A_idx_me2m6A,
                                    fread(paste0("~/Downstream_analysis/Cross_tis_consistency/", i, ".m6A_idx.me2m6A.txt"), sep = "\t", header = FALSE,
                                          col.names = c("m6A_idx", "m6A_peak")) %>% mutate(Tissue = i))
}

MR_results_me2m6A <- data.frame()
for(i in c("Lung", "Muscle")){
        MR_results_me2m6A <- bind_rows(MR_results_me2m6A,
                                       fread(paste0("~/2SMR/Results/", i, ".me2m6A.mr.fdr.res"), sep = "\t", header = TRUE) %>%
                                               select(exposure, outcome, b, se, pval, fdr) %>% mutate(Tissue = i))
}

MR_results_me2m6A_Tis1 <- MR_results_me2m6A %>%
        rename(Exposure_Tis1 = exposure, Outcome_Tis1 = outcome, b_Tis1 = b, se_Tis1 = se, pval_Tis1 = pval, fdr_Tis1 = fdr, Tissue1 = Tissue) %>%
        inner_join(y = m6A_idx_me2m6A, by = c("Outcome_Tis1" = "m6A_peak", "Tissue1" = "Tissue")) 

MR_results_me2m6A_Tis2 <- MR_results_me2m6A %>%
        rename(Exposure_Tis2 = exposure, Outcome_Tis2 = outcome, b_Tis2 = b, se_Tis2 = se, pval_Tis2 = pval, fdr_Tis2 = fdr, Tissue2 = Tissue) %>%
        inner_join(y = m6A_idx_me2m6A, by = c("Outcome_Tis2" = "m6A_peak", "Tissue2" = "Tissue")) 

MR_results_me2m6A_combine <- MR_results_me2m6A_Tis1 %>%
        inner_join(y = MR_results_me2m6A_Tis2, by = c("Exposure_Tis1" = "Exposure_Tis2", "m6A_idx")) %>%
        filter(!(Tissue1 == Tissue2 & pval_Tis1 != pval_Tis2))

##Calculate weighted correlations using standard errors of MR results in Tissue2
library(weights)

##m6A-H3K27ac across four tissues
wtd.cor(MR_results_m6A2ha_combine[MR_results_m6A2ha_combine$fdr_Tis1 < 0.1,]$b_Tis1, MR_results_m6A2ha_combine[MR_results_m6A2ha_combine$fdr_Tis1 < 0.1,]$b_Tis2, 1/(MR_results_m6A2ha_combine[MR_results_m6A2ha_combine$fdr_Tis1 < 0.1,]$se_Tis2)^2)
wtd.cor(MR_results_ha2m6A_combine[MR_results_ha2m6A_combine$fdr_Tis1 < 0.1,]$b_Tis1, MR_results_ha2m6A_combine[MR_results_ha2m6A_combine$fdr_Tis1 < 0.1,]$b_Tis2, 1/(MR_results_ha2m6A_combine[MR_results_ha2m6A_combine$fdr_Tis1 < 0.1,]$se_Tis2)^2)
##m6A-DNAme across lung and muscle
wtd.cor(MR_results_m6A2me_combine[MR_results_m6A2me_combine$Tissue1 == "Lung" & MR_results_m6A2me_combine$fdr_Tis1 < 0.1,]$b_Tis1, MR_results_m6A2me_combine[MR_results_m6A2me_combine$Tissue1 == "Lung" & MR_results_m6A2me_combine$fdr_Tis1 < 0.1,]$b_Tis2, 1/(MR_results_m6A2me_combine[MR_results_m6A2me_combine$Tissue1 == "Lung" & MR_results_m6A2me_combine$fdr_Tis1 < 0.1,]$se_Tis2)^2)
wtd.cor(MR_results_m6A2me_combine[MR_results_m6A2me_combine$Tissue1 == "Muscle" & MR_results_m6A2me_combine$fdr_Tis1 < 0.1,]$b_Tis1, MR_results_m6A2me_combine[MR_results_m6A2me_combine$Tissue1 == "Muscle" & MR_results_m6A2me_combine$fdr_Tis1 < 0.1,]$b_Tis2, 1/(MR_results_m6A2me_combine[MR_results_m6A2me_combine$Tissue1 == "Muscle" & MR_results_m6A2me_combine$fdr_Tis1 < 0.1,]$se_Tis2)^2)
wtd.cor(MR_results_me2m6A_combine[MR_results_me2m6A_combine$Tissue1 == "Lung" & MR_results_me2m6A_combine$fdr_Tis1 < 0.1,]$b_Tis1, MR_results_me2m6A_combine[MR_results_me2m6A_combine$Tissue1 == "Lung" & MR_results_me2m6A_combine$fdr_Tis1 < 0.1,]$b_Tis2, 1/(MR_results_me2m6A_combine[MR_results_me2m6A_combine$Tissue1 == "Lung" & MR_results_me2m6A_combine$fdr_Tis1 < 0.1,]$se_Tis2)^2)
wtd.cor(MR_results_me2m6A_combine[MR_results_me2m6A_combine$Tissue1 == "Muscle" & MR_results_me2m6A_combine$fdr_Tis1 < 0.1,]$b_Tis1, MR_results_me2m6A_combine[MR_results_me2m6A_combine$Tissue1 == "Muscle" & MR_results_me2m6A_combine$fdr_Tis1 < 0.1,]$b_Tis2, 1/(MR_results_me2m6A_combine[MR_results_me2m6A_combine$Tissue1 == "Muscle" & MR_results_me2m6A_combine$fdr_Tis1 < 0.1,]$se_Tis2)^2)
