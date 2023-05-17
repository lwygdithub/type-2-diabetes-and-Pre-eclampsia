
library(TwoSampleMR)


library(data.table)


setwd("C:/Users/lwy/Desktop/糖尿病")


a<-fread(input = "finngen_R8_ALLERG_RHINITIS.tsv", sep = "\t", header = T)
View(a)
head(a)

b<-subset(a,p<5e-08)
write.csv(b, file="exposure.csv")
head(b)

bmi_exp_dat_clumped<-read_exposure_data(filename = "exposure.csv",sep = ",",snp_col = "SNP",beta_col = "b",se_col = "se",effect_allele_col = "A1",other_allele_col = "A2",eaf_col = "Freq1.Hapmap", clump = FALSE)
c_clumped <- clump_data(c,clump_kb = 10000,clump_r2 = 0.01,clump_p1 = 1,clump_p2 = 1,pop = "EUR")


setwd("D:/5-孟德尔随机化GWAS数据/Alzheimer IGAP/IGAP_summary_statistics")



outcome_dat<-read_outcome_data(filename="IGAP_summary_statistics.txt", snps = bmi_exp_dat_clumped$SNP, sep = "\t",snp_col = "SNP",beta_col = "beta",se_col = "se",effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")


dat<-harmonise_data(exposure_dat = bmi_exp_dat_clumped,outcome_dat = outcome_dat)


mr(dat)
generate_odds_ratios(mr_res = mr(dat))
mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
mr_scatter_plot(mr_results = mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),dat)
mr_heterogeneity(dat)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
mr_pleiotropy_test(dat)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
run_mr_presso(dat,NbDistribution = 3000)

for(i in 1:length(BMI$eaf.exposure)){
if(BMI$eaf.exposure[i] > 0.5){BMI$maf[i] = 1 - BMI$eaf.exposure[i]}else{BMI$maf[i] = BMI$eaf.exposure[i]}
}