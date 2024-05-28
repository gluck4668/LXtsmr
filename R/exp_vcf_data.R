
# VCF data were downloaed from https://gwas.mrcieu.ac.uk/datasets/
exp_vcf_data <- function(exp_dat,exp_p){

exp_df <- VariantAnnotation::readVcf(exp_dat)  # 读取vcf数据
exp_df <- gwasglue::gwasvcf_to_TwoSampleMR(vcf =exp_df,type = "exposure") #格式化数据
exp_df <- distinct(exp_df,SNP,.keep_all = T) %>%
            subset(pval.exposure<exp_p) %>%
            subset(grepl("rs",exp_df$SNP))

if(nrow(exp_df)>0){
  print(paste0("Note: At the condition of exp_p<",exp_p,", there are ",nrow(exp_df), " SNPs."))
  } else
    stop(paste0("Note: At the condition of exp_p<",exp_p,", there is no SNP."))

#-------Linkage Disequilibrium (LD) test--------
#write.xlsx(exp_df,"analysis results/exp_df.xlsx")
#exp_df <- read.xlsx("analysis results/exp_df.xlsx")

ld_clum <- exp_ld_clum(exp_df,clump_kb,clump_r2)
exp_clum <- ld_clum$exp_clum

#------------------------------
exp_data <- exp_clum
exp_file <- paste0("analysis results/",stringr::str_extract(exp_dat,".*(?=[.])"),".csv")
write.csv(exp_data,exp_file)

exp_snp <- data.frame(exp_data$SNP)
names(exp_snp) <- "SNP"

result <- list(exp_data=exp_data,exp_snp=exp_snp)
return(result)

}
