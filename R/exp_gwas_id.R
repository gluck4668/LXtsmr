
# The gwas data ID were obtained from https://gwas.mrcieu.ac.uk/datasets/

exp_gwas_id <- function(exp_dat,exp_p,clump_r2,clump_kb){

exp_df <- extract_instruments(outcomes=exp_dat,p1 = exp_p,
                                clump=F, r2=clump_r2, kb=clump_kb)

exp_df <- subset(exp_df,grepl("rs",exp_df$SNP))

if(nrow(exp_df)>0)
  print(paste0("Note: At the condition of exp_p<",exp_p,", there are ",nrow(exp_df), " SNPs.")) else
  stop(paste0("Note: At the condition of exp_p<",exp_p,", there is no SNP."))

#-------Linkage Disequilibrium (LD) test--------
#write.xlsx(exp_df,"analysis results/exp_df.xlsx")
#exp_df <- read.xlsx("analysis results/exp_df.xlsx")

ld_clum <- exp_ld_clum(exp_df,clump_kb,clump_r2)
exp_clum <- ld_clum$exp_clum

#----------------------------------
exp_data <- exp_clum
exp_file <- paste0("analysis results/",exp_dat,".csv")
write.csv(exp_data,exp_file)

exp_snp <- data.frame(exp_data$SNP)
names(exp_snp) <- "SNP"

result <- list(exp_data=exp_data,exp_snp=exp_snp)
return(result)

}
