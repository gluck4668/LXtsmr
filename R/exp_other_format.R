
exp_other_format <- function(exp_dat,snp_exp,pval_exp,beta_exp,
                             se_exp,effect_allele_exp,other_allele_exp,
                             eaf_exp,clump_kb,clump_r2){

#-----read data file-----
exp_file_type <- str_extract(exp_dat,"(?<=\\.)[^\\.]+$")

if(tolower(exp_file_type)=="gz")
  exp_df <- fread(exp_dat) else
      {if(tolower(exp_file_type)=="txt")
        exp_df <- read_table(exp_dat) else
          exp_df <- eval(str2expression(paste0("read.",exp_file_type,"(exp_dat)")))
        }

#---------------
exp_10 <- head(exp_df)

 for (i in 1:ncol(exp_10)) {
      col_val <- eval(str2expression(paste0("exp_10[,",i,"]")))
      if(any(grepl("rs",col_val))){
        snp_exp_sit <- i
        snp_exp <- names(exp_10)[i]}
    }

exp_df <- eval(str2expression(paste0("distinct(exp_df,",snp_exp,",.keep_all = T)")))
exp_df <-eval(str2expression(paste0("subset(exp_df,",pval_exp,"<",exp_p,")")))

if(nrow(exp_df)>0)
  print(paste0("Note: At the condition of exp_p<",exp_p,", there are ",nrow(exp_df), " SNPs.")) else
  stop(paste0("Note: At the condition of exp_p<",exp_p,", there is no SNP."))

exp_file <- paste0("analysis results/",stringr::str_extract(exp_dat,".*(?=[.])"),".csv")
write.csv(exp_df,exp_file)

colnames(exp_df)[snp_exp_sit] <- "SNP"
exp_file02 <- paste0("analysis results/",stringr::str_extract(exp_dat,".*(?=[.])"),"02.csv")
write.csv(exp_df,exp_file02)

exp_df <- read_exposure_data(filename = exp_file02,
                                   sep = ",",
                                   snp_col = "SNP",
                                   beta_col = beta_exp,
                                   se_col = se_exp,
                                   effect_allele_col = effect_allele_exp,
                                   other_allele_col = other_allele_exp,
                                   eaf_col = eaf_exp,
                                   pval_col = pval_exp)

exp_df <- subset(exp_df,grepl("rs",exp_df$SNP))

if(file.exists(exp_file02))
  file.remove(exp_file02)

if(trimws(toupper(beta_exp)) == "OR")
exp_df$beta.exposure <- log(exp_df$beta.exposure)

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
