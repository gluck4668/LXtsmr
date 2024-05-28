
mr_data <- function(exp_dat,ukbb_anno,clump_kb,clump_r2,exposure_name,
                       exp_p,beta_exp,se_exp,effect_allele_exp,
                       other_allele_exp,eaf_exp,pval_exp,
                       outcome_name,
                       out_dat,beta_out,se_out,effect_allele_out,
                       other_allele_out,eaf_out,pval_out){


#---- Exposure data ---
if(!grepl("[.]",exp_dat,ignore.case = T)){
cat("Begin to read the exposure data......, it may take a long time......")
exp_data <- extract_instruments(outcomes=exp_dat,p1 = exp_p,
                                 clump=TRUE, r2=clump_r2, kb=clump_kb)
cat("Reading the exposure data is done......")
if(nrow(exp_data)>0)
  print(paste0("Note: At the condition of exp_p<",exp_p,", there are ",nrow(exp_data), " SNPs."))
if(nrow(exp_data)==0)
  stop(paste0("Note: At the condition of exp_p<",exp_p,", there is no SNP."))


exp_file <- paste0("analysis results/",exp_dat,".csv")
write.csv(exp_data,exp_file)

exp_snp <- data.frame(exp_data$SNP)
names(exp_snp) <- "SNP"

} else {

if(grepl("vcf.gz",exp_dat,ignore.case = T)) {
#---exposure data---
  cat("Begin to read the exposure data......, it may take a long time......")
  exp_data <- VariantAnnotation::readVcf(exp_dat)  # 读取vcf数据
  exp_data <- gwasglue::gwasvcf_to_TwoSampleMR(vcf =exp_data,type = "exposure") #格式化数据
  exp_data <- distinct(exp_data,SNP,.keep_all = T) %>%
              subset(pval.exposure<exp_p)
  cat("Reading the exposure data is done......")
  if(nrow(exp_data)==0){
    stop(paste0("Note: At the condition of exp_p<",exp_p,", there is no SNP."))
    }
  if(nrow(exp_data)>0){
    print(paste0("Note: At the condition of exp_p<",exp_p,", there are ",nrow(exp_data), " SNPs."))
    }
  exp_data <- clump_data(exp_data,clump_kb = clump_kb,clump_r2 =clump_r2)
  exp_file <- paste0("analysis results/",stringr::str_extract(exp_dat,".*(?=[.])"),".csv")
  write.csv(exp_data,exp_file)

  exp_snp <- data.frame(exp_data$SNP)
  names(exp_snp) <- "SNP"

 } else {

if(grepl(".bgz",exp_dat,ignore.case = T)){
  cat("Begin to read the exposure data......, it may take a long time......")
  ukbb <- fread(exp_dat)
  ukbb_anno <- fread(ukbb_anno)
  exp_df <- inner_join(ukbb_anno,ukbb,by = "variant")
  rm(ukbb, ukbb_anno)
  cat("Reading the exposure data is done......")
  exp_10 <- head(exp_df)

  for (i in 1:ncol(exp_10)) {
    col_val <- eval(str2expression(paste0("exp_10[,",i,"]")))
    if(any(grepl("rs",col_val))){
      snp_exp_sit <- i
      snp_exp <- names(exp_10)[i]}
  }

  exp_df <- eval(str2expression(paste0("distinct(exp_df,",snp_exp,",.keep_all = T)")))
  exp_df <-eval(str2expression(paste0("subset(exp_df,",pval_exp,"<",exp_p,")")))
  if(nrow(exp_df)>0){
    print(paste0("Note: At the condition of exp_p<",exp_p,", there are ",nrow(exp_df), " SNPs."))
  }
  if(nrow(exp_df)==0){
    stop(paste0("Note: At the condition of exp_p<",exp_p,", there is no SNP."))
  }

  exp_file <- paste0("analysis results/",stringr::str_extract(exp_dat,".*(?=[.])"),".csv")
  write.csv(exp_df,exp_file)

  colnames(exp_df)[snp_exp_sit] <- "SNP"
  exp_file02 <- paste0("analysis results/",stringr::str_extract(exp_dat,".*(?=[.])"),"02.csv")
  write.csv(exp_df,exp_file02)

  exp_data <- read_exposure_data(filename = exp_file02,
                                 sep = ",",
                                 snp_col = "SNP",
                                 beta_col = beta_exp,
                                 se_col = se_exp,
                                 effect_allele_col = effect_allele_exp,
                                 other_allele_col = other_allele_exp,
                                 eaf_col = eaf_exp,
                                 pval_col = pval_exp)

  if(trimws(toupper(beta_exp)) == "OR")
    exp_data$beta.exposure <- log(exp_data$beta.exposure)

#--------------------
tryCatch(exp_clum <- TwoSampleMR::clump_data(exp_data,
                                      clump_kb = clump_kb,
                                      clump_r2 = clump_r2),
error=function(e){print(paste("Linkage Disequilibrium (LD) test can't be conducted online",
                             "due to lack of acces to https://api.opengwas.io/api/ld/clump"))}
            )

if(!exists("exp_clum")){
#devtools::install_github("explodecomputer/plinkbinr")
library(plinkbinr)
get_plink_exe()

clum_df <- dplyr::tibble(rsid=exp_data$SNP,
                           pval=exp_data$pval.exposure,
                           id=exp_data$id.exposure)

exp_clum <- ieugwasr::ld_clump(dat =clum_df,
                                       clump_kb = clump_kb,
                                       clump_r2 = clump_r2,
                                       clump_p = 1,
                                       bfile = "1kg.v3/EUR",
                                       plink_bin = get_plink_exe())


exp_clum <- exp_data %>% dplyr::filter(exp_data$SNP %in% exp_clum$rsid)
}

#--------
exp_data <- exp_clum
exp_snp <- data.frame(exp_data$SNP)
names(exp_snp) <- "SNP"

#-------------------
} else

{
  cat("Begin to read the exposure data......, it may take a long time......")
  exp_file_type <- str_extract(exp_dat,"(?<=\\.)[^\\.]+$")
  if(tolower(exp_file_type)=="gz")
    exp_df <- fread(exp_dat) else
    {if(tolower(exp_file_type)=="txt")
      exp_df <- read_table(exp_dat) else
        exp_df <- eval(str2expression(paste0("read.",exp_file_type,"(exp_dat)")))
    }

 cat("Reading the exposure data is done......")

 exp_10 <- head(exp_df)

  for (i in 1:ncol(exp_10)) {
    col_val <- eval(str2expression(paste0("exp_10[,",i,"]")))
    if(any(grepl("rs",col_val))){
      snp_exp_sit <- i
      snp_exp <- names(exp_10)[i]}
  }

exp_df <- eval(str2expression(paste0("distinct(exp_df,",snp_exp,",.keep_all = T)")))
exp_df <-eval(str2expression(paste0("subset(exp_df,",pval_exp,"<",exp_p,")")))
if(nrow(exp_df)>0){
  print(paste0("Note: At the condition of exp_p<",exp_p,", there are ",nrow(exp_df), " SNPs."))
  }
if(nrow(exp_df)==0){
  stop(paste0("Note: At the condition of exp_p<",exp_p,", there is no SNP."))
  }

exp_file <- paste0("analysis results/",stringr::str_extract(exp_dat,".*(?=[.])"),".csv")
write.csv(exp_df,exp_file)

colnames(exp_df)[snp_exp_sit] <- "SNP"
exp_file02 <- paste0("analysis results/",stringr::str_extract(exp_dat,".*(?=[.])"),"02.csv")
write.csv(exp_df,exp_file02)

exp_data <- read_exposure_data(filename = exp_file02,
                                 sep = ",",
                                 snp_col = "SNP",
                                 beta_col = beta_exp,
                                 se_col = se_exp,
                                 effect_allele_col = effect_allele_exp,
                                 other_allele_col = other_allele_exp,
                                 eaf_col = eaf_exp,
                                 pval_col = pval_exp)

if(trimws(toupper(beta_exp)) == "OR")
exp_data$beta.exposure <- log(exp_data$beta.exposure)

if(file.exists(exp_file02))
file.remove(exp_file02)

#--------------------
tryCatch(exp_clum <- TwoSampleMR::clump_data(exp_data,
                                             clump_kb = clump_kb,
                                             clump_r2 = clump_r2),
         error=function(e){print(paste("Linkage Disequilibrium (LD) test can't be conducted online",
                                       "due to lack of acces to https://api.opengwas.io/api/ld/clump"))}
)

if(!exists("exp_clum")){
  #devtools::install_github("explodecomputer/plinkbinr")
  library(plinkbinr)
  get_plink_exe()

  clum_df <- dplyr::tibble(rsid=exp_data$SNP,
                           pval=exp_data$pval.exposure,
                           id=exp_data$id.exposure)

  exp_clum <- ieugwasr::ld_clump(dat =clum_df,
                                 clump_kb = clump_kb,
                                 clump_r2 = clump_r2,
                                 clump_p = 1,
                                 bfile = "1kg.v3/EUR",
                                 plink_bin = get_plink_exe())


exp_clum <- exp_data %>% dplyr::filter(exp_data$SNP %in% exp_clum$rsid)
}

#--------
exp_data <- exp_clum

exp_snp <- data.frame(exp_data$SNP)
names(exp_snp) <- "SNP"
}

}

}

#---- Outcome data ------------------------------------------------------------
if(!grepl("[.]",out_dat,ignore.case = T)){
cat("Begin to read the outcome data, it may take a long time......")
out_data <- extract_outcome_data(snps=exp_data$SNP,
                                 outcomes=out_dat,
                                 proxies = FALSE,
                                 maf_threshold =0.01)
cat("Reading the outcome data is done......")
out_file <- paste0("analysis results/",out_dat,".csv")
write.csv(out_data,out_file)


} else{

if(grepl("vcf.gz",out_dat,ignore.case = T)){
  cat("Begin to read the outcome data, it may take a long time......")
  out_data <- VariantAnnotation::readVcf(out_dat) %>%
              gwasglue::gwasvcf_to_TwoSampleMR(type = "outcome") #格式化数据
  cat("Reading the outcome data is done......")
  out_data <- merge(exp_snp,out_data,by="SNP")

  out_file <- paste0("analysis results/",stringr::str_extract(out_dat,".*(?=[.])"),".csv")
  write.csv(out_data,out_file)

} else

{
cat("Begin to read the outcome data, it may take a long time......")
out_file_type <- str_extract(out_dat,"(?<=\\.)[^\\.]+$")
if(tolower(out_file_type)=="gz")
    out_df <- fread(out_dat) else
    {if(tolower(out_file_type)=="txt")
      out_df <- read_table(out_dat) else
        out_df <- eval(str2expression(paste0("read.",out_file_type,"(out_dat)")))
    }
cat("Reading the outcome data is done......")

out_10 <- head(out_df)

for (x in 1:ncol(out_10)) {
  col_val <- eval(str2expression(paste0("out_10[,",x,"]")))
  if(any(grepl("rs",col_val))){
    snp_out_sit <- x
    snp_out <- names(out_10)[x]
  }
}

colnames(out_df)[snp_out_sit] <- "SNP"
head(out_df)
out_df <- merge(exp_snp,out_df,by="SNP")
out_df <- distinct(out_df,SNP,.keep_all = T)
out_file <- "analysis results/out_data.csv"
write.csv(out_df,out_file)

out_data <- read_outcome_data(
                             snps=exp_data$SNP,
                             filename=out_file,
                             sep = ",",
                             snp_col = "SNP",
                             beta_col = beta_out,
                             se_col = se_out,
                             effect_allele_col = effect_allele_out,
                             other_allele_col = other_allele_out,
                             eaf_col = eaf_out,
                             pval_col = pval_out)

if(trimws(toupper(beta_out)) == "OR")
  out_data$beta.outcome <- log(out_data$beta.outcome)

out_df02 <- eval(str2expression(paste0("dplyr::rename(out_df,",snp_out,"='SNP')")))
out_file02 <- paste0("analysis results/",stringr::str_extract(out_dat,".*(?=[.])"),".csv")
write.csv(out_df02,out_file02)

if(file.exists(out_file))
file.remove(out_file)

}
}
#----harmonise data------------------------------------------------------------
mr_df_local <- harmonise_data(exposure_dat=data.frame(exp_data),
                        outcome_dat=data.frame(out_data),
                        action= 2)

if(exists("exposure_name") & exists("outcome_name")){
   mr_df_local$exposure <- exposure_name
   mr_df_local$outcome <- outcome_name}

write.xlsx(mr_df_local,"analysis results/mr_data.xlsx")

local_data <- list(mr_data=mr_df_local)

return(local_data)

}
