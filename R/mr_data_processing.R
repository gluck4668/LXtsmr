

mr_data_processing <- function(exp_dat,clump_kb,clump_r2,
                    exposure_name, exp_p,beta_exp,se_exp,
                    effect_allele_exp, other_allele_exp,eaf_exp,pval_exp,
                    outcome_name,out_dat,beta_out,se_out,effect_allele_out,
                    other_allele_out,eaf_out,pval_out){

#------- Exposure data -------------------------------------------------------
#------ gwas id format------
if(!grepl("[.]",exp_dat,ignore.case = T)){
 print("Begin to read the exposure data, it may take a long time......")
 exp_set <- exp_gwas_id (exp_dat,exp_p,clump_r2,clump_kb)
 exp_data <- exp_set$exp_data
 exp_snp <- exp_set$exp_snp
 print("Processing the exposure data is done......")} else{

#------VCF.gz format------
 if(grepl("vcf.gz",exp_dat,ignore.case = T)){
  print("Begin to read the exposure data, it may take a long time......")
  exp_set <- exp_vcf_data(exp_dat,exp_p)
  exp_data <- exp_set$exp_data
  exp_snp <- exp_set$exp_snp
  print("Processing the exposure data is done......")
  } else {

#-------bgz format------
if(grepl(".bgz",exp_dat,ignore.case = T)){
    print("Begin to read the exposure data, it may take a long time......")

    exp_set <- exp_bgz_data (exp_dat,snp_exp,pval_exp,beta_exp,
                             se_exp,effect_allele_exp,other_allele_exp,
                             eaf_exp,clump_kb,clump_r2)

    exp_data <- exp_set$exp_data
    exp_snp <- exp_set$exp_snp
    print("Processing the exposure data is done......")
  } else


#-------the other format------
{
  print("Begin to read the exposure data, it may take a long time......")
    exp_set <- exp_other_format(exp_dat,snp_exp,pval_exp,beta_exp,
                             se_exp,effect_allele_exp,other_allele_exp,
                             eaf_exp,clump_kb,clump_r2)

    exp_data <- exp_set$exp_data
    exp_snp <- exp_set$exp_snp
    print("Processing the exposure data is done......")
  }

}
}


#---- Outcome data ------------------------------------------------------------
if(grepl(".bgz",out_dat,ignore.case = T)){
print("Begain to read the outcome data, it may take a long time....")
out_df <- fread(out_dat)
if(grepl("variant",names(out_df)) %>% any()){
    if(!file.exists("variants.tsv.bgz")){
    ukbb_url <- "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz"
    download.file(ukbb_url,destfile = "variants.tsv.bgz",method = "libcurl" )}
    ukbb_anno <- fread("variants.tsv.bgz")
    out_df <- inner_join(ukbb_anno,out_df,by = "variant")
    rm(ukbb_anno)
  }

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


} else{
#-------

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







