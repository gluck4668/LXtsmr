
LXtsmr <- function(){

#-------R packages ----
inst_pack <- function(){

packs <- c("BiocManager","devtools","remote","openxlsx","dplyr","readr",
           "purrr","ggplot2","tidyr","forestploter","tidyverse","stringr",
           "shiny","remotes","MendelianRandomization","simex",
           "data.table","VariantAnnotation")

pack_uninst <- packs[!packs %in% installed.packages()[,1]]

if(length(pack_uninst)>0){
      tryCatch(install.packages(pack_uninst),error=function(e){e})
      tryCatch(BiocManager::install(pack_uninst),error=function(e){e})}
#------
mr_packs <- c("MRInstruments","TwoSampleMR","gwasvcf","gwasglue",
              "RMVMR","MRPRESSO","ieugwasr","plinkbinr")

if(!"ieugwasr" %in% installed.packages()[,1])
  install.packages("ieugwasr-master/ieugwasr_1.0.0.tar.gz",
                   repos = NULL, type = "source")

if(!"MRInstruments" %in% installed.packages()[,1])
  remotes::install_github("MRCIEU/MRInstruments")

if(!"TwoSampleMR" %in% installed.packages()[,1])
  remotes::install_github("MRCIEU/TwoSampleMR")

if(!"MRPRESSO" %in% installed.packages()[,1])
  remotes::install_github("rondolab/MR-PRESSO")

if(!"RMVMR" %in% installed.packages()[,1])
  remotes::install_github("WSpiller/RMVMR")

if(!"gwasvcf" %in% installed.packages()[,1])
  remotes::install_github("mrcieu/gwasvcf")

if(!"gwasglue" %in% installed.packages()[,1])
  remotes::install_github("mrcieu/gwasglue")

if(!"plinkbinr" %in% installed.packages()[,1])
devtools::install_github("explodecomputer/plinkbinr")

#------


for(i in c(packs,mr_packs)){
  library(i,character.only = T)}

 }

inst_pack()

#---create dir-----------------------------------------------------------------
if(!dir.exists("analysis results"))
    dir.create("analysis results")

#---MR data (gz,vcf.gz,csv,txt format)-----------------------------------------
mr_df <- mr_data_processing (exp_dat,clump_kb,clump_r2,
                             exposure_name, exp_p,beta_exp,se_exp,
                             effect_allele_exp, other_allele_exp,eaf_exp,pval_exp,
                             outcome_name,out_dat,beta_out,se_out,effect_allele_out,
                             other_allele_out,eaf_out,pval_out)

mr_data <- mr_df$mr_data

#---MR analysis----------------------------------------------------------------
mr_analysis <- mr_analysis(mr_data)

print("The results can be found in the folder of <analysis results> ")


#---BWMR analysis--------------------------------------------------------------
print("begin the BWMR analysis......")
bwmr_analysis <- bwmr_analysis(mr_data,exposure_name,outcome_name)
print("The BWMR analysis is done......")

#--- show one of the plot of BWMR analysis----
bwmr_analysis$result

}
