
exp_ld_clum <- function(exp_df,clump_kb,clump_r2){

#-------Linkage Disequilibrium (LD) test--------
tryCatch(exp_clum <- TwoSampleMR::clump_data(exp_df,clump_kb = clump_kb,clump_r2 = clump_r2),
         error=function(e){cat(paste("Note: The LD test cannot be conducted online due to ",
                                     "the unavailability of https://api.opengwas.io/api/ld/clump.",
                                     "Attempting a local LD test using the ieugwasr::ld_clump function.",sep="\n" ))}
        )

if(exists("exp_clum")){
  result <- list(exp_clum=exp_clum)
  return(result)} else
{ print("-------------------------------------------------")
  print("1. Please create a folder named 'local_clum_data'")
  print("2. Download the '1kg.v3.tgz' from http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz")
  print("3. Unzip '1kg.v3.tgz' and copy all the 15 files to 'local_clum_data'. ")

  if(!file.exists("local_clum_data/EUR.bed")){
    url_ld <- "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz"
    download.file(url_ld,"1kg.v3.tgz")}

  #devtools::install_github("explodecomputer/plinkbinr")
  library(plinkbinr)
  get_plink_exe()

  clum_df <- dplyr::tibble(rsid=exp_df$SNP,
                           pval=exp_df$pval.exposure,
                           id=exp_df$id.exposure)

  tryCatch(exp_clum <- ieugwasr::ld_clump(dat =clum_df,
                                          clump_kb = clump_kb,
                                          clump_r2 = clump_r2,
                                          clump_p = 1,
                                          bfile = "local_clum_data/EUR",
                                          plink_bin = get_plink_exe()),
           error=function(e){stop(paste0("Note: At the condition of exp_p<",exp_p,
                                         ", there is no suitable SNPs for LD test.Please try to alter exp_p value."))}
  )

  exp_clum <- exp_df %>% dplyr::filter(exp_df$SNP %in% exp_clum$rsid)

  result <- list(exp_clum=exp_clum)

  return(result)

}

}
