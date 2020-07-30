# explanation -------------------------------------------------------------

# Original codes outputs only final result contains gene and its p-value
# We need regularized betas, SNP class and gene assignment.


# load library ------------------------------------------------------------

library(genee)
library(dplyr)


# gabage collection before run --------------------------------------------

# this code use huge memory, so clear memory before run
rm(list=ls())
gc()
gc()

# to avoid error of big.matrix
file.rename("ld_mat.bin","ld_mat.txt")

# set path ----------------------------------------------------------------

# working directory decralation is needed for genee package
setwd("genee/")

path_head <- ""
path_hg19 <- "genee/"

load("data/glist.hg19.rda")


# loop --------------------------------------------------------------------

result_list <- lapply(seq(1,22),function(x){
  
  cat(paste("start CHR",formatC(x,width = 2,flag="0")),"\n\n")
  
  
  #########################################################################
  # data setting part
  
  # set filename and file_path
  gwas_name <- paste("snps_or_chr",formatC(x,width = 2,flag="0"),".txt",sep = "")
  ld_name <- paste("corr_mat_chr",formatC(x,width = 2,flag="0"),".csv",sep = "")
  path_gwas <- file.path(path_head,gwas_name)
  path_ld <- file.path(path_head,ld_name)
  
  cat("load data...")
  
  # load data
  data_gwas <- read.table(path_gwas,stringsAsFactors = F, header = T, sep="\t")
  data_ld <- read.table(path_ld,
                        stringsAsFactors = F,
                        header = T,
                        row.names = "SNP_B",
                        sep="")
  ld_mat <- as.matrix(data_ld)
  write.table(ld_mat,file="ld_mat", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
  ld <- setupX("ld_mat", sep = ' ')
  # prepare_LD("ld_mat")
  
  # data convert for input
  input_gwas <- data_gwas %>%
    select("CHR","SNP","BP","OR") %>%
    mutate_at(vars("OR"),log)
  
  # SNP assign to Gene
  # This list is output to result
  # Now GENEE DOES NOT RUN ON SEX CHROMOSOME
  gene_assign <- genee_list(glist.hg19 = glist.hg19,
                            all_chr = input_gwas$CHR,
                            all_pos = input_gwas$BP,
                            upper = 50000,
                            lower = 50000)
  gene_info <- gene_assign[[1]]
  gene_list <- gene_assign[[2]]
  
  # prepare prior weight
  prior_weight <- rep(1, length(input_gwas[,4]))
  
  cat("finish!","\n")

  
  #########################################################################
  # beta regularization part
  
  cat("calc regularization...")
  
  # calculate regularized beta
  test_beta <- genee_regularize_regression_bigdata(betas = input_gwas$OR,
                                                   ld_path = "ld_mat",
                                                   alpha = 0.5,
                                                   lambda.min = 0.0001,
                                                   nlambda = 100,
                                                   nfolds = 10)
  
  # create data for output
  compare_beta <- data.frame(input_gwas,reg_beta=test_beta[-1])
  
  cat("finish!","\n")
  
  #########################################################################
  # SNP classification and epsilon effect extraction part
  
  cat("calc classification...")
  
  # SNP classification by EM algorithm
  # result_em_class contains original data and its class
  # this class data will be returned with final results
  betas_nonzero <- test_beta[which(test_beta!=0)]
  EM_fit <- inference_result_lasso_mclust<-Mclust(betas_nonzero, G = 1:9, modelNames = "V")
  result_em_class <- data.frame(data=EM_fit$data,class=EM_fit$classification)
  
  # extract epsilon effect
  # this part was copied from original codes
  if(length(EM_fit$parameters$variance$sigmasq)>1){
    if(EM_fit$parameters$pro[which(EM_fit$parameters$variance$sigmasq == max(EM_fit$parameters$variance$sigmasq))]>=0.5){
      epsilon_effect = sort(EM_fit$parameters$variance$sigmasq, decreasing = TRUE)[1]
    }else{
      epsilon_effect = sort(EM_fit$parameters$variance$sigmasq, decreasing = TRUE)[2]
    }
  }else{
    epsilon_effect = sort(EM_fit$parameters$variance$sigmasq, decreasing = TRUE)[1]
  }
  
  cat("finish!","\n")
  
  
  #########################################################################
  # Statistical test part
  
  cat("statistical test...")
  
  #run test
  tempresults <- genee_loop(betas = test_beta, 
                            ld = data_ld,
                            epsilon_effect = epsilon_effect,
                            prior_weight = prior_weight, gene_list = gene_list)
  
  #test statistics
  test_statistics <- tempresults[[1]]
  
  #pvals for genes
  tvar <- tempresults[[2]]
  
  #pvals for genes
  pvals <- tempresults[[3]]
  
  
  #########################################################################
  # Final results creation part
  
  chr <- as.numeric(gene_info[,1])
  gene <- as.character(gene_info[,2])
  start <- as.numeric(gene_info[,3])
  end <- as.numeric(gene_info[,4])
  nsnp <- as.numeric(gene_info[,5])
  test_q <- test_statistics
  q_var <- tvar
  pval <- pvals
  final_results <- data.frame(chr, gene, start, end, nsnp, test_q, q_var, pval)
  
  final_results_list <- list(gene_assign = gene_assign,
                             betas_result = compare_beta,
                             SNP_classification = result_em_class,
                             final_results = final_results)
  
  cat("finish!","\n\n")
  
  # remove temporal data
  rm(list=c("data_gwas","data_ld","input_gwas","gene_assign","ld_mat",
            "gene_info","gene_list","prior_weight","test_beta","compare_beta",
            "betas_nonzero","EM_fit","result_em_class","epsilon_effect",
            "tempresults","test_statistics","tvar","pvals",
            "chr","gene","start","end","nsnp","test_q","q_var","pval","final_results"))
  
  # gabage collection
  gc()
  gc()
  
  # remove ld_matrix file
  file.remove("ld_mat")
  file.remove("ld_mat.desc")
  file.rename("ld_mat.bin","ld_mat.txt") # this file remove manually after calculation
  
  # save result for each chromosome for safe
  result_name <- paste("genee_result_CHR",formatC(x,width = 2,flag="0"),".rds",sep = "")
  saveRDS(final_results_list,result_name)
  
  # return final result
  return(final_results_list)
})


# save result -------------------------------------------------------------

saveRDS(result_list,"final_results_list.rds")
