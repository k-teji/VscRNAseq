##' Add WFs to scRNA-seq
##' A function to add cell type specific weight factors to scRNAseq, resulting weighted mean and variance scRNAseq count
##' 
##' Argument 
##' scRNAseq_raw_count: A data-frame of normalised count matrix in scRNA-seq
##' WF_summary: A data-frame of the output from WF calculation
##' 
##' Argument format
##' * scRNAseq_raw_count: A data-frame with 'gene' column
##' * WF_summary: A data-frame from the output in WF calculation (data from '..._WF_summary.csv')
##' 
##' Return 
##' A list with:
##' * Weighted_mean: the mean of count data with weighting factors for each cell-type
##' * Weighted_variance: the variance of count data with weighting factors for each cell-type
##' 
##' 
##' 

library(tidyverse)
library(magrittr)

Add_WF_to_scRNAseq <- function(scRNAseq_raw_count, WF_summary){
  
  # list cell types
  cell_type_list <- unique(WF_summary$cell_type)
  cell_type_list <- cell_type_list[order(cell_type_list)]
  
  # celltype loop for addition of WF to scRNAseq
  for(ct in 1:length(cell_type_list)){
    
    # filter data by cell type
    cur_WF_summary <- WF_summary %>%
      dplyr::filter(., cell_type == cell_type_list[ct]) %>%
      dplyr::arrange(sample_label)
    cur_scRNAseq_count <- scRNAseq_raw_count %>%
      dplyr::select(., c('gene', all_of(cur_WF_summary$sample_label))) %>%
      tibble::column_to_rownames('gene')
    
    # variable name conversion for simplicity
    mu <- cur_WF_summary %>%
      dplyr::select(., c(sample_label, mean)) %>%
      tibble::column_to_rownames('sample_label') %>%
      as.matrix()
    sigma <- cur_WF_summary %>%
      dplyr::select(., c(sample_label, variance)) %>%
      tibble::column_to_rownames('sample_label') %>%
      as.matrix()
    x <- cur_scRNAseq_count %>%
      as.matrix()
    
    # Add WF 
    stopifnot(all(colnames(x) == rownames(mu)))
    WF_mu <- (x %*% mu)/ncol(x)
    colnames(WF_mu) <- cell_type_list[ct]
    
    stopifnot(all(colnames(x) == rownames(sigma)))
    WF_var <- ((x^2)%*%sigma)/ncol(x)
    colnames(WF_var) <- cell_type_list[ct]
    
    # data cumulating
    WF_mu <- WF_mu %>%
      data.frame(., stringsAsFactors = F, check.names = F, check.rows = F) %>%
      tibble::rownames_to_column('gene') 
    WF_mu_all <- if(ct == 1L){
      WF_mu
    } else {
      merge(WF_mu_all, WF_mu, by = 'gene')
    }
    
    WF_var <- WF_var %>%
      data.frame(., stringsAsFactors = F, check.names = F, check.rows = F) %>%
      tibble::rownames_to_column('gene')
    WF_var_all <- if(ct == 1L){
      WF_var
    } else {
      merge(WF_var_all, WF_var, by = 'gene')
    }
    
  }# ct
  
  # output
  output <- list()
  
  output[['Weighted_mean']] <- WF_mu_all
  output[['Weighted_variance']] <- WF_var_all
  
  return(output)
}# Add_WF ... function