#' Sampling cells on ratio for calculating weighting factors
#' Sampling cells on ratio for calculating weighting factors
#' 
#' Argument 
#' scRNAseq_annotation_data: A data-frame of sample_labels and cell-types labels
#' celltype_reference_data: A data-frame of cell_types and the corresponding ratio data
#' sample_scale: A integer value to sample as total
#' seed: seed
#' 
#' 

library(dplyr)
library(magrittr)

##=== function to sample cells on ratio ===##
sample_cells_for_WF <- function(scRNAseq_annotation_data, 
                                celltype_reference_data, 
                                sample_scale = sample_scale, 
                                seed){
  
  # decide sample numbers for each cell type
  sampling_num_list <- celltype_reference_data %>%
    dplyr::mutate(., sampling_num = round(sample_scale*Ratio))
  if(any(sampling_num_list$sampling_num == 0L)){
    sampling_num_list$sampling_num[sampling_num_list$sampling_num == 0] <- 1L
  }
  
  # sampling cells
  sample_cells <- c()
  
  for(sc in 1:nrow(sampling_num_list)){
    cur_celltype_sampled <- scRNAseq_annotation_data %>%
      dplyr::filter(., cell_type == sampling_num_list$cell_type[sc])
    set.seed(seed)
    cur_sample_cells <- sample(cur_celltype_sampled$sample_label, size = sampling_num_list$sampling_num[sc], replace = F)
    
    sample_cells <- if(sc == 1L){
      cur_sample_cells
    } else {
      c(sample_cells, cur_sample_cells)
    }
  }# sc
  
  return(sample_cells)
}# sample_cells_for_WF


##############################################################################

#' Calculate Weighting factors
#' Calculate Weighting factors
#' 
#' Argument 
#' woRNAseq_data: A data-frame of whole-organ RNA-seq data
#' scRNAseq_count_data: A data-frame of scRNA-seq data
#' scRNAseq_annotation_data: A data-frame of sample_labels and cell-types' labels
#' signature_genes: A vector of signature genes
#' celltype_reference_data: A data-frame of cell-types labels and the ratio
#' save_names: A save header name
#' constrain_value: A value for constraints of a quadratic problem
#' sample_scale: A integer value to downsample for WF calculation
#' verbose: logical
#' 
#' Argument format
#' * woRNAseq_data: A data-frame with 'gene' column
#' * scRNAseq_count_data: A data-frame with 'gene' column
#' * scRNAseq_annotation_data: A data-frame with 'sample_labels' and 'cell_type' columns
#' * signature_genes: A vector of characters
#' * celltype_reference_data: A data-frame with 'cell_type' and 'Ratio' columns
#' * save_name: A characters for header of saving name
#' * constrain_value: Constraint value (default: 0)
#' * sample_scale: A integer value to downsample for WF calculation (default:100)
#' * verbose: logical
#' 
#' Return csv-files as follows:
#' 
#' * ..._WF_raw.csv: the results of raw WF value
#' * ..._WF_summary.csv: the results of WF summary, which is used for the downstream analysis
#' * ..._pearson_comparison_results.csv: the results of pearson correlation between 'synthetic' and 'real' whole-organ RNA-seq. Ctrl is no WF, while WF is with WF.
#' 
#' 
#' 

library(tidyverse)
library(data.table)
library(magrittr)
library(osqp)

Weighting_factors_calculation <- function(woRNAseq_data, scRNAseq_count_data,
                                          scRNAseq_annotation_data, 
                                          signature_genes,
                                          celltype_reference_data, 
                                          save_names, constrain_value = 0, 
                                          sample_scale = 100L, 
                                          verbose = FALSE){

  ##%%%% data filter by signature genes %%%%##
  filter_scRNAseq_count_data <- scRNAseq_count_data %>%
    dplyr::filter(., gene %in% signature_genes) %>%
    dplyr::arrange(., gene)
  filter_woRNAseq_count_data <- woRNAseq_data %>%
    dplyr::filter(., gene %in% signature_genes) %>%
    dplyr::arrange(., gene) %>%
    tibble::column_to_rownames('gene')
  
  ##%%%% variable conversion for QP solver %%%%##
  y_num <- ncol(filter_woRNAseq_count_data)
  
  for(y_i in 1:length(y_num)){
    cur_y <- filter_woRNAseq_count_data[,y_i] %>%
      as.matrix()
    
    y_all <- if(y_i == 1L){
      cur_y
    } else {
      y_all + cur_y
    }
  }# y_i
  
  y1 <- filter_woRNAseq_count_data[,1] %>%
    as.matrix()
  
  ##%%%% loop for efficiently calculate WFs %%%%##
  # loop
  last_trials_num <- 0L
  
  for(wt_num in 1:10){
    
    cat('\nCalculate WFs in ', length(signature_genes), ' signature genes... ', wt_num, ' round in 10 rounds :: ', as.character(Sys.time()), '\n')
    
    # initialisation i
    i <- 1L
    min_count <- 0L
    
    repeat{
      if(verbose){
        cat(i, '-th :: ')
      }
      
      if(min_count > 10L){
        break
      } else {
        # select count
        seed <- i + 10*wt_num
        sampling_cells_list <- sample_cells_for_WF(scRNAseq_annotation_data = scRNAseq_annotation_data, celltype_reference_data = celltype_reference_data, sample_scale = sample_scale, seed = seed)
        X <- filter_scRNAseq_count_data %>%
          dplyr::select(., c('gene', all_of(sampling_cells_list))) %>%
          dplyr::arrange(., gene) %>%
          tibble::column_to_rownames('gene') %>%
          as.matrix()
        
        stopifnot(all(rownames(y_all) == rownames(X)))
        
        ## prepare q, P, A for QP
        q <- -2*t(X)%*%(dim(X)[2]*y_all)
        P <- 2*y_num*t(X)%*%X
        A <- diag(1, dim(X)[2], dim(X)[2])
        
        set.seed(seed)
        settings <- osqp::osqpSettings(verbose = FALSE)
        WF_result <- solve_osqp(P = P, q = q, A = A, l = rep(constrain_value, dim(X)[2]), u = rep(Inf, dim(X)[2]), settings)
        WF <- WF_result$x %>%
          data.frame(., stringsAsFactors = F, check.names = F, check.rows = F) %>%
          dplyr::rename(., proxy = 1) %>%
          dplyr::mutate(., proxy = ifelse(proxy < 0, constrain_value, proxy))
        colnames(WF) <- paste0('trial_', i + last_trials_num)
        rownames(WF) <- colnames(X)
        
        # reconstruct y1
        reconst_y1_ctrl <- rowSums(X)/dim(X)[2]
        reconst_y1_WF <- X %*% as.matrix(WF)/dim(X)[2]
        
        # calculate pearson coeff
        pearson_res_mtx <- matrix(NA, 1, 3) %>%
          data.frame(., stringsAsFactors = F, check.names = F, check.rows = F)
        colnames(pearson_res_mtx) <- c('trial', 'ctrl', 'WF')
        pearson_res_mtx$trial <- paste0('trial_', i + last_trials_num)
        pearson_res_mtx$ctrl <- cor(y1, reconst_y1_ctrl, method = 'p')
        pearson_res_mtx$WF <- cor(y1, reconst_y1_WF, method = 'p')
        
        pearson_res_cum <- if(i == 1L){
          pearson_res_mtx
        } else {
          bind_rows(pearson_res_cum, pearson_res_mtx)
        }
        
        # WF result cum
        WF <- WF %>%
          data.frame(., stringsAsFactors = F, check.names = F, check.rows = F) %>%
          tibble::rownames_to_column('sample_name') %>%
          tidyr::gather(., trials, WF_value, -sample_name)
        WF_cum <- if(i == 1L){
          WF
        } else {
          bind_rows(WF_cum, WF)
        }
        
        if(length(unique(WF_cum$sample_name)) != nrow(scRNAseq_annotation_data)){
          min_count <- 0L
        } else {
          WF_summary <- WF_cum %>%
            dplyr::group_by(sample_name) %>%
            dplyr::summarise(., count = n())
          min_count <- min(WF_summary$count)
        }

      }# if min_count ...
      
      i <- i + 1L
    }# repeat
    
    # save temporarily data
    if(wt_num == 1L){
      write.csv(pearson_res_cum, paste0(save_names, '_pearson_comparison_results.csv'), row.names = F)
      write.csv(WF_cum, paste0(save_names, '_WF_raw.csv'), row.names = F)
    } else {
      pearson_res_fin <- data.table::fread(paste0(save_names, '_pearson_comparison_results.csv'), header = T, stringsAsFactors = F, check.names = F) %>%
        rbind(., pearson_res_cum)
      write.csv(pearson_res_fin, paste0(save_names, '_pearson_comparison_results.csv'), row.names = F)
      WF_cum_fin <- data.table::fread(paste0(save_names, '_WF_raw.csv'), header = T, stringsAsFactors = F, check.names = F) %>%
        rbind(., WF_cum)
      write.csv(WF_cum_fin, paste0(save_names, '_WF_raw.csv'), row.names = F)
      rm(pearson_res_fin, WF_cum_fin)
      gc()
    }
    
    # renew last trials number
    last_trials_num <- last_trials_num + i - 1L
    
  }# wt_num
  
  # data summarise
  WF_cum_fin <- data.table::fread(paste0(save_names, '_WF_raw.csv'), header = T, stringsAsFactors = F, check.names = F)
  WF_summarised <- WF_cum_fin %>%
    dplyr::group_by(sample_name) %>%
    dplyr::summarise(., numbers = n(), mean = mean(WF_value), variance = var(WF_value)) %>%
    merge(scRNAseq_annotation_data, ., by.x = 'sample_label', by.y = 'sample_name')
  
  write.csv(WF_summarised, paste0(save_names, '_WF_summary.csv'), row.names = F)
}# function

