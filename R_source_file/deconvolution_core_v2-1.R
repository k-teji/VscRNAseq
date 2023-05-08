##' V-scRNAseq deconvolution
##' A function to deconvolve whole-organ RNA-seq and estimate virtual scRNA-seq
##' 
##' Argument
##' wscRNAseq_mean: A data-frame of the mean of weighted scRNA-seq count
##' wscRNAseq_var: A data-frame of the variance of weighted scRNA-seq count
##' reference_celltype_ratio: A vector of cell-type reference ratio data
##' Signature_genes: A vector of signature genes
##' whole_RNAseq: A data-frame of one normalised count in whole-organ RNAseq
##' stan_iter_num: The number of sampling iterations (default:5000)
##' stan_thin_num: The number of thin for accepted sampling (default:1)
##' stan_parameter_list: A list of stan parameters such as adapt_delta, max_treedepth, etc. (default:list(adapt_delta = 0.8, max_treedepth = 15)) 
##' 
##' Argument format
##' * wscRNAseq_mean: A data-frame with 'gene' column
##' * wscRNAseq_var: A data-frame with 'gene' column
##' * reference_celltype_ratio: A vector of 'ratio' ordered by 'cell_type'
##' * Signature_genes: A vector of top N signature gene sets
##' * whole_RNAseq: A data-frame with 'gene' column and only one count data
##' * stan_iter_num: A integer for sampling iteration
##' * stan_thin_num: A integer for sampling thinning
##' * stan_parameter_list: A list of stan parameters such as adapt_delta, max_treedepth, etc.
##' 
##' Return 
##' A list with:
##' * estimated_scRNAseq_mean: estimated weighted scRNA-seq mean count
##' * estimated_scRNAseq_variance: estimated weighted scRNA-seq variance count
##' * estimated_celltype_ratio_mean: estimated cellular mean ratio
##' * estimated_celltype_ratio_variance: estimated cellular variance ratio
##' * other_estimated_scRNAseq_all-results: A list of estimation results in non-signature genes expression patterns
##' * other_estimation_rmse_y-to-gamma_list: the table for gamma selection by RMSE
##' 
##' 

library(tidyverse)
library(rstan)
library(magrittr)
library(doParallel)
library(tibble)
library(dplyr)

## main code

VscRNAseq_deconvolute <- function(wscRNAseq_mean, wscRNAseq_var, reference_celltype_ratio, 
                                  Signature_genes, whole_RNAseq, stan_iter_num = 5000, stan_thin_num = 1, 
                                  stan_parameter_list = list(adapt_delta = 0.8, max_treedepth = 15)){
 
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  ## make result list
  result_list <- list()
  
  ## import rstan model
  rstan_model_input <- '
    data {
      int N; // gene numbers
      int M; // cell type numbers
      vector[N] y; // woRNAseq
      matrix[N,M] x_mu;
      matrix[N,M] x_sigma;
    }
    
    parameters {
      real<lower=0> alpha;
      real<lower=0> beta;
      matrix[N,M] x_raw;
      simplex[M] r; // cell type ratio
    }
    
    transformed parameters{
      matrix<lower=0>[N,M] x; // scRNAseq
      vector[N] y_error;
    
      for(m in 1:M){
        x[,m] = x_mu[,m] + x_sigma[,m] .* x_raw[,m];
        
        for(n in 1:N){
          if(x[n,m] < 0){
            x[n,m] = 0;
          }
        }
      }
      
      y_error = (y - x*r)/beta;
    }
    
    model {
      target += std_normal_lpdf(alpha);
      target += std_normal_lpdf(beta);
      target += std_normal_lpdf(to_vector(x_raw));
      target += dirichlet_lpdf(r | rep_vector(alpha, M));
      target += std_normal_lpdf(y_error);
    }
  '
  rstan_model <- rstan::stan_model(rstan_model_input)

  ## make input data for stan sampling
  filter_wscRNAseq_mean <- wscRNAseq_mean %>%
    dplyr::filter(., gene %in% Signature_genes) %>%
    dplyr::arrange(., gene) %>%
    tibble::column_to_rownames('gene') %>%
    as.matrix()
  filter_wscRNAseq_mean[filter_wscRNAseq_mean == 0] <- 1e-20
  
  filter_wscRNAseq_var <- wscRNAseq_var %>%
    dplyr::filter(., gene %in% Signature_genes) %>%
    dplyr::arrange(., gene) %>%
    tibble::column_to_rownames('gene') %>%
    as.matrix() %>%
    sqrt()
  filter_wscRNAseq_var[filter_wscRNAseq_var == 0] <- 1e-20
  
  stopifnot(all(rownames(filter_wscRNAseq_mean) == rownames(filter_wscRNAseq_var)))
  stopifnot(all(colnames(filter_wscRNAseq_mean) == rownames(reference_celltype_ratio)))
  stopifnot(all(colnames(filter_wscRNAseq_mean) == colnames(filter_wscRNAseq_var)))
  
  filter_woRNAseq <- whole_RNAseq %>%
    dplyr::rename(., y = 2) %>%
    dplyr::filter(., gene %in% Signature_genes) %>%
    dplyr::arrange(., gene) %>%
    tibble::column_to_rownames('gene') %>%
    as.vector(.)
  stopifnot(all(rownames(filter_wscRNAseq_mean) == rownames(filter_woRNAseq)))
  
  ## prep stan running
  input_data <- list(
    N = nrow(filter_wscRNAseq_mean),
    M = ncol(filter_wscRNAseq_mean),
    y = filter_woRNAseq$y,
    x_mu = filter_wscRNAseq_mean,
    x_sigma = filter_wscRNAseq_var
  )
  
  init_data <- function(...){
    list(x_raw = matrix(0, nrow = nrow(filter_wscRNAseq_mean), ncol = ncol(filter_wscRNAseq_mean)),
         r = reference_celltype_ratio$Ratio,
         alpha = 1, beta = 1)
  }
  
  ## stan sampling
  message('stan Sampling for signature genes count and ratio calculation')
  nuts_results <- rstan::sampling(rstan_model, data = input_data, init = init_data, warmup = round(stan_iter_num*0.5), refresh = 100,
                                  iter = stan_iter_num, chains = 3, thin = stan_thin_num, verbose = TRUE, seed = 29,
                                  control = stan_parameter_list)
  message('stan sampling done')
  
  result_list[['stan_raw_results']] <- nuts_results
  
  ## extract estimated results
  est_nuts_results <- rstan::extract(nuts_results)
  
  est_wscRNAseq_mean <- apply(est_nuts_results$x, c(2,3), mean)
  rownames(est_wscRNAseq_mean) <- rownames(filter_wscRNAseq_mean)
  colnames(est_wscRNAseq_mean) <- colnames(filter_wscRNAseq_mean)
  
  est_wscRNAseq_var <- apply(est_nuts_results$x, c(2,3), var)
  rownames(est_wscRNAseq_var) <- rownames(filter_wscRNAseq_mean)
  colnames(est_wscRNAseq_var) <- colnames(filter_wscRNAseq_mean)
  
  est_celltype_ratio_mean <- apply(est_nuts_results$r, 2, mean) %>%
    as.matrix()
  rownames(est_celltype_ratio_mean) <- colnames(filter_wscRNAseq_mean)
  
  est_celltype_ratio_var <- apply(est_nuts_results$r, 2, var) %>%
    as.matrix()
  rownames(est_celltype_ratio_var) <- colnames(filter_wscRNAseq_mean)
  
  ## calculate other
  other_woRNAseq <- whole_RNAseq %>%
    dplyr::filter(., !gene %in% Signature_genes) %>%
    dplyr::arrange(., gene) %>%
    tibble::column_to_rownames('gene')
  other_wscRNAseq_mean <- wscRNAseq_mean %>%
    dplyr::filter(., gene %in% rownames(other_woRNAseq)) %>%
    dplyr::arrange(., gene) %>%
    tibble::column_to_rownames('gene') %>%
    as.matrix()
  other_wscRNAseq_var <- wscRNAseq_var %>%
    dplyr::filter(., gene %in% rownames(other_woRNAseq)) %>%
    dplyr::arrange(., gene) %>%
    tibble::column_to_rownames('gene') %>%
    as.matrix()
  
  stopifnot(all(colnames(other_wscRNAseq_mean) == colnames(other_wscRNAseq_var)))
  stopifnot(all(rownames(other_wscRNAseq_mean) == rownames(other_wscRNAseq_var)))
  stopifnot(all(rownames(other_woRNAseq) == rownames(other_wscRNAseq_mean)))
  stopifnot(all(colnames(other_wscRNAseq_mean) == rownames(est_celltype_ratio_mean)))
  
  gamma_list <- c(10^seq(-5,5,1))
  
  other_estimate_all <- list()
  rmse_y2gamma_list <- NULL
  
  for(g in 1:length(gamma_list)){
    gamma <- gamma_list[g]
    other_estimate <- estimate_scRNAseq(other_woRNAseq = other_woRNAseq, other_wscRNAseq_mean = other_wscRNAseq_mean, other_wscRNAseq_var = other_wscRNAseq_var, est_ratio = est_celltype_ratio_mean, gamma = gamma, init_wscRNAseq_mean = other_wscRNAseq_mean)
    
    other_estimate_all[[g]] <- other_estimate
    
    rmse_y2gamma_list <- if(g == 1L){
      other_estimate[['rmse_y']]
    } else {
      c(rmse_y2gamma_list, other_estimate[['rmse_y']])
    }
  }# gamma_list
  
  best_other_est <- other_estimate_all[[which.min(rmse_y2gamma_list)]]
  
  ## extract other wscRNAseq data
  est_other_wscRNAseq_mean <- best_other_est[['estimated_X_mean']]
  est_other_wscRNAseq_var <- best_other_est[['estimated_X_variance']]
  
  ## combine all genes
  est_wscRNAseq_all_mean <- rbind(est_wscRNAseq_mean, est_other_wscRNAseq_mean)
  est_wscRNAseq_all_var <- rbind(est_wscRNAseq_var, est_other_wscRNAseq_var)
  
  ## data input result list
  result_list[['estimated_scRNAseq_mean']] <- est_wscRNAseq_all_mean
  result_list[['estimated_scRNAseq_variance']] <- est_wscRNAseq_all_var
  result_list[['estimated_celltype_ratio_mean']] <- est_celltype_ratio_mean
  result_list[['estimated_celltype_ratio_variance']] <- est_celltype_ratio_var
  result_list[['other_estimated_scRNAseq_all-results']] <- other_estimate_all
  result_list[['other_estimation_rmse_y-to-gamma_list']] <- rmse_y2gamma_list
  
  return(result_list)
}# VscRNAseq_deconvolute

###################################################################################################################
##' X prediction
##' function: X prediction (Predict gene expression pattern)
##' 
##' Argument
##' y: A vector of whole organ RNA-seq
##' mu: A matrix of scRNAseq mean count
##' gamma: A hyperparameter
##' sigma: A matrix of scRNAseq variance count
##' r: the estimated celltype ratio
##' X: A matrix of scRNAseq mean count predicted in the prior iteration
##' est_wscRNAseq_var: A matrix of scRNAseq variance count predicted in the prior iteration
##' 
##' Return X_all: estimated weighted scRNAseq mean count
##' 
##' 
##' 

cal_est_wscRNAseq_mean <- function(y, mu, gamma, sigma, r, X, est_wscRNAseq_var){
  
  for(j in 1:ncol(X)){
    X_j <- est_wscRNAseq_var[,j]*( (r[j,]/(gamma^2))*y - (r[j,]/(gamma^2))*X[,-j]%*%r[-j,] + (1/sigma[,j])*mu[,j] )
    colnames(X_j) <- colnames(X)[j]
    
    X_all <- if(j == 1L){
      X_j
    } else {
      cbind(X_all, X_j)
    }
    
  }# j
  
  X_all[X_all < 0] <- 0
  
  return(X_all)
  
}# cal_est_wscRNAseq_mean

#######################################################################################
##' Estimate non-signature genes expression pattern
##' Estimate non-signature genes expression pattern and select the best gamma
##' 
##' Argument
##' other_woRNAseq: A data-frame of whole-organ RNA-seq in non-signature genes
##' other_wscRNAseq_mean: A data-frame of weighted scRNA-seq mean count in non-signature genes
##' other_wscRNAseq_var: A data-frame of weighted scRNA-seq variance count in non-signature genes
##' est_ratio: A vector of estimated cellular ratio
##' gamma: A hyperparameter of gamma
##' init_wscRNAseq_mean: A data-frame of initial weighted scRNA-seq mean count in non-signature genes
##' 

library(ggplot2)

estimate_scRNAseq <- function(other_woRNAseq, other_wscRNAseq_mean, 
                              other_wscRNAseq_var, est_ratio, gamma, 
                              init_wscRNAseq_mean){
  
  # data conversion to symbolic name
  y <- other_woRNAseq %>% as.matrix(.)
  mu <- other_wscRNAseq_mean %>% as.matrix(.)  
  sigma <- other_wscRNAseq_var %>% as.matrix(.)
  sigma[sigma == 0] <- 1e-20
  r <- est_ratio %>% as.matrix(.)
  X_init <- init_wscRNAseq_mean %>% as.matrix(.)
  
  # calculate estimated wscRNAseq variance
  for(i in 1:ncol(sigma)){
    sigma_i <- ( gamma^2*sigma[,i] )/( sigma[,i]*(r[i,]^2) + gamma^2 ) %>%
      as.matrix(.)
    colnames(sigma_i) <- colnames(sigma)[i]
    
    est_wscRNAseq_var <- if(i == 1L){
      sigma_i
    } else {
      cbind(est_wscRNAseq_var, sigma_i)
    }
    
  }# i
  
  # calculate X_j by loop
  
  X_old <- X_init
  
  rmse_y_list <- NULL
  
  message('\ngamma : ', gamma)
  for(h in 1:100){
    
    cat(h, '/100 ::')
    
    X_new <- cal_est_wscRNAseq_mean(y = y, mu = mu, gamma = gamma, sigma = sigma, r = r, X = X_old, est_wscRNAseq_var = est_wscRNAseq_var)
    
    vir_y <- X_new %*% r
    
    rmse_y <- sqrt(mean((y - vir_y)^2))
    
    if(h == 1L){
      rmse_y_list <- rmse_y
      names(rmse_y_list) <- c(h)
    } else {
      rmse_y_list <- c(rmse_y_list, rmse_y)
      names(rmse_y_list) <- c(1:h)
    }
    
    if(rmse_y < 1){
      break
    }
    
    X_old <- X_new
    
  }# h
  
  rmse_y_list_df <- data.frame(value = rmse_y_list) %>%
    dplyr::mutate(., No = as.integer(names(rmse_y_list)))
  
  rmse_plot <- ggplot(rmse_y_list_df, aes(x = No, y = value)) +
    geom_line() +
    theme_bw()
  
  ## output data list
  output <- list()
  output[['rmse_y']] <- rmse_y
  output[['estimated_X_mean']] <- X_new
  output[['estimated_X_variance']] <- est_wscRNAseq_var
  output[['rmse_y_plot']] <- rmse_plot
  output[['gamma']] <- gamma
  
  return(output)
}# estimate_scRNAseq