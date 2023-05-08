#' Extract Signature genes
#' Extract signature genes by Random Forest
#' 
#' Argument
#' scRNAseq_count: A data-frame of scRNA-seq normalised count data
#' scRNAseq_annotation: A data-frame of sample_labels and cell-types labels
#' celltype_ref_ratio: A data-frame of cell-types' label and the reference ratio
#' seed: seed
#' 
#' Argument format
#' * scRNAseq_count: A data-frame with 'gene' column
#' * scRNAseq_annotation: A data-frame with 'sample_label' and 'cell_type' columns
#' * celltype_ref_ratio: A data-frame with 'cell_type' and 'Ratio' columns
#' * seed: A integer value
#' 
#' Return 
#' A list with values:
#' * seed: The used value for seed
#' * training_data: The data used for training
#' * test_data: The data used for testing
#' * randomForest_tuning: The list object of the tuning result
#' * randomForest_model: The list object of the generated ranodomForest model
#' * F1_score: The value of F1-score in the generated model
#' * signature_genes: A data-frame of signature genes list with feature importance (MeanDecreaseGini)
#' 
#' 

library(tidyverse)
library(magrittr)
library(tictoc)
library(randomForest)
library(MLmetrics)

##--- function define ---##
extract_signature_genes <- function(scRNAseq_count, scRNAseq_annotation, celltype_ref_ratio = NULL, seed = 29){
  
  ## calculate cell type reference ratio
  if(is.null(celltype_ref_ratio)){
    celltype_ref_ratio <- scRNAseq_annotation %>%
      dplyr::group_by(cell_type) %>%
      dplyr::summarise(., count = n()) %>%
      dplyr::mutate(., Ratio = count/sum(count))
  }
  
  ## Make sample-name labeled scRNAseq count data
  annotated_scRNAseq_count <- scRNAseq_count %>%
    tibble::column_to_rownames('gene') %>%
    t( ) %>%
    data.frame(., stringsAsFactors = F, check.names = F, check.rows = F) %>%
    tibble::rownames_to_column('sample_label') %>%
    base::merge(dplyr::select(scRNAseq_annotation, c(sample_label, cell_type)), ., by = 'sample_label')
  
  # List cell types based on scRNAseq annotation
  celltype.list <- unique(annotated_scRNAseq_count$cell_type)
  
  # Divide scRNAseq data into two: train data and test data
  # Perform division step at every cell type
  for(i in seq_len(length(celltype.list))){
    # Set current cell type
    cur.celltype <- celltype.list[i]
    
    # Filter annotation by current cell type
    cur_scRNAseq <- annotated_scRNAseq_count %>%
      dplyr::filter(., annotated_scRNAseq_count$cell_type == cur.celltype)
    
    # Split training and test data at random at 8:2
    set.seed(seed)
    cur_train <- cur_scRNAseq %>%
      dplyr::sample_frac(size = 0.8, replace = F)
    cur_test <- cur_scRNAseq %>%
      dplyr::filter(., !sample_label %in% cur_train$sample_label)

    # Deposit train data to finally integrate all cell types
    train.data <- if(i == 1){
      cur_train
    } else {
      rbind(train.data, cur_train)
    }
    
    # Deposit test data to finally integrate all cell types
    test.data <- if(i == 1){
      cur_test
    } else {
      rbind(test.data, cur_test)
    }
    
  }#i
  
  ## data tidy up
  train.data <- train.data %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames('sample_label') %>%
    dplyr::mutate(., cell_type = factor(cell_type, levels = celltype.list))
  
  test.data <- test.data %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames('sample_label') %>%
    dplyr::mutate(., cell_type = factor(cell_type, celltype.list))
  
  ## RF parameter tuning
  message('RF tuning')
  tic()
  set.seed(seed = seed)
  tuning_RF <- tuneRF(x = train.data[,-1], y = train.data[,1], stepFactor = 2)
  toc()
  
  tuning_RF <- data.frame(tuning_RF, stringsAsFactors = F, check.names = F, check.rows = F)
  mtry_best <- tuning_RF$mtry[which.min(tuning_RF$OOBError)]
  
  ## make RF model
  message('Make RF model')
  tic()
  set.seed(seed = seed)
  rf_result <- randomForest(x = train.data[,-1], y = train.data[,1], mtry = mtry_best, ntree = 500, importance = T)
  toc()

  ## check test data
  test_rf <- predict(rf_result, test.data)
  test.classifiction <- test_rf %>%
    data.frame(stringsAsFactors = F, check.names = F)
  colnames(test.classifiction) <- c("Rf_class")
  test.classifiction %<>% tibble::rownames_to_column("sample_label")
  test.real.class <- test.data %>%
    tibble::rownames_to_column("sample_label") %>%
    dplyr::select(., c(sample_label, cell_type))
  colnames(test.real.class)[2] <- c("Real_class") 
  combine.test.data <- merge(test.real.class, test.classifiction, by = "sample_label")
  
  ## calculate F1 score
  F1_score <- MLmetrics::F1_Score(y_true = combine.test.data$Real_class, y_pred = combine.test.data$Rf_class)
  
  ## Extract feature importance on MeanDecreaseGini
  signature_genes <- importance(rf_result) %>%
    data.frame(., stringsAsFactors = F, check.names = F, check.rows = F) %>%
    tibble::rownames_to_column('gene') %>%
    dplyr::select(., c(gene, MeanDecreaseGini)) %>%
    dplyr::filter(., MeanDecreaseGini != 0) %>%
    dplyr::arrange(., desc(MeanDecreaseGini))
  
  ## output data
  output_result <- list()
  output_result[['seed']] <- seed
  output_result[['training_data']] <- train.data
  output_result[['test_data']] <- test.data
  output_result[['randomForest_tuning']] <- tuning_RF
  output_result[['randomForest_model']] <- rf_result
  output_result[['F1_score']] <- F1_score
  output_result[['signature_genes']] <- signature_genes
  
  return(output_result)
}# pickup_signature_genes
