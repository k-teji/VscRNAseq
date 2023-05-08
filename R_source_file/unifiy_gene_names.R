#' Gene names unification between whole-organ and scRNA-seq
#' 
#' Unify gene names by gene alias and normalise count data by the designated scale 
#' 
#' Argument
#' woRNAseq: A data-frame of whole-organ RNA-seq or bulk RNA-seq count data
#' scRNAseq: A data-frame of single cell RNA-seq count data
#' gene_converter: A data-frame of a converter of gene names from old name to new name by gene aliases
#' remove_genes: A vector of genes removed for the analysis
#' normalise_scale: A value of noramlisation scale for count data
#' 
#' Argument format 
#' * woRNAseq : A data-frame with 'gene' column
#' * scRNAseq : A data-frame with 'gene' column
#' * gene_converter : A data-frame of gene symbols and the corresponding entrezID (See the example in data)
#' * remove_genes : NULL or A vector of genes (e.g., c('Rs45s', 'Akap5', 'Lrrc17'))
#' * normalise_scale : A numeric value
#' 
#' Return 
#' A list with values:
#' * scRNAseq_raw: Single cell RNA-seq input count data
#' * scRNAseq_normalised: Single cell RNA-seq normalised count data
#' * woRNAseq_raw: Whole-organ RNA-seq input coun data
#' * woRNAseq_normalised: Whole-organ RNA-seq normalised count data
#' 

library(tidyverse)
library(magrittr)

unify_gene_names <- function(woRNAseq, scRNAseq, gene_converter, 
                             remove_genes = NULL, normalise_scale = 1e+6){
  
  ## make converter from scRNAseq gene names to woRNAseq gene names
  woRNAseq_gene <- woRNAseq %>%
    dplyr::select(., all_of(c('gene'))) %>%
    merge(., gene_converter, by = 'gene') %>%
    dplyr::rename(., woRNAseq = gene)
  scRNAseq_gene <- scRNAseq %>%
    dplyr::select(., all_of(c('gene'))) %>%
    merge(., gene_converter, by = 'gene') %>%
    dplyr::rename(., scRNAseq = gene)
  merge_gene_conv <- merge(scRNAseq_gene, woRNAseq_gene, by = 'id_master') %>%
    dplyr::select(., -all_of(c('id_master')))
  
  ## gene filtering
  scRNAseq_filter <- scRNAseq %>%
    merge(merge_gene_conv, ., by.x = 'scRNAseq', by.y = 'gene') %>%
    dplyr::select(., -all_of(c('scRNAseq'))) %>%
    dplyr::rename(., gene = woRNAseq) %>%
    dplyr::filter(., !gene %in% remove_genes)
  woRNAseq_filter <- woRNAseq %>%
    dplyr::filter(., gene %in% scRNAseq_filter$gene) %>%
    dplyr::filter(., !gene %in% remove_genes)
  
  ## sum up counts which the duplicated genes have
  scDup_genes <- scRNAseq_filter %>%
    dplyr::select(., all_of(c('gene'))) %>%
    dplyr::group_by(gene) %>%
    dplyr::filter(., n() > 1)
  if(nrow(scDup_genes) > 0){
    scDup_genes_data <- scRNAseq_filter %>%
      dplyr::group_by(gene) %>%
      dplyr::filter(., n() > 1) %>%
      dplyr::ungroup() %>%
      tidyr::gather(., sample, count, -gene) %>%
      dplyr::group_by(gene, sample) %>%
      dplyr::summarise(., value = sum(count)) %>%
      tidyr::spread(., sample, value)
    
    scRNAseq_filter <- scRNAseq_filter %>%
      dplyr::group_by(gene) %>%
      dplyr::filter(., n() == 1) %>%
      dplyr::bind_rows(., scDup_genes_data) %>%
      dplyr::arrange(., gene)
  }
  
  woRNAseq_genes <- woRNAseq_filter %>%
    dplyr::select(., all_of(c('gene'))) %>%
    dplyr::group_by(gene) %>%
    dplyr::filter(., n() > 1)
  
  if(nrow(woRNAseq_genes) > 0){
    woDup_data <- woRNAseq_filter %>%
      dplyr::group_by(gene) %>%
      dplyr::filter(., n() > 1) %>%
      dplyr::ungroup() %>%
      tidyr::gather(., sample, count, -gene) %>%
      dplyr::group_by(gene, sample) %>%
      dplyr::summarise(., value = sum(count)) %>%
      tidyr::spread(., sample, value)
    woRNAseq_filter <- woRNAseq_filter %>%
      dplyr::group_by(gene) %>%
      dplyr::filter(., n() == 1) %>%
      dplyr::bind_rows(., woDup_data) %>%
      dplyr::arrange(., gene)
  }
  
  ## normalisation
  scRNAseq_norm <- scRNAseq_filter %>%
    tibble::column_to_rownames('gene') %>%
    apply(., 2, function(x){
      normalise_scale*x/sum(x)
    }) %>%
    data.frame(., stringsAsFactors = F, check.names = F, check.rows = F) %>%
    tibble::rownames_to_column('gene')
  
  woRNAseq_norm <- woRNAseq_filter %>%
    tibble::column_to_rownames('gene') %>%
    apply(., 2, function(x){
      normalise_scale*x/sum(x)
    }) %>%
    data.frame(., stringsAsFactors = F, check.names = F, check.rows = F) %>%
    tibble::rownames_to_column('gene')
  
  result_list <- list()
  
  result_list[['scRNAseq_raw']] <- scRNAseq_filter
  result_list[['scRNAseq_normalised']] <- scRNAseq_norm
  result_list[['woRNAseq_raw']] <- woRNAseq_filter
  result_list[['woRNAseq_normalised']] <- woRNAseq_norm
  
  return(result_list)
}# unify_gene_names
