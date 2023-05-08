# VscRNAseq

# 'Virtual single-cell RNA-seq'

The VscRNAseq is a cellular deconvolution with cell-type specific weighting factors, resulting in cellular ratio and gene expression patterns.

### Package Requirement

You require the below packages for processing VscRNAseq.

-   tidyverse
-   tibble
-   magrittr
-   tictoc
-   randomForest
-   MLmetrics
-   osqp
-   rstan
-   doParallel

### Instruction for VscRNAseq

### 5 steps for VscRNAseq

1.  Preprocessing
2.  Extract Signature gene sets
3.  Calculate Weighting factors
4.  Add Weighting factors to scRNA-seq reference data
5.  Estimate cellular ratio and gene expression patterns by deconvolution

### Step0 Prepare data etc.

-   Import the VscRNAseq R codes

``` r
source('unify_gene_names.R')
source('extract_signature_genes_by_RF.R')
source('calculate_WFs.R')
source('Add_WF_to_scRNAseq.R')
source('deconvolution_core_v2-1.R')
```

-   Prepare whole-organ RNA-seq count data whose gene column-name is 'gene'

-   Prepare scRNA-seq count data whose gene column-name is 'gene'

-   Prepare scRNA-seq annotation table at least containing 'sample_label' and 'cell_type' columns

    ('Sample_label' is cell sample ID and 'cell_type' is as it is.)

-   Prepare cell-type reference ratio with 'cell_type' and 'Ratio' columns

### Step1 Preprocessing

-   Match gene symbols between whole-organ RNA-seq and scRNA-seq count data

    (If your data is from mice, you can use gene_converter deposited in data.)

-   Normalize count data by a certain scale (e.g. 1e+6)

-   You can use a function 'unify_gene_names' (See below)

    (Sample data: whole-organ Heart RNA-seq from i-organs atlas and Heart scRNA-seq from Tabula Muris)

``` r
## Heart health WTs in i-organs atlas
load('data/woRNAseq.rda')

## Heart SmartSeq2 data in Tabula Muris
load('data/scRNAseq.rda')

## gene converteing table from old to new names for mouse
load('data/gene_converter.rda')

## Unify gene names between whole-organ and single cell RNA-seq
unify_result <- unify_gene_names(woRNAseq = woRNAseq, scRNAseq = scRNAseq,
                                 gene_converter = gene_converter, 
                                 remove_genes = c('Rs45s', 'Akap5', 'Lrrc17'),
                                 normalise_scale = 1e+6)
```

### Step2 Extract Signature gene sets

``` r
## scRNAseq annotation
load('scRNAseq_annotation.rda')

## Cellular reference ratio
load('celltype_ref_ratio.rda')

## Extract Signature genes
sigG_result <- extract_signature_genes(scRNAseq_count = unify_result[['scRNAseq_normalised']], 
                                       scRNAseq_annotation = scRNAseq_annotation, 
                                       celltype_ref_ratio = celltype_ref_ratio)
```

### Step3 Calculate Weighting factors

-   Select the number of signature genes extracted at Step2 (In the example below, select 300)
-   Calculate Weighting factors
-   The results are automatically outputted to the designated directory as csv-formatted files

``` r
sigGenes <- sigG_result[['signature_genes']] %>%
  dplyr::top_n(x = ., n = 300, wt = MeanDecreaseGini)
sigGenes <- sigGenes$gene
celltype_ref_ratio <- celltype_ref_ratio %>%
  dplyr::select(., all_of(c('cell_type', 'Ratio')))
Weighting_factors_calculation(woRNAseq_data = unify_result[['woRNAseq_normalised']],
                             scRNAseq_count_data = unify_result[['scRNAseq_normalised']],
                             scRNAseq_annotation_data = scRNAseq_annotation,
                             signature_genes = sigGenes,
                             celltype_reference_data = celltype_ref_ratio,
                             save_name = 'xxx',
                             constraint_value = 0,
                             sample_scale = 100L)
```

### Step4 Add Weighting factors to scRNA-seq normalised count data

``` r
WF_summary <- read.csv('xxx_WF_summry.csv', header = T, stringsAsFactors = F, check.names = F)
WF_scRNAseq <- Add_WF_to_scRNAseq(scRNAseq_raw_count = unify_result[['scRNAseq_normalised]],
                                  WF_summary = WF_summary)
```

### Step5 Estimate cellular ratio and gene expression patterns by deconvolution

-   If the deconvolution sampling is not converged well, change the stan parameters, such as stan_iter_num, adapt_delta, max_treedepth, etc.

``` r
 wscRNAseq_mean <- WF_scRNAseq[['Weighted_mean']] ## WF_scRNAseq: Add_WF_to_scRNAseq results
 wscRNAseq_var <- WF_scRNAseq[['Weighted_variance']] ## WF_scRNAseq: Add_WF_to_scRNAseq results
 reference_celltype_ratio <- celltype_ref_ratio %>%
                             tibble::column_to_rownames('cell_type') %>%
                             as.vector()
 whole_RNAseq <- unify_result[['woRNAseq_normalised']] %>%
   dplyr::select(., all_of(c('gene', 'Normalised_WT2')))
 deconvolution_res <- VscRNAseq_deconvolute(wscRNAseq_mean = wscRNAseq_mean, 
                                            wscRNAseq_var = wscRNAseq_var,
                                            reference_celltype_ratio = reference_celltype_ratio, 
                                            Signature_genes = SigGenes,
                                            whole_RNAseq = whole_RNAseq,
                                            stan_iter_num = 3000,
                                            stan_parameter_list = list(adapt_delta = 0.8, max_treedepth = 15))
```

-   Cellular ratio and gene expression patterns are outputted in the list.
