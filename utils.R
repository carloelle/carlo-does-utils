library(DropletUtils)
library(Seurat)
library(Giotto)
library(Matrix)

#using DropUtils, we create a .h5 files from Read10X folder outs. 
#Expected: barcodes.tsv.gz, features.tsv.gz matrix.mtx.gz

filter_matrix<-Read10X('/yourdirectory/filtered_feature_bc_matrix/')
setwd('/yourdirectory/filtered_feature_bc_matrix/')
write10xCounts("./filtered_feature_bc_matrix.h5", filter_matrix, type = "HDF5",
               genome = "hm38", version = "3", overwrite = TRUE,
               gene.id = rownames(filter_matrix),
               gene.symbol = rownames(filter_matrix))


#very simply calculate percentages for C2L results
process_C2L<-function(seurat_object){
  mt <- seurat_object@meta.data
  cell_type_columns <- mt[, grep('*_C2L', colnames(mt))]
  row_totals <- rowSums(cell_type_columns)
  cell_type_percentages <- sweep(cell_type_columns, 1, row_totals, "/") * 100
  colnames(cell_type_percentages) <- paste0(colnames(cell_type_percentages), "_percentage")
  mt <- cbind(mt, cell_type_percentages)
  seurat_object@meta.data <- mt
  return(seurat_object)
}

#Do a PAGE Enrichment analysis using the Giotto package, calculate p-value and store enrichment and a T/F column inside the metadata
#internally it select a sample-specific negative control

process_PAGEanalysis <- function(seurat_object, signatures_all, selected_signatures,only_fibro=T) {
  if(only_fibro==T){signatures_all <- signatures_all[!names(signatures_all) %in% c("CCM", "EMRM")]}
  raw_exprs <- seurat_object@assays$Spatial@layers$counts
  rownames(raw_exprs) <- rownames(seurat_object@assays$Spatial@features@.Data)
  colnames(raw_exprs) <- rownames(seurat_object@meta.data)
  spatial_locs <- as.data.table(GetTissueCoordinates(seurat_object)[,1:2])
  colnames(spatial_locs) <- c("x", "y")
  
  myGiottoObj <- createGiottoObject(raw_exprs = raw_exprs, spatial_locs = spatial_locs)
  myGiottoObj <- normalizeGiotto(gobject = myGiottoObj)
  
  # Create signature matrix for initial PAGE analysis
  all_signatures <- names(signatures_all)
  signature_matrix_complete <- makeSignMatrixPAGE(
    sign_names = all_signatures,
    sign_list = lapply(all_signatures, function(sign) signatures_all[[sign]])
  )
  
  # Run initial PAGE enrichment analysis on all signatures
  myGiottoObj_initial <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_complete,
    min_overlap_genes = 2,
    output_enrichment = c("original", "zscore")
  )
  
  
  # Convert initial PAGE z-score results to data frame
  initial_zscore_df <- as.data.frame(myGiottoObj_initial@spatial_enrichment$PAGE)
  
  
  # Select a negative control signature (most different from selected signatures)
  cosine_similarity <- function(a, b) {
    sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
  }
  
  
  selected_zscores <- initial_zscore_df[, selected_signatures]
  similarities <- sapply(all_signatures, function(sign) {
    if (sign %in% selected_signatures) return(NA)
    sign_zscores <- initial_zscore_df[, sign]
    mean(apply(selected_zscores, 1, function(row) cosine_similarity(row, sign_zscores)))
  })
  
  
  negative_control <- all_signatures[which.min(similarities)]
  
  
  # Integrate the negative control into the signature matrix
  final_signatures <- c(selected_signatures, negative_control)
  signature_matrix_with_control <- makeSignMatrixPAGE(
    sign_names = final_signatures,
    sign_list = lapply(final_signatures, function(sign) signatures_all[[sign]])
  )
  
  # Run PAGE enrichment analysis with the negative control
  myGiottoObj_zscore_with_control <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_with_control,
    min_overlap_genes = 2,
    output_enrichment = c("original", "zscore")
  )
  
  # Convert PAGE z-score results with control to data frame and store in Seurat metadata
  zscore_with_control_df <- as.data.frame(myGiottoObj_zscore_with_control@spatial_enrichment$PAGE)
  colnames(zscore_with_control_df) <- paste0(colnames(zscore_with_control_df), "PAGE_w/NegCtrl")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, zscore_with_control_df)
  
  # Run PAGE enrichment analysis with p-values
  myGiottoObj_pval <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_with_control,
    min_overlap_genes = 2,
    p_value = TRUE,
    reverse_log_scale = FALSE
  )
  
  # Extract p-value results
  pval_df <- as.data.frame(myGiottoObj_pval@spatial_enrichment$PAGE)
  
  # Identify significantly enriched spots (p-value < 0.05)
  significant_spots <- pval_df < 0.05
  
  # Store significant enrichment information in Seurat metadata
  significant_df <- as.data.frame(significant_spots)
  colnames(significant_df) <- paste0(colnames(significant_df), "PAGE_significant")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, significant_df)
  
  return(seurat_object)
  
}





