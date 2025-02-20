##### Use harmony to integrate gene expression data for computing GEX dissimilarity
##### between across people clones with same TCR sequence

library(Seurat)
library(stringr)
library(tibble)
folder_path = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis"
file_list = list.files(folder_path)
file_GSE = file_list[str_detect(file_list,"GSE|EMTAB9357_pbmc|PRJCA002413_pbmc|10X|Zhang_2020")]
file_matrix = list.files(paste0(folder_path,"/", file_GSE, "/matrix"), full.names = TRUE)

across_sam_pair <- readRDS("~/projects/Research/TCR and Gene Expression/PC distance analysis/across_sam_pair.rds")
sam1_count = sam2_count = c()
for (i in 1:nrow(across_sam_pair)) {
  sam1_count[i] = file_matrix[grep(paste0("(^|/|_)", across_sam_pair$sam1[i], ".rds"), file_matrix)]
  sam2_count[i] = file_matrix[grep(paste0("(^|/|_)", across_sam_pair$sam2[i], ".rds"), file_matrix)]
}
across_sam_pair$sam1_count = sam1_count
across_sam_pair$sam2_count = sam2_count

#nucleotide
across_sam_pair <- readRDS("~/projects/Research/TCR and Gene Expression/PC distance analysis/across_sam_pair_nt.rds")
sam1_count = sam2_count = c()
for (i in 1:nrow(across_sam_pair)) {
  sam1_count[i] = file_matrix[grep(paste0("(^|/|_)", across_sam_pair$sam1[i], ".rds"), file_matrix)]
  sam2_count[i] = file_matrix[grep(paste0("(^|/|_)", across_sam_pair$sam2[i], ".rds"), file_matrix)]
}
across_sam_pair$sam1_count = sam1_count
across_sam_pair$sam2_count = sam2_count

i=1
for (i in 1:nrow(across_sam_pair)){
  data1 = CreateSeuratObject(counts = readRDS(file = across_sam_pair$sam1_count[i]), 
                             project=across_sam_pair$sam1[i], min.cells = 3, min.features = 200)
  data1@meta.data$nCount_RNA = colSums(x = data1, slot = "counts")  # nCount_RNA
  data1@meta.data$nFeature_RNA = colSums(x = GetAssayData(object = data1, slot = "counts") > 0)
  data2 = CreateSeuratObject(counts = readRDS(file = across_sam_pair$sam2_count[i]), 
                             project=across_sam_pair$sam2[i], min.cells = 3, min.features = 200)
  data2@meta.data$nCount_RNA = colSums(x = data2, slot = "counts")  # nCount_RNA
  data2@meta.data$nFeature_RNA = colSums(x = GetAssayData(object = data2, slot = "counts") > 0)
  
  combined <- merge(data1, y = data2, add.cell.ids = c(across_sam_pair$sam1[i], across_sam_pair$sam2[i]))
  combined$sample = combined@meta.data[["orig.ident"]]
  combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
  combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
  combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
  combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(combined)
  combined <- ScaleData(combined, features = all.genes)
  combined <- RunPCA(combined, features = VariableFeatures(object = combined))
  
  combined <- IntegrateLayers(
    object = combined, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    verbose = T
  )
  harmony = combined@reductions[["harmony"]]@cell.embeddings[,1:10] %>%
    as.data.frame() %>%
    rownames_to_column(var="sample_barcode")
  
  saveRDS(harmony, 
          paste0("~/projects/Research/TCR and Gene Expression/PC distance analysis/harmony/har_",
                 across_sam_pair$sam1[i], "_", across_sam_pair$sam2[i], ".rds"))
  remove(combined)
  remove(data1)
  remove(data2)
  remove(harmony)
  gc()
}
