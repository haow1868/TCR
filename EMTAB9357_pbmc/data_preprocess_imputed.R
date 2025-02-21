
countoriname <- readRDS("~/projects/Research/TCR and Gene Expression/PC distance analysis/EMTAB9357_pbmc/proc/countoriname.rds")
barcode = colnames(countoriname)
library(stringr)
sample = sort(unique(str_split(barcode, ":", simplify = TRUE)[,3]))
proc <- readRDS("~/projects/Research/TCR and Gene Expression/3.24/proc.rds")
sample_tcr = sort(unique(proc$sample))

# get the gene expression matrix; do not run
for (i in 1:length(sample)) {
  bc = barcode[endsWith(barcode, paste0(":",sample[i]))]
  if(unique(str_split(bc, ":", simplify = TRUE)[,3]) != sample[i]){
    stop("multiple samples")
  }
  mt = countoriname[,colnames(countoriname) %in% bc]
  saveRDS(mt, paste0("~/projects/Research/TCR and Gene Expression/PC distance analysis/EMTAB9357_pbmc/matrix/mtx_",
                     sample[i],".rds"))
  remove(bc)
  remove(mt)
  gc()
}

library(dplyr)
library(tidyverse)
proc <- readRDS("~/projects/Research/TCR and Gene Expression/3.24/proc.rds")
sample = sort(unique(proc$sample))
proc_single = proc %>%
  filter(chain_pairing == "Single pair") %>%
  select(sample, TRA_1_cdr3, TRB_1_cdr3, TRA_1_cdr3_nt, TRB_1_cdr3_nt) 

# gene expression matrix
file_matrix = list.files("~/projects/Research/TCR and Gene Expression/PC distance analysis/EMTAB9357_pbmc/matrix",
                         full.names = TRUE)
library(Seurat)
library(R.utils)
library(SAVER)
df = vector("list", length = length(file_matrix))
for (i in seq_along(file_matrix)) {
  seurat.data = readRDS(file_matrix[i])
  # Initialize the Seurat object with the raw (non-normalized data).
  seurat <- CreateSeuratObject(counts = seurat.data)
  
  # Standard pre-processing workflow
  # QC and selecting cells for further analysis
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  
  # Remove nfeature < 200 and > 2500, remove percent.mt > 20
  seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
  seurat = subset(seurat, cells = rownames(proc_single))
  
  if(ncol(seurat) < 500){
    next
  }
  
  # Normalize the data
  seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Identification of highly variable features (feature selection), select 2000 features
  seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
  
  # Scaling the data
  all.genes <- rownames(seurat)
  seurat <- ScaleData(seurat, features = all.genes)
  
  # Perform linear dimensional reduction
  seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
  
  # Extract the PC values and add predicted cell subtypes
  cell_PC = as.data.frame(seurat@reductions[["pca"]]@cell.embeddings)[1:10] 
  
  data = merge(proc_single, cell_PC, by="row.names") %>% 
    column_to_rownames(var="Row.names")
  
  exp = t(as.data.frame(seurat@assays[["RNA"]]@data[which(seurat@assays[["RNA"]]@data@Dimnames[[1]] %in% c("CD4", "CD8A","CD8B")), ]))
  data = merge(data, exp, by="row.names") %>% 
    column_to_rownames(var="Row.names")
  
  # imputed
  genes.idx = which(rownames(seurat@assays[["RNA"]]@counts) %in% c("CD4", "CD8A","CD8B"))
  seurat_imputed = saver(seurat@assays[["RNA"]]@counts, estimates.only = TRUE, ncores = 23, 
                       pred.genes = genes.idx, pred.genes.only = T)
  exp_imp = t(seurat_imputed)
  colnames(exp_imp) = paste0(colnames(exp_imp),"_imp")
  data = merge(data, exp_imp, by="row.names")
  
  colnames(data)[which(colnames(data) == "Row.names")] = "barcode"
  
  # connect with TCR sequence
  df[[i]] = data
  
  remove(seurat)
  remove(cell_PC)
  gc()
  print("done")
}
names(df) = sample

df_immu = bind_rows(df)
colnames(df_immu)[which(colnames(df_immu) == "TRA_1_cdr3")] = "cdr3_TRA"
colnames(df_immu)[which(colnames(df_immu) == "TRB_1_cdr3")] = "cdr3_TRB"
colnames(df_immu)[which(colnames(df_immu) == "TRA_1_cdr3_nt")] = "cdr3_nt_TRA"
colnames(df_immu)[which(colnames(df_immu) == "TRB_1_cdr3_nt")] = "cdr3_nt_TRB"

df_immu = df_immu %>%
  group_by(cdr3_TRA, cdr3_TRB) %>%
  mutate(ID = cur_group_id()) %>%
  ungroup()
df_immu$sample_ID = paste0(df_immu$sample, "-", df_immu$ID)

saveRDS(df_immu, file = "~/projects/Research/TCR and Gene Expression/PC distance analysis/EMTAB9357_pbmc/df_immu_imputed.rds")


subtype_annotation = function(data){
  CD4 = data$CD4
  CD8A = data$CD8A
  CD8B = data$CD8B
  
  CD4_imp = data$CD4_imp
  CD8A_imp = data$CD8A_imp
  CD8B_imp = data$CD8B_imp
  
  subtype = c()
  for (i in seq_along(CD4)) {
    if (CD4[i] == 0 & CD8A[i] == 0 & CD8B[i] == 0){
      CD4[i] = CD4_imp[i]
      CD8A[i] = CD8A_imp[i]
      CD8B[i] = CD8B_imp[i]
    }
    subtype[i] = if_else( CD4[i] > max(CD8A[i], CD8B[i]), "CD4+ T cells", "CD8+ T cells")
  }
  
  return(subtype)
}

subtype_annotation_imp = function(data){
  CD4_imp = data$CD4_imp
  CD8A_imp = data$CD8A_imp
  CD8B_imp = data$CD8B_imp
  
  subtype = c()
  for (i in seq_along(CD4_imp)) {
    subtype[i] = if_else( CD4_imp[i] > max(CD8A_imp[i], CD8B_imp[i]), "CD4+ T cells", "CD8+ T cells")
  }
  return(subtype)
}



df_immu_imputed$subtype = subtype_annotation(df_immu_imputed)
df_immu_imputed$subtype_imp = subtype_annotation_imp(df_immu_imputed)



df_immu_imputed$CD8 = if_else(df_immu_imputed$CD8A>df_immu_imputed$CD8B, df_immu_imputed$CD8A, df_immu_imputed$CD8B)


# obtain sample information
sample_information = proc %>% 
  select(sample, age, sex, tissue, disease_status, disease_severity) %>%
  unique()
rownames(sample_information) = 1:nrow(sample_information)
write_csv(sample_information, file = "~/projects/Research/TCR and Gene Expression/PC distance analysis/EMTAB9357_pbmc/sample_information_EMTAB9357_pbmc.csv")






















