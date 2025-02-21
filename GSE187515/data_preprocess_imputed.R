

setwd("/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis/GSE187515")
files = list.files("target")
mtx_files_tar = files[str_detect(files, "matrix.tar")]

for (i in seq_along(mtx_files_tar)) {
  untar(paste0("target/", mtx_files_tar[i]), exdir = "target")
}

files = list.files("target")
tcr_files = files[str_detect(files, "filtered_contig_igblast_db-pass_parse-select.tab")]
file_matrix = files[endsWith(files, "matrix")]
sample = str_split(file_matrix, "matrix", simplify = T)[,1]



for (i in seq_along(file_matrix)) {
  seurat.data = Read10X(data.dir = paste0("target/", file_matrix[i]))
  
  if (length(seurat.data) == 2){
    seurat.data = seurat.data[["Gene Expression"]]
  }
  saveRDS(seurat.data, paste0("/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis/GSE187515/matrix/", sample[i], ".rds"))
  remove(seurat.data)
  gc()
  print("done")
}


for (i in seq_along(file_matrix)) {
  contig_anno = read.table(paste0("target/", tcr_files[i]), header = T, sep = "\t", fill = T)
  
  saveRDS(contig_anno, paste0("~/projects/Research/TCR and Gene Expression/PC distance analysis/GSE187515/tcr/", sample[i],".rds"))
  print("done")
  gc()
}


library(Seurat)
library(R.utils)
library(tidyverse)
library(SAVER)
df = vector("list", length = length(file_matrix))

for (i in seq_along(file_matrix)) {
  seurat.data = Read10X(data.dir = paste0("target/", file_matrix[i]))
  
  if (length(seurat.data) == 2){
    seurat.data = seurat.data[["Gene Expression"]]
  }
  
  seurat <- CreateSeuratObject(counts = seurat.data, 
                               min.cells = 3, min.features = 200)
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
  
  contig_anno = read.table(paste0("target/", tcr_files[i]), header = T, sep = "\t", fill = T)
  contig_anno = contig_anno[contig_anno$FUNCTIONAL == "TRUE" &
                              (contig_anno$LOCUS=="TRA" | contig_anno$LOCUS=="TRB"),] %>%
    select(CELL, LOCUS, JUNCTION_10X, JUNCTION_10X_AA) %>%
    filter(JUNCTION_10X != "", JUNCTION_10X_AA != "")
  colnames(contig_anno) = c("barcode", "chain", "cdr3_nt", "cdr3")
  
  # get the single-pair cells
  single_pair_idx = contig_anno %>% 
    dplyr::count(barcode, chain) %>% 
    filter(n==1) %>% 
    dplyr::count(barcode) %>% 
    filter(n==2) %>% 
    pull(barcode) %>% unique()
  contig_anno = contig_anno %>% filter(barcode %in% single_pair_idx)
  
  contig_anno_wide = contig_anno %>% 
    pivot_wider(
      names_from = chain,
      values_from = c(cdr3, cdr3_nt)) %>% 
    column_to_rownames(var="barcode")
  
  seurat = subset(seurat, cells = rownames(contig_anno_wide))
  
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
  
  data = merge(contig_anno_wide, cell_PC, by="row.names") %>% 
    column_to_rownames(var="Row.names")
  
  exp = t(as.data.frame(seurat@assays[["RNA"]]@data[which(seurat@assays[["RNA"]]@data@Dimnames[[1]] %in% c("CD4", "CD8A", "CD8B")), ]))
  
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
  
  df[[i]] = data
  remove(seurat)
  remove(cell_PC)
  remove(data)
  gc()
  print("done")
}

names(df) = sample

# Targeted dataframe
df_immu = bind_rows(df, .id = "sample") 

# filter sample with # of cells < 500
sam_id = df_immu %>% 
  dplyr::count(sample) %>% 
  filter(n>=500) %>% arrange(desc(n)) %>% pull(sample)

df_immu = df_immu %>%
  filter(sample %in% sam_id) %>%
  group_by(cdr3_TRA, cdr3_TRB) %>%
  mutate(ID = cur_group_id()) %>%
  ungroup() %>%
  mutate(subtype = "CD4+CD8-")
df_immu$sample_ID = paste0(df_immu$sample, "-", df_immu$ID)
saveRDS(df_immu, file ="df_immu_imputed_anno.rds")



df_immu_imputed_anno <- readRDS("~/projects/Research/TCR and Gene Expression/PC distance analysis/GSE187515/df_immu_imputed_anno.rds")
df_immu_imputed_anno = 
  df_immu_imputed_anno %>%
  group_by(cdr3_nt_TRA, cdr3_nt_TRB) %>%
  mutate(ID_nt = cur_group_id()) %>%
  ungroup()
df_immu_imputed_anno$sample_ID_nt = paste0(df_immu_imputed_anno$sample, "-",
                                           df_immu_imputed_anno$ID_nt)

saveRDS(df_immu_imputed_anno, "~/projects/Research/TCR and Gene Expression/PC distance analysis/GSE187515/df_immu.rds")


# replace M1, M2, and M3 since we have same sample names in Zhang_2020
sam = df_immu$sample
sam[which(sam == "M1")] = "2023-1"
sam[which(sam == "M2")] = "2023-2"
sam[which(sam == "M3")] = "2024"

df_immu$sample = sam
df_immu$sample_ID = paste0(df_immu$sample, "-", df_immu$ID)

# add CD4 and CD8 level
df_immu$CD4_level = rep(1, nrow(df_immu))
df_immu$CD8_level = rep(0, nrow(df_immu))

sample_information_GSE187515 <- read_csv("PC distance analysis/GSE187515/sample_information_GSE187515.csv")

df_immu = df_immu %>% merge(sample_information_GSE187515)

saveRDS(df_immu, file ="df_immu_imputed_anno.rds")

