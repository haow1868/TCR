library(Seurat)
library(stringr)
folder_path = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis/GSE175450"

cell_metadata = read.table(paste0(folder_path, "/GSE175450_cell_metadata_table.tsv.gz"),
                           header = TRUE)
sample = str_split(cell_metadata$sample_name, "_PBMC|_Tcells", n = 2 ,simplify = TRUE)[,1]
cell_metadata$sample = sample
barcode = str_split(cell_metadata$cell_id, "_", n = 2, simplify = TRUE)[,1]
cell_metadata$barcode = barcode
cell_metadata$id = str_split(cell_metadata$cell_id, "_", n = 2, simplify = TRUE)[,2]
cell_metadata_list = split(cell_metadata, cell_metadata$orig.ident)

h5_file = list.files(paste0(folder_path, "/h5_file"), full.names = TRUE)
file_anno = list.files(paste0(folder_path, "/filtered_contig_annotations"), full.names = TRUE)

ident_id = cell_metadata %>% select(orig.ident, id) %>% unique() %>% arrange(orig.ident)
id = ident_id$id
countmtx = vector("list", length(h5_file))
for (i in seq_along(h5_file)) {
  seurat = Read10X_h5(filename = h5_file[i])[["Gene Expression"]]
  col = colnames(seurat)
  new_col = paste0(col, "_", id[i])
  colnames(seurat) = new_col
  countmtx[[i]] = seurat
  remove(seurat)
  gc()
}

mtx = RowMergeSparseMatrices(countmtx[[1]], countmtx[[2]])
for (i in 3:length(countmtx)) {
  mtx = RowMergeSparseMatrices(mtx, countmtx[[i]])
  print("done")
}

sample_cell_id = cell_metadata %>% select(sample, cell_id)
sample_cell_id_list = split(sample_cell_id, sample_cell_id$sample)

for (i in seq_along(sample_cell_id_list)){
  barcode = sample_cell_id_list[[i]]$cell_id
  sample_mtx = mtx[,colnames(mtx) %in% barcode]
  saveRDS(sample_mtx,
          file = paste0(folder_path, "/matrix/", names(sample_cell_id_list)[i], ".rds"))
  print("done")
}

file_anno = list.files(paste0(folder_path, "/filtered_contig_annotations"), full.names = TRUE)
contig_anno_list = vector("list", length(file_anno))
for (i in seq_along(file_anno)) {
  contig_anno = read.csv(file_anno[i])
  contig_anno$barcode = paste0(contig_anno$barcode, "_", id[i])
  contig_anno_list[[i]] = contig_anno
  write_csv(contig_anno, file = file_anno[i])
}

contig_anno = bind_rows(contig_anno_list)


# the above is to split two data files
folder_path = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis/GSE175450"
file_anno = list.files(paste0(folder_path, "/filtered_contig_annotations"), full.names = TRUE)
library(purrr)
contig_anno = file_anno %>% map_dfr(read.csv)
contig_anno = contig_anno[contig_anno$is_cell %in% c("TRUE","true") &
                            contig_anno$high_confidence %in% c("TRUE","true") &
                            (contig_anno$chain=="TRA" | contig_anno$chain=="TRB") &
                            contig_anno$productive %in% c("TRUE","true"),] %>%
  select(barcode, chain, v_gene, d_gene, j_gene, cdr3, cdr3_nt) %>%
  filter(cdr3 != "None", cdr3_nt != "None")

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
    values_from = c(v_gene, d_gene, j_gene, cdr3, cdr3_nt)) %>%
  column_to_rownames(var="barcode")
saveRDS(contig_anno_wide, paste0(folder_path, "/tcr.rds"))
library(R.utils)
library(tidyverse)
library(Seurat)
library(dplyr)
library(SAVER)
matrix_file = list.files(paste0(folder_path, "/matrix"), full.names = TRUE)
df = vector("list", length = length(matrix_file))
for (i in seq_along(matrix_file)) {
  
  # Initialize the Seurat object with the raw (non-normalized data).
  seurat <- CreateSeuratObject(counts = readRDS(matrix_file[i]), 
                               min.cells = 3, min.features = 200)
  
  # Standard pre-processing workflow
  # QC and selecting cells for further analysis
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  
  # Remove nfeature < 200 and > 2500, remove percent.mt > 20
  seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
  
  seurat = subset(seurat, cells = contig_anno_wide$barcode)
  
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
  data = merge(contig_anno_wide, cell_PC, by="row.names") %>% 
    column_to_rownames(var="Row.names")
  
  exp = t(as.data.frame(seurat@assays[["RNA"]]@data[which(seurat@assays[["RNA"]]@data@Dimnames[[1]] %in% 
                                                            c("CD4", "CD8A", "CD8B")), ]))
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
sample = str_split(list.files(paste0(folder_path, "/matrix")), pattern = ".rds", simplify = TRUE)[,1]
names(df) = sample

# targeted dataframe
df_immu = bind_rows(df, .id = "sample")

# filter sample with # of cells < 500
sam_id = df_immu %>% 
  dplyr::count(sample) %>% 
  filter(n>=500) %>% arrange(desc(n)) %>% pull(sample)

df_immu = df_immu %>%
  filter(sample %in% sam_id) %>%
  group_by(cdr3_TRA, cdr3_TRB) %>%
  mutate(ID = cur_group_id()) %>%
  ungroup()
df_immu$sample_ID = paste0(df_immu$sample, "-", df_immu$ID)

saveRDS(df_immu, file = paste0(folder_path, "/df_immu_imputed.rds"))


sample_information = 
  cell_metadata %>% 
  select(sample, disease_group, disease_severity, disease_phase) %>%
  unique()
write_csv(sample_information, 
          file = paste0(folder_path, "/sample_information.csv"))




