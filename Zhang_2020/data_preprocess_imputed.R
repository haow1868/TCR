
folder_path = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis"


library(Seurat)
seurat = read.csv(paste0(folder_path, "/Zhang_2020/scRNA_gene_expression.csv.gz"))

library(data.table)
count = fread(paste0(folder_path, "/Zhang_2020/scRNA_gene_expression.csv"), header=TRUE)
count = count %>% column_to_rownames(var = "V1")

gbm <- as.matrix(count)

# save sparse matrix
library(Matrix)
sparse.gbm <- Matrix(gbm , sparse = T)
head(sparse.gbm)
## Market Exchange Format (MEX) format
writeMM(obj = sparse.gbm, file="matrix.mtx")

# save genes and cells names
write(x = rownames(gbm), file = "genes.tsv")
write(x = colnames(gbm), file = "barcodes.tsv")

saveRDS(sparse.gbm,
        file = paste0(folder_path, "/Zhang_2020/matrix_all_cells.rds"))

seurat = CreateSeuratObject(counts = sparse.gbm)


untar(paste0(folder_path, "/Zhang_2020/SARS_COV_2_immune_processed_data.tar.gz"),
      exdir = paste0(folder_path, "/Zhang_2020/immune_data"))

barcode = colnames(sparse.gbm)
id = str_split(barcode, "_", n = 2, simplify = TRUE)[,2] %>% unique()
sample = list.files(paste0(folder_path, "/Zhang_2020/vdj_annotation"))
for (i in 2:14) {
  cell = barcode[endsWith(barcode, paste0("_",i))]
  count = sparse.gbm[, colnames(sparse.gbm) %in% cell]
  colnames(count) = str_split(cell, "_", n = 2, simplify = TRUE)[,1]
  saveRDS(count,
          paste0(folder_path, "/Zhang_2020/matrix/", sample[i], ".rds"))
  remove(count)
  gc()
  print("done")
}

sample = list.files(paste0(folder_path, "/Zhang_2020/vdj_annotation"))
file_anno = paste0(folder_path, "/Zhang_2020/vdj_annotation/", sample, "/filtered_contig_annotations.csv")
file_matrix = paste0(folder_path, "/Zhang_2020/matrix/", sample, ".rds")

library(Seurat)
library(R.utils)
library(tidyverse)
df = vector("list", length = length(file_matrix))
for (i in seq_along(file_matrix)) {
  
  # Initialize the Seurat object with the raw (non-normalized data).
  seurat <- CreateSeuratObject(counts = readRDS(file_matrix[i]), 
                               min.cells = 3, min.features = 200)
  
  # Standard pre-processing workflow
  # QC and selecting cells for further analysis
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  
  # Remove nfeature < 200 and > 2500, remove percent.mt > 20
  seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
  # Read annotation information
  contig_anno = read_csv(file_anno[i])
  
  contig_anno = contig_anno[contig_anno$is_cell %in% c("TRUE", "True") &
                              contig_anno$high_confidence %in% c("TRUE", "True") &
                              (contig_anno$chain=="TRA" | contig_anno$chain=="TRB") &
                              contig_anno$productive %in% c("TRUE","True"),] %>%
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
  
  seurat = subset(seurat, cells = rownames(contig_anno_wide))
  
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
  ungroup()
df_immu$sample_ID = paste0(df_immu$sample, "-", df_immu$ID)

saveRDS(df_immu, 
        file = "~/projects/Research/TCR and Gene Expression/PC distance analysis/Zhang_2020/df_immu_imputed.rds")



