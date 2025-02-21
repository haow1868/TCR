
folder_path <- "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis/GSE168732/untar"

file_list = list.files(folder_path)

file_matrix = file_list[endsWith(file_list, "matrix.mtx.gz")]
file_genes = file_list[endsWith(file_list, "features.tsv.gz")]
file_barcodes = file_list[endsWith(file_list, "barcodes.tsv.gz")]
file_anno = file_list[endsWith(file_list, "TCR_filtered_contig_annotations.csv.gz")]
library(stringr)
sample = paste0(str_split(file_matrix, "_", simplify = TRUE)[,2],
                "_",
                str_split(file_matrix, "_", simplify = TRUE)[,3])

for (i in seq_along(file_matrix)) {
  
  # Load the seurat dataset
  seurat.data <- ReadMtx(mtx = file.path(folder_path, file_matrix[i]),
                         cells = file.path(folder_path, file_barcodes[i]),
                         features = file.path(folder_path, file_genes[i]))
  
  if (length(seurat.data) == 2){
    seurat.data = seurat.data[["Gene Expression"]]
  }
  saveRDS(seurat.data, paste0("/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis/GSE168732/matrix/", sample[i], ".rds"))
  remove(seurat.data)
  gc()
  print("done")
}


library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(patchwork)
library(R.utils)
library(SAVER)
df = list()

for (i in seq_along(file_matrix)) {
  
  # Load the seurat dataset
  seurat.data <- ReadMtx(mtx = file.path(folder_path, file_matrix[i]),
                         cells = file.path(folder_path, file_barcodes[i]),
                         features = file.path(folder_path, file_genes[i]))
  
  if (length(seurat.data) == 2){
    seurat.data = seurat.data[["Gene Expression"]]
  }
  
  # Initialize the Seurat object with the raw (non-normalized data).
  seurat <- CreateSeuratObject(counts = seurat.data, min.cells = 3, min.features = 200)
  
  # Standard pre-processing workflow
  # QC and selecting cells for further analysis
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  
  # Remove nfeature < 200 and > 2500, remove percent.mt >= 5
  seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
  
  # Read annotation information
  contig_anno = read.table(file.path(folder_path, file_anno[i]), header=TRUE, sep = ",")
  
  contig_anno = contig_anno[contig_anno$is_cell=="True" &
                              contig_anno$high_confidence=="True" &
                              (contig_anno$chain=="TRA" | contig_anno$chain=="TRB") &
                              contig_anno$productive=="True",] %>%
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
  
  # Normalize the data
  seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Identification of highly variable features (feature selection), select 2000 features
  seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
  
  # Scaling the data
  all.genes <- rownames(seurat)
  seurat <- ScaleData(seurat, features = all.genes)
  
  # Perform linear dimensional reduction
  seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
  
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
  ungroup()
df_immu$sample_ID = paste0(df_immu$sample, "-", df_immu$ID)

saveRDS(df_immu, file = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis/GSE168732/df_immu_imputed.rds")
