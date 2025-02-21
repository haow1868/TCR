
# set the path to the folder containing the tar archives
folder_path <- "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis/GSE145370/untar"

# get a list of all files in the folder
file_list <- list.files(path=folder_path)

# loop through the list and extract each tar archive
for (file in file_list) {
  if (endsWith(file, ".tar.gz")) {
    tar_file <- file.path(folder_path, file)
    tar_folder = gsub(".tar.gz$", "", basename(file))
    untar(tar_file, exdir = file.path(folder_path, tar_folder))
  }
}

file_list <- list.files(folder_path)
file_matrix = file_list[endsWith(file_list, "_matrix")]
file_anno = file_list[endsWith(file_list, ".csv.gz")]
library(stringr)
sample = str_split(file_matrix, "_", simplify = TRUE)[,2]



for (i in seq_along(file_matrix)) {
  # Load the PBMC dataset
  seurat.data <- Read10X(data.dir = 
                           paste0(file.path(folder_path, file_matrix[i]),
                                  "/", "filtered_feature_bc_matrix"))
  saveRDS(seurat.data, paste0("/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis/GSE145370/matrix/", sample[i], ".rds"))
  remove(seurat.data)
  gc()
  print("done")
}



# standard pipeline for calculating Principal component distance
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(patchwork)
library(R.utils)
library(SAVER)
df = list()

for (i in seq_along(file_matrix)) {
  
  # Load the PBMC dataset
  pbmc.data <- Read10X(data.dir = 
                         paste0(file.path(folder_path, file_matrix[i]),
                                "/", "filtered_feature_bc_matrix"))
  
  if (length(pbmc.data) == 2){
    pbmc.data = pbmc.data[["Gene Expression"]]
  }
  
  # Initialize the Seurat object with the raw (non-normalized data).
  pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200)
  
  # Standard pre-processing workflow
  # QC and selecting cells for further analysis
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  
  # Remove nfeature < 200 and > 2500, remove percent.mt >= 5
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
  
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
  
  pbmc = subset(pbmc, cells = rownames(contig_anno_wide))
  
  # Normalize the data
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Identification of highly variable features (feature selection), select 2000 features
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  
  # Scaling the data
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  
  # Perform linear dimensional reduction
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  
  cell_PC = as.data.frame(pbmc@reductions[["pca"]]@cell.embeddings)[1:10]
  data = merge(contig_anno_wide, cell_PC, by="row.names") %>% 
    column_to_rownames(var="Row.names")
  
  exp = t(as.data.frame(pbmc@assays[["RNA"]]@data[which(pbmc@assays[["RNA"]]@data@Dimnames[[1]] %in% c("CD4", "CD8A", "CD8B")), ]))
  data = merge(data, exp, by="row.names") %>% 
    column_to_rownames(var="Row.names")
  
  # imputed
  genes.idx = which(rownames(pbmc@assays[["RNA"]]@counts) %in% c("CD4", "CD8A","CD8B"))
  seurat_imputed = saver(pbmc@assays[["RNA"]]@counts, estimates.only = TRUE, ncores = 23, 
                         pred.genes = genes.idx, pred.genes.only = T)
  exp_imp = t(seurat_imputed)
  colnames(exp_imp) = paste0(colnames(exp_imp),"_imp")
  data = merge(data, exp_imp, by="row.names")
  
  colnames(data)[which(colnames(data) == "Row.names")] = "barcode"
  
  df[[i]] = data
  remove(pbmc)
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

saveRDS(df_immu, file = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis/GSE145370/df_immu_imputed.rds")



