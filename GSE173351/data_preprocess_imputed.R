

folder_path <- "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis/GSE173351/untar"

test = read.delim(file.path(folder_path,file_list[1]))

file_list <- list.files(folder_path)
file_txt = file_list[endsWith(file_list, "txt.gz")]

file_targz = file_list[endsWith(file_list, "tar.gz")]
file_anno = file_list[endsWith(file_list, "vdj.tar.gz")]
file_matrix = setdiff(file_targz, file_anno)
file_matrix = file_matrix[!str_detect(file_matrix, "vdt|tdj")]
library(stringr)
file_anno_sam = str_split(file_anno, "\\.vdj", n=2, simplify=TRUE)[,1]
file_matrix_sam = str_split(file_matrix, "\\.tar", n=2, simplify = TRUE)[,1]
file_matrix = file_matrix[is.element(file_matrix_sam, file_anno_sam)]
diff = setdiff(file_matrix_sam, file_anno_sam)
sample = str_split(file_anno_sam, "_", n=2, simplify = TRUE)[,2]

# loop through the list and extract each tar archive
for (file in file_list) {
  if (endsWith(file, ".vdj.tar.gz")) {
    tar_file <- file.path(folder_path, file)
    untar(gunzip(tar_file))
  }
}

for (file in file_list) {
  if (endsWith(file, ".tar.gz")) {
    tar_file <- file.path(folder_path, file)
    untar(tar_file, exdir = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis/GSE173351/count")
  }
}

sample = file_matrix = list.files("count")
file_anno = list.files("vdj_annotation")

df_immu_list = split(df_immu, df_immu$sample)
patientID = str_split(sample, "_", n=2, simplify = TRUE)[,1]
# sample source
tissue = str_split(sample, "_", n=2, simplify = TRUE)[,2]
tissue[str_detect(tissue, "mettumor")] = "mettumor"
tissue[str_detect(tissue, "Flu|MANA|tumor_|W2|W4")] = "tumor"
tissue[str_detect(tissue, "LN")] = "LN"
tissue[str_detect(tissue, "normal")] = "normal"
tissue[str_detect(tissue, "M3")] = "blood"

#time
time = str_split(sample, "_", n=2, simplify = TRUE)[,2]
time[!str_detect(time, "W2|W4")] = "post-treatment"
time[str_detect(time, "W2")] = "W2"
time[str_detect(time, "W4")] = "W4"

demo_info = read.csv("~/projects/Research/TCR and Gene Expression/PC distance analysis/GSE173351/demo_inifo.csv")
sample_information = data.frame(sample = sample,
                                patientID = patientID,
                                tissue = tissue,
                                time = time,
                                disease = "lung cancers")
library(dplyr)
sample_information = left_join(sample_information, demo_info, by = "patientID")
saveRDS(sample_information, file = "sample_information.rds")



# standard pipeline for calculating Principal component distance

sample = file_matrix = list.files("count")
file_anno = list.files("vdj_annotation")
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(patchwork)
library(R.utils)
library(SAVER)
library(factoextra)

df = vector("list", length = length(file_matrix))

for (i in seq_along(file_matrix)) {
  
  # Load the PBMC dataset
  pbmc.data <- Read10X(data.dir = paste0("count", "/", file_matrix[i]))
  
  
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
  contig_anno = read.csv(paste0("vdj_annotation", "/",
                                file_anno[i], "/filtered_contig_annotations.csv"), 
                         header=TRUE, sep = ",")
  
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
  
  if (ncol(pbmc) < 500){
    next
  }
  
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
  
  exp = t(as.data.frame(pbmc@assays[["RNA"]]@data[which(pbmc@assays[["RNA"]]@data@Dimnames[[1]] %in% c("CD4", "CD8A","CD8B")), ]))
  
  data = merge(data, exp, by="row.names") %>% 
    column_to_rownames(var="Row.names")
  
  genes.idx = which(rownames(pbmc@assays[["RNA"]]@counts) %in% c("CD4", "CD8A", "CD8B"))
  pbmc_imputed = saver(pbmc@assays[["RNA"]]@counts, estimates.only = TRUE, ncores = 23, 
                       pred.genes = genes.idx, pred.genes.only = T)
  exp_imp = t(pbmc_imputed)
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

saveRDS(df_immu, file = "df_immu_imputed.rds")

