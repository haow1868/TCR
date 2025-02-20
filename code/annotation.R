library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
folder_path = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis"
file_list = list.files(folder_path)
file_GSE = file_list[str_detect(file_list, "GSE|EMTAB9357_pbmc|PRJCA002413_pbmc|10X|Zhang_2020")]
file_GSE = file_GSE[which(!file_GSE %in% c("GSE187515"))]

df_immu_imputed_file = paste0(folder_path, "/", file_GSE, "/df_immu_imputed.rds")
df_immu_imputed_all = df_immu_imputed_file %>% map_dfr(readRDS)
df_immu_imputed_all$CD4_imp_log2 = log2(df_immu_imputed_all$CD4_imp + 1)
df_immu_imputed_all$CD8A_imp_log2 = log2(df_immu_imputed_all$CD8A_imp + 1)

sample_information_file = paste0(folder_path, "/", file_GSE,"/sample_information_", file_GSE, ".csv")
sample_information = sample_information_file %>% map_dfr(read.csv)

df_immu_imputed_all = df_immu_imputed_all %>% merge(sample_information)
sam = df_immu_imputed_all$sample %>% unique()

derivative = function(x,y){
  der = c()
  for (i in 1:(length(x)-1)) {
    der[i] = (y[i+1] - y[i])/(x[i+1] - x[i])
  }
  return(der)
}

trough = function(x, der1, der2){
  tro = c()
  for (i in 1:length(der2)) {
    if (der1[i]<0  & der1[i+1]>0 & der2[i]>0){
      tro[i] = x[i]
    } else {
      tro[i] = 0
    }
  }
  tro = tro[!tro==0][1]
  return(tro)
}

study = df_immu_imputed_all %>% pull(GEO) %>% unique() %>% sort()
trough_CD4 = c()
trough_CD8A = c()
for (i in seq_along(study)){
  data = df_immu_imputed_all %>% filter(GEO %in% study[i])
  
  CD4_den = density(data$CD4_imp_log2)
  CD4_x = CD4_den[["x"]]
  CD4_y = CD4_den[["y"]]
  CD4_1_der = derivative(x = CD4_x, y = CD4_y)
  CD4_2_der = derivative(x = CD4_x[1:(length(CD4_x)-1)],y = CD4_1_der)
  trough_CD4[i] = trough(x = CD4_x, der1 = CD4_1_der, der2 = CD4_2_der)
  
  CD8A_den = density(data$CD8A_imp_log2)
  CD8A_x = CD8A_den[["x"]]
  CD8A_y = CD8A_den[["y"]]
  CD8A_1_der = derivative(x = CD8A_x, y = CD8A_y)
  CD8A_2_der = derivative(x = CD8A_x[1:(length(CD8A_x)-1)],y = CD8A_1_der)
  trough_CD8A[i] = trough(x = CD8A_x, der1 = CD8A_1_der, der2 = CD8A_2_der)
}

cutoff =  data.frame(GEO = study,
                     trough_CD4 = trough_CD4,
                     trough_CD8A = trough_CD8A)
saveRDS(cutoff, file = "PC distance analysis/data_allStudy/cutoff_all.rds")


for (i in seq_along(study)) {
  df_immu_imputed = df_immu_imputed_all %>% filter(GEO %in% study[i])
  df_immu_imputed = df_immu_imputed %>% merge(cutoff)
  
  CD4_imp_log2 = df_immu_imputed$CD4_imp_log2
  CD8A_imp_log2 = df_immu_imputed$CD8A_imp_log2
  cutoff_CD4 = df_immu_imputed$trough_CD4
  cutoff_CD8A = df_immu_imputed$trough_CD8A
  
  CD4_level = if_else(CD4_imp_log2 > cutoff_CD4, "CD4+", "CD4-")
  CD8_level = if_else(CD8A_imp_log2 > cutoff_CD8A, "CD8+", "CD8-")
  subtype = paste0(CD4_level, CD8_level)
  df_immu_imputed$subtype = subtype
  
  df_immu_imputed = 
    df_immu_imputed %>%
    mutate(CD4_level = if_else(CD4_level == "CD4+", 1, 0),
           CD8_level = if_else(CD8_level == "CD8+", 1, 0)) %>%
    select(!c("v_gene_TRA", "v_gene_TRB", "d_gene_TRA", "d_gene_TRB", "j_gene_TRA", "j_gene_TRB"))

  df_immu_imputed = df_immu_imputed[ , colSums(is.na(df_immu_imputed))==0]
  
  df_immu_imputed = df_immu_imputed %>% filter(subtype %in% c("CD4+CD8-", "CD4-CD8+"))
  saveRDS(df_immu_imputed, file = paste0(folder_path, "/", file_GSE[i], "/df_immu_imputed_anno.rds"))
  print("done")
}



