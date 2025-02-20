library(parallel)
library(dplyr)
library(tidyr)
library(tibble)
library(R.utils)
library(stringr)

folder_path = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis"
file_list = list.files(folder_path)
file_GSE = file_list[str_detect(file_list, "GSE|EMTAB9357_pbmc|PRJCA002413_pbmc|10X|Zhang_2020")]
df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_imputed_anno.rds")

for (i in 1:length(df_immu_file)) {
  df_immu = readRDS(df_immu_file[i])
  df_immu = df_immu %>%
    group_by(cdr3_TRA) %>%
    mutate(ID_TRA = cur_group_id()) %>%
    ungroup()
  df_immu$sample_ID_TRA = paste0(df_immu$sample, "-", df_immu$ID_TRA)
  df_immu = df_immu %>%
    group_by(cdr3_TRB) %>%
    mutate(ID_TRB = cur_group_id()) %>%
    ungroup()
  df_immu$sample_ID_TRB = paste0(df_immu$sample, "-", df_immu$ID_TRB)
  df_immu = df_immu %>%
    group_by(cdr3_nt_TRA) %>%
    mutate(ID_nt_TRA = cur_group_id()) %>%
    ungroup()
  
  df_immu$sample_ID_nt_TRA = paste0(df_immu$sample, "-", df_immu$ID_nt_TRA)
  df_immu = df_immu %>%
    group_by(cdr3_nt_TRB) %>%
    mutate(ID_nt_TRB = cur_group_id()) %>%
    ungroup()
  df_immu$sample_ID_nt_TRB = paste0(df_immu$sample, "-", df_immu$ID_nt_TRB)
  df_immu = df_immu %>%
    group_by(cdr3_nt_TRA, cdr3_nt_TRB) %>%
    mutate(ID_nt = cur_group_id()) %>%
    ungroup()
  df_immu$sample_ID_nt = paste0(df_immu$sample, "-", df_immu$ID_nt)
  
  saveRDS(df_immu, file = paste0(folder_path, "/", file_GSE[i], "/df_immu_", file_GSE[i], ".rds"))
  print("done")
}





