##### This file is for computing GEX dissimilarity between clones, including 
##### "mixed clones with purity < 100%" and "GEX dissimilarity between CD4+
##### and CD8+ cell clones". 

library(purrr)
folder_path = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis"
source(paste0(folder_path, "/code/GEX_dissim_fun.R"))
file_list = list.files(folder_path)
file_GSE = file_list[str_detect(file_list, "GSE|EMTAB9357_pbmc|PRJCA002413_pbmc|10X|Zhang_2020")]
df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_", file_GSE, ".rds")

# amino acid (aa)
# TRA
for (j in 1:length(file_GSE)) {
  df_immu = readRDS(df_immu_file[j])
  df_immu_list = split(df_immu, df_immu$sample)
  sample_list = sort(unique(df_immu$sample))
  
  # compute GEX dissimilarity within each sample
  for (i in 1:length(df_immu_list)) {
    df_immu = df_immu_list[[i]]
    dist_df_TRA = lv_dist(df_immu_list[[i]], ID = "ID_TRA", sample_ID = "sample_ID_TRA",
                          nucleotide = F, TRA = T, TRB = F, cdr3 = "cdr3_TRA")
    PC_distance_lv = PC_dis_lv(data = df_immu, lv_dis_data = dist_df_TRA, ID="ID_TRA")
    saveRDS(PC_distance_lv,
            paste0(folder_path, "/", file_GSE[j], "/output_CD4_CD8/PC_lv_est/PC_distance_lv_", sample_list[i], ".rds"))
    print("done")
    # clear space for saving memory
    remove("dist_df_TRA")
    remove("PC_distance_lv")
    gc()
    }
  PC_distance_lv_file = list.files(paste0(folder_path, "/", file_GSE[j], "/output_CD4_CD8/PC_lv_est"), full.names = TRUE)
  PC_distance_lv = PC_distance_lv_file %>% map_dfr(readRDS)
  saveRDS(PC_distance_lv, paste0(folder_path, "/", file_GSE[j], "/output_CD4_CD8/PC_distance_lv_TRA.rds"))
  remove(PC_distance_lv)
  gc()
}

# TRB
for (j in 1:length(file_GSE)) {
  
  df_immu = readRDS(df_immu_file[j])
  df_immu_list = split(df_immu, df_immu$sample)
  sample_list = sort(unique(df_immu$sample))
  for (i in 1:length(df_immu_list)) {
    df_immu = df_immu_list[[i]]
    dist_df_TRB = lv_dist(df_immu_list[[i]], ID = "ID_TRB", sample_ID = "sample_ID_TRB",
                          nucleotide = F, TRA = F, TRB = T, cdr3 = "cdr3_TRB")
    PC_distance_lv = PC_dis_lv(data = df_immu, lv_dis_data = dist_df_TRB, ID="ID_TRB")
    saveRDS(PC_distance_lv,
          paste0(folder_path, "/", file_GSE[j], "/output_CD4_CD8/PC_lv_est_TRB/PC_distance_lv_", sample_list[i], ".rds"))
    print("done")
    remove("dist_df_TRB")
    remove("PC_distance_lv")
    gc()
    }
  PC_distance_lv_file = list.files(paste0(folder_path, "/", file_GSE[j], "/output_CD4_CD8/PC_lv_est"), full.names = TRUE)
  PC_distance_lv = PC_distance_lv_file %>% map_dfr(readRDS)
  saveRDS(PC_distance_lv, paste0(folder_path, "/", file_GSE[j], "/output_CD4_CD8/PC_distance_lv_TRB.rds"))
  remove(PC_distance_lv)
  gc()
  print(paste0(j, "done"))
}


# TRA+TRB
for (j in 1:length(file_GSE)) {
  
  df_immu = readRDS(df_immu_file[j])
  df_immu_list = split(df_immu, df_immu$sample)
  sample_list = sort(unique(df_immu$sample))
  
  for (i in 1:length(df_immu_list)) {
    df_immu = df_immu_list[[i]]
    dist_df = lv_dist(df_immu_list[[i]], ID = "ID", sample_ID = "sample_ID",
                      nucleotide = F, TRA = T, TRB = T, cdr3 = "none")
    PC_distance_lv = PC_dis_lv(data = df_immu, lv_dis_data = dist_df, ID="ID")
    saveRDS(PC_distance_lv, 
            paste0(folder_path, "/", file_GSE[j], "/output_CD4_CD8/PC_lv_est/PC_distance_lv_", sample_list[i], ".rds"))
    print("done")
    remove("dist_df")
    #remove("PC_distance_lv")
    gc()
  }
  PC_distance_lv_file = list.files(paste0(folder_path, "/", file_GSE[j], "/output_CD4_CD8/PC_lv_est"), full.names = TRUE)
  PC_distance_lv = PC_distance_lv_file %>% map_dfr(readRDS) %>% bind_rows()
  saveRDS(PC_distance_lv, paste0(folder_path, "/", file_GSE[j], "/output_CD4_CD8/PC_distance_lv.rds"))
  remove(PC_distance_lv)
  gc()
  print(paste0(j, "done"))
}
  
# nucleotides(nt)
# TRA
for (j in 1:length(file_GSE)) {
  
  df_immu = readRDS(df_immu_file[j])
  df_immu_list = split(df_immu, df_immu$sample)
  sample_list = sort(unique(df_immu$sample))
  
  for (i in 1:length(df_immu_list)) {
    df_immu = df_immu_list[[i]]
    dist_df_TRA = lv_dist(df_immu_list[[i]], ID = "ID_nt_TRA", sample_ID = "sample_ID_nt_TRA",
                          nucleotide = T, TRA = T, TRB = F, cdr3 = "cdr3_nt_TRA")
    PC_distance_lv = PC_dis_lv(data = df_immu, lv_dis_data = dist_df_TRA, ID="ID_nt_TRA")
    saveRDS(PC_distance_lv,
            paste0(folder_path, "/", file_GSE[j], "/nt/output_CD4_CD8/PC_lv_est/PC_distance_lv_", sample_list[i], ".rds"))
    print("done")
    remove("dist_df_TRA")
    remove("PC_distance_lv")
    gc()
  }
  PC_distance_lv_file = list.files(paste0(folder_path, "/", file_GSE[j], "/nt/output_CD4_CD8/PC_lv_est"), full.names = TRUE)
  PC_distance_lv = PC_distance_lv_file %>% map_dfr(readRDS)
  saveRDS(PC_distance_lv, paste0(folder_path, "/", file_GSE[j], "/nt/output_CD4_CD8/PC_distance_lv_TRA.rds"))
  remove(PC_distance_lv)
  gc()
  print(paste0(j, "done"))
}


# TRB
for (j in 1:length(file_GSE)) {
  
  df_immu = readRDS(df_immu_file[j])
  df_immu_list = split(df_immu, df_immu$sample)
  sample_list = sort(unique(df_immu$sample))
  
  for (i in 1:length(df_immu_list)) {
    df_immu = df_immu_list[[i]]
    dist_df_TRB = lv_dist(df_immu_list[[i]], ID = "ID_nt_TRB", sample_ID = "sample_ID_nt_TRB",
                          nucleotide = T, TRA = F, TRB = T, cdr3 = "cdr3_nt_TRB")
    PC_distance_lv = PC_dis_lv(data = df_immu, lv_dis_data = dist_df_TRB, ID="ID_nt_TRB")
    saveRDS(PC_distance_lv,
            paste0(folder_path, "/", file_GSE[j], "/nt/output_CD4_CD8/PC_lv_est/PC_distance_lv_", sample_list[i], ".rds"))
    print("done")
    remove("dist_df_TRB")
    remove("PC_distance_lv")
    gc()
  }
  PC_distance_lv_file = list.files(paste0(folder_path, "/", file_GSE[j], "/nt/output_CD4_CD8/PC_lv_est"), full.names = TRUE)
  PC_distance_lv = PC_distance_lv_file %>% map_dfr(readRDS) %>% bind_rows()
  saveRDS(PC_distance_lv, paste0(folder_path, "/", file_GSE[j], "/nt/output_CD4_CD8/PC_distance_lv_TRB.rds"))
  remove(PC_distance_lv)
  gc()
  print(paste0(j, "done"))
}

# TRA+TRB
for (j in 1:length(file_GSE)) {
  
  df_immu = readRDS(df_immu_file[j])
  df_immu_list = split(df_immu, df_immu$sample)
  sample_list = sort(unique(df_immu$sample))
  
  for (i in 1:length(df_immu_list)) {
    df_immu = df_immu_list[[i]]
    dist_df = lv_dist(df_immu_list[[i]], ID = "ID_nt", sample_ID = "sample_ID_nt",
                      nucleotide = T, TRA = T, TRB = T, cdr3 = "none")
    PC_distance_lv = PC_dis_lv(data = df_immu, lv_dis_data = dist_df, ID="ID_nt")
    saveRDS(PC_distance_lv, 
            paste0(folder_path, "/", file_GSE[j], "/nt/output_CD4_CD8/PC_lv_est/PC_distance_lv_", sample_list[i], ".rds"))
    print("done")
    remove("dist_df")
    #remove("PC_distance_lv")
    gc()
  }
  PC_distance_lv_file = list.files(paste0(folder_path, "/", file_GSE[j], "/nt/output_CD4_CD8/PC_lv_est"), full.names = TRUE)
  PC_distance_lv = PC_distance_lv_file %>% map_dfr(readRDS) %>% bind_rows()
  saveRDS(PC_distance_lv, paste0(folder_path, "/", file_GSE[j], "/nt/output_CD4_CD8/PC_distance_lv.rds"))
  remove(PC_distance_lv)
  gc()
  print(paste0(j, "done"))
}



