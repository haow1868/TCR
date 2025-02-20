##### This file is for getting data of CD4+ cell clones and CD8+ cell clones
##### separately since the results from PC_distance_lv_pipeline and 
##### PC_distance_within_clones_pipeline consider both "mixed clones with < 100% purity"
##### and "GEX disimilarity between CD4+ and CD8+ cell clones".

library(stringr)
library(purrr)
folder_path = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis"
file_list = list.files(folder_path)
file_GSE = file_list[str_detect(file_list, "GSE|EMTAB9357_pbmc|PRJCA002413_pbmc|10X|Zhang_2020")]

# get the CD4 and CD8 clones and corresponding normalizing constant
#df_immu_file = paste0(folder_path, "/", file_GSE,"/df_immu_imputed_anno.rds")
df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_", file_GSE, ".rds")

# clone ID with 100% purity
pure_clone_ID = function(data, file_GSE, ID, sample_ID) {
  
  # Convert string arguments to symbols
  ID_sym = sym(ID)
  sample_ID_sym = sym(sample_ID)
  
  total_num = df_immu %>% 
    dplyr::count(sample, !!ID_sym, !!sample_ID_sym, name="total_num")
  
  sub_num = df_immu %>% 
    dplyr::count(sample, !!ID_sym, !!sample_ID_sym, subtype) %>%
    merge(total_num) %>%
    mutate(per = n/total_num) %>%
    filter(per == 1)
  
  ID_CD4 = sub_num[sub_num$subtype=="CD4+CD8-",] %>% select(sample, !!ID_sym, !!sample_ID_sym, subtype, total_num)
  saveRDS(ID_CD4, file = paste0(folder_path, "/", file_GSE, "/pure_clone_ID/", ID, "_CD4.rds"))
  ID_CD8 = sub_num[sub_num$subtype=="CD4-CD8+",] %>% select(sample, !!ID_sym, !!sample_ID_sym, subtype, total_num)
  saveRDS(ID_CD8, file = paste0(folder_path, "/", file_GSE, "/pure_clone_ID/", ID, "_CD8.rds"))
  
}

# save clone with 100% purity
for (i in 1:length(file_GSE)) {
  df_immu = readRDS(df_immu_file[i])
  pure_clone_ID(data = df_immu, file_GSE[i], "ID", "sample_ID")
  pure_clone_ID(data = df_immu, file_GSE[i], "ID_TRA", "sample_ID_TRA")
  pure_clone_ID(data = df_immu, file_GSE[i], "ID_TRB", "sample_ID_TRB")
  pure_clone_ID(data = df_immu, file_GSE[i], "ID_nt", "sample_ID_nt")
  pure_clone_ID(data = df_immu, file_GSE[i], "ID_nt_TRA", "sample_ID_nt_TRA")
  pure_clone_ID(data = df_immu, file_GSE[i], "ID_nt_TRB", "sample_ID_nt_TRB")
}


distance_full_separation = function(data, file_GSE, ID, sample_ID, nucleotide, subtype, clone_type){
  if (subtype == "CD4") {
    pure_clone_ID = paste0(folder_path, "/", file_GSE, "/pure_clone_ID/", ID, "_CD4.rds") %>% map_dfr(readRDS)
  } else if (subtype == "CD8") {
    pure_clone_ID = paste0(folder_path, "/", file_GSE, "/pure_clone_ID/", ID, "_CD8.rds") %>% map_dfr(readRDS)
  }
  
  if (nucleotide == T) {
    folder = paste0(file_GSE, "/nt/output_", subtype, "/", clone_type)
  } else if (nucleotide == F) {
    folder = paste0(file_GSE, "/output_", subtype, "/", clone_type)
  }
  pure_clone_sample_ID = pure_clone_ID[[sample_ID]]
  distance_full_filtered = data %>% filter(!!sym(sample_ID) %in% pure_clone_sample_ID) 
  saveRDS(distance_full_filtered, file = paste0(folder_path, "/", folder, "/distance_full.rds"))
}

# for aa
for (i in 1:length(file_GSE)) {
  distance_full_file = paste0(folder_path, "/", file_GSE, "/output_CD4_CD8/distance_full.rds")
  distance_full = readRDS(distance_full_file[i])
  distance_full_separation(data = distance_full, file_GSE = file_GSE[i], ID = "ID", 
                           sample_ID = "sample_ID", nucleotide = F, subtype = "CD4", clone_type = "TRA+TRB")
  distance_full_separation(data = distance_full, file_GSE = file_GSE[i], ID = "ID",
                           sample_ID = "sample_ID", nucleotide = F, subtype = "CD8", clone_type = "TRA+TRB")
  
  distance_full_TRA_file = paste0(folder_path, "/", file_GSE, "/output_CD4_CD8/distance_full_TRA.rds")
  distance_full_TRA = readRDS(distance_full_TRA_file[i])
  distance_full_separation(data = distance_full_TRA, file_GSE = file_GSE[i], ID = "ID_TRA",
                           sample_ID = "sample_ID_TRA", nucleotide = F, subtype = "CD4", clone_type = "TRA")
  distance_full_separation(data = distance_full_TRA, file_GSE = file_GSE[i], ID = "ID_TRA",
                           sample_ID = "sample_ID_TRA", nucleotide = F, subtype = "CD8", clone_type = "TRA")
  
  distance_full_TRB_file = paste0(folder_path, "/", file_GSE, "/output_CD4_CD8/distance_full_TRB.rds")
  distance_full_TRB = readRDS(distance_full_TRB_file[i])
  distance_full_separation(data = distance_full_TRB, file_GSE = file_GSE[i], ID = "ID_TRB",
                           sample_ID = "sample_ID_TRB", nucleotide = F, subtype = "CD4", clone_type = "TRB")
  distance_full_separation(data = distance_full_TRB, file_GSE = file_GSE[i], ID = "ID_TRB",
                           sample_ID = "sample_ID_TRB", nucleotide = F, subtype = "CD8", clone_type = "TRB")
  print(paste0(i,"done"))
}

# for nt
for (i in 1:length(file_GSE)) {
  distance_full_file = paste0(folder_path, "/", file_GSE, "/nt/output_CD4_CD8/distance_full.rds")
  distance_full = readRDS(distance_full_file[i])
  distance_full_separation(data = distance_full, file_GSE = file_GSE[i], ID = "ID_nt", 
                           sample_ID = "sample_ID_nt", nucleotide = T, subtype = "CD4", clone_type = "TRA+TRB")
  distance_full_separation(data = distance_full, file_GSE = file_GSE[i], ID = "ID_nt",
                           sample_ID = "sample_ID_nt", nucleotide = T, subtype = "CD8", clone_type = "TRA+TRB")
  
  distance_full_TRA_file = paste0(folder_path, "/", file_GSE, "/nt/output_CD4_CD8/distance_full_TRA.rds")
  distance_full_TRA = readRDS(distance_full_TRA_file[i])
  distance_full_separation(data = distance_full_TRA, file_GSE = file_GSE[i], ID = "ID_nt_TRA",
                           sample_ID = "sample_ID_nt_TRA", nucleotide = T, subtype = "CD4", clone_type = "TRA")
  distance_full_separation(data = distance_full_TRA, file_GSE = file_GSE[i], ID = "ID_nt_TRA",
                           sample_ID = "sample_ID_nt_TRA", nucleotide = T, subtype = "CD8", clone_type = "TRA")
  
  distance_full_TRB_file = paste0(folder_path, "/", file_GSE, "/nt/output_CD4_CD8/distance_full_TRB.rds")
  distance_full_TRB = readRDS(distance_full_TRB_file[i])
  distance_full_separation(data = distance_full_TRB, file_GSE = file_GSE[i], ID = "ID_nt_TRB",
                           sample_ID = "sample_ID_nt_TRB", nucleotide = T, subtype = "CD4", clone_type = "TRB")
  distance_full_separation(data = distance_full_TRB, file_GSE = file_GSE[i], ID = "ID_nt_TRB",
                           sample_ID = "sample_ID_nt_TRB", nucleotide = T, subtype = "CD8", clone_type = "TRB")
  print(paste0(i,"done"))
}

# a new normalization way: average over GEX dissimilarity when levenshtein distance > 5
PC_distance_lv_file = paste0(folder_path, "/", file_GSE, "/output_CD4_CD8/PC_distance_lv.rds")
PC_distance_lv_file = paste0(folder_path, "/", file_GSE, "/output_CD4_CD8/PC_distance_lv_TRA.rds")
PC_distance_lv_file = paste0(folder_path, "/", file_GSE, "/output_CD4_CD8/PC_distance_lv_TRB.rds")

clone_separation = function(data, ID, sample_ID, threshold, CD4, folder){
  
  if (CD4 == T) {
    pure_clone_ID = paste0(folder_path, "/", file_GSE, "/pure_clone_ID/", ID, "_CD4.rds") %>% map_dfr(readRDS)
  } else if (CD4 == F){
    pure_clone_ID = paste0(folder_path, "/", file_GSE, "/pure_clone_ID/", ID, "_CD8.rds") %>% map_dfr(readRDS)
  }
  
  data[[paste0(sample_ID, ".x")]] = paste0(data$sample, "-", data[[paste0(ID, ".x")]])
  data[[paste0(sample_ID, ".y")]] = paste0(data$sample, "-", data[[paste0(ID, ".y")]])
  
  PC_distance_lv = data %>%
    filter(!!sym(paste0(sample_ID, ".x")) %in% pure_clone_ID[[sample_ID]] & 
             !!sym(paste0(sample_ID, ".y")) %in% pure_clone_ID[[sample_ID]])
  saveRDS(PC_distance_lv, file = paste0(folder_path, "/", folder, "/PC_distance_lv.rds"))
  
  sample_median = 
    PC_distance_lv %>% group_by(sample, distance) %>%
    dplyr::summarise(median = median(PC_dis), .groups = "drop")
  saveRDS(sample_median, file = paste0(folder_path, "/", folder, "/sample_median.rds"))
  
  # normalizing constant
  med_dist_bet_clones = 
    PC_distance_lv %>% filter(distance > threshold) %>%
    group_by(sample) %>% dplyr::summarise(median_distance_two_clones = median(PC_dis))
  saveRDS(med_dist_bet_clones, file = paste0(folder_path, "/", folder, "/med_dist_bet_clones.rds"))
}

# for aa
for(i in 1:length(file_GSE)){
  # TRA+TRB
  PC_distance_lv_file = paste0(folder_path, "/", file_GSE, "/output_CD4_CD8/PC_distance_lv.rds")
  PC_distance_lv = readRDS(PC_distance_lv_file[i])
  # CD4
  clone_separation(data = PC_distance_lv, ID = "ID", sample_ID = "sample_ID",
                   threshold = 10, CD4 = T, folder = paste0(file_GSE[i], "/output_CD4/TRA+TRB"))
  # CD8
  clone_separation(data = PC_distance_lv, ID = "ID", sample_ID = "sample_ID",
                   threshold = 10, CD4 = F, folder = paste0(file_GSE[i], "/output_CD8/TRA+TRB"))
  remove(PC_distance_lv)
  gc()
  
  # TRA
  PC_distance_lv_file = paste0(folder_path, "/", file_GSE, "/output_CD4_CD8/PC_distance_lv_TRA.rds")
  PC_distance_lv = readRDS(PC_distance_lv_file[i])
  # CD4
  clone_separation(data = PC_distance_lv, ID = "ID_TRA", sample_ID = "sample_ID_TRA",
                   threshold = 5, CD4 = T, folder = paste0(file_GSE[i], "/output_CD4/TRA"))
  # CD8
  clone_separation(data = PC_distance_lv, ID = "ID_TRA", sample_ID = "sample_ID_TRA",
                   threshold = 5, CD4 = F, folder = paste0(file_GSE[i], "/output_CD8/TRA"))
  remove(PC_distance_lv)
  gc()
  
  # TRB
  PC_distance_lv_file = paste0(folder_path, "/", file_GSE, "/output_CD4_CD8/PC_distance_lv_TRB.rds")
  PC_distance_lv = readRDS(PC_distance_lv_file[i])
  # CD4
  clone_separation(data = PC_distance_lv, ID = "ID_TRB", sample_ID = "sample_ID_TRB",
                   threshold = 5, CD4 = T, folder = paste0(file_GSE[i], "/output_CD4/TRB"))
  # CD8
  clone_separation(data = PC_distance_lv, ID = "ID_TRB", sample_ID = "sample_ID_TRB",
                   threshold = 5, CD4 = F, folder = paste0(file_GSE[i], "/output_CD8/TRB"))
  remove(PC_distance_lv)
  gc()
  print(paste0(i, "done"))
}


# for nt
for(i in 1:length(file_GSE)){
  # TRA+TRB
  PC_distance_lv_file = paste0(folder_path, "/", file_GSE, "/nt/output_CD4_CD8/PC_distance_lv.rds")
  PC_distance_lv = readRDS(PC_distance_lv_file[i])
  # CD4
  clone_separation(data = PC_distance_lv, ID = "ID_nt", sample_ID = "sample_ID_nt",
                   threshold = 10, CD4 = T, folder = paste0(file_GSE[i], "/nt/output_CD4/TRA+TRB"))
  # CD8
  clone_separation(data = PC_distance_lv, ID = "ID_nt", sample_ID = "sample_ID_nt",
                   threshold = 10, CD4 = F, folder = paste0(file_GSE[i], "/nt/output_CD8/TRA+TRB"))
  remove(PC_distance_lv)
  gc()
  
  # TRA
  PC_distance_lv_file = paste0(folder_path, "/", file_GSE, "/nt/output_CD4_CD8/PC_distance_lv_TRA.rds")
  PC_distance_lv = readRDS(PC_distance_lv_file[i])
  # CD4
  clone_separation(data = PC_distance_lv, ID = "ID_nt_TRA", sample_ID = "sample_ID_nt_TRA",
                   threshold = 5, CD4 = T, folder = paste0(file_GSE[i], "/nt/output_CD4/TRA"))
  # CD8
  clone_separation(data = PC_distance_lv, ID = "ID_nt_TRA", sample_ID = "sample_ID_nt_TRA",
                   threshold = 5, CD4 = F, folder = paste0(file_GSE[i], "/nt/output_CD8/TRA"))
  remove(PC_distance_lv)
  gc()
  
  # TRB
  PC_distance_lv_file = paste0(folder_path, "/", file_GSE, "/nt/output_CD4_CD8/PC_distance_lv_TRB.rds")
  PC_distance_lv = readRDS(PC_distance_lv_file[i])
  # CD4
  clone_separation(data = PC_distance_lv, ID = "ID_nt_TRB", sample_ID = "sample_ID_nt_TRB",
                   threshold = 5, CD4 = T, folder = paste0(file_GSE[i], "/nt/output_CD4/TRB"))
  # CD8
  clone_separation(data = PC_distance_lv, ID = "ID_nt_TRB", sample_ID = "sample_ID_nt_TRB",
                   threshold = 5, CD4 = F, folder = paste0(file_GSE[i], "/nt/output_CD8/TRB"))
  remove(PC_distance_lv)
  gc()
  print(paste0(i, "done"))
}










