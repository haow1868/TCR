##### This file is for computing GEX dissimilarity within the same T cell clones
##### i.e., TCR dissimilarity = 0.

library(parallel)
library(dplyr)
library(tidyr)
library(tibble)
library(R.utils)
library(stringr)

folder_path = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis"
file_list = list.files(folder_path)
file_GSE = file_list[str_detect(file_list, "GSE|EMTAB9357_pbmc|PRJCA002413_pbmc|10X|Zhang_2020")]

# select clones with 100% purity
  same_subtype = function(data, ID, sample_ID){
    clone_ID = data %>% 
      count(sample, !!sym(ID), !!sym(sample_ID), subtype) %>% 
      count(sample, !!sym(ID), !!sym(sample_ID)) %>%
      filter(n == 1) %>%
      pull(!!sym(sample_ID))
    
    data_filtered = data %>% filter(!!sym(sample_ID) %in% clone_ID)
    return(data_filtered)
  }
  
  # PC distance calculation
  PC_dis = function(data){
    # euclidean distance based on 10 PCs
    dis = as.vector(dist(data %>% select(matches("PC")) ))
    dis_clone = mean(dis)
  }
  
  PC_dis_clone = function(data, ID, sample_ID) {
    
    # Convert string arguments to symbols
    ID_sym = sym(ID)
    sample_ID_sym = sym(sample_ID)
    
    # the number of cells
    tot_num = data %>%
      group_by(sample = !!sym("sample"), !!ID_sym, !!sample_ID_sym) %>%
      summarise(total_num = n(), .groups = "drop")
    
    # Join data with tot_num using dynamic column names
    df = left_join(data, tot_num, by=c("sample", as.character(ID_sym), as.character(sample_ID_sym))) %>% 
      filter(total_num != 1)
    
    # Split the data frame into a list by sample_ID
    df_list = split(df, df[[as.character(sample_ID_sym)]])
    
    # Assume PC_dis is a function you have defined elsewhere that calculates distances
    distance = data.frame(PC_dis = unlist(lapply(df_list, PC_dis))) %>%
      tibble::rownames_to_column(var = as.character(sample_ID_sym))
    
    # Join distance with tot_num and add a new column
    distance_clone = left_join(distance, tot_num, by=c(as.character(sample_ID_sym))) %>% 
      mutate(type = "estimated")
  }
  
  PC_permuted = function(data, seed){
    set.seed(seed = seed)
    idx = sample(1:nrow(data)) # reassign cell idx
    data[, grep("PC", colnames(data))] = data[idx, grep("PC", colnames(data))]
    data
  }
  
  # permutation within CD4/CD8 separately for avoiding confounding on cell subtypes
  df_permuted = function(data){
    data_CD4 = 
      data %>% filter(subtype %in% c("CD4+CD8-")) %>%
      PC_permuted(seed=123)
    data_CD8 = 
      data %>% filter(subtype %in%  c("CD4-CD8+")) %>%
      PC_permuted(seed=123)
    
    data_per = rbind(data_CD4, data_CD8)
    return(data_per)
  }
  
  PC_dist_per = function(data, ID, sample_ID, n){
    data_list = split(data, data$sample)
    
    distance_per <- mclapply(1:n, function(i) {
      data_permuted <- bind_rows(mclapply(data_list, df_permuted,
                                          mc.cores = detectCores()-1))
      PC_dis_clone(data = data_permuted, ID, sample_ID)
    }, mc.cores = detectCores()-1)
    
    return(bind_rows(distance_per) %>% mutate(type="permuted"))
  }
  
  
  distance_full_fun = function(data, ID, sample_ID, n){
    
    df_immu_filtered = same_subtype(data, ID, sample_ID)
    distance_clone = PC_dis_clone(df_immu_filtered, ID, sample_ID)
    distance_per = PC_dist_per(data = df_immu_filtered, ID, sample_ID, n)
    # average over 100 permutations
    mean_distance_per = distance_per %>% group_by(sample, !!sym(ID) , !!sym(sample_ID)) %>%  
      mutate(PC_dis = mean(PC_dis)) %>% ungroup() %>% distinct()
    
    distance_full = rbind(distance_clone, mean_distance_per)
    return(distance_full)
  }

# clone is defined by cdr3_TRA and cdr3_TRB
df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_", file_GSE, ".rds")
for (i in 1:length(file_GSE)) {
  df_immu = readRDS(df_immu_file[i])
  distance_full = distance_full_fun(data = df_immu, ID="ID", sample_ID = "sample_ID", n=100)
  saveRDS(distance_full, file = paste0(folder_path, "/", file_GSE[i], "/output_CD4_CD8/distance_full.rds"))
  print(paste0(i, "done"))
  remove(distance_full)
  remove(df_immu)
  gc()
}

# aa_TRA
df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_", file_GSE, ".rds")
for (i in 1:length(file_GSE)) {
  df_immu = readRDS(df_immu_file[i])
  distance_full = distance_full_fun(data = df_immu, ID="ID_TRA", sample_ID = "sample_ID_TRA", n=100)
  saveRDS(distance_full, file = paste0(folder_path, "/", file_GSE[i], "/output_CD4_CD8/distance_full_TRA.rds"))
  print(paste0(i, "done"))
  remove(distance_full)
  remove(df_immu)
  gc()
}  

# aa_TRB
df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_", file_GSE, ".rds")
for (i in 1:length(file_GSE)) {
  df_immu = readRDS(df_immu_file[i])
  distance_full = distance_full_fun(data = df_immu, ID="ID_TRB", sample_ID = "sample_ID_TRB", n=100)
  saveRDS(distance_full, file = paste0(folder_path, "/", file_GSE[i], "/output_CD4_CD8/distance_full_TRB.rds"))
  print(paste0(i, "done"))
  remove(distance_full)
  remove(df_immu)
  gc()
}

# nt_TRA_TRB
df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_", file_GSE, ".rds")
for (i in 1:length(file_GSE)) {
  df_immu = readRDS(df_immu_file[i])
  distance_full = distance_full_fun(data = df_immu, ID="ID_nt", sample_ID = "sample_ID_nt", n=100)
  saveRDS(distance_full, file = paste0(folder_path, "/", file_GSE[i], "/nt/output_CD4_CD8/distance_full.rds"))
  print(paste0(i, "done"))
  remove(distance_full)
  remove(df_immu)
  gc()
}

# nt_TRA
df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_", file_GSE, ".rds")
for (i in 1:length(file_GSE)) {
  df_immu = readRDS(df_immu_file[i])
  distance_full = distance_full_fun(data = df_immu, ID="ID_TRA", sample_ID = "sample_ID_TRA", n=100)
  saveRDS(distance_full, file = paste0(folder_path, "/", file_GSE[i], "/nt/output_CD4_CD8/distance_full_TRA.rds"))
  print(paste0(i, "done"))
  remove(distance_full)
  remove(df_immu)
  gc()
}

# nt_TRB
df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_", file_GSE, ".rds")
for (i in 1:length(file_GSE)) {
  df_immu = readRDS(df_immu_file[i])
  distance_full = distance_full_fun(data = df_immu, ID="ID_nt_TRB", sample_ID = "sample_ID_nt_TRB", n=100)
  saveRDS(distance_full, file = paste0(folder_path, "/", file_GSE[i], "/nt/output_CD4_CD8/distance_full_TRB.rds"))
  print(paste0(i, "done"))
  remove(distance_full)
  remove(df_immu)
  gc()
}

