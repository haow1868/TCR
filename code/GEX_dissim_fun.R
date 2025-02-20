library(parallel)
library(dplyr)
library(tidyr)
library(tibble)
library(R.utils)
library(stringr)
library(data.table)
library(Rfast)
library(stringdist)

combination = function(vec){
  combs <- combn(vec, 2)
  self_combs <- combn(vec, 1, FUN=function(x) c(x, x))
  final_combs <- t(cbind(combs, self_combs)) %>% as.data.frame() %>% 
    mutate(across(where(is.numeric), as.character))
  return(final_combs)
}

# Levenshtein distance
lv_dist = function(data, ID, sample_ID, nucleotide, TRA, TRB, cdr3){
  
  # Convert string arguments to symbols
  ID_sym = sym(ID)
  sample_ID_sym = sym(sample_ID)
  cdr3_sym = sym(cdr3)
  
  ID_multiple = data %>% 
    group_by(sample, !!ID_sym, !!sample_ID_sym) %>%
    summarise(total_num = n(),
              .groups = "drop") %>%
    filter(total_num != 1) %>%
    select(!!ID_sym) %>% distinct() %>% pull()
  
  if (nucleotide == T & TRA == T & TRB == T) {
    df_target = distinct(data, cdr3_nt_TRA, cdr3_nt_TRB, !!ID_sym, sample) 
  } else if (nucleotide == F & TRA == T & TRB == T) {
    df_target = distinct(data, cdr3_TRA, cdr3_TRB, !!ID_sym, sample)
  } else {
    df_target = distinct(data, !!cdr3_sym, !!ID_sym, sample)
  }
  df_target = df_target %>% mutate(across(where(is.numeric), as.character))
  
  # all pairwise clones
  ID_comb = combination(df_target %>% pull(!!ID_sym))
  colnames(ID_comb) = c(paste0(as.character(ID_sym), ".x"), paste0(as.character(ID_sym), ".y"))
  
  dist_df <- ID_comb %>%
    left_join(df_target, by = setNames(ID, paste0(ID, ".x"))) %>%
    left_join(df_target, by = setNames(c(ID, "sample"), c(paste0(ID, ".y"), "sample")))
  
  # drop clones with single cell when distance = 0 (cannot calculate dissimilarity)
  dist_df <- dist_df[!(dist_df[[paste0(ID, ".x")]] == dist_df[[paste0(ID, ".y")]] & 
                         !(dist_df[[paste0(ID, ".x")]] %in% ID_multiple)), ]
  # LV distance
  if (nucleotide == T & TRA == T & TRB == T) {
    dist_df$distance_TRA = stringdist(dist_df[["cdr3_nt_TRA.x"]], dist_df[["cdr3_nt_TRA.y"]], method = "lv")
    dist_df$distance_TRB = stringdist(dist_df[["cdr3_nt_TRB.x"]], dist_df[["cdr3_nt_TRB.y"]], method = "lv")
    dist_df$distance = dist_df$distance_TRA + dist_df$distance_TRB
  } else if (nucleotide == F & TRA == T & TRB == T) {
    dist_df$distance_TRA = stringdist(dist_df[["cdr3_TRA.x"]], dist_df[["cdr3_TRA.y"]], method = "lv")
    dist_df$distance_TRB = stringdist(dist_df[["cdr3_TRB.x"]], dist_df[["cdr3_TRB.y"]], method = "lv")
    dist_df$distance = dist_df$distance_TRA + dist_df$distance_TRB
  } else {
    dist_df$distance = stringdist(dist_df[[paste0(cdr3, ".x")]], dist_df[[paste0(cdr3, ".y")]], method = "lv")
  }
  return(dist_df)
}


# GEX dissimilarity
PC_dis_lv = function(data, lv_dis_data, ID){
  
  data_list = split(data, data[[ID]])
  data_list = lapply(data_list, function(df){df %>% select(starts_with("PC"))})
  n_cores = detectCores() - 4
  
  # Chunk size depending on the number of available cores
  chunk_size = ceiling(nrow(lv_dis_data)/n_cores)
  # Split data into chunks
  chunks = split(lv_dis_data, ceiling(seq_len(nrow(lv_dis_data))/chunk_size))
  
  results_list = mclapply(chunks, function(chunk) {
    PC_distance_list = vector("list", nrow(chunk))
    
    for (i in seq_len(nrow(chunk))) {
      id_1 = chunk[i,][[paste0(ID, ".x")]]
      id_2 = chunk[i,][[paste0(ID, ".y")]]
      
      df_lv_ID1 = data_list[[id_1]] %>% as.matrix()
      df_lv_ID2 = data_list[[id_2]] %>% as.matrix()
      
      PC_distance_pairwise = as.vector(Rfast::dista(df_lv_ID1, df_lv_ID2))
      PC_distance_list[[i]] = mean(PC_distance_pairwise[PC_distance_pairwise != 0])
    }
    
    chunk$PC_dis = unlist(PC_distance_list)
    chunk$type = "estimated"
    return(chunk)
  }, mc.cores = n_cores)
  
  # Combine results
  lv_dis_data = do.call(rbind, results_list)
  
  return(lv_dis_data)
}


