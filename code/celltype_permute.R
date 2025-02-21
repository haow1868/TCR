##### Permutation via cell subtype (for calculating clone purity)

library(stringr)
library(purrr)
library(dplyr)
library(viridis)
library(ggsci)
library(ggplot2)
library(ggpattern)
library(ggpubr)
library(parallel)

folder_path = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis"
file_list = list.files(folder_path)
file_GSE = file_list[str_detect(file_list, "GSE|EMTAB9357_pbmc|PRJCA002413_pbmc|10X|Zhang_2020")]

sample_information_file = paste0(folder_path, "/", file_GSE,"/sample_information_", file_GSE, ".csv")
sample_information = sample_information_file %>% map_dfr(read.csv)

df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_", file_GSE, ".rds")
df_immu_all = df_immu_file %>% map_dfr(readRDS)
df_immu_all_list = split(df_immu_all, df_immu_all$sample)
# permute celltype
df_permuted = function(data){
  idx = sample(1:nrow(data))
  data[, "subtype"] = data[idx, "subtype"]
  data
}

df_immu_all_per = list()
for (i in 1:100) {
  df_immu_all_per[[i]] = lapply(df_immu_all_list, df_permuted) %>% bind_rows()
  print(paste0(i,"done"))
}
saveRDS(df_immu_all_per, file = paste0(folder_path, "/data_allStudy/df_immu_all_per.rds"))

# Calculating purity for each clone
clone_prop = function(df, ID, sample_ID){
  
  # Convert string arguments to symbols
  ID_sym = sym(ID)
  sample_ID_sym = sym(sample_ID)
  
  total_num = df_immu_all %>% 
    dplyr::count(sample, !!sym(ID), !!sym(sample_ID), name="total_num")
  
  sub_num = df %>% 
    dplyr::count(sample, !!ID_sym, !!sample_ID_sym, subtype) %>%
    merge(total_num) %>%
    filter(total_num > 1) %>%
    mutate(per = n/total_num) %>%
    mutate(diff = abs(per-0.5)) 
  
  sub_num = sub_num %>%
    dplyr::mutate(category = case_when(
      total_num == 1 ~ "Single",
      total_num >= 2  & total_num < 11  ~ "Small",
      total_num >= 11 & total_num < 100 ~ "Medium",
      total_num >= 100               ~ "Large"
    ))
  
  sub_num = sub_num %>% 
    dplyr::mutate(clone_subtype = ifelse(per == 1, "Same celltype", "Mixed celltype")) %>%
    dplyr::mutate(clone_subtype = case_when(
      clone_subtype == "Same celltype" ~ subtype,
      clone_subtype == "Mixed celltype" ~ "Mixture",
    )) 
  
  sub_num = sub_num %>%
    group_by(!!sample_ID_sym) %>%
    mutate(max_prop = max(per)) %>%
    distinct(sample, !!ID_sym, !!sample_ID_sym, category, clone_subtype, total_num, max_prop) %>%
    ungroup()
  return(sub_num)
}


df_immu_all_per = readRDS(paste0(folder_path, "/data_allStudy/df_immu_all_per.rds"))
# aa
# TRA
purity_perm_list = mclapply(df_immu_all_per, clone_prop, 
                             ID="ID_TRA", sample_ID="sample_ID_TRA", mc.cores = 23)
purity_perm = purity_perm_list %>% bind_rows()
purity_perm = purity_perm %>%
  group_by(sample, ID_TRA, sample_ID_TRA, category, total_num) %>%
  summarise(max_prop = mean(max_prop), .groups = "drop") %>%
  mutate(type = "Permuted")
saveRDS(purity_perm, paste0(folder_path, "/data_allStudy/aa/purity_perm_TRA.rds"))

# TRB
purity_perm_list = mclapply(df_immu_all_per, clone_prop, 
                             ID="ID_TRB", sample_ID="sample_ID_TRB", mc.cores = 23)
purity_perm = purity_perm_list %>% bind_rows()
purity_perm = purity_perm %>%
  group_by(sample, ID_TRB, sample_ID_TRB, category, total_num) %>%
  summarise(max_prop = mean(max_prop), .groups = "drop") %>%
  mutate(type = "Permuted")
saveRDS(purity_perm, paste0(folder_path, "/data_allStudy/aa/purity_perm_TRB.rds"))

# TRA+TRB
purity_perm_list = mclapply(df_immu_all_per, clone_prop, 
                             ID="ID", sample_ID="sample_ID", mc.cores = 23)
purity_perm = purity_perm_list %>% bind_rows()
purity_perm = purity_perm %>%
  group_by(sample, ID, sample_ID, category, total_num) %>%
  summarise(max_prop = mean(max_prop), .groups = "drop") %>%
  mutate(type = "Permuted")
saveRDS(purity_perm, paste0(folder_path, "/data_allStudy/aa/purity_perm.rds"))

# nucleotide
# TRA
purity_perm_list = mclapply(df_immu_all_per, clone_prop, 
                             ID="ID_nt_TRA", sample_ID="sample_ID_nt_TRA", mc.cores = 23)
purity_perm = purity_perm_list %>% bind_rows()
purity_perm = purity_perm %>%
  group_by(sample, ID_nt_TRA, sample_ID_nt_TRA, category, total_num) %>%
  summarise(max_prop = mean(max_prop), .groups = "drop") %>%
  mutate(type = "Permuted")
saveRDS(purity_perm, paste0(folder_path, "/data_allStudy/nt/purity_perm_TRA.rds"))

# TRB
purity_perm_list = mclapply(df_immu_all_per, clone_prop, 
                             ID="ID_nt_TRB", sample_ID="sample_ID_nt_TRB", mc.cores = 23)
purity_perm = purity_perm_list %>% bind_rows()
purity_perm = purity_perm %>%
  group_by(sample, ID_nt_TRB, sample_ID_nt_TRB, category, total_num) %>%
  summarise(max_prop = mean(max_prop), .groups = "drop") %>%
  mutate(type = "Permuted")
saveRDS(purity_perm, paste0(folder_path, "/data_allStudy/nt/purity_perm_TRB.rds"))

# TRA+TRB
purity_perm_list = mclapply(df_immu_all_per, clone_prop, 
                             ID="ID_nt", sample_ID="sample_ID_nt", mc.cores = 23)
purity_perm = purity_perm_list %>% bind_rows()
purity_perm = purity_perm %>%
  group_by(sample, ID_nt, sample_ID_nt, category, total_num) %>%
  summarise(max_prop = mean(max_prop), .groups = "drop") %>%
  mutate(type = "Permuted")
saveRDS(purity_perm, paste0(folder_path, "/data_allStudy/nt/purity_perm.rds"))


# Real purity
# aa
# TRA
purity = clone_prop(df_immu_all, ID = "ID_TRA", sample_ID = "sample_ID_TRA") %>% 
  mutate(type = "Estimated")
saveRDS(purity, paste0(folder_path, "/data_allStudy/aa/purity_TRA.rds"))

# TRB
purity = clone_prop(df_immu_all, ID = "ID_TRB", sample_ID = "sample_ID_TRB") %>% 
  mutate(type = "Estimated")
saveRDS(purity, paste0(folder_path, "/data_allStudy/aa/purity_TRB.rds"))

# TRA+TRB
purity = clone_prop(df_immu_all, ID = "ID", sample_ID = "sample_ID") %>% 
  mutate(type = "Estimated")
saveRDS(purity, paste0(folder_path, "/data_allStudy/aa/purity.rds"))

# nucleotide
# TRA
purity = clone_prop(df_immu_all, ID = "ID_nt_TRA", sample_ID = "sample_ID_nt_TRA") %>% 
  mutate(type = "Estimated")
saveRDS(purity, paste0(folder_path, "/data_allStudy/nt/purity_TRA.rds"))

# TRB
purity = clone_prop(df_immu_all, ID = "ID_nt_TRB", sample_ID = "sample_ID_nt_TRB") %>% 
  mutate(type = "Estimated")
saveRDS(purity, paste0(folder_path, "/data_allStudy/nt/purity_TRB.rds"))

# TRA+TRB
purity = clone_prop(df_immu_all, ID = "ID_nt", sample_ID = "sample_ID_nt") %>% 
  mutate(type = "Estimated")
saveRDS(purity, paste0(folder_path, "/data_allStudy/nt/purity.rds"))





