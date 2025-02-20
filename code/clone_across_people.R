##### Get across people GEX dissimilarity for clones with same TCR sequence via
##### integrative results from harmony.

library(stringr)
library(readr)
library(readxl)
library(purrr)
library(dplyr)
library(tidyr)
library(tibble)
library(viridis)
library(ggsci)
library(ggplot2)
library(ggpattern)
library(ggpubr)
library(parallel)
library(ggpval)

folder_path = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis"
file_list = list.files(folder_path)
file_GSE = file_list[str_detect(file_list, "GSE|EMTAB9357_pbmc|PRJCA002413_pbmc|10X|Zhang_2020")]

#df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_imputed_anno.rds")
df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_", file_GSE, ".rds")
df_immu_all = df_immu_file %>% map_dfr(readRDS)

sample_information_file = paste0(folder_path, "/", file_GSE,"/sample_information_", file_GSE, ".xlsx")
sample_information = sample_information_file %>% 
  map_dfr(~read_excel(.) %>% mutate(age = as.character(age))) %>%
  select(sample, ind)

df_immu_all = df_immu_all %>% merge(sample_information)

#clone with more than one person 
clone = df_immu_all %>%
  count(cdr3_TRA, cdr3_TRB, ind) %>%
  count(cdr3_TRA, cdr3_TRB, name = "clone_count") %>%
  filter(clone_count > 1) %>%
  arrange(desc(clone_count))

df_immu_all_clone = df_immu_all %>%
  inner_join(clone, by = c("cdr3_TRA", "cdr3_TRB"))

# celltype percentage for each individual and clone
subtype_prop = df_immu_all_clone %>%
  group_by(cdr3_TRA, cdr3_TRB, ind) %>%
  count(subtype) %>% 
  mutate(prop = n/sum(n)) %>%
  ungroup()
  
mix_clone = subtype_prop %>%
  filter(prop != 1) %>%
  distinct(cdr3_TRA, cdr3_TRB) %>%
  inner_join(subtype_prop, by=c("cdr3_TRA", "cdr3_TRB"))

pure_clone = anti_join(subtype_prop, mix_clone)
TRA = pure_clone %>% distinct(cdr3_TRA, cdr3_TRB) %>% pull(cdr3_TRA)
TRB = pure_clone %>% distinct(cdr3_TRA, cdr3_TRB) %>% pull(cdr3_TRB)
subtype_num = c()
for (i in 1:length(TRA)) {
  subtype_num[i] = pure_clone %>% 
    filter(cdr3_TRA == TRA[i] &
             cdr3_TRB == TRB[i]) %>%
    pull(subtype) %>%
    unique() %>%
    length()
}
mixtype_idx = which(subtype_num == 2)
pure_clone_mixtype = pure_clone %>% 
  inner_join(data.frame(cdr3_TRA = TRA[mixtype_idx],
                        cdr3_TRB = TRB[mixtype_idx]),
             by = c("cdr3_TRA", "cdr3_TRB"))

sample_pair = function(i){
  df = df_immu_all_clone %>%
    filter(cdr3_TRA == TRA[i] & cdr3_TRB == TRB[i]) %>%
    select(sample, ind)
  
  sam = df$sample
  ind = df$ind

  sam_comb = t(combn(sam, 2)) %>% 
    as.data.frame() %>%
    cbind(t(combn(ind, 2)))
  colnames(sam_comb) = c("sam1", "sam2", "ind1", "ind2")
  sam_comb = sam_comb %>% filter(ind1 != ind2) %>% distinct()
  
  return(sam_comb)
}

#sam_pair = mclapply(1:length(TRA), sample_pair, mc.cores = 23) %>% bind_rows() %>% distinct()
#saveRDS(sam_pair, "~/projects/Research/TCR and Gene Expression/PC distance analysis/across_sam_pair.rds")
across_sam_pair <- readRDS("~/projects/Research/TCR and Gene Expression/PC distance analysis/across_sam_pair.rds")

##### across_sam_pair is used in integrative_analysis.R for avoiding redundant 
##### calculation of harmony since we only considers people who have "across people clones"
---------------------------------------------------------------

harmony_files = list.files("~/projects/Research/TCR and Gene Expression/PC distance analysis/harmony")

df_permuted = function(data){
  idx = sample(1:nrow(data))
  data[, grep("harmony", colnames(data))] = data[idx, grep("harmony", colnames(data))]
  data
}

harmony_distance = function(i, permute=FALSE){
  df = df_immu_all_clone %>% 
    filter(cdr3_TRA == TRA[i] &
             cdr3_TRB == TRB[i])
  
  df_sample_pair = sample_pair(i)
  df_dis_mean = data.frame()
  for (j in 1:nrow(df_sample_pair)) {
    
    if (permute == FALSE){
      harmony = readRDS(paste0(folder_path,"/harmony/har_", df_sample_pair$sam1[j], 
                               "_", df_sample_pair$sam2[j], ".rds")) %>%
        select(sample_barcode, paste0("harmony_", 1:10))
    } else if (permute == TRUE) {
      harmony = readRDS(paste0(folder_path,"/harmony/har_", df_sample_pair$sam1[j], 
                               "_", df_sample_pair$sam2[j], ".rds")) %>%
        select(sample_barcode, paste0("harmony_", 1:10)) %>%
        df_permuted()
    }
    
    df_harmony = df %>% 
      mutate(sample_barcode = paste0(sample, "_", barcode)) %>%
      merge(harmony, by="sample_barcode")
    
    df_harmony_sam1 = df_harmony %>% filter(sample == df_sample_pair$sam1[j]) %>% select(paste0("harmony_", 1:10))
    df_harmony_sam2 = df_harmony %>% filter(sample == df_sample_pair$sam2[j]) %>% select(paste0("harmony_", 1:10))
    dis = as.vector(proxy::dist(df_harmony_sam1, df_harmony_sam2))
    dis_mean = mean(dis[dis!=0])
    df_dis_mean = rbind(df_dis_mean,
                        data.frame(cdr3_TRA = TRA[i],
                             cdr3_TRB = TRB[i],
                             sam1 = df_sample_pair$sam1[j],
                             sam2 = df_sample_pair$sam2[j],
                             harmony_dis = dis_mean))
  }
  return(df_dis_mean)
}

harmony_dis = mclapply(1:length(TRA), harmony_distance, mc.cores = 23) %>% 
  bind_rows() %>% mutate(type = "real data")
harmony_dis_per = mclapply(1:length(TRA), harmony_distance, permute=TRUE, mc.cores = 23) %>% 
  bind_rows() %>% mutate(type = "permuted data")

harmony = rbind(harmony_dis, harmony_dis_per)
harmony$type = factor(harmony$type, levels = c("real data", "permuted data"))

harmony_pure_same = harmony %>% 
  anti_join(data.frame(cdr3_TRA = TRA[mixtype_idx],
                        cdr3_TRB = TRB[mixtype_idx]),
             by = c("cdr3_TRA", "cdr3_TRB")) %>%
  left_join(pure_clone %>%
              select(cdr3_TRA, cdr3_TRB, subtype) %>%
              distinct())
harmony_pure_same$subtype = factor(harmony_pure_same$subtype,
                                   levels = c("CD4+CD8-", "CD4-CD8+"))


pdf(paste0(folder_path, "/plots_paper/figure1/across_people_figure1e.pdf"), width = 7, height = 3)
ggplot(data = harmony_pure_same) +
  geom_density(aes(harmony_dis, fill = type)) +
  facet_grid(.~subtype) +
  scale_fill_manual(values = alpha(c('pink', 'skyblue'), alpha = 0.5),
                    name = "")  +
  theme_bw() +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.box.spacing = unit(0.05, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        strip.background =element_rect(fill="white", colour = "black", linewidth = rel(2)),
        panel.border = element_rect(colour = "black")) +
  labs(x="GEX dissimilarity", y="Density")
dev.off()



#nucleotide

folder_path = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis"
file_list = list.files(folder_path)
file_GSE = file_list[str_detect(file_list, "GSE|EMTAB9357_pbmc|PRJCA002413_pbmc|10X|Zhang_2020")]

df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_", file_GSE, ".rds")
df_immu_all = df_immu_file %>% map_dfr(readRDS)

sample_information_file = paste0(folder_path, "/", file_GSE,"/sample_information_", file_GSE, ".xlsx")
sample_information = sample_information_file %>% 
  map_dfr(~read_excel(.) %>% mutate(age = as.character(age))) %>%
  select(sample, ind)

df_immu_all = df_immu_all %>% merge(sample_information)

#clone with more than one person 
clone = df_immu_all %>%
  count(cdr3_nt_TRA, cdr3_nt_TRB, ind) %>%
  count(cdr3_nt_TRA, cdr3_nt_TRB, name = "clone_count") %>%
  filter(clone_count > 1) %>%
  arrange(desc(clone_count))

df_immu_all_clone = df_immu_all %>%
  inner_join(clone, by = c("cdr3_nt_TRA", "cdr3_nt_TRB"))

# celltype percentage for each individual and clone
subtype_prop = df_immu_all_clone %>%
  group_by(cdr3_nt_TRA, cdr3_nt_TRB, ind) %>%
  count(subtype) %>% 
  mutate(prop = n/sum(n)) %>%
  ungroup()

mix_clone = subtype_prop %>%
  filter(prop != 1) %>%
  distinct(cdr3_nt_TRA, cdr3_nt_TRB) %>%
  inner_join(subtype_prop, by=c("cdr3_nt_TRA", "cdr3_nt_TRB"))

pure_clone = anti_join(subtype_prop, mix_clone)
TRA = pure_clone %>% distinct(cdr3_nt_TRA, cdr3_nt_TRB) %>% pull(cdr3_nt_TRA)
TRB = pure_clone %>% distinct(cdr3_nt_TRA, cdr3_nt_TRB) %>% pull(cdr3_nt_TRB)
subtype_num = c()
for (i in 1:length(TRA)) {
  subtype_num[i] = pure_clone %>% 
    filter(cdr3_nt_TRA == TRA[i] &
             cdr3_nt_TRB == TRB[i]) %>%
    pull(subtype) %>%
    unique() %>%
    length()
}
mixtype_idx = which(subtype_num == 2)
pure_clone_mixtype = pure_clone %>% 
  inner_join(data.frame(cdr3_nt_TRA = TRA[mixtype_idx],
                        cdr3_nt_TRB = TRB[mixtype_idx]),
             by = c("cdr3_nt_TRA", "cdr3_nt_TRB"))

sample_pair = function(i){
  df = df_immu_all_clone %>%
    filter(cdr3_nt_TRA == TRA[i] & cdr3_nt_TRB == TRB[i]) %>%
    select(sample, ind)
  
  sam = df$sample
  ind = df$ind
  
  sam_comb = t(combn(sam, 2)) %>% 
    as.data.frame() %>%
    cbind(t(combn(ind, 2)))
  colnames(sam_comb) = c("sam1", "sam2", "ind1", "ind2")
  sam_comb = sam_comb %>% filter(ind1 != ind2) %>% distinct()
  
  return(sam_comb)
}

#sam_pair = mclapply(1:length(TRA), sample_pair, mc.cores = 23) %>% bind_rows() %>% distinct()
#saveRDS(sam_pair, "~/projects/Research/TCR and Gene Expression/PC distance analysis/across_sam_pair_nt.rds")

across_sam_pair <- readRDS("~/projects/Research/TCR and Gene Expression/PC distance analysis/across_sam_pair_nt.rds")
harmony_files = list.files("~/projects/Research/TCR and Gene Expression/PC distance analysis/harmony")

df_permuted = function(data){
  idx = sample(1:nrow(data))
  data[, grep("harmony", colnames(data))] = data[idx, grep("harmony", colnames(data))]
  data
}

harmony_distance = function(i, permute=FALSE){
  df = df_immu_all_clone %>% 
    filter(cdr3_nt_TRA == TRA[i] &
             cdr3_nt_TRB == TRB[i])
  
  df_sample_pair = sample_pair(i)
  df_dis_mean = data.frame()
  for (j in 1:nrow(df_sample_pair)) {
    
    if (permute == FALSE){
      harmony = readRDS(paste0(folder_path,"/harmony/har_", df_sample_pair$sam1[j], 
                               "_", df_sample_pair$sam2[j], ".rds")) %>%
        select(sample_barcode, paste0("harmony_", 1:10))
    } else if (permute == TRUE) {
      harmony = readRDS(paste0(folder_path,"/harmony/har_", df_sample_pair$sam1[j], 
                               "_", df_sample_pair$sam2[j], ".rds")) %>%
        select(sample_barcode, paste0("harmony_", 1:10)) %>%
        df_permuted()
    }
    
    df_harmony = df %>% 
      mutate(sample_barcode = paste0(sample, "_", barcode)) %>%
      merge(harmony, by="sample_barcode")
    
    df_harmony_sam1 = df_harmony %>% filter(sample == df_sample_pair$sam1[j]) %>% select(paste0("harmony_", 1:10))
    df_harmony_sam2 = df_harmony %>% filter(sample == df_sample_pair$sam2[j]) %>% select(paste0("harmony_", 1:10))
    dis = as.vector(proxy::dist(df_harmony_sam1, df_harmony_sam2))
    dis_mean = mean(dis[dis!=0])
    df_dis_mean = rbind(df_dis_mean,
                        data.frame(cdr3_nt_TRA = TRA[i],
                                   cdr3_nt_TRB = TRB[i],
                                   sam1 = df_sample_pair$sam1[j],
                                   sam2 = df_sample_pair$sam2[j],
                                   harmony_dis = dis_mean))
  }
  return(df_dis_mean)
}

harmony_dis = mclapply(1:length(TRA), harmony_distance, mc.cores = 23) %>% 
  bind_rows() %>% mutate(type = "real data")
harmony_dis_per = mclapply(1:length(TRA), harmony_distance, permute=TRUE, mc.cores = 23) %>% 
  bind_rows() %>% mutate(type = "permuted data")

harmony = rbind(harmony_dis, harmony_dis_per)
harmony$type = factor(harmony$type, levels = c("real data", "permuted data"))

harmony_pure_same = harmony %>% 
  anti_join(data.frame(cdr3_nt_TRA = TRA[mixtype_idx],
                       cdr3_nt_TRB = TRB[mixtype_idx]),
            by = c("cdr3_nt_TRA", "cdr3_nt_TRB")) %>%
  left_join(pure_clone %>%
              select(cdr3_nt_TRA, cdr3_nt_TRB, subtype) %>%
              distinct())


pdf(paste0(folder_path, "/plots_paper/figure1/across_people_figure1e_nt.pdf"), width = 2.5, height = 5)
ggplot(data = harmony_pure_same) +
  geom_boxplot(aes(x=type, y=harmony_dis, fill = type),outlier.shape = NA) +
  labs(x="", y="GEX dissimilarity") +
  scale_fill_manual(values = alpha(c('pink', 'skyblue'), alpha = 0.5)) +
  scale_y_continuous(limits = c(0,40)) +
  theme_bw() +
  facet_grid(subtype ~.) +
  theme(panel.border = element_rect(colour = "black", linewidth = 0.8),
        strip.background =element_rect(fill="white"),
        legend.position="none")
dev.off()

