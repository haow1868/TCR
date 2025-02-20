#Permutation between clones with the same size (comparison between clones)
library(dplyr)
library(stringr)
library(purrr)
library(data.table)
folder_path = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis"
file_list = list.files(folder_path)
file_GSE = file_list[str_detect(file_list, "GSE|EMTAB9357_pbmc|PRJCA002413_pbmc|10X|Zhang_2020")]

df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_", file_GSE, ".rds")
df_immu_all = df_immu_file %>% map_dfr(readRDS)

Random_permuted = function(data, file_GSE, ID, sample_ID, nucleotide, subtype, clone_type){
  
  if(file_GSE == "GSE187515" & subtype ==  "CD8") {
    #GSE187515 does not have CD8+ T cells
    folder = paste0(file_GSE, "/output_", subtype, "/", clone_type)
    saveRDS(data.frame(), file = paste0(folder_path, "/", folder, "/Ran_sample_median.rds"))
  } else {
    total_num = df_immu_all %>% 
      dplyr::count(sample, !!sym(ID), !!sym(sample_ID), name="total_num")
    setDT(data)
    setDT(total_num)
    
    # First merge (on sample_ID.x)
    PC_distance_lv <- merge(data, total_num, by.x = paste0(sample_ID, ".x"), by.y = sample_ID)
    setnames(PC_distance_lv, "total_num", "total_num.x")
    
    # Second merge (on sample_ID.y)
    PC_distance_lv <- merge(PC_distance_lv, total_num, by.x = paste0(sample_ID, ".y"), by.y = sample_ID)
    setnames(PC_distance_lv, "total_num", "total_num.y")
    
    Ran_PC_distance_lv = PC_distance_lv[, c(.SD, .(PC_dis_ran = as.numeric(sample(PC_dis, .N)))), 
                                        by=.(total_num.x, total_num.y, sample)] %>% as.data.frame()
    Ran_PC_distance_lv = Ran_PC_distance_lv[!duplicated(as.list(Ran_PC_distance_lv))]
    
    # sample median
    Ran_sample_median = 
      Ran_PC_distance_lv %>% group_by(sample, distance) %>%
      dplyr::summarise(median = median(PC_dis_ran),
                       .groups = "drop")
    if(nucleotide == F) {
      folder = paste0(file_GSE, "/output_", subtype, "/", clone_type)
    } else if (nucleotide == T) {
      folder = paste0(file_GSE, "/nt/output_", subtype, "/", clone_type)
    }
    
    saveRDS(Ran_sample_median, file = paste0(folder_path, "/", folder, "/Ran_sample_median.rds"))
  }
  
}


#Amino acid
for (i in 1:length(file_GSE)) {
  #TRA
  PC_distance_lv_CD4_file = paste0(folder_path, "/", file_GSE, "/output_CD4/TRA/PC_distance_lv.rds")
  PC_distance_lv_CD8_file = paste0(folder_path, "/", file_GSE, "/output_CD8/TRA/PC_distance_lv.rds")
  
  PC_distance_lv_CD4 = readRDS(PC_distance_lv_CD4_file[i])
  Random_permuted(data = PC_distance_lv_CD4, file_GSE[i], ID="ID_TRA",
                  sample_ID = "sample_ID_TRA", nucleotide = F, subtype="CD4", clone_type = "TRA")

  PC_distance_lv_CD8 = readRDS(PC_distance_lv_CD8_file[i])
  Random_permuted(data = PC_distance_lv_CD8, file_GSE[i], ID="ID_TRA",
                    sample_ID = "sample_ID_TRA", nucleotide = F, subtype="CD8", clone_type = "TRA")
  #TRB
  PC_distance_lv_CD4_file = paste0(folder_path, "/", file_GSE, "/output_CD4/TRB/PC_distance_lv.rds")
  PC_distance_lv_CD8_file = paste0(folder_path, "/", file_GSE, "/output_CD8/TRB/PC_distance_lv.rds")
  
  PC_distance_lv_CD4 = readRDS(PC_distance_lv_CD4_file[i])
  Random_permuted(data = PC_distance_lv_CD4, file_GSE[i], ID="ID_TRB",
                  sample_ID = "sample_ID_TRB", nucleotide = F, subtype="CD4", clone_type = "TRB")
  
  PC_distance_lv_CD8 = readRDS(PC_distance_lv_CD8_file[i])
  Random_permuted(data = PC_distance_lv_CD8, file_GSE[i], ID="ID_TRB",
                  sample_ID = "sample_ID_TRB", nucleotide = F, subtype="CD8", clone_type = "TRB")
  #TRA+TRB
  PC_distance_lv_CD4_file = paste0(folder_path, "/", file_GSE, "/output_CD4/TRA+TRB/PC_distance_lv.rds")
  PC_distance_lv_CD8_file = paste0(folder_path, "/", file_GSE, "/output_CD8/TRA+TRB/PC_distance_lv.rds")
  
  PC_distance_lv_CD4 = readRDS(PC_distance_lv_CD4_file[i])
  Random_permuted(data = PC_distance_lv_CD4, file_GSE[i], ID="ID",
                  sample_ID = "sample_ID", nucleotide = F, subtype="CD4", clone_type = "TRA+TRB")
  
  PC_distance_lv_CD8 = readRDS(PC_distance_lv_CD8_file[i])
  Random_permuted(data = PC_distance_lv_CD8, file_GSE[i], ID="ID",
                  sample_ID = "sample_ID", nucleotide = F, subtype="CD8", clone_type = "TRA+TRB")
  
  remove(PC_distance_lv_CD4)
  remove(PC_distance_lv_CD8)
  gc()
  print(paste0(i,"done"))
}

#nucleotide
for (i in 1:length(file_GSE)) {
  #TRA
  PC_distance_lv_CD4_file = paste0(folder_path, "/", file_GSE, "/nt/output_CD4/TRA/PC_distance_lv.rds")
  PC_distance_lv_CD8_file = paste0(folder_path, "/", file_GSE, "/nt/output_CD8/TRA/PC_distance_lv.rds")
  
  PC_distance_lv_CD4 = readRDS(PC_distance_lv_CD4_file[i])
  Random_permuted(data = PC_distance_lv_CD4, file_GSE[i], ID="ID_nt_TRA",
                  sample_ID = "sample_ID_nt_TRA", nucleotide = T, subtype="CD4", clone_type = "TRA")
  
  PC_distance_lv_CD8 = readRDS(PC_distance_lv_CD8_file[i])
  Random_permuted(data = PC_distance_lv_CD8, file_GSE[i], ID="ID_nt_TRA",
                  sample_ID = "sample_ID_nt_TRA", nucleotide = T, subtype="CD8", clone_type = "TRA")
  #TRB
  PC_distance_lv_CD4_file = paste0(folder_path, "/", file_GSE, "/nt/output_CD4/TRB/PC_distance_lv.rds")
  PC_distance_lv_CD8_file = paste0(folder_path, "/", file_GSE, "/nt/output_CD8/TRB/PC_distance_lv.rds")
  
  PC_distance_lv_CD4 = readRDS(PC_distance_lv_CD4_file[i])
  Random_permuted(data = PC_distance_lv_CD4, file_GSE[i], ID="ID_nt_TRB",
                  sample_ID = "sample_ID_nt_TRB", nucleotide = T, subtype="CD4", clone_type = "TRB")
  
  PC_distance_lv_CD8 = readRDS(PC_distance_lv_CD8_file[i])
  Random_permuted(data = PC_distance_lv_CD8, file_GSE[i], ID="ID_nt_TRB",
                  sample_ID = "sample_ID_nt_TRB", nucleotide = T, subtype="CD8", clone_type = "TRB")
  #TRA+TRB
  PC_distance_lv_CD4_file = paste0(folder_path, "/", file_GSE, "/nt/output_CD4/TRA+TRB/PC_distance_lv.rds")
  PC_distance_lv_CD8_file = paste0(folder_path, "/", file_GSE, "/nt/output_CD8/TRA+TRB/PC_distance_lv.rds")
  
  PC_distance_lv_CD4 = readRDS(PC_distance_lv_CD4_file[i])
  Random_permuted(data = PC_distance_lv_CD4, file_GSE[i], ID="ID_nt",
                  sample_ID = "sample_ID_nt", nucleotide = T, subtype="CD4", clone_type = "TRA+TRB")
  
  PC_distance_lv_CD8 = readRDS(PC_distance_lv_CD8_file[i])
  Random_permuted(data = PC_distance_lv_CD8, file_GSE[i], ID="ID_nt",
                  sample_ID = "sample_ID_nt", nucleotide = T, subtype="CD8", clone_type = "TRA+TRB")
  
  remove(PC_distance_lv_CD4)
  remove(PC_distance_lv_CD8)
  gc()
  print(paste0(i,"done"))
}










