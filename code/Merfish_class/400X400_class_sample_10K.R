setwd('/picb/neurosys/chenrenrui/Mouse_ST/A molecularly defined and spatially resolved cell atlas of the whole mouse brain/')
.libPaths('/data/neurosys-svr2/chenrenrui/lib/')

library(vroom)
library(stringr)
library(data.table)
library(doParallel)
library(foreach)
library(dplyr)




start_time.t <- Sys.time()
file_names <- list.files('./data_filter/neocortex/')

for (file_name in file_names) {
  Merfish_mouse_neocortex_layer_region <- read.delim('./Merfish_mouse_neocortex_layer_region.txt',header = T)
  Merfish_brain_cell_type <- read.delim('./Merfish_brain_cell_type.txt',header = T)
  
  ##Neocortex
  file_chose_names <- str_split(file_name,".txt",simplify = T)[,1]
  
  load(paste(paste('./class_1K/neocortex/400X400_20/class/cell_window/',file_chose_names,sep = ''),'RData',sep = '.'))
  
  
  start_time <- Sys.time()
  all_all_result <- list()
  layer_number = 1
  for (layer.1 in sort(unique(file_tmp$layer))) {
    file_tmp_sample <- file_tmp[layer ==  layer.1]
    cell_windows_layer_d_sample <- cell_windows_layer_d[layer ==layer.1]
    #setDT(cell_windows_layer_d_sample)
    all_all_result_s <- list()
    
    temp1 <- file_tmp_sample[,c('cell_label','cell_Neuron_type')]
    temp1$cell_label <- temp1$cell_label
    # 替换原有的合并操作
    setkey(temp1, cell_label)
    setkey(cell_windows_layer_d_sample, cell_label)
    temp2 <- merge(temp1, cell_windows_layer_d_sample)
    temp2$num <- 1
    
    # 使用数据表的聚合操作进行分组计算
    temp3 <- temp2[, .(sum_value = sum(num)), by = .(cell_Neuron_type, loc)]
    
    all_all_result_s[[1]] <- temp3
    all_all_result[[layer_number]] <- all_all_result_s
    names(all_all_result)[layer_number] <- layer.1
    layer_number = layer_number + 1
    
  }
  
  # 设置并行计算的核心数
  cl <- makeCluster(40)
  registerDoParallel(cl)
  
  sample_num =10000
  
  all_all_result_sample <- list()
  layer_number = 1
  for (layer.1 in sort(unique(file_tmp$layer))) {
    file_tmp_sample <- file_tmp[layer == layer.1]
    cell_windows_layer_d_sample <- cell_windows_layer_d[layer == layer.1]
    #setDT(cell_windows_layer_d_sample)
    
    all_all_result_s  <- foreach(i = 1:sample_num, .combine = c,.packages = c('data.table','stringr')) %dopar% {
      file_tmp_sample$cell_Neuron_type <- sample(file_tmp_sample$cell_Neuron_type,nrow(file_tmp_sample))
      temp1 <- file_tmp_sample[,c('cell_label','cell_Neuron_type')]
      temp1$cell_label <- temp1$cell_label
      # 替换原有的合并操作
      setkey(temp1, cell_label)
      setkey(cell_windows_layer_d_sample, cell_label)
      temp2 <- merge(temp1, cell_windows_layer_d_sample)
      temp2$num <- 1
      
      # 使用数据表的聚合操作进行分组计算
      temp3 <- temp2[, .(sum_value = sum(num)), by = .(cell_Neuron_type, loc)]
      
      return(list(temp3))
    }
    all_all_result_sample[[layer_number]] <- all_all_result_s
    names(all_all_result_sample)[layer_number] <- layer.1
    layer_number = layer_number + 1
  }
  
  end_time <- Sys.time()
  execution_time <- difftime(end_time, start_time, units = "mins")
  execution_time  
  
  save(all_all_result,all_all_result_sample,file=paste(paste('./class_1K/neocortex/400X400_20/class/cell_sample/ ',file_chose_names,sep = ''),'RData',sep = '.'))
  
  stopCluster(cl)
  
  print(paste(paste('./class_1K/neocortex/400X400_20/class/cell_sample/',file_chose_names,sep = ''),'done',execution_time,sep = ' ')) 
  
  
  rm(temp1)
  rm(temp2)
  rm(temp3)
  gc()
  
  
  
  
  
  start_time <- Sys.time()
  cell_Sliding_window_result_p <- list()
  a = 1
  for (i in names(all_all_result)) {
    temp.t <- rbindlist(all_all_result[[i]])
    temp.s <- rbindlist(all_all_result_sample[[i]])
    temp.t$type <- 'TRUE'
    temp.s$type <- 'sample'
    
    temp1 <- rbind(temp.t,temp.s)
    temp1 <- temp1[, rank := frank(-sum_value,ties.method = "last"), by = .(cell_Neuron_type, loc)]
    temp.1.rank <- temp1[type=='TRUE']
    temp.1.rank$p <- temp.1.rank$rank/(sample_num+1)
    
    temp.1.rank$layer <- i
    cell_Sliding_window_result_p[[a]] <- temp.1.rank
    names(cell_Sliding_window_result_p)[a] <- i
    a = a + 1
    #print(i)
  }
  end_time <- Sys.time()
  execution_time <- difftime(end_time, start_time, units = "mins")
  execution_time 
  
  print(paste(paste('./class_1K/neocortex/400X400_20/class/cell_pvalue/',file_chose_names,sep = ''),'done',execution_time,sep = ' ')) 
  
  cell_Sliding_window_result_p_d <- rbindlist(cell_Sliding_window_result_p)
  
  write.csv(cell_Sliding_window_result_p_d,file = paste('./class_1K/neocortex/400X400_20/class/cell_pvalue/',file_chose_names,'_','cell_sliding_window_result_p.csv',sep=''))
  
  cell_Sliding_window_result_p_d_sig_high <- cell_Sliding_window_result_p_d[cell_Sliding_window_result_p_d$p<0.05,]
  
  write.csv(cell_Sliding_window_result_p_d_sig_high,
            file = paste('./class_1K/neocortex/400X400_20/class/cell_pvalue/',file_chose_names,'_',
                         'cell_sliding_window_result_p_sig_high.csv',sep=''))
  
  rm(all_all_result)
  rm(all_all_result_s)
  rm(all_all_result_sample)
  rm(cell_Sliding_window_result_p)
  rm(cell_Sliding_window_result_p_d)
  rm(cell_Sliding_window_result_p_d_sig_high)
  rm(cell_windows_layer_d)
  rm(cell_windows_layer_d_sample)
  rm(temp.1.rank)
  rm(temp.s)
  rm(temp.t)
  rm(temp1)
  rm(file_tmp)
  rm(file_tmp_sample)
  gc()
}
end_time.t <- Sys.time()
execution_time <- difftime(end_time.t, start_time.t, units = "mins")
execution_time  
print(paste('total',execution_time,sep = ' '))


