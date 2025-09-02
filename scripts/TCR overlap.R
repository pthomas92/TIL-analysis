
setup_session <- function(wd){
  
  library(tidyverse)
  library(patchwork)
  
  setwd(wd) 
  
}

setup_session('~/OneDrive - University College London/_Leo Post Doc/_Ari-TCR-analysis/')
chain_type = 'alpha'

readData = function(summary_path, base_path){
  cluster_summary = read.csv(summary_path) %>% 
    mutate(mice_shared_count = sapply(strsplit(mice_shared, ','), length))
  base_clusters = read.csv(base_path)
  return(list(cluster_summary = cluster_summary, base_clusters = base_clusters))
}
heatmapSummary = function(input, split_on, title){
  
  cluster_members <- input %>%
    separate_rows(!!split_on, sep = ",") %>%
    mutate(value = T) %>%   
    pivot_wider(names_from = !!split_on, 
                values_from = value, 
                values_fill = F) %>% 
    pivot_longer(cols = c(BALBC, F0, N2, CBA),
                 names_to = 'mouse',
                 values_to = 'present')
  
  plt = ggplot(cluster_members, aes(x = factor(mouse, levels = c('F0', 'N2', 'BALBC', 'CBA')),
                                   y = factor(cluster_id), fill = factor(present)))+
    geom_tile()+
    labs(x = "Mouse", y = "Cluster",
         fill = 'Present',
         title = title)+
    scale_fill_manual(values = c('grey', 'firebrick'))+
    theme_bw()+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  return(list(data = cluster_members,
              plot = plt))
  
}
jaccard = function(input, col_idx){
  
  ji = list()
  counter_outer = 1
  
  for(i in col_idx){
    holder = list()
    counter_inner = 1
    
    for(j in col_idx){
      
      test  = as.numeric(unlist(input[,i]))
      query = as.numeric(unlist(input[,j]))
      
      intersection_size = sum(test == 1 & query == 1)
      union_size = sum(test == 1 | query == 1)
      
      result = data.frame(
        source = colnames(input)[i],
        dest   = colnames(input)[j],
        jaccard_index = ifelse(i != j, intersection_size / union_size, NA)
      )
      holder[[counter_inner]] = result
      counter_inner = counter_inner + 1
      
    }
    ji[[counter_outer]] = do.call('rbind', holder)
    counter_outer = counter_outer + 1
  }
  
  df <- do.call('rbind', ji)
  return(df)
}
calcJaccard = function(input, split_on, title){
  
  data = input %>%
    separate_rows(!!split_on, sep = ",") %>%
    mutate(value = 1) %>%   
    pivot_wider(names_from = !!split_on, 
                values_from = value, 
                values_fill = 0)
  
  cols_ = ncol(input):ncol(data)
  jaccard_ = jaccard(data, cols_)
  
  p1 = ggplot(data = jaccard_, aes(x = factor(source, levels = c('F0', 'N2', 'BALBC', 'CBA')),
                        y = factor(dest, levels = c('F0', 'N2', 'BALBC', 'CBA')),
                        fill = jaccard_index))+
    geom_tile(colour = 'black', linewidth = 1)+
    theme_bw()+
    scale_fill_gradient(low = 'white', high = 'forestgreen')+
    labs(x = '', y = '', fill = 'Jaccard Index')
  
  p2 = ggplot(data = jaccard_, aes(x = factor(source, levels = c('F0', 'N2', 'BALBC', 'CBA')),
                                 y = jaccard_index))+
    geom_boxplot()+
    geom_point(aes(colour = dest), size = 3, position = position_jitter(width = 0.2))+
    theme_bw()+
    labs(x = 'Mouse', y = 'Jaccard Index', colour = 'Mouse overlap')
  
  plt = patchwork::wrap_plots(list(p1, p2))+
    plot_annotation(title = title) &
    theme(plot.title = element_text(size = 15, hjust = 0.5))
  
  return(list(data = jaccard_,
              plot = plt))
  
}

summary_path = paste('outputs', chain_type, 'tcrdist/cluster_sizes/0_complete-data_cluster-summary-stats.csv', sep = '/')
base_path = paste('outputs', chain_type, 'tcrdist/0_complete-data_pgens_clustered.csv', sep = '/')

data = readData(summary_path, base_path)
cluster_summary = data$cluster_summary %>% filter(!is.na(cluster_id))
base_clusters = data$base_clusters
rm(data)

heatmaps = list()
heatmaps[[1]] = heatmapSummary(cluster_summary, split_on = 'mice_shared', title = 'TCRdist cluster presence')
heatmaps[[2]] = heatmapSummary(base_clusters %>% select(-mouse), split_on = 'aa_present_in', title = 'AA clonotype presence')
names(heatmaps) = c('TCRdist', 'AA_clonotype')

jaccard_data = list()
jaccard_data[[1]] = calcJaccard(input = cluster_summary, split_on = 'mice_shared', title = 'TCRdist cluster Jaccard Index')
jaccard_data[[2]] = calcJaccard(input = base_clusters, split_on = 'aa_present_in', title = 'AA clonotype Jaccard Index')
names(jaccard_data) = c('TCRdist', 'AA_clonotype')

output_path = paste('outputs', chain_type, 'tcrdist/cluster_overlaps', sep = '/')
source_name = sapply(strsplit(basename(summary_path), '_'), function(x) paste(x[1:2], collapse = '_'))

for(i in 1:length(heatmaps)){
  write.csv(heatmaps[[i]]$data,
            paste(output_path, paste(source_name, names(heatmaps)[i], 'presence-absence_data.csv', sep = '_'), sep = '/'))
}

for(i in 1:length(jaccard_data)){
  write.csv(jaccard_data[[i]]$data,
            paste(output_path, paste(source_name, names(jaccard_data)[i], 'jaccard-index_data.csv', sep = '_'), sep = '/'))
  
  pdf(file = paste(output_path, paste(source_name, names(jaccard_data)[i], 'jaccard-index_plots.pdf', sep = '_'), sep = '/'),
      height = 5, width = 7)
    print(jaccard_data[[i]]$plot)
  dev.off()
  
}

pdf(file = paste(output_path, paste(source_name, names(heatmaps)[i], 'presence-absence_heatmaps.pdf', sep = '_'), sep = '/'),
    height = 5, width = 7)
  print(patchwork::wrap_plots(list(heatmaps[[1]]$plot,
                                   heatmaps[[2]]$plot)))
dev.off()
