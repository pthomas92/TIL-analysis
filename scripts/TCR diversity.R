
setup_session <- function(wd){
  
  library(data.table)
  library(tidyverse)
  library(patchwork)
  library(msa)
  library(igraph)
  library(edgeR)
  library(rlang)
  library(patchwork)
  
  setwd(wd) 
  
}
plotSharing <- function(df, yvar, title, rotate_x = F, .scales = 'fixed'){
  
  if(rotate_x == T){
    x_axis = theme(axis.text.x = element_text(hjust = 1, angle = 45))
  } else {
    x_axis = theme(axis.text.x = element_text())
  }
  
  ggplot(data = df,
         aes(x = !!sym(yvar)))+
    geom_bar(colour = 'black', fill = 'steelblue')+
    facet_wrap(~mouse, scales = .scales)+
    scale_y_log10()+
    labs(title = title,
         x = 'Mice shared by',
         y = 'Count')+
    theme_bw()+
    x_axis
}
plotDiversity <- function(df, yvar, title, ylabel) {
  ggplot(df, aes(x = factor(mouse, levels = c("F0","N2","BALBC","CBA")), y = !!sym(yvar))) +
    geom_boxplot(colour = "steelblue") +
    geom_point(size = 2, position = position_jitter(width = 0.04)) +
    theme_bw() +
    labs(x = "Mouse", y = ylabel, title = title)
}

setup_session('~/OneDrive - University College London/_Leo Post Doc/_Ari-TCR-analysis/')
chain_type = 'alpha'

plots_ = list()
plot_idx = 1

files = read.csv(paste('outputs', chain_type, 'tcrdist', '0_complete-data_pgens_clustered.csv', sep = '/'))
write_loc = paste('outputs', chain_type, 'tcrdist/diversity', sep = '/')

##### plot identical clone sharing between and within mice #####

### Plot identical clone (nt and aa) sharing within mice of the same breed
plots_[[plot_idx]] = patchwork::wrap_plots(list((plotSharing(files, yvar = 'intramouse_score_nt', title = 'Nucleotide', rotate_x = F)),
                           (plotSharing(files, yvar = 'intramouse_score_aa', title = 'Amino Acid', rotate_x = F))
                           ))+
  plot_annotation(title = paste('Exact', chain_type, 'chain sharing within mice')) &
  theme(plot.title = element_text(hjust = 0.5))
plot_idx = plot_idx + 1

### Plot identical clone (nt and aa) sharing within mice longitudinal tumour passages
plots_[[plot_idx]] = patchwork::wrap_plots(list((plotSharing(files, 'nt_present_in', 'Nucleotide', rotate_x = T, .scales = 'free_x')),
                           (plotSharing(files, 'aa_present_in', 'Amino Acid', rotate_x = T, .scales = 'free_x'))
))+
  plot_annotation(title = paste('Exact', chain_type, 'chain sharing within mice')) &
  theme(plot.title = element_text(hjust = 0.5))
plot_idx = plot_idx + 1

#####

##### pivot the file into long format for diversity calculation #####

files = files %>% 
  pivot_longer(cols = c(B1, B2, B3),
               names_to = 'replicate',
               values_to = 'reads') %>% 
  filter(reads > 0)

#####

##### calculate diversity, and scale by the maximum possible diversity (for shannon) #####

max_shannon = sapply(files %>%
                       group_by(mouse, replicate) %>%
                       summarise(count = n()) %>%
                       pull(count),
                     function(x){vegan::diversity(rep(1, times = x))})
max_chao1 = sapply(files %>%
                     group_by(mouse, replicate) %>%
                     summarise(count = n()) %>%
                     pull(count),
                   function(x){vegan::estimateR(rep(1, times = x))[2]})

### TCRdist cluster read proportions (normalised read frequency per cluster; per mouse, read sum per cluster / total)

#### make data
tcrdist_reads = files %>%
  group_by(mouse, replicate, cluster_id) %>%
  summarise(read_sums = sum(reads)) %>% 
  mutate(read_norm = read_sums / sum(read_sums)) %>% 
  summarise(shannon = vegan::diversity(read_norm),
            simpson = vegan::diversity(read_norm, index = 'simpson')) %>%
  ungroup() %>% 
  mutate(normalised_shannon = shannon / max_shannon,
         grouping = 'tcrdist',
         measure = 'reads') %>% 
  dplyr::select(grouping, measure, mouse, replicate, everything())

#### plot
plots_[[plot_idx]] = patchwork::wrap_plots(list((plotDiversity(tcrdist_reads, yvar = 'shannon', title = 'Shannon', ylabel = 'Shannon')),
                           (plotDiversity(tcrdist_reads, yvar = 'normalised_shannon', title = 'Normalised Shannon', ylabel = 'Normalised Shannon')),
                           (plotDiversity(tcrdist_reads, yvar = 'simpson', title = 'Simpson', ylabel = 'Simpson'))
))+
  plot_annotation(title = 'Diversity measurement',
                  subtitle = 'TCRdist clone TMM normalised read sums') &
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
plot_idx = plot_idx + 1

### TCRdist cluster size

#### make data
tcrdist_size = files %>%
  group_by(mouse, replicate, cluster_id) %>%
  summarise(cluster_size = n()) %>% 
  ungroup() %>% 
  group_by(mouse, replicate) %>% 
  summarise(chao1 = vegan::estimateR(cluster_size)[2],
            shannon = vegan::diversity(cluster_size),
            simpson = vegan::diversity(cluster_size, index = 'simpson')) %>%
  ungroup() %>% 
  mutate(normalised_shannon = shannon / max_shannon,
         grouping = 'tcrdist',
         measure = 'size') %>% 
  dplyr::select(grouping, measure, mouse, replicate, everything())

#### plot
plots_[[plot_idx]] = patchwork::wrap_plots(list((plotDiversity(tcrdist_size, yvar = 'shannon', title = 'Shannon', ylabel = 'Shannon')),
                           (plotDiversity(tcrdist_size, yvar = 'normalised_shannon', title = 'Normalised Shannon', ylabel = 'Normalised Shannon')),
                           (plotDiversity(tcrdist_size, yvar = 'simpson', title = 'Simpson', ylabel = 'Simpson'))
))+
  plot_annotation(title = 'Diversity measurement',
                  subtitle = 'TCRdist clone size (no weighting)') &
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
plot_idx = plot_idx + 1

### AA clonotype read proportions

#### make data
clonotype_reads = files %>%
  group_by(mouse, replicate, id_aa) %>%
  summarise(read_sums = sum(reads)) %>% 
  mutate(read_norm = read_sums / sum(read_sums)) %>% 
  summarise(shannon = vegan::diversity(read_norm),
            simpson = vegan::diversity(read_norm, index = 'simpson')) %>%
  ungroup() %>% 
  mutate(normalised_shannon = shannon / max_shannon,
         grouping = 'aa_clonotype',
         measure = 'reads') %>% 
  dplyr::select(grouping, measure, mouse, replicate, everything())

#### plot
plots_[[plot_idx]] = patchwork::wrap_plots(list((plotDiversity(clonotype_reads, yvar = 'shannon', title = 'Shannon', ylabel = 'Shannon')),
                           (plotDiversity(clonotype_reads, yvar = 'normalised_shannon', title = 'Normalised Shannon', ylabel = 'Normalised Shannon')),
                           (plotDiversity(clonotype_reads, yvar = 'simpson', title = 'Simpson', ylabel = 'Simpson'))
))+
  plot_annotation(title = 'Diversity measurement',
                  subtitle = 'Amino acid clone TMM normalised read sums') &
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
plot_idx = plot_idx + 1

### AA clonotype sizes
clonotype_size = files %>%
  group_by(mouse, replicate, id_aa) %>%
  summarise(cluster_size = n()) %>% 
  ungroup() %>% 
  group_by(mouse, replicate) %>% 
  summarise(chao1 = vegan::estimateR(cluster_size)[2],
            shannon = vegan::diversity(cluster_size),
            simpson = vegan::diversity(cluster_size, index = 'simpson')) %>%
  ungroup() %>% 
  mutate(normalised_shannon = shannon / max_shannon,
         grouping = 'aa_clonotype',
         measure = 'size') %>% 
  dplyr::select(grouping, measure, mouse, replicate, everything())

#### plot
plots_[[plot_idx]] = patchwork::wrap_plots(list((plotDiversity(clonotype_size, yvar = 'shannon', title = 'Shannon', ylabel = 'Shannon')),
                           (plotDiversity(clonotype_size, yvar = 'normalised_shannon', title = 'Normalised Shannon', ylabel = 'Normalised Shannon')),
                           (plotDiversity(clonotype_size, yvar = 'simpson', title = 'Simpson', ylabel = 'Simpson'))
))+
  plot_annotation(title = 'Diversity measurement',
                  subtitle = 'Amino Acid clone frequency (no weighting)') &
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

#####

f_ = ls(pattern = '(reads|size)')

for(i in f_){
  
  write.csv(get(i), paste(write_loc, paste(i, '_diversity.csv', sep = ''), sep = '/'), row.names = F)
  
}

pdf(paste(write_loc, 'diversity_plots.pdf', sep = '/'), onefile = T, height = 4)
print(plots_)
dev.off()
