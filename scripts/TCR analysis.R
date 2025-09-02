
setup_session <- function(wd){
  
  library(data.table)
  library(tidyverse)
  library(patchwork)
  library(msa)
  library(igraph)
  library(edgeR)
  library(rlang)
  
  setwd(wd) 
  
}

setup_session('~/OneDrive - University College London/_Leo Post Doc/_Ari-TCR-analysis/')
chain_type = 'beta'
write_files = T

readData = function(read_params){
  read_params = read.delim(read_params)
  files = list()
  for(i in 1:nrow(read_params)){
    files[[i]] = read.delim(read_params$path[i]) %>% 
      mutate(
        mouse = read_params$mouse[i],
        tumour = read_params$tumour[i],
        chain = read_params$chain[i]
      )
  }
  return(files)
}
filterNormaliseData = function(x, chain_filter){
  normalised_data = do.call('rbind', x) %>% 
    filter(chain == chain_filter) %>% 
    mutate(mouse = ifelse(mouse == 'BALB_C', 'BALBC', mouse)) %>% 
    group_by(mouse, tumour, chain, v_call, j_call, junction, junction_aa) %>% 
    summarise(duplicate_count = round(mean(duplicate_count), digits = 0)) %>% 
    ungroup() %>% 
    group_by(mouse, tumour, chain, v_call, j_call, junction, junction_aa) %>% 
    summarise(reads = sum(duplicate_count)) %>% 
    mutate(tcrs_million = 1e6 * (reads / sum(.$reads))) %>% 
    TMM_norm() %>% 
    ungroup() %>% 
    filter(if_any(everything(), ~ .x > 0)) %>% 
    filter(B1 > 0 | B2 > 0 | B3 > 0)
  return(normalised_data)
}
TMM_norm = function(x){
  
  counts = x %>% 
    pivot_wider(id_cols = c(
      chain, v_call, j_call, junction, junction_aa
    ),
    names_from = c(mouse, tumour),   # still separates replicates
    values_from = reads,
    names_sep = '_',
    values_fill = 0
    )
  
  count_ids = counts[, 1:6]
  counts = counts[,7:ncol(counts)]
  
  y = edgeR::DGEList(counts = counts)
  y = edgeR::calcNormFactors(y, method = 'TMM')
  
  cpm_norm <- edgeR::cpm(y, normalized.lib.sizes = TRUE)
  y = cbind(count_ids, cpm_norm)
  
  # keep replicate-level data, but extract mouse (timepoint) correctly
  y = y %>% 
    filter(junction != '') %>% 
    pivot_longer(cols = grep('_B\\d$', colnames(y), value = T),
                 names_to = 'mouse_tumour',
                 values_to = 'reads') %>% 
    mutate(mouse = sapply(strsplit(mouse_tumour, '_'), function(x) x[1]), # first part = timepoint
           replicate = sapply(strsplit(mouse_tumour, '_'), function(x) x[2]), # B1/B2/B3
           id_nt = paste(chain, v_call, j_call, junction, sep = '--'),
           id_aa = paste(chain, v_call, j_call, junction_aa, sep = '--')) %>% 
    select(-mouse_tumour) %>% 
    pivot_wider(id_cols = c(mouse, chain, v_call, j_call, junction, junction_aa, id_nt, id_aa),
                names_from = replicate, values_from = reads)
  
  return(y)
  
}
findConvergence = function(df, type){
  if(type == 'intra'){
    # find intramouse (type) convergence
    y.1 = df %>% 
      mutate(B1_bin = ifelse(B1 > 0, 1, 0),
             B2_bin = ifelse(B2 > 0, 1, 0),
             B3_bin = ifelse(B3 > 0, 1, 0),
             intramouse_score_nt = B1_bin + B2_bin + B3_bin) %>% 
      select(-ends_with('bin')) %>% 
      arrange(desc(intramouse_score_nt))
    
    y.2 = df %>% 
      group_by(id_aa, mouse) %>% 
      summarise(B1 = sum(B1),
                B2 = sum(B2),
                B3 = sum(B3)) %>% 
      mutate(B1_bin = ifelse(B1 > 0, 1, 0),
             B2_bin = ifelse(B2 > 0, 1, 0),
             B3_bin = ifelse(B3 > 0, 1, 0),
             intramouse_score_aa = B1_bin + B2_bin + B3_bin) %>% 
      select(-ends_with('bin')) %>% 
      arrange(desc(intramouse_score_aa)) %>% 
      select(mouse, id_aa, intramouse_score_aa)
    
    df = left_join(y.1, y.2) %>% arrange(desc(intramouse_score_nt), desc(intramouse_score_nt), id_aa)
  } else if(type == 'inter'){
    
    cross_mouse_tcrs = df %>% 
      group_by(id_nt, id_aa) %>% 
      summarise(count = n()) %>% 
      arrange(desc(count)) %>% 
      filter(count > 1) %>% 
      pull(id_nt)
    
    y.1 = df %>%
      group_by(id_nt, id_aa) %>%
      summarise(
        nt_present_in = paste(sort(unique(mouse)), collapse = ","),
        intermouse_score_nt = n_distinct(mouse),
        .groups = "drop"
      )
    
    y.2 = df %>%
      group_by(id_aa) %>%
      summarise(
        aa_present_in = paste(sort(unique(mouse)), collapse = ","),
        intermouse_score_aa = n_distinct(mouse),
        .groups = "drop"
      )
    
    df = left_join(df, left_join(y.1, y.2))
  }
  return(df)
}
summariseAllData = function(input_data, summary_plot = 'aa'){
  df = findConvergence(input_data, type = 'intra')
  df = findConvergence(df, type = 'inter')
  
  if(summary_plot == 'aa'){
    cols = c('aa_present_in', 'intermouse_score_aa')
  } else if(summary_plot == 'nt'){
    cols = c('nt_present_in', 'intermouse_score_nt')
  }
  
  plot_data = df %>% 
    group_by_at(cols) %>% 
    summarise(count = n()) %>%
    arrange(desc(count))
  colnames(plot_data)[1:2] = c('x', 'colour')
  
  plt = ggplot(data = plot_data,
               aes(x = factor(x, levels = plot_data$x),
                   y = count,
                   fill = factor(colour)))+
    geom_col(colour = 'black')+
    theme_bw()+
    labs(y = 'Count',
         x = 'Mouse group(s)',
         fill = 'Overlap size',
         title = ifelse(grepl('aa', cols), 'AA sequence overlap', 'Nucleotide sequence overlap'))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
    geom_hline(yintercept = 1)
  
  return(
    list(
      data = df,
      plot = plt
    )
  )
}


files = readData('repertoires.txt')
normalised_data = filterNormaliseData(x = files, chain_filter = chain_type)
normalised_data = summariseAllData(normalised_data)

if(write_files == T){

  normalised_data$data %>% 
    write_csv(paste('outputs', chain_type, '0_complete-data.csv', sep = '/'))
  
  normalised_data$data %>% 
    filter(aa_present_in == 'F0') %>% 
    write_csv(paste('outputs', chain_type, '1.0_F0_only.csv', sep = '/'))
  
  normalised_data$data %>% 
    filter(aa_present_in == 'N2') %>% 
    write_csv(paste('outputs', chain_type, '1.1_N2_only.csv', sep = '/'))
  
  normalised_data$data %>% 
    filter(aa_present_in == 'BALBC') %>% 
    write_csv(paste('outputs', chain_type, '1.2_BALBC_only.csv', sep = '/'))
  
  normalised_data$data %>% 
    filter(aa_present_in == 'CBA') %>% 
    write_csv(paste('outputs', chain_type, '1.3_CBA_only.csv', sep = '/'))
  
  normalised_data$data %>% 
    filter(aa_present_in == 'F0,N2') %>% 
    write_csv(paste('outputs', chain_type, '2.0_F0_N2.csv', sep = '/'))
  
  normalised_data$data %>% 
    filter(aa_present_in == 'N2,BALBC') %>% 
    write_csv(paste('outputs', chain_type, '2.1_BALBC_N2.csv', sep = '/'))
  
  normalised_data$data %>% 
    filter(aa_present_in == 'BALBC,CBA') %>% 
    write_csv(paste('outputs', chain_type, '2.2_BALBC_CBA.csv', sep = '/'))
  
  normalised_data$data %>% 
    filter(aa_present_in == 'BALBC,F0,N2') %>% 
    write_csv(paste('outputs', chain_type, '3.1_BALBC_F0_N2.csv', sep = '/'))
  
  normalised_data$data %>% 
    filter(aa_present_in == 'BALBC,CBA,N2') %>% 
    write_csv(paste('outputs', chain_type, '3.2_BALBC_CBA_N2.csv', sep = '/'))
  
  normalised_data$data %>% 
    filter(aa_present_in == 'BALBC,CBA,F0,N2') %>% 
    write_csv(paste('outputs', chain_type, '4.0_BALBC_CBA_F0_N2.csv', sep = '/'))

  pdf(paste('outputs', chain_type, 'inter-mouse-overlap.pdf', sep = '/'))
    print(normalised_data$plot)
  dev.off()
    
}  
