
setup_session <- function(wd){
  
  library(data.table)
  library(tidyverse)
  library(patchwork)
  library(msa)
  library(igraph)
  library(edgeR)
  library(rlang)
  library(fitdistrplus)
  library(e1071)
  
  setwd(wd) 
  
}
setup_session('~/OneDrive - University College London/_Leo Post Doc/_Ari-TCR-analysis/')

getClusterSizes = function(dm, cut_at){
  
  cluster_sizes = dm %>% group_by(cluster_id) %>% summarise(count = n()) %>% arrange(desc(count))
  
  var = try(expr = {
    fitdist(cluster_sizes$count, distr = "gamma", method = "mle")
  })
  if(as.character(summary(var)[2]) == 'try-error'){
    
    plt = ggplot(cluster_sizes, aes(x = count)) +
      geom_histogram(binwidth = 1, 
                     fill = "lightgrey", color = "black")+
      theme_bw()+
      theme(text = element_text(size = 12))+
      labs(x = 'Cluster size',
           y = 'Counts',
           title = paste('Cluster size distribution (could not fit line)'),
           subtitle = paste('TCRdist threshold', cut_at, sep = ': '))
    
  } else {
    fit_ <- fitdist(cluster_sizes$count, distr = "gamma", method = "mle")
    shape <- fit_$estimate["shape"]
    rate  <- fit_$estimate["rate"]
    binwidth = 1
    x_vals <- seq(min(cluster_sizes$count),
                  max(cluster_sizes$count),
                  length.out = 500)
    n <- length(cluster_sizes$count)
    y_vals <- dgamma(x_vals, shape = shape, rate = rate) * n * binwidth
    
    fit_df <- data.frame(x = x_vals, y = y_vals)
  
  plt = ggplot(cluster_sizes, aes(x = count)) +
    geom_histogram(binwidth = binwidth, 
                   fill = "lightgrey", color = "black") +
    geom_line(data = fit_df, aes(x = x, y = y), 
              color = "steelblue", size = 1)+
    theme_bw()+
    theme(text = element_text(size = 12))+
    labs(x = 'Cluster size',
         y = 'Counts',
         title = paste('Cluster size distribution'),
         subtitle = paste('TCRdist threshold', cut_at, sep = ': '))
  }
  
  return(list(cluster_sizes = cluster_sizes,
              plot = plt))
  
}
makeClusterLogo = function(data, cluster_name, cdr3_col, trim5 = T, trim3 = T){
  
  tmp_cluster = data %>% filter(cluster_id == cluster_name)
  seqs = as.character(unlist(tmp_cluster[,cdr3_col]))
  seqs = sapply(strsplit(seqs, ''),
                function(y){
                  paste(y[ifelse(trim5 == T, 2, 1):ifelse(trim3 == T, length(y)-1, length(y))], collapse = '')
                  })
  
  msa_ = msa(inputSeqs = seqs,
             method = 'ClustalOmega', type = 'protein')
  
  aligned_seqs = paste(ifelse(trim5 == T, 'C', ''),
                       as.character(unmasked(msa_)),
                       ifelse(trim3 == T, 'F', ''),
                       sep = '')  # extract aligned strings
  plt = ggseqlogo::ggseqlogo(aligned_seqs)+ggtitle(paste('cluster', cluster_name))
  
  fasta = paste(paste('>', tmp_cluster$id_nt, sep = ''), aligned_seqs, sep = '\n')
  
  return(list(fasta = fasta,
              plot = plt))
  
}

chain = 'alpha'
write_aln = T

basepath = paste('outputs', chain, 'tcrdist', sep = '/')
files = basename(list.files(basepath, pattern = '.csv'))
filesets = unique(gsub('_.*?$', '', files))

for(i in 1:length(filesets)){

  dynamic_path = paste(basepath, grep(paste(paste('^', filesets[i], sep = ''), 'pgens', sep = '.*?'), files, value = T), sep = '/')
  dist_path = paste(basepath, grep(paste(paste('^', filesets[i], sep = ''), 'tcrdistmatrix', sep = '.*?'), files, value = T), sep = '/')
  fasta_out = gsub('^(.*?tcrdist)/(.*?)_pgens.csv', '\\1/alignments/\\2.fasta', dynamic_path)
  
  dynamic_data = read.csv(dynamic_path)
  dynamic_distmat = read.csv(dist_path)
  dynamic_distmat = dynamic_distmat[,2:ncol(dynamic_distmat)]
  rownames(dynamic_distmat) = colnames(dynamic_distmat) = gsub('\\.', '-', colnames(dynamic_distmat))
  
  cluster_id = hclust(as.dist(dynamic_distmat))
  cluster_id = cutree(cluster_id, h = 55)
  dynamic_data = left_join(dynamic_data, as.data.frame(cluster_id) %>% rownames_to_column('id_nt'))
  cluster_sizes = getClusterSizes(dynamic_data, cut_at = 55)
  
  write.csv(dynamic_data, file = gsub('.csv', '_clustered.csv', dynamic_path))
  
  pdf(paste(paste(basepath, 'cluster_sizes', sep = '/'),
            gsub('pgens.csv', 'cluster-size-histogram.pdf', grep(paste(filesets[i], 'pgens', sep = '.*?'), files, value = T)),
            sep = '/'))
    print(cluster_sizes$plot)
  dev.off()
  
  cluster_sizes = cluster_sizes$cluster_sizes
  
  cluster_sizes = left_join(cluster_sizes, 
    lapply(cluster_sizes$cluster_id, function(x){
      d_ = dynamic_data %>% filter(cluster_id == x)
      return(
        data.frame(cluster_id = x,
                   mice_shared = paste(unique(d_$mouse), collapse = ','),
                   v_genes = paste(unique(d_[,grep('^v_\\w_gene', colnames(d_))]), collapse = ','),
                   j_genes = paste(unique(d_[,grep('^j_\\w_gene', colnames(d_))]), collapse = ','),
                   unique_cdr3s = length(unique(d_[,grep('^cdr3_\\w_aa', colnames(d_))])),
                   mean_pgen = mean(d_[,grep('^pgen', colnames(d_))]),
                   pgen_skew = tryCatch(
                     skewness(d_[,grep('^pgen', colnames(d_))], type = 2),
                     error = function(e) NA_real_ )
                   )
        )
  }) %>% bind_rows()) %>% dplyr::select(cluster_id, count, everything()) %>% dplyr::rename('cluster_count' = count)
  
  write.csv(cluster_sizes,
            paste(paste(basepath, 'cluster_sizes', sep = '/'),
                  gsub('pgens.csv', 'cluster-summary-stats.csv',
                       grep(paste(paste('^', filesets[i], sep = ''), 'pgens', sep = '.*?'), files, value = T)),
                  sep = '/'))
  
  if(length(which(cluster_sizes$cluster_count >= 10)) > 0){
    clusters_select = as.numeric(cluster_sizes %>% filter(cluster_count >= 10) %>% pull(cluster_id))
  } else {
    clusters_select = as.numeric(cluster_sizes %>% slice_head(n = 10) %>% pull(cluster_id))
  }
  for(cluster_test in clusters_select[!is.na(clusters_select)]){
    
    aln = makeClusterLogo(data = dynamic_data,
                          cluster_name = cluster_test,
                          cdr3_col = ifelse(chain == 'beta', 'cdr3_b_aa', 'cdr3_a_aa'),
                          trim5 = F, trim3 = F)
    
    if(write_aln == T){
      aln_name = gsub('.fasta', paste('-cluster', cluster_test, '.fasta', sep = ''), fasta_out)
      img_name = gsub('.fasta', '.pdf', gsub('alignments', 'alignments/img', aln_name))
      write_lines(aln$fasta, aln_name, sep = '\n')
      
      pdf(img_name, height = 4, width = 8)
        print(aln$plot)
      dev.off()
      
    }
  }
}
