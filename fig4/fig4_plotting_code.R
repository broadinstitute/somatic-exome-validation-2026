library(ggplot2)
library(ggpubr)
library(ggbreak)
library(dplyr)


gbm_variant_overlap_table = read.delim("gbm_variant_overlap_table.txt")
gbm_plotting_data_split_var_type = as.data.frame(gbm_variant_overlap_table[which(gbm_variant_overlap_table$validation_power_wex>0.95),] %>% 
                                                   mutate(is_snp = variant_type == "SNP") %>% group_by(sample, status, is_snp) %>% summarize(n=n()))

for(sample in unique(gbm_plotting_data_split_var_type$sample)){
  for(var_status in c("v2_only", "v6_only", "overlap")){
    if(length(which(gbm_plotting_data_split_var_type$sample == sample &
                    gbm_plotting_data_split_var_type$status == var_status &
                    gbm_plotting_data_split_var_type$is_snp)) == 0){
      gbm_plotting_data_split_var_type[nrow(gbm_plotting_data_split_var_type)+1,] = c(sample, var_status, TRUE, 0)
      gbm_plotting_data_split_var_type$is_snp = as.logical(gbm_plotting_data_split_var_type$is_snp)
    }
    if(length(which(gbm_plotting_data_split_var_type$sample == sample &
                    gbm_plotting_data_split_var_type$status == var_status &
                    !gbm_plotting_data_split_var_type$is_snp)) == 0){
      gbm_plotting_data_split_var_type[nrow(gbm_plotting_data_split_var_type)+1,] = c(sample, var_status, FALSE, 0)
      gbm_plotting_data_split_var_type$is_snp = as.logical(gbm_plotting_data_split_var_type$is_snp)
    }
  }
}

gbm_plotting_data_split_var_type$status = factor(gbm_plotting_data_split_var_type$status, levels=c("v6_only", "v2_only", "overlap"))
gbm_plotting_data_split_var_type$n = as.numeric(gbm_plotting_data_split_var_type$n)

## Panel A Plotting
GBM_Snp_Overlap_Plot = ggplot(gbm_plotting_data_split_var_type[which(gbm_plotting_data_split_var_type$is_snp),], 
                                      aes(x=sample, y=n, fill=status))+geom_col()+labs(x="", y="Powered SNP Count")+
  scale_fill_discrete(name="Status",labels=c("WES_V6_Only", "WES_V2_Only","Shared"))+scale_y_break(c(110,1100), scale=0.5)+
  theme_bw()+theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())+
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20), plot.title = element_text(size=20), 
        legend.text = element_text(size=16), legend.title = element_text(size=20), strip.text = element_text(size=18), 
        axis.title.x = element_text(hjust=0.4))
GBM_InDel_Overlap_Plot = ggplot(gbm_plotting_data_split_var_type[which(!gbm_plotting_data_split_var_type$is_snp),],
                                        aes(x=sample, y=n, fill=status))+geom_col()+labs(x="", y="Powered InDel Count")+
  scale_fill_discrete(name="Status",labels=c("WES_V6_Only", "WES_V2_Only","Shared"))+scale_y_break(c(30,100), scale=0.4)+
  theme_bw()+theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1), axis.ticks.x=element_blank())+
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20), plot.title = element_text(size=20), 
        legend.text = element_text(size=16), legend.title = element_text(size=20), strip.text = element_text(size=18), 
        axis.title.x = element_text(hjust=0.4))

## Panel B Plotting
GBM_variant_overlap_scatterplot = ggplot(gbm_variant_overlap_table[which(gbm_variant_overlap_table$validation_power_wex>=0.95),],
                                  aes(x=v2_alt_frac, y=v6_alt_frac, color=status, shape=status,alpha=status))+
  geom_point()+scale_alpha_manual(values=c(0.7,0.7,0.2),name="Status", labels=c("WES_V6_Only", "WES_V2_Only","Shared"))+
  scale_shape_manual(values=c(16,16,1),name="Status", labels=c("WES_V6_Only", "WES_V2_Only","Shared"))+
  geom_abline(slope=1,color="grey",linetype="dashed")+scale_color_discrete(name="Status", labels=c("WES_V6_Only", "WES_V2_Only","Shared"))+
  theme_bw()+facet_wrap(~var_categ, nrow=1)+labs(x="WES V2 Alt Fraction", y="WES V6 Alt Fraction", title="FFPE GBM Powered Variant Alt Fractions")+
  theme(axis.text = element_text(size=18), axis.title = element_text(size=22), plot.title = element_text(size=22), 
        legend.text = element_text(size=18), legend.title = element_text(size=22), strip.text = element_text(size=20), 
        axis.text.x = element_text(angle=60, vjust=1, hjust=1))+ylim(0,1)+xlim(0,1)

## Panel C Plotting
PCR_scatterplot = ggplot(somatic_wes_pcr_results_plotting, aes(x=v2_alt_frac, y=v6_alt_frac, color=Status, shape=interaction(as.logical(detected), powered)))
+geom_point(size=3)+geom_abline(slope=1, linetype="dashed", color="grey")+labs(x="WES V2 Alt Frac", y="WES V6 Alt Frac")
+scale_shape_manual(name="PCR Result", labels=c("FP (Unpowered)", "TP (Unpowered)", "TP (Powered)"), values=c(4,1,16))
+theme_bw()+theme(axis.text = element_text(size=18), axis.title = element_text(size=22), plot.title = element_text(size=22),
                  legend.text = element_text(size=18), legend.title = element_text(size=22), strip.text = element_text(size=20), 
                  axis.text.x = element_text(angle=60, vjust=1, hjust=1))

## Panel D Plotting
mutlist_to_96_contexts <- function(mutlist, genomeFile) {
  samples <- unique(mutlist$SampleID)
  genomeFile <- genomeFile
  trinuc_mut_mat <- matrix(0, ncol = 96, nrow = length(samples))
  for (n in 1:length(samples)) {
    s <- samples[n]
    mutations <- as.data.frame(mutlist[mutlist$SampleID == s, c("Chr", "Pos", "Ref", "Alt")])
    colnames(mutations) <- c("chr", "pos", "ref", "mut")
    mutations$pos <- as.numeric(mutations$pos)
    mutations$chr <- as.character(mutations$chr)
    mutations <- mutations[(mutations$ref %in% c("A", "C", "G", "T")) & (mutations$mut %in% c("A", "C", "G", "T")) & mutations$chr %in% c(as.character(1:22), "X", "Y"),]
    mutations$trinuc_ref <- as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(as.numeric(mutations$pos) - 1, as.numeric(mutations$pos) + 1))))
    ntcomp <- c(T = "A", G = "C", C = "G", A = "T")
    mutations$sub <- paste(mutations$ref, mutations$mut, sep = ">")
    mutations$trinuc_ref_py <- mutations$trinuc_ref
    for (j in 1:nrow(mutations)) {
      if (mutations$ref[j] %in% c("A", "G")) { # Purine base
        mutations$sub[j] <- paste(ntcomp[mutations$ref[j]], ntcomp[mutations$mut[j]], sep = ">")
        mutations$trinuc_ref_py[j] <- paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j], split = "")[[1]])], collapse = "")
      }
    }
    freqs <- table(paste(mutations$sub, paste(substr(mutations$trinuc_ref_py, 1, 1), substr(mutations$trinuc_ref_py, 3, 3), sep = "-"), sep = ","))
    sub_vec <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    ctx_vec <- paste(rep(c("A", "C", "G", "T"), each = 4), rep(c("A", "C", "G", "T"), times = 4), sep = "-")
    full_vec <- paste(rep(sub_vec, each = 16), rep(ctx_vec, times = 6), sep = ",")
    freqs_full <- freqs[full_vec]
    freqs_full[is.na(freqs_full)] <- 0
    names(freqs_full) <- full_vec
    trinuc_mut_mat[n, ] <- freqs_full
    print(s)
  }
  colnames(trinuc_mut_mat) <- full_vec
  rownames(trinuc_mut_mat) <- samples
  return(trinuc_mut_mat)
}

plot_spectrum_short <- function(freqs_full, sample_id, add_to_title = "") {
  sub_vec <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ctx_vec <- paste(rep(c("A", "C", "G", "T"), each = 4), rep(c("A", "C", "G", "T"), times = 4), sep = "-")
  full_vec <- paste(rep(sub_vec, each = 16), rep(ctx_vec, times = 6), sep = ",")
  freqs_full <- as.numeric(freqs_full)
  names(freqs_full) <- full_vec
  
  xstr <- paste(substr(full_vec, 5, 5), substr(full_vec, 1, 1), substr(full_vec, 7, 7), sep = "")
  
  colvec <- rep(c("dodgerblue", "black", "red", "grey70", "olivedrab3", "plum2"), each = 16)
  y <- freqs_full
  maxy <- max(y)
  h <- barplot(y, las = 2, col = colvec, border = NA, ylim = c(0, maxy * 1.5), space = 1, cex.names = 0.6, names.arg = xstr, ylab = "Number of mutations", main = sample_id)
  
  for (j in 1:length(sub_vec)) {
    xpos <- h[c((j - 1) * 16 + 1, j * 16)]
    rect(xpos[1] - 0.5, maxy * 1.2, xpos[2] + 0.5, maxy * 1.3, border = NA, col = colvec[j * 16])
    text(x = mean(xpos), y = maxy * 1.3, pos = 3, label = sub_vec[j])
  }
  
  # Add the number of mutations to the upper left corner
  mtext(paste0("Number of mutations: ", sum(freqs_full), add_to_title), side = 3, adj = 0, line = 0.5, cex = 1)
  
  return(h)
}

plot_mut_context = function(v2_muts_path, v6_muts_path, shared_muts_path, ref_fasta){
  gbm_contexts_v2 = mutlist_to_96_contexts(read.delim(v2_muts_path), ref_fasta)
  gbm_contexts_v6 = mutlist_to_96_contexts(read.delim(v6_muts_path), ref_fasta)
  gbm_contexts_agg = mutlist_to_96_contexts(read.delim(shared_muts_path), ref_fasta)
  
  pdf("gbm_spectrum_barplot_wide.pdf", width=8, height=10)
  nf <- layout( matrix(c(1,2,3), ncol=1, byrow = TRUE))
  plot_spectrum_short(t(as.matrix(apply(gbm_contexts_agg, 2, sum))), "WES Shared Variants\n")
  plot_spectrum_short(t(as.matrix(apply(gbm_contexts_v6, 2, sum))), "WES V6-Only Variants\n")
  plot_spectrum_short(t(as.matrix(apply(gbm_contexts_v2, 2, sum))), "WES V2-Only Variants\n")
  dev.off()
  
  return()
}
