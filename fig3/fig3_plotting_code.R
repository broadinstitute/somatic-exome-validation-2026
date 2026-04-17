library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)

tableau20 = c("#1F77B4","#AEC7E8","#FF7F0E","#FFBB78","#2CA02C",
              "#98DF8A","#D62728","#FF9896","#9467BD","#C5B0D5",
              "#8C564B","#C49C94","#E377C2","#F7B6D2","#7F7F7F",
              "#C7C7C7","#BCBD22","#DBDB8D","#17BECF","#9EDAE5")
## Pooled Sample Variant Recall Boxplot
read_summary_table = function(v2_summary, v6_summary){
  v2_summary_df = read.delim(v2_summary, skip=1, blank.lines.skip = TRUE)[c(1:6,9:14,17:22),]
  v2_summary_df$replicate = rep(c("A","A","B","B","C","C"),3)
  v2_summary_df$plex = rep(c("10plex", "20plex", "5plex"), each=6)
  v2_summary_df$pipeline = "WES V2"
  v6_summary_df = read.delim(v6_summary, skip=1, blank.lines.skip = TRUE)[c(1:6,9:14,17:22),]
  v6_summary_df$replicate = rep(c("A","A","B","B","C","C"),3)
  v6_summary_df$plex = rep(c("10plex", "20plex", "5plex"), each=6)
  v6_summary_df$pipeline = "WES V6"
  
  sensitivity_summary = rbind(v2_summary_df, v6_summary_df)
  sensitivity_summary$TP = as.numeric(sensitivity_summary$TP) 
  sensitivity_summary$FP = as.numeric(sensitivity_summary$FP)
  sensitivity_summary$FN = as.numeric(sensitivity_summary$FN)
  sensitivity_summary$RECALL = as.numeric(sensitivity_summary$RECALL)
  sensitivity_summary$PRECISION = as.numeric(sensitivity_summary$PRECISION)
  
  all_plex_summary = as.data.frame(sensitivity_summary %>% group_by(replicate, pipeline, type) %>% summarize(TP=sum(TP), FP=sum(FP), FN=sum(FN)))
  all_plex_summary$plex = "all_plex"
  all_plex_summary$RECALL = all_plex_summary$TP/(all_plex_summary$TP+all_plex_summary$FN)
  all_plex_summary$PRECISION =  all_plex_summary$TP/(all_plex_summary$TP+all_plex_summary$FP)
  
  sensitivity_summary = rbind(sensitivity_summary, all_plex_summary[colnames(sensitivity_summary)])
  sensitivity_summary$pipeline = factor(sensitivity_summary$pipeline, levels=c("WES V6", "WES V2"))
  sensitivity_summary$type = factor(sensitivity_summary$type, levels=c("SNP", "INDEL"))
  return(sensitivity_summary)
}

label_snv = c("SNV", "InDel")
names(label_snv) = c("SNP", "INDEL")

ggplot(sensitivity_summary) + geom_boxplot(aes(x=plex, y=RECALL, fill=pipeline, color=pipeline), alpha=0.7, linewidth=0.2) +
  geom_point(aes(x=plex, y=RECALL, color=pipeline), size=2, position=position_dodge(0.75), alpha=0.7) +
  labs(x="") + theme_bw() + facet_wrap(~type, labeller = labeller(type=label_snv), ncol=1)+
  scale_color_manual(name="Workflow", labels=c("WES V6", "WES V2"), values=c("#F8766D","#00BA38"))+
  scale_fill_manual(name="Workflow", labels=c("WES V6", "WES V2"), values=c("#F8766D","#00BA38"))+
  labs(x="HapMap Sample", y="Recall", title="HapMap Pooled Sample\nVariant Recall")+
  theme(axis.text = element_text(size=16),axis.text.x = element_text(angle = 60, vjust=1, hjust=1), 
        axis.title = element_text(size=20), plot.title = element_text(size=20), legend.text = element_text(size=16), 
        legend.title = element_text(size=20), strip.text = element_text(size=18), legend.position = "left")

## Pooled Sample Detection Sensitivity

autovalRawToMelt = function(autoval_raw_table, var_type, af_override=NULL){
  stopifnot(toupper(var_type) %in% c("SNP","SNV","INDEL"))
  if(toupper(var_type) %in% c("SNP","SNV")){
    af_cols = paste0('X', c(0.05,0.1,0.2,0.4,0.8,'1.0'))
  }
  else{
    af_cols = paste0('X', c(0.1,0.2,'1.0'))
  }
  if(!is.null(af_override)){
    af_cols = af_override
  }
  lci_cols = paste0(af_cols, "_LCI")
  lci_vals = c()
  uci_cols = paste0(af_cols, "_UCI")
  uci_vals = c()
  for(i in 1:length(af_cols)){
    lci_vals = c(lci_vals, autoval_raw_table[,lci_cols[i]])
    uci_vals = c(uci_vals, autoval_raw_table[,uci_cols[i]])
  }
  
  autoval_melt = melt(autoval_raw_table, id.vars=c("depth_bin"), measure.vars=af_cols)
  autoval_melt$LCI = lci_vals
  autoval_melt$UCI = uci_vals
  return(autoval_melt)
}


plot_sensitivity = function(raw_table_v2, raw_table_v6){
  v2_snv = read.delim(raw_table_v2, skip=1, nrows = 17)
  v2_snv_plotting = autovalRawToMelt(v2_snv, "SNV")
  v2_indel = read.delim(raw_table_v2, skip=21, nrows = 17)
  v2_indel_plotting = autovalRawToMelt(v2_indel, "INDEL")
  v6_snv = read.delim(raw_table_v6, skip=1, nrows = 17)
  v6_snv_plotting = autovalRawToMelt(v6_snv, "SNV")
  v6_indel = read.delim(raw_table_v6, skip=21, nrows = 17)
  v6_indel_plotting = autovalRawToMelt(v6_indel, "INDEL")
  
  v2_autoval_snv_plot = ggplot(v2_snv_plotting[v2_snv_plotting$variable %in% c("X0.05","X0.1","X0.2","X1.0"),],
                               aes(x=depth_bin, y=value, color=variable))+geom_line(linewidth=1.2)+
    geom_point(alpha=0.5)+geom_errorbar(aes(ymax=UCI, ymin=LCI), width=0, linewidth=1)+
    scale_color_manual(name="VAF Bin", labels=c("5%", "10%", "20%","100%"), values=tableau20)+
    labs(x="Read Depth Bin", y="Sensitivity", title="WES V2 HapMap SNV Sensitivity")+theme_bw()+
    theme(axis.text = element_text(size=16), axis.title = element_text(size=20), 
          plot.title = element_text(size=20), legend.text = element_text(size=16), 
          legend.title = element_text(size=20))+ylim(0,1)
  
  v2_autoval_indel_plot = ggplot(v2_indel_plotting, aes(x=depth_bin, y=value, color=variable))+
    geom_line(linewidth=1.2)+geom_point(size=2)+geom_errorbar(aes(ymax=UCI, ymin=LCI), width=0, linewidth=1)+
    scale_color_manual(name="VAF Bin", labels=c("10%", "20%", "100%"), values=tableau20[2:9])+
    labs(x="Read Depth Bin", y="Sensitivity", title="WES V2 HapMap InDel Sensitivity")+
    theme_bw()+theme(axis.text = element_text(size=16), axis.title = element_text(size=20), 
                     plot.title = element_text(size=20), legend.text = element_text(size=16), 
                     legend.title = element_text(size=20))
  
  v6_autoval_snv_plot = ggplot(v6_snv_plotting[v6_snv_plotting$variable %in% c("X0.05","X0.1","X0.2","X1.0"),], 
                               aes(x=depth_bin, y=value, color=variable))+geom_line(linewidth=1.2)+
    geom_point(alpha=0.5)+geom_errorbar(aes(ymax=UCI, ymin=LCI), width=0, linewidth=1)+
    scale_color_manual(name="VAF Bin", labels=c("5%", "10%", "20%","100%"), values=tableau20)+
    labs(x="Read Depth Bin", y="Sensitivity", title="WES V6 HapMap SNV Sensitivity")+theme_bw()+
    theme(axis.text = element_text(size=16), axis.title = element_text(size=20), 
          plot.title = element_text(size=20), legend.text = element_text(size=16), 
          legend.title = element_text(size=20))+ylim(0,1)
  
  v6_autoval_indel_plot = ggplot(v6_indel_plotting, aes(x=depth_bin, y=value, color=variable))+
    geom_line(linewidth=1.2)+geom_point(size=2)+geom_errorbar(aes(ymax=UCI, ymin=LCI), width=0, linewidth=1)+
    scale_color_manual(name="VAF Bin", labels=c("10%", "20%", "100%"), values=tableau20[2:9])+
    labs(x="Read Depth Bin", y="Sensitivity", title="WES V6 HapMap InDel Sensitivity")+
    theme_bw()+theme(axis.text = element_text(size=16), axis.title = element_text(size=20), 
                     plot.title = element_text(size=20), legend.text = element_text(size=16), 
                     legend.title = element_text(size=20))
  
  ggarrange(v6_autoval_snv_plot+
              xlim(0,800)+labs(x="")+geom_vline(xintercept=125, linetype="dotted",
                                                color="red3", linewidth=1),
            v2_autoval_snv_plot+
              labs(y="", x="")+xlim(0,800)+geom_vline(xintercept=125, linetype="dotted", 
                                     color="red3", linewidth=1),
            v6_autoval_indel_plot+
              xlim(0,800)+geom_vline(xintercept=125, linetype="dotted", 
                                     color="red3", linewidth=1), 
            v2_autoval_indel_plot+
              labs(y="")+xlim(0,800)+geom_vline(xintercept=125, linetype="dotted", 
                                                color="red3", linewidth=1), 
            nrow = 2, ncol=2, common.legend = TRUE, legend = "right", align = "hv")
}

## False Positive Rate

combine_fp_summary = function(summary_v2, summary_v6){
  summary_v2_df = read.delim(summary_v2)
  summary_v2_df$pipeline_label = "WES V2"
  summary_v6_df = read.delim(summary_v6)
  summary_v6_df$pipeline_label = "WES V6"
  summary_df = rbind(summary_v2_df, summary_v6_df)
  summary_df$pipeline_label = factor(summary_df$pipeline_label, levels=c("WES V2", "WES V6"))
  
  fp_summary = melt(summary_df, id.vars=c(1,7), measure.vars=c(4,5))
  fp_summary$var_label = "SNV"
  fp_summary$var_label[grep("indel", fp_summary$variable)] = "InDel"
  fp_summary$var_label = factor(fp_summary$var_label, levels = c("SNV", "InDel"))
  return(fp_summary)
}

ggplot(fp_summary, aes(x=pipeline_label, fill=pipeline_label, color=pipeline_label, y=value))+
  geom_boxplot(alpha=0.7, width=0.5)+
  geom_point(inherit.aes = FALSE, aes(x=pipeline_label, y=value, color=pipeline_label), 
             size=2, position = position_jitterdodge(jitter.width = 0.2), alpha=0.7)+
  labs(x="", y="FP Rate")+theme_bw()+ facet_wrap(~var_label, ncol=1)+
  scale_color_manual(name="Workflow", labels=c("WES V2", "WES V6"), 
                     values=c("#00BA38","#F8766D"))+
  scale_fill_manual(name="Workflow", labels=c("WES V2", "WES V6"), 
                    values=c("#00BA38","#F8766D"))+
  labs(x="", y="FP Calls per Million Bases", 
       title="Normal-Normal Analysis False\nPositive Rate")+
  theme(axis.text = element_text(size=16), axis.title = element_text(size=20), 
        plot.title = element_text(size=20), legend.text = element_text(size=16),
        legend.title = element_text(size=20), strip.text = element_text(size=18))


