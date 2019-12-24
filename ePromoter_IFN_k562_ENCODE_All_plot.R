#http://www.gastonsanchez.com/gotplot/plots/heatmap-ggplot.html

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggrepel)

##################### WRITING PLOT DATA ##################################################################
write_data <- function(dataframe,outputfile){
  write.table(dataframe,outputfile,sep="\t",row.names = F)
}

##################### DENSITY PLOT FUNCTIONS #############################################################
read_file_density_plot<-function(withInh_file,withoutInh_file){
  inh <- read.csv(withInh_file,sep='\t',header=T)
  inh$category <- 'With Inhibitor'
  no_inh <- read.csv(withoutInh_file,sep="\t",header=T)
  no_inh$category <- 'No Inhibitor'
  df <- rbind(inh[,11:12],no_inh[,11:12])
  return(df)
}
#simple density plot
density_plot<-function(df,x_axis_max_limit){
  write_data(df,'plot_data/density_plot_data.tsv')

  p<- ggplot()+
    geom_density(data=df,aes(tss_closest.peak_distance,colour=category))+
    xlim(0,x_axis_max_limit)+
    #theme(aspect.ratio=4/3)+
    #scale_x_continuous(breaks=df$tss_closest.peak_distance,labels = scientific_10())+
    #scale_x_log10(breaks = trans_breaks('log10',function(x) 10^x),
    #              labels = trans_format('log10', math_format(10^.x)))+
    #yScale=c('log10')+
    scale_fill_discrete(yScale=c('log10'),breaks = rev(levels(df$category)))+ 
    theme_bw()+
    theme(text = element_text(size = 8),
          plot.margin=unit(c(0.09,0.09,0.09,0.09),"cm"),
          panel.grid.major = element_line(size=0.1), panel.grid.minor = element_line(size = 0.1),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position=c(0.82,0.88),
          legend.title = element_text(size = 6), 
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.2, "cm"))
  ggsave("plot_data/density_plot.png",p, width = 2.9, height = 2.4)
  return(p)
}
density_plot_5kb = density_plot(df,5000)
#log10 transformed density plot
log10_density_plot<-function(df,x_axis_max_limit){
  write_data(df,'plot_data/log10_density_plot_data.tsv')
  p<-ggplot()+
    geom_density(data=df,aes(tss_closest.peak_distance,colour=category))+
    xlim(0,x_axis_max_limit)+
    scale_x_log10(breaks = trans_breaks('log10',function(x) 10^x),
                  labels = trans_format('log10', math_format(10^.x)))+
    theme(aspect.ratio=4/3)+
    theme(text = element_text(size = 14),legend.position = 'right')+
    scale_fill_discrete(breaks = rev(levels(df$category)))+
    theme(plot.margin=unit(c(0.09,0.09,0.09,0.09),"cm"))+
    theme_bw()+
    theme(panel.grid.major = element_line(size=0.1), 
          panel.grid.minor = element_line(size = 0.1),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position=c(0.22,0.88),
          legend.title = element_text(size = 6), 
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.2, "cm"))
  ggsave("plot_data/log10_density_plot.png",p, width = 2.9, height = 2.4)
  return(p)
}

##################### HEATMAP PLOT FUNCTIONS #############################################################
read_file_heatmap<-function(filename){
  data <- read.csv(filename,sep='\t',header=T)
  df1 = data.frame(data$TF,data$log10Pval_Distal_Muerdter2017inh,'Distal (With Inhibitor)')
  colnames(df1) <- c('TF','log10Pvalue','Category')
  df2 = data.frame(data$TF,data$log10Pval_Epromoter_Muerdter2017inh,'Epromoter (With Inhibitor)')
  colnames(df2) <- c('TF','log10Pvalue','Category')
  
  df3 = data.frame(data$TF,data$log10Pval_Distal_Muerdter2017noinh,'Distal (No Inhibitor)')
  colnames(df3) <- c('TF','log10Pvalue','Category')
  df4 = data.frame(data$TF,data$log10Pval_Epromoter_Muerdter2017noinh,'Epromoter (No Inhibitor)')
  colnames(df4) <- c('TF','log10Pvalue','Category')
  
  df = rbind(df1,df2,df3,df4)
  return(df)
}
heatmap<-function(filename){
  t = strsplit(filename,'_')
  title = paste0(gsub('.tsv','',t[[1]][6]),' motifs ',t[[1]][3],'\n(',t[[1]][4],' ',t[[1]][5],')')
  data = read_file_heatmap(filename)
  write_data(data,paste0('plot_data/',t[[1]][3],'_',t[[1]][4],'_',t[[1]][5],'_heatmap_plot_data.tsv'))
  p<-ggplot(data=data,aes(x=Category,y=TF))+
    geom_tile(aes(fill=log10Pvalue),color='white',size=0.3)+
    scale_fill_gradient(low='grey90',high = 'darkred',limits=c(0,325))+
    theme(text = element_text(size = 7),axis.text.x = element_text(angle = 40,hjust=1,vjust=1))+
    xlab("")+
    ylab("")+
    ggtitle(title)+
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          legend.position=c(1.215,0.5),
          #aspect.ratio = 4/3,
          #plot.title = element_blank())
          plot.title = element_text(hjust=0.5,size = 5, colour = "black"))+
    theme(plot.margin=unit(c(0.09,1.09,0.09,0.09),"cm"))+
    guides(fill= guide_colorbar(barheight=4.0,barwidth = 0.6))
  
  ggsave(paste0('plot_data/',t[[1]][3],'_',t[[1]][4],'_',t[[1]][5],'_heatmap.png'),p, width = 2.5, height = 2.1)
  return(p)
}
##################### VOLCANO PLOT FUNCTIONS #############################################################
read_file_volcano_plot <-function(filename){
  genes <- read.table(filename, header = TRUE)[,c('gene','DESeq2_log2FoldChange','DESeq2_pvalue','DESeq2_padj')]
  colnames(genes) <- c('Gene','log2FoldChange','pValue','pAdj')
  genes$Significant <- ifelse(((genes$pAdj < 0.05) & (genes$log2FoldChange>1)) | ((genes$pAdj < 0.05) & (genes$log2FoldChange < -1)), "pAdj < 0.001", "Not Sig")
  genes <- na.omit(genes)
  
  return(genes)
}
volcano_plot<-function(filename){
  genes <- read_file_volcano_plot(filename)
  genes = genes[order(rev(genes$Significant)),]
  write_data(genes,'plot_data/volcano_plot_data.tsv')
  p<-ggplot(genes,aes(x=log2FoldChange,y=-log10(pAdj)))+
    geom_point(aes(color = Significant),size=0.001)+
    guides(fill = guide_legend(reverse = F))+
    scale_color_manual(values = c('grey','darkmagenta'),name='Significance',labels = c('Not Sig','log2FC >  1 (pAdj<0.001)\nlog2FC < -1 (pAdj<0.001)'))+
    #scale_y_log10(breaks = trans_breaks('log10',function(x) 10^x),
    #              labels = trans_format('log10', math_format(10^.x)))+
    theme_bw()+
    theme(axis.text=element_text(size=8),
          axis.ticks = element_blank(),
          #plot.title = element_text(hjust=0.5,size = 50, colour = "black"),
          legend.title = element_text(size=6),
          plot.margin=unit(c(0.09,3.09,0.09,0.09),"cm"),
          panel.grid.major = element_line(size=0.1), panel.grid.minor = element_line(size = 0.1),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.text=element_text(margin = margin(r=10,unit='pt'),size=5),
          legend.position=c(1.45,0.6),
          legend.key.size = unit(6, "pt"),
          legend.spacing.x = unit(1, 'pt'))
          #aspect.ratio = 4/3,
          #plot.title = element_blank())
  ggsave('plot_data/2volcano_plot.png',p, width = 3.2, height = 2.1)
  return(p)
}
volcano = volcano_plot('data/RNAseq_K562_DESeq2_EdgeR_pval_and_norm_count_log2.txt')
##################### SCATTER PLOT IFN VS NS FUNCTIONS ###################################################
read_file_starseq_ifn_vs_ns<-function(filename){
  data = read.csv(filename,sep="\t",header=T)[,c('gene','NS1_fold_change','NS2_fold_change','IFN1_fold_change','IFN2_fold_change','IFN_vs_NS.plot_condition')]
  data = na.omit(data)
  data$NS.log2FC <- log((data$NS1_fold_change+data$NS2_fold_change)/2,2)
  data$IFN.log2FC <- log((data$IFN1_fold_change+data$IFN2_fold_change)/2,2)
  data$IFN_vs_NS.log2FC <- data$IFN.log2FC - data$NS.log2FC
  data$Epromoters <-   ifelse(data$IFN_vs_NS.plot_condition == 'Never Epromoter','No Epromoter',
                              ifelse(data$IFN_vs_NS.plot_condition != 'Never Epromoter' & data$IFN_vs_NS.log2FC>1,'Induced Epromoter',
                                     ifelse(data$IFN_vs_NS.plot_condition != 'Never Epromoter' & data$IFN_vs_NS.log2FC< -1,'Repressed Epromoter','Non-induced Epromoter')))
  return(data)
}
scatter_plot_starseq_ifn_vs_ns<-function(filename){
  data = read_file_starseq_ifn_vs_ns(filename)
  data = na.omit(data)
  write_data(data,'plot_data/scatter_plot_starseq_ifn_vs_ns_plot_data.tsv')
  
  p<-ggplot(data,aes(NS.log2FC,IFN.log2FC,col=Epromoters))+
    geom_point(pch=20,size=0.001)+
    scale_colour_manual(values = c('darkred','gray60','orange','cyan3'),
                        breaks=c('Induced Epromoter','Non-induced Epromoter','Repressed Epromoter','No Epromoter'),
                        labels = function(x) format(x, scientific = TRUE))+
    theme(axis.text=element_text(size=6),
          axis.ticks = element_blank(),
          #plot.title = element_text(hjust=0.5,size = 50, colour = "black"),
          legend.title = element_text(size=6),
          plot.margin=unit(c(0.09,2.59,0.09,0.09),"cm"),
          panel.grid.major = element_line(size=0.1), panel.grid.minor = element_line(size = 0.1),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.text=element_text(margin = margin(r=10,unit='pt'),size=5),
          legend.position=c(1.35,0.6),
          legend.key.size = unit(8, "pt"),
          legend.spacing.x = unit(1, 'pt'))
  ggsave('plot_data/scatter_starseq_ifn_vs_ns.png',p, width = 3.2, height = 2.1)
  return(p)
}

##################### SCATTER PLOT STAR SEQ vs RNA SEQ FUNCTIONS #########################################
read_file_starseq_rnaseq<-function(filename){
  data = read.csv(filename,sep="\t",header=TRUE)[c('gene','DESeq2_log2FoldChange','DESeq2_padj','IFN_vs_NS_fold_change.log2FC','IFN_vs_NS.plot_condition')]
  data = na.omit(data)
  
  data$Condition <- ifelse((data$IFN_vs_NS.plot_condition != 'Never Epromoter'),'Epromoter',
                           ifelse(((data$DESeq2_log2FoldChange>1) | (data$DESeq2_log2FoldChange< -1) & (data$DESeq2_padj < 0.001)),'Diff Gene (No Epromoter)','No Change')
  )
  data <- data[data$Condition != 'No Change',]
  return(data)
}
scatter_plot_starseq_rnaseq<-function(filename){
  data = read_file_starseq_rnaseq(filename)
  write_data(data,'scatter_plot_starseq_rnaseq_plot_data.tsv')
  p<-ggplot(data,aes(IFN_vs_NS_fold_change.log2FC,DESeq2_log2FoldChange,col=Condition))+
    geom_point(pch=20,size=0.001)+
    scale_color_manual(values = c("red","mediumblue"))+
    theme(axis.text=element_text(size=6),
          axis.ticks = element_blank(),
          #plot.title = element_text(hjust=0.5,size = 50, colour = "black"),
          legend.title = element_text(size=6),
          plot.margin=unit(c(0.09,2.59,0.09,0.09),"cm"),
          panel.grid.major = element_line(size=0.1), panel.grid.minor = element_line(size = 0.1),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.text=element_text(margin = margin(r=10,unit='pt'),size=5),
          legend.position=c(1.35,0.6),
          legend.key.size = unit(8, "pt"),
          legend.spacing.x = unit(1, 'pt'))
  ggsave('plot_data/scatter_starseq_rnaseq.png',p, width = 3.2, height = 2.1)
  return(p)
}
scatter_plot_starseq_rnaseq('data/plot_data_K562_RNAseq_Capstarseq.hg19.tsv')
##########################################################################################################


##################### DENSITY PLOT DATA ##################################################################
df = read_file_density_plot('data/Muerdter2017_annotated_peakfile_with_inhibitor.hg19.tsv',
                            'data/Muerdter2017_annotated_peakfile_without_inhibitor.hg19.tsv')
density_plot_5kb = density_plot(df,5000)
density_plot_log10_all = log10_density_plot(df,max(df$tss_closest.peak_distance))

##################### HEATMAP PLOT DATA ##################################################################
heatmap_epromoter_noinh <- heatmap('data/homerResults/heatmap_data/knownResults.Muerdter2017_Epromoter_No_Inhibitor_top10.tsv')
heatmap_epromoter_inh <- heatmap('data/homerResults/heatmap_data/knownResults.Muerdter2017_Epromoter_With_Inhibitor_top10.tsv')
heatmap_distal_noinh <- heatmap('data/homerResults/heatmap_data/knownResults.Muerdter2017_Distal_No_Inhibitor_top10.tsv')
heatmap_distal_inh <- heatmap('data/homerResults/heatmap_data/knownResults.Muerdter2017_Distal_With_Inhibitor_top10.tsv')
heatmap_distal_inh
##################### VOLCANO PLOT DATA ##################################################################
volcano = volcano_plot('data/RNAseq_K562_DESeq2_EdgeR_pval_and_norm_count_log2.txt')
##################### STARSEQ IFN vs NS PLOT DATA ########################################################
star_seq_ifn_vs_ns <- scatter_plot_starseq_ifn_vs_ns('data/plot_data_K562_RNAseq_Capstarseq.hg19.tsv')

###################### GROUP ALL PLOTS ###################################################################
p1<-grid.arrange(density_plot_5kb,
                 density_plot_log10_all,
                 arrangeGrob(heatmap_epromoter_noinh,heatmap_epromoter_inh,
                             heatmap_distal_noinh,heatmap_distal_inh,
                             nrow=2,ncol = 2),
                 nrow=1)

p2<-grid.arrange(volcano,
                 star_seq_ifn_vs_ns,
                 volcano,
                 nrow=1,ncol=3)

p<-grid.arrange(p1,p2,
                nrow=2,ncol=1)

ggsave("1testing_density_plot_with_various_x-axis_windows.png",p, width = 15, height = 9)
##########################################################################################################

