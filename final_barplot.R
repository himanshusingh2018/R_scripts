library(ggplot2)

##############not_refseq_promoter/ologram_IRF1STAT1STAT2_IRF9STAT1STAT2/00_ologram_stats.tsv###################################
barplot_from_ologram_output_file<-function(filename){
  df1 <- read.csv(filename,sep="\t",header=T)[c('summed_bp_overlaps_expectation_shuffled','summed_bp_overlaps_variance_shuffled','summed_bp_overlaps_true','summed_bp_overlaps_pvalue')]
  label <- c('summed_bp_overlaps_expectation_shuffled','summed_bp_overlaps_true')
  y <- c(df1$summed_bp_overlaps_expectation_shuffled,df$summed_bp_overlaps_true)
  cell <- c(1,2)
  df <- data.frame(stringsAsFactors = FALSE,
                   value = c(df1$summed_bp_overlaps_expectation_shuffled,df1$summed_bp_overlaps_true),
                   label = c('Shuffled','True'))
  var = df1$summed_bp_overlaps_expectation_shuffled
  
  p<-ggplot(df, aes(x=label,y=value,fill=label))+
    geom_bar(stat='identity')+
    geom_errorbar(data=subset(df,label == 'Shuffled'),
                  aes(ymin=df$value[1]-sqrt(var),ymax=df$value[1]+sqrt(var)),width=0.5)+
    
    ylab('Nb. of Overlapping base pairs')+
    ggtitle('Total overlap length per region type')+
    scale_fill_manual(values=c('grey65','skyblue3'),name = "", labels = c("Shuffled", "True"))+
    scale_y_continuous(expand = c(0,0)) + #PLOT WITHOUT P VALUE
    theme_bw()+  
    theme(plot.title = element_text(hjust=0.5,size=6.5,face='bold'),
          #axis.title.x=element_text(size=9), axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=9),
          axis.ticks.x=element_blank(),
          legend.position = 'none',
          legend.text = element_text(size=6.5),
          legend.key.size = unit(0.3,'cm'),
          legend.margin = margin(-3,-3,-3,-3),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
    return(p)
}

p<-barplot_from_ologram_output_file('not_refseq_promoter/ologram_IRF1STAT1STAT2_IRF9STAT1STAT2/00_ologram_stats.tsv')+
  geom_text(aes(x=2.55,y=60000,label='p = 4.6e-86'),angle=-90,color='red',size = 2.0)
ggsave('not_refseq_promoter_IRF1STAT1STAT2_IRF9STAT1STAT2.png',p,width=2.2,height=2.8,dpi=300)
#############################################################################################################################################################################

#not_refseq_promoter/ologram_IRF9STAT1STAT2/00_ologram_stats.tsv

barplot_from_ologram_output_file<-function(filename){
  df1 <- read.csv(filename,sep="\t",header=T)[c('summed_bp_overlaps_expectation_shuffled','summed_bp_overlaps_variance_shuffled','summed_bp_overlaps_true','summed_bp_overlaps_pvalue')]
  label <- c('summed_bp_overlaps_expectation_shuffled','summed_bp_overlaps_true')
  y <- c(df1$summed_bp_overlaps_expectation_shuffled,df$summed_bp_overlaps_true)
  cell <- c(1,2)
  df <- data.frame(stringsAsFactors = FALSE,
                   value = c(df1$summed_bp_overlaps_expectation_shuffled,df1$summed_bp_overlaps_true),
                   label = c('Shuffled','True'))
  var = df1$summed_bp_overlaps_variance_shuffled
  
  p<-ggplot(df, aes(x=label,y=value,fill=label))+
    geom_bar(stat='identity')+
    geom_errorbar(data=subset(df,label == 'Shuffled'),
                  aes(ymin=df$value[1]-sqrt(var),ymax=df$value[1]+sqrt(var)),width=0.5)+
    
    ylab('Nb. of Overlapping base pairs')+
    ggtitle('Total overlap length per region type')+
    scale_fill_manual(values=c('grey65','skyblue3'),name = "", labels = c("Shuffled", "True"))+
    scale_y_continuous(expand = c(0,0)) + #PLOT WITHOUT P VALUE
    theme_bw()+  
    theme(plot.title = element_text(hjust=0.5,size=6.5,face='bold'),
          #axis.title.x=element_text(size=9), axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=9),
          axis.ticks.x=element_blank(),
          legend.position = 'none',
          legend.text = element_text(size=6.5),
          legend.key.size = unit(0.3,'cm'),
          legend.margin = margin(-3,-3,-3,-3),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
  return(p)
}

p1<-barplot_from_ologram_output_file('not_refseq_promoter/ologram_IRF9STAT1STAT2/00_ologram_stats.tsv')+
  geom_text(aes(x=2.55,y=45000,label='p = 1.8e-84'),angle=-90,color='red',size = 2.0)

ggsave('not_refseq_promoter_IRF9STAT1STAT2.png',p1,width=2.2,height=2.8,dpi=300)

#############################################################################################################################################################################

#refseq_promoter/ologram_IRF9STAT1STAT2/00_ologram_stats.tsv

barplot_from_ologram_output_file<-function(filename){
  df1 <- read.csv(filename,sep="\t",header=T)[c('summed_bp_overlaps_expectation_shuffled','summed_bp_overlaps_variance_shuffled','summed_bp_overlaps_true','summed_bp_overlaps_pvalue')]
  label <- c('summed_bp_overlaps_expectation_shuffled','summed_bp_overlaps_true')
  y <- c(df1$summed_bp_overlaps_expectation_shuffled,df$summed_bp_overlaps_true)
  cell <- c(1,2)
  df <- data.frame(stringsAsFactors = FALSE,
                   value = c(df1$summed_bp_overlaps_expectation_shuffled,df1$summed_bp_overlaps_true),
                   label = c('Shuffled','True'))
  var = df1$summed_bp_overlaps_variance_shuffled
  
  p<-ggplot(df, aes(x=label,y=value,fill=label))+
    geom_bar(stat='identity')+
    geom_errorbar(data=subset(df,label == 'Shuffled'),
                  aes(ymin=df$value[1]-sqrt(var),ymax=df$value[1]+sqrt(var)),width=0.5)+
    
    ylab('Nb. of Overlapping base pairs')+
    ggtitle('Total overlap length per region type')+
    scale_fill_manual(values=c('grey65','skyblue3'),name = "", labels = c("Shuffled", "True"))+
    scale_y_continuous(expand = c(0,0)) + #PLOT WITHOUT P VALUE
    theme_bw()+  
    theme(plot.title = element_text(hjust=0.5,size=6.5,face='bold'),
          #axis.title.x=element_text(size=9), axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=9),
          axis.ticks.x=element_blank(),
          legend.position = 'none',
          legend.text = element_text(size=6.5),
          legend.key.size = unit(0.3,'cm'),
          legend.margin = margin(-3,-3,-3,-3),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
  return(p)
}

p2<-barplot_from_ologram_output_file('refseq_promoter/ologram_IRF9STAT1STAT2/00_ologram_stats.tsv')+
  geom_text(aes(x=2.55,y=12000,label='p = 3.1e-13'),angle=-90,color='red',size = 2.0)
p2
ggsave('refseq_promoter_IRF9STAT1STAT2.png',p2,width=2.2,height=2.8,dpi=300)

#############################################################################################################################################################################

#refseq_promoter/ologram_IRF1STAT1STAT2_IRF9STAT1STAT2/00_ologram_stats.tsv

barplot_from_ologram_output_file<-function(filename){
  df1 <- read.csv(filename,sep="\t",header=T)[c('summed_bp_overlaps_expectation_shuffled','summed_bp_overlaps_variance_shuffled','summed_bp_overlaps_true','summed_bp_overlaps_pvalue')]
  label <- c('summed_bp_overlaps_expectation_shuffled','summed_bp_overlaps_true')
  y <- c(df1$summed_bp_overlaps_expectation_shuffled,df$summed_bp_overlaps_true)
  cell <- c(1,2)
  df <- data.frame(stringsAsFactors = FALSE,
                   value = c(df1$summed_bp_overlaps_expectation_shuffled,df1$summed_bp_overlaps_true),
                   label = c('Shuffled','True'))
  print(df)
  var = df1$summed_bp_overlaps_variance_shuffled
  print(var)
  p<-ggplot(df, aes(x=label,y=value,fill=label))+
    geom_bar(stat='identity')+
    geom_errorbar(data=subset(df,label == 'Shuffled'),
                  aes(ymin=df$value[1]-sqrt(var),ymax=df$value[1]+sqrt(var)),width=0.5)+
    
    ylab('Nb. of Overlapping base pairs')+
    ggtitle('Total overlap length per region type')+
    scale_fill_manual(values=c('grey65','skyblue3'),name = "", labels = c("Shuffled", "True"))+
    scale_y_continuous(expand = c(0,0)) + #PLOT WITHOUT P VALUE
    theme_bw()+  
    theme(plot.title = element_text(hjust=0.5,size=6.5,face='bold'),
          #axis.title.x=element_text(size=9), axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=9),
          axis.ticks.x=element_blank(),
          legend.position = 'none',
          legend.text = element_text(size=6.5),
          legend.key.size = unit(0.3,'cm'),
          legend.margin = margin(-3,-3,-3,-3),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
  return(p)
}

p3<-barplot_from_ologram_output_file('refseq_promoter/ologram_IRF1STAT1STAT2_IRF9STAT1STAT2/00_ologram_stats.tsv')+
  geom_text(aes(x=2.55,y=13000,label='p = 2.2e-16'),angle=-90,color='red',size = 2.0)
p3
ggsave('refseq_promoter_ologram_IRF1STAT1STAT2_IRF9STAT1STAT2.png',p3,width=2.2,height=2.8,dpi=300)
#############################################################################################################################################################
