library(ggplot2)

ologram_output <- read.csv('STARRseq_FAIREseq_Overlap_Enrichment/00_ologram_stats.tsv',sep="\t",
                           header=T)[c('feature_type',
                                       'summed_bp_overlaps_expectation_shuffled',
                                       'summed_bp_overlaps_variance_shuffled','summed_bp_overlaps_true',
                                       'summed_bp_overlaps_pvalue')]

df_gene = ologram_output[(ologram_output$feature_type == 'gene'),]
df_promoters = ologram_output[(ologram_output$feature_type == 'Promoters'),]
df_intergenic = ologram_output[(ologram_output$feature_type == 'Intergenic'),]
df_gene
df <- data.frame(stringsAsFactors = FALSE,
                 value = c(df_gene$summed_bp_overlaps_expectation_shuffled,df_gene$summed_bp_overlaps_true,
                           df_promoters$summed_bp_overlaps_expectation_shuffled,df_promoters$summed_bp_overlaps_true,
                           df_intergenic$summed_bp_overlaps_expectation_shuffled,df_intergenic$summed_bp_overlaps_true),
                 sd = c(sqrt(df_gene$summed_bp_overlaps_variance_shuffled)[1],NA,
                        sqrt(df_promoters$summed_bp_overlaps_variance_shuffled)[1],NA,
                        sqrt(df_intergenic$summed_bp_overlaps_variance_shuffled)[1],NA),
                 label = rep(c('Shuffled','True'),3),
                 location = c('gene','gene','promoters','promoters','intergenic','intergenic'))
df
df$location <- factor(df$location,levels = c('promoters','gene','intergenic'))

dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = value + sd,
              ymin = value - sd)


p<-ggplot(data = df, aes(x = location, 
                         y = value,
                         location=factor(location, levels=c('intergenic','promoter','gene')),
                         fill = label))+
  geom_bar(stat = "identity", position = dodge)+
  geom_errorbar(limits, position = dodge, width = 0.25)+
  ylab('Nb. of Overlapping base pairs')+
  #ggtitle('Total overlap length per region type')+
  scale_fill_manual(values=c('grey65','skyblue3'),name = "", labels = c("Shuffled", "True"))+
  scale_y_continuous(expand = c(0,0)) + #PLOT WITHOUT P VALUE
  #coord_cartesian(ylim = c(0, max(df$value+1000)))+
  theme_bw()+
  theme(#plot.title = element_text(hjust=0.5,size=6.5,face='bold'),
    axis.title=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=8),
    axis.ticks.x=element_blank(),
    axis.text.x = element_text(hjust=1,size=6,angle=45), 
    axis.text.y = element_text(hjust=1,size=6),#,angle=45), 
    #legend.position = 'top',
    legend.position = c(0.2, 0.990),
    legend.text = element_text(size=6),
    legend.key.size = unit(0.3,'cm'),
    legend.margin = margin(-3,-3,-3,-3),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"))
ggsave('STARRseq_FAIREseq_Overlap_Enrichment/nb.of.overlapping.basepairs.png',p,width=2.2,height=2.2,dpi=200)

