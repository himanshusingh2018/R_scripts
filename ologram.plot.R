library(data.table)
library(reshape2)
library(ggplot2)
library(ggpubr)


#READ OLOGRAM OUTPUT FILES FROM THE OUTPUT FOLDER
read_ologram_output<-function(filename){
  df <- read.csv(filename,sep="\t",
                 header=T)[c('nb_intersections_expectation_shuffled','nb_intersections_variance_shuffled',
                             'nb_intersections_true','nb_intersections_pvalue',
                             'summed_bp_overlaps_expectation_shuffled','summed_bp_overlaps_variance_shuffled',
                             'summed_bp_overlaps_true','summed_bp_overlaps_pvalue')]
  
  return(df)
}

#PLOT THE INTERSECTION REGIONS
intersection_regions_bar_plot<-function(dataframe){
  data = reshape2::melt(dataframe)#column to row conversion
  intersect = data[(data$variable %like% 'intersection'),]#extract all intersection data rows
  #calculate sd of shuffled regions
  distal_nb_intersections_sd_shuffled = sqrt(intersect$value[(intersect$variable == 'nb_intersections_variance_shuffled')&(intersect$feature == 'distal')][1])
  proximal_nb_intersections_sd_shuffled = sqrt(intersect$value[(intersect$variable == 'nb_intersections_variance_shuffled')&(intersect$feature == 'proximal')][1])
  #delete extra colums
  intersect = intersect[(intersect$variable != 'nb_intersections_variance_shuffled' & intersect$variable != 'nb_intersections_pvalue'),]
  #reorder the column values in proximal and distal
  intersect$feature = factor(intersect$feature, levels=c('proximal','distal'))
  #Add sd column in the dataframe
  intersect$sd = c(distal_nb_intersections_sd_shuffled,proximal_nb_intersections_sd_shuffled,NA,NA)
  dodge <- position_dodge(width = 0.9)#Dodge overlapping objects side-to-side
  limits <- aes(ymax = value + sd,ymin = value - sd)
  #plot generation
  p <- ggplot(intersect,aes(x=feature,
                            y=value,
                            fill=variable))+
    geom_bar(position = 'dodge', stat = 'identity')+
    geom_errorbar(limits,position = dodge, width = 0.25)+
    scale_y_continuous(expand = c(0, 0))+
    ylab('Nb. of intersections')+
    #ggtitle('Total overlap length per region type')+
    scale_fill_manual(values=c('grey65','skyblue3'),name = "", labels = c("Shuffled", "True"))+
    theme_bw()+  
    theme(#plot.title = element_text(hjust=0.5,size=6.5,face='bold'),
      axis.title=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=11),
      axis.ticks.x=element_blank(),
      axis.text.x = element_text(hjust=0.5,size=10),#angle=90), 
      axis.text.y = element_text(hjust=0.5,size=10),#angle=90), 
      legend.position = 'none',
      legend.text = element_text(size=0.5),
      legend.key.size = unit(0.2,'cm'),
      legend.margin = margin(-3,-3,-3,-3),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"))
  return(p)
}

#PLOT THE OVERLAPPING REGIONS
overlapping_regions_bar_plot<-function(dataframe){
  data = reshape2::melt(dataframe)#column to row conversion
  overlap = data[(data$variable %like% 'overlaps'),]#extract all intersection data rows
  
  #calculate sd of shuffled regions
  distal_summed_bp_overlaps_sd_shuffled = sqrt(overlap$value[(overlap$variable == 'summed_bp_overlaps_variance_shuffled')&(overlap$feature == 'distal')][1])
  proximal_summed_bp_overlaps_sd_shuffled = sqrt(overlap$value[(overlap$variable == 'summed_bp_overlaps_variance_shuffled')&(overlap$feature == 'proximal')][1])
  
  #delete extra colums
  overlap = overlap[(overlap$variable != 'summed_bp_overlaps_variance_shuffled' & overlap$variable != 'summed_bp_overlaps_pvalue'),]
  #reorder the column values in proximal and distal
  overlap$feature = factor(overlap$feature, levels=c('proximal','distal'))
  #Add sd column in the dataframe
  overlap$sd = c(distal_summed_bp_overlaps_sd_shuffled,proximal_summed_bp_overlaps_sd_shuffled,NA,NA)
  print(overlap)
  dodge <- position_dodge(width = 0.9)#Dodge overlapping objects side-to-side
  limits <- aes(ymax = value + sd,ymin = value - sd)
  
  p <- ggplot(overlap,aes(x=feature,
                          y=value,
                          fill=variable))+
    geom_bar(position = 'dodge', stat = 'identity')+
    geom_errorbar(limits,position = dodge, width = 0.25)+
    scale_y_continuous(expand = c(0, 0))+
    ylab('Nb. of overlapping bps.')+
    #ggtitle('Total overlap length per region type')+
    scale_fill_manual(values=c('grey65','skyblue3'),name = "", labels = c("Shuffled", "True"))+
    theme_bw()+  
    theme(#plot.title = element_text(hjust=0.5,size=6.5,face='bold'),
      axis.title=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=11),
      axis.ticks.x=element_blank(),
      axis.text.x = element_text(hjust=0.5,size=10),#angle=90), 
      axis.text.y = element_text(hjust=0.5,size=10),#angle=90), 
      legend.position = 'none',
      legend.text = element_text(size=0.5),
      legend.key.size = unit(0.2,'cm'),
      legend.margin = margin(-3,-3,-3,-3),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"))
  return(p)
}


#READ THE DISTAL DRIVER ENRICHMENT OLOGRAM OUTPUT
distal_driver = read_ologram_output('distal_driverelements_ologram/00_ologram_stats.tsv')

#READ THE DISTAL DRIVER ENRICHMENT OLOGRAM OUTPUT
proximal_driver = read_ologram_output ('proximal_driverelements_ologram/00_ologram_stats.tsv' ) 

#COMBINE ALL DRIVER ENRICHMENT OLOGRAM OUTPUT
driver = rbind (distal_driver, proximal_driver) 
driver $ feature <- c ('distal', 'proximal') 

#READ DISTAL ENHANCER OLOGRAM OUTPUT
distal_enhancer = read_ologram_output ('distal_enhancer_ologram/00_ologram_stats.tsv') 
#READ PROXIMAL OLOGRAM OUTPUT
proximal_enhancer = read_ologram_output('proximal_enhancer_ologram/00_ologram_stats.tsv') 

#COMBINE ALL ENHANCER OLOGRAM OUTPUT
enhancer = rbind (distal_enhancer, proximal_enhancer) 
enhancer$feature <- c ('distal','proximal') 

#ENHANCER PROXIMAL AND DISTAL INTERSECTION REGION BAR PLOT
enhancer_intersection_barplot= intersection_regions_bar_plot(enhancer)
#DRIVER PROXIMAL AND DISTAL OVERLAPPING REGION BAR PLOT
driver_intersection_barplot= intersection_regions_bar_plot(driver)

#ENHANCER PROXIMAL AND DISTAL OVERLAPPING REGION BAR PLOT
enhancer_overlap_barplot= overlapping_regions_bar_plot(enhancer)
#DRIVER PROXIMAL AND DISTAL OVERLAPPING REGION BAR PLOT
driver_overlap_barplot= overlapping_regions_bar_plot(driver)


ggsave(enhancer_intersection_barplot, file="Enhancer_intersection_regions.png", width = 2.0, height = 1.9, dpi = 150)
ggsave(enhancer_overlap_barplot, file="Enhancer_overlap_regions.png", width = 2.0, height = 1.9, dpi = 150)

ggsave(driver_intersection_barplot, file="Driver_intersection_regions.png", width = 2.0, height = 1.9, dpi = 150)
ggsave(driver_overlap_barplot, file="Driver_overlap_regions.png", width = 2.0, height = 1.9, dpi = 150)



#p<-ggarrange(enhancer_intersection_barplot,enhancer_overlap_barplot,
#             driver_intersection_barplot,driver_overlap_barplot,
#             nrow=2, ncol=2,
#             labels = c('Enhancer','',
#                        'Driver'),
#             hjust = 1.5, #vjust = 1.5,
#             font.label = list(size=9),common.legend = FALSE)
#p
#ggsave(p, file="Gualdrini_SRF_binding_gene_promoter_frequency.table.png", width = 4, height = 3, dpi = 150)


