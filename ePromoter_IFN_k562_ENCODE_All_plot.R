library(ggplot2)
library(ggpubr)
library(gridExtra)

#################################################################################################
############# MUERDTER DATA PLOT
#READ FILE
read_file_density_plot<-function(withInh_file,withoutInh_file){
  inh <- read.csv(withInh_file,sep='\t',header=T)
  inh$category <- 'With Inhibitor'
  no_inh <- read.csv(withoutInh_file,sep="\t",header=T)
  no_inh$category <- 'No Inhibitor'
  df <- rbind(inh[,11:12],no_inh[,11:12])
  return(df)
}

#DENSITY PLOT FUNCTION
#simple density plot
density_plot<-function(df,x_axis_max_limit){
  ggplot()+
    geom_density(data=df,aes(tss_closest.peak_distance,colour=category))+
    xlim(0,x_axis_max_limit)+
    theme(legend.position = 'top')+
    scale_fill_discrete(breaks = rev(levels(df$category)))
}

#log10 transformed density plot
log10_density_plot<-function(df,x_axis_max_limit){
  ggplot()+
    geom_density(data=df,aes(tss_closest.peak_distance,colour=category))+
    xlim(0,x_axis_max_limit)+
    scale_x_log10()+
    theme(legend.position = 'top')+
    scale_fill_discrete(breaks = rev(levels(df$category)))
}


df = read_file_density_plot('data/Muerdter2017_annotated_peakfile_with_inhibitor.hg19.tsv',
          'data/Muerdter2017_annotated_peakfile_without_inhibitor.hg19.tsv')

density_plot_5kb = density_plot(df,5000)
density_plot_log10_all = log10_density_plot(df,max(df$tss_closest.peak_distance))

ggarrange(density_plot_5kb,density_plot_all,density_plot_log10_all,
          ncol = 3, nrow = 3)

p<-grid.arrange(density_plot_5kb,
             arrangeGrob(density_plot_5kb,density_plot_5kb, nrow=2),
             arrangeGrob(density_plot_5kb, nrow = 2),
             nrow=1)

ggsave("density_plot_with_various_x-axis_windows.png",p, width = 18, height = 18)

'
density_plot_5kb = density_plot(df,5000)
density_plot_10kb = density_plot(df,10000)
density_plot_15kb = density_plot(df,15000)
density_plot_20kb = density_plot(df,20000)
density_plot_25kb = density_plot(df,25000)
density_plot_50kb = density_plot(df,50000)
density_plot_100kb = density_plot(df,100000)
density_plot_all = density_plot(df,max(df$tss_closest.peak_distance))
density_plot_log10_all = log10_density_plot(df,max(df$tss_closest.peak_distance))


ggarrange(density_plot_5kb, density_plot_10kb, density_plot_15kb,
          density_plot_20kb, density_plot_25kb, density_plot_50kb,
          density_plot_100kb,density_plot_all,density_plot_log10_all,
          ncol = 3, nrow = 3)

ggsave("density_plot_with_various_x-axis_windows.png", width = 18, height = 18)
'
#################################################################################################

#################################################################################################
############# MUERDTER DATA PLOT

