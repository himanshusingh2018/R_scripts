library(ggplot2)
library(ggpubr)


###########################################################################################################################################
#Bubble Plot Function
###########################################################################################################################################
bubble_plot<- function(data)
{
  ggplot(data=data)+
    geom_count(mapping=aes(x=number.of.promoters.per.cluster, y= number.of.promoter.binding.TFs.in.cluster, size=Frequency.Percentage))+
    theme(text = element_text(size=10))+ scale_y_continuous(limits = c(0,2),breaks = seq(0,2,1))+ xlim(2,6) +
    scale_size_continuous(limits = c(0,100))+
    theme_gray(base_size = 8)+
    labs(x = "No. of Promoters/Cluster",y = "No. of Promoters binding TFs in Cluster",size = "Frequency (%)")
}

#blank plot
p00 <- ggplot() + theme_void()

x11<-read.table("output/tf_analysis/bubble_chart_gene_promoter_cluster_tf_binding_frequency.table",header=T,sep="\t")

p11<-bubble_plot(x11)

p<-ggarrange(p11,
               nrow=1, ncol=1,
               labels = c('SRF_TPA_30min'),
               hjust = -1.5, #vjust = 1.5,
               font.label = list(size=9),common.legend = FALSE)
p
ggsave(p, file="Gualdrini_SRF_binding_gene_promoter_frequency.table.png", width = 4, height = 3, dpi = 150)



