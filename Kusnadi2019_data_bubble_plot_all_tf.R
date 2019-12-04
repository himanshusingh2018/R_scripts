library(ggplot2)
library(ggpubr)


###########################################################################################################################################
#Bubble Plot Function
###########################################################################################################################################
bubble_plot<- function(data)
{
  ggplot(data=data)+
    geom_count(mapping=aes(x=number.of.promoters.per.cluster, y= number.of.promoter.binding.TFs.in.cluster, size=Frequency.Percentage))+
    theme(text = element_text(size=10))+ scale_y_continuous(limits = c(0,6),breaks = seq(0,8,1))+ #xlim(c(2,11),break=1) +
    scale_x_continuous(limits = c(2, 11),breaks=seq(2,11,1))+
    scale_size_continuous(limits = c(0,100))+
    theme_gray(base_size = 8)+
    labs(x = "No. of Promoters/Cluster",y = "No. of Promoters binding TFs in Cluster",size = "Frequency (%)")
}

#blank plot
#p00 <- ggplot() + theme_void()
x11<-read.table("output/tf_analysis_GSM3702335_ChIPseq_SREBF2_1_Peaks/bubble_chart_gene_promoter_cluster_tf_binding_frequency.table",header=T, sep="\t")
x12<-read.table("output/tf_analysis_GSM3702336_ChIPseq_SREBF2_2_Peaks/bubble_chart_gene_promoter_cluster_tf_binding_frequency.table",header=T, sep="\t")
x21<-read.table("output/tf_analysis_GSM3702337_ChIPseq_TNF_SREBF2_1_Peaks/bubble_chart_gene_promoter_cluster_tf_binding_frequency.table",header=T, sep="\t")
x22<-read.table("output/tf_analysis_GSM3702338_ChIPseq_TNF_SREBF2_2_Peaks/bubble_chart_gene_promoter_cluster_tf_binding_frequency.table",header=T, sep="\t")


p11<-bubble_plot(x11)
p12<-bubble_plot(x12)
p21<-bubble_plot(x21)
p22<-bubble_plot(x12)

p<-ggarrange(p11,p12,
             p21,p22,
             nrow=2, ncol=2,
             labels = c('GSM3702335_ChIPseq_SREBF2_1_Peaks','GSM3702336_ChIPseq_SREBF2_2_Peaks',
                          'GSM3702337_ChIPseq_TNF_SREBF2_1_Peaks','GSM3702338_ChIPseq_TNF_SREBF2_2_Peaks'),
               #hjust = -3.5, #vjust = 1.5,
               font.label = list(size=9),common.legend = FALSE)
p
ggsave(p, file="Kusnadi2019_all_tf_binding_gene_promoter_frequency.table.png", width = 10, height = 12, dpi = 150)



