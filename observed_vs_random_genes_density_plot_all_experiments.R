library(ggplot2)
library(ggpubr)
library(cowplot)
setwd('/home/singh/Desktop/Analysis/stress_data_plot/least_distant_observed_vs_randomlySelected/25Aug2019')
#Read least distant genes files
theme_set(theme_cowplot())
read_tss_least_distant_genes <- function(observed_filename){
  observed = read.table(observed_filename,header=T,sep="\t")
  #randomly_selected = read.table(randomly_selected_filename,header=T,sep="\t")
  return(observed)
}

#density plot
density_plot <- function(data1,data2){
  #Observed data
  data1<-data.frame(data1$distance)
  data1['experiment'] <- 'observed'
  colnames(data1) <- c('distance','Experiment')
  #Randomly selected
    data2 <- data.frame(data2$distance)
  data2['experiment'] <- 'randomly_selected'
  colnames(data2) <- c('distance','Experiment')
  df <- rbind(data1,data2)
  #pval = paste0(c('pvalue =',2.4))
  ggplot(df)+
    geom_density(aes(x=df$distance, colour=Experiment),show.legend=FALSE)+
    stat_density(aes(x=df$distance, y = ..scaled.., colour=Experiment),
                 geom="line",position="identity")+
    xlim(0,10000000)+
    scale_y_continuous()+
    #scale_y_discrete(limits = rev(levels(df$distance)))+
    scale_colour_manual( values = c("red","blue"))+
    labs(y="Density", x="Least distance between TSS of two genes")+
    #geom_text(aes(x=7000000,y=0.7,family="myFont3",label=paste('ks.test pvalue = ',round(ks.test(data1$distance,data2$distance)$p.value,4))))
    geom_text(aes(x=6000000,y=0.8,family="myFont3",label=paste('ks.test pvalue = ',formatC(ks.test(data1$distance,data2$distance,alternative='two.sided',
                                                                                                   exact = NULL)$p.value,format="e"))))
}


#density plot
density_plot <- function(data1,data2){
  #Observed data
  data1<-data.frame(data1$distance)
  data1['experiment'] <- 'observed'
  colnames(data1) <- c('distance','Experiment')
  #Randomly selected
  data2 <- data.frame(data2$distance)
  data2['experiment'] <- 'randomly_selected'
  colnames(data2) <- c('distance','Experiment')
  df <- rbind(data1,data2)
  #pval = paste0(c('pvalue =',2.4))
  ggplot(df)+
    geom_density(aes(x=df$distance, colour=Experiment),show.legend=FALSE)+
    stat_density(aes(x=df$distance, y = ..scaled.., colour=Experiment),
                 geom="line",position="identity")+
    xlim(0,10000000)+
    scale_y_continuous()+
    #scale_y_discrete(limits = rev(levels(df$distance)))+
    scale_colour_manual( values = c("red","blue"))+
    labs(y="Density", x="Least distance between TSS of two genes")+
    geom_text(aes(x=6000000,y=0.8,family="myFont3",
                  label=ifelse(ks.test(data1$distance,data2$distance,alternative = 'two.sided',exact = NULL)$p.value == 0,'ks.test pvalue = 1.0000e-322',
                               paste('ks.test pValue = ', formatC(ks.test(data1$distance,data2$distance,alternative = 'two.sided',exact=NULL)$p.value,format = "e")))))
    #geom_text(aes(x=7000000,y=0.7,family="myFont3",label=paste('ks.test pvalue = ',round(ks.test(data1$distance,data2$distance)$p.value,4))))
    #geom_text(aes(x=6000000,y=0.8,family="myFont3",
    #              label=ifelse(ks.test(data1$distance,data$distance,alternative = 'two.sided',exact = NULL)$p.value == 0,'ks.test pvalue = 1.0000e-322',
    #                           paste('ks.test pValue = ', formatC(ks.test(data1$distance,data2$distance,alternative = 'two.sided',exact=NULL)$p.value,format = "e")))))
}




#############################################################################################################################################################################
#############################################################################################################################################################################
#############################################################################################################################################################################

#Human Data
#random_selection_of_genes_in_cluster.genes_random_selection(hg38_genes, hg38_refgene_data,778, 98, 1)#Viervaahra HS human
Viervaara_HS<-density_plot(read_tss_least_distant_genes('observed/Viervaara_HS_least_distant_induced_genes.distance'),
                           read_tss_least_distant_genes('randomly_selected/Viervaahra_HS_hg_randomly_selected_1_times.798_induced_genes.least_distance'))



#random_selection_of_genes_in_cluster.genes_random_selection(hg38_genes, hg38_refgene_data,1155,138, 1)#Porter DSB human
Porter_DSB<-density_plot(read_tss_least_distant_genes('observed/Porter_DSB_least_distant_induced_genes.distance'),
                         read_tss_least_distant_genes('randomly_selected/Porter_DSB_hg_randomly_selected_1_times.1155_induced_genes.least_distance'))

#random_selection_of_genes_in_cluster.genes_random_selection(hg38_genes, hg38_refgene_data,1155,285, 1)#Lyu HS human
Lyu_HS<-density_plot(read_tss_least_distant_genes('observed/Lyu_HS_least_distant_induced_genes.distance'),
                     read_tss_least_distant_genes('randomly_selected/Lyu_HS_hg_randomly_selected_1_times.1155_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(hg38_genes, hg38_refgene_data,613,77, 1)#Jin TNFa human
Jin_TNFa<-density_plot(read_tss_least_distant_genes('observed/Jin_TNFa_least_distant_induced_genes.distance'),
                       read_tss_least_distant_genes('randomly_selected/Jin_TNFa_hg_randomly_selected_1_times.613_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(hg38_genes, hg38_refgene_data,674, 115,1)#Hogan TNFa Human
Hogan_TNFa<-density_plot(read_tss_least_distant_genes('observed/Hogan_TNFa_least_distant_induced_genes.distance'),
                         read_tss_least_distant_genes('randomly_selected/Hogan_TNFa_hg_randomly_selected_1_times.674_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(hg38_genes, hg38_refgene_data,820,167, 1)#Hogan IL1b human
Hogan_IL1B<-density_plot(read_tss_least_distant_genes('observed/Hogan_IL1B_least_distant_induced_genes.distance'),
                         read_tss_least_distant_genes('randomly_selected/Hogan_IL1b_hg_randomly_selected_1_times.820_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(hg38_genes, hg38_refgene_data,409, 78,1)#Park TNF + IFN
Park_IFN_TNF<-density_plot(read_tss_least_distant_genes('observed/Park_IFN_TNF_least_distant_induced_genes.distance'),
                           read_tss_least_distant_genes('randomly_selected/Park_TNF_IFN_hg_randomly_selected_1_times.409_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(hg38_genes, hg38_refgene_data,997,230, 1)#Park TNF
Park_TNF<-density_plot(read_tss_least_distant_genes('observed/Park_TNF_least_distant_induced_genes.distance'),
                       read_tss_least_distant_genes('randomly_selected/Park_TNF_hg_randomly_selected_1_times.997_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(hg38_genes, hg38_refgene_data,591,98, 1)#Franco E2 human
Franco_E2<-density_plot(read_tss_least_distant_genes('observed/Franco_E2_least_distant_induced_genes.distance'),
                        read_tss_least_distant_genes('randomly_selected/Franco_E2_hg_randomly_selected_1_times.591_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(hg38_genes, hg38_refgene_data,397, 55,1)#Franco TNFa human
Franco_TNFa<-density_plot(read_tss_least_distant_genes('observed/Franco_TNFa_least_distant_induced_genes.distance'),
                          read_tss_least_distant_genes('randomly_selected/Franco_TNF_hg_randomly_selected_1_times.397_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(hg38_genes, hg38_refgene_data,777, 147,1)#Franco E2+TNFa human
Franco_E2_TNFa<-density_plot(read_tss_least_distant_genes('observed/Franco_E2_TNFa_least_distant_induced_genes.distance'),
                             read_tss_least_distant_genes('randomly_selected/Franco_E2_TNF_hg_randomly_selected_1_times.777_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(hg38_genes, hg38_refgene_data,482, 107,1)#k562 ENCODE human
k562_IRF1<-density_plot(read_tss_least_distant_genes('observed/k562_IRF1_least_distant_induced_genes.distance'),
                        read_tss_least_distant_genes('randomly_selected/k562_IRF1_hg_randomly_selected_1_times.482_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(hg38_genes, hg38_refgene_data,160, 19,1)#Jubb Dex hMDM
Jubb_hMDM<-density_plot(read_tss_least_distant_genes('observed/Jubb_GR_hMDM_least_distant_induced_genes.distance'),
                        read_tss_least_distant_genes('randomly_selected/Jubb_DEX_hMDM_hg_randomly_selected_1_times.160_induced_genes.least_distance'))
#Mouse Data
#random_selection_of_genes_in_cluster.genes_random_selection(mm10_genes, mm10_refgene_data,252,60, 1)#Mancino LPS Mouse
Mancino_LPS<-density_plot(read_tss_least_distant_genes('observed/Mancino_LPS_least_distant_induced_genes.distance'),
                          read_tss_least_distant_genes('randomly_selected/Mancino_LPS_mm_randomly_selected_1_times.252_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(mm10_genes, mm10_refgene_data,896, 103,1)#Picolo IFNg Mouse
Piccolo_IFNg<-density_plot(read_tss_least_distant_genes('observed/Piccolo_IFNg_least_distant_induced_genes.distance'),
                           read_tss_least_distant_genes('randomly_selected/Picolo_IFN_mm_randomly_selected_1_times.896_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(mm10_genes, mm10_refgene_data,370,34, 1)#Picolo IL4 Mouse
Piccolo_IL4<-density_plot(read_tss_least_distant_genes('observed/Piccolo_IL4_least_distant_induced_genes.distance'),
                          read_tss_least_distant_genes('randomly_selected/Picolo_IL4_mm_randomly_selected_1_times.370_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(mm10_genes, mm10_refgene_data,280,33, 1)#Gualdrini SRF mouse
Gualdrini_SRF<-density_plot(read_tss_least_distant_genes('observed/Gualdrini_SRF_least_distant_induced_genes.distance'),
                            read_tss_least_distant_genes('randomly_selected/Gualdrini_SRF_mm_randomly_selected_1_times.280_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(mm10_genes, mm10_refgene_data,1439,354, 1)#Ensault Serum  mouse
Esnault_Serum<-density_plot(read_tss_least_distant_genes('observed/Esnault_Serum_least_distant_induced_genes.distance'),
                            read_tss_least_distant_genes('randomly_selected/Ensault_Serum_mm_randomly_selected_1_times.1439_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(mm10_genes, mm10_refgene_data,438, 54,1)#Cardamone Mito Infection 
Cardamone_mito_infect<-density_plot(read_tss_least_distant_genes('observed/Cardamone_mito_infect_least_distant_induced_genes.distance'),
                                    read_tss_least_distant_genes('randomly_selected/Cardamone_MitoInfect_mm_randomly_selected_1_times.438_induced_genes.least_distance'))
#random_selection_of_genes_in_cluster.genes_random_selection(mm10_genes, mm10_refgene_data,224, 19,1)#Jubb Dex mBMDM
Jubb_mBMDM<-density_plot(read_tss_least_distant_genes('observed/Jubb_GR_mBMDM_least_distant_induced_genes.distance'),
                         read_tss_least_distant_genes('randomly_selected/Jubb_Dex_mBMDM_mm_randomly_selected_1_times.224_induced_genes.least_distance'))

#random_selection_of_genes_in_cluster.genes_random_selection('Kusnadi_TNF_SREBP_hg',hg38_genes, hg38_refgene_data,1115, 297,1)#Kusnadi TNF SREBP human
Kusnadi_TNF_SREBP<-density_plot(read_tss_least_distant_genes('observed/Kusnadi_TNF_SREBP_least_distant_induced_genes.distance'),
                         read_tss_least_distant_genes('randomly_selected/Kusnadi_TNF_SREBP_hg_randomly_selected_1_times.1115_induced_genes.least_distance'))
Jubb_mBMDM
#blank plot
g0 <- ggplot() + theme_void()

p<-ggarrange(Viervaara_HS,Porter_DSB,Lyu_HS,
             Jin_TNFa,Hogan_TNFa,Hogan_IL1B,
             Park_IFN_TNF,Park_TNF,Franco_E2,
             Franco_TNFa,Franco_E2_TNFa,k562_IRF1,
             Jubb_hMDM,Mancino_LPS,Piccolo_IFNg,
             Piccolo_IL4,Gualdrini_SRF,Esnault_Serum,
             Cardamone_mito_infect,Jubb_mBMDM,Kusnadi_TNF_SREBP,
             nrow=7,ncol=3,
             labels = c('Viervaara_HS','Porter_DSB','Lyu_HS',
                        'Jin_TNFa','Hogan_TNFa','Hogan_IL1B',
                        'Park_IFN_TNF','Park_TNF','Franco_E2',
                        'Franco_TNFa','Franco_E2_TNFa','k562_IRF1',
                        'Jubb_hMDM','Mancino_LPS','Piccolo_IFNg',
                        'Piccolo_IL4','Gualdrini_SRF','Esnault_Serum',
                        'Cardamone_mito_infect','Jubb_mBMDM','Kusnadi_TNF_SREBP'),
             hjust = -2.5, vjust =2,
             font.label = list(size=11),common.legend = T)#,legend="top")
p
ggsave(p, file="testing_observed_vs_expected.png", width = 18, height = 21, dpi = 150)
