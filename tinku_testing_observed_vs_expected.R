library(ggplot2)
library(ggpubr)
library(cowplot)

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
    geom_text(aes(x=6000000,y=0.8,family="myFont3",
                  label=ifelse(ks.test(data1$distance,data2$distance,alternative = 'two.sided',exact = NULL)$p.value == 0,'ks.test pvalue = 1.0000e-322',
                               paste('ks.test pValue = ', formatC(ks.test(data1$distance,data2$distance,alternative = 'two.sided',exact=NULL)$p.value,format = "e")))))
  #geom_text(aes(x=7000000,y=0.7,family="myFont3",label=paste('ks.test pvalue = ',round(ks.test(data1$distance,data2$distance)$p.value,4))))
  #geom_text(aes(x=6000000,y=0.8,family="myFont3",
  #              label=ifelse(ks.test(data1$distance,data$distance,alternative = 'two.sided',exact = NULL)$p.value == 0,'ks.test pvalue = 1.0000e-322',
  #                           paste('ks.test pValue = ', formatC(ks.test(data1$distance,data2$distance,alternative = 'two.sided',exact=NULL)$p.value,format = "e")))))
}

Biddie_Dex <- density_plot(read_tss_least_distant_genes('observed/Biddie_Dex_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Biddie_Dex_randomly_selected_1_times.72_induced_genes.least_distance'))
Biddie_Dex_Tetracycline <- density_plot(read_tss_least_distant_genes('observed/Biddie_Dex_Tetracycline_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Biddie_Dex_Tetracycline_randomly_selected_1_times.269_induced_genes.least_distance'))
Brown_TNFa <- density_plot(read_tss_least_distant_genes('observed/Brown_TNFa_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Brown_TNFa_randomly_selected_1_times.260_induced_genes.least_distance'))
Camps_Hypoxia <- density_plot(read_tss_least_distant_genes('observed/Camps_Hypoxia_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Camps_Hypoxia_randomly_selected_1_times.212_induced_genes.least_distance'))
Cardamone_mito_infect <- density_plot(read_tss_least_distant_genes('observed/Cardamone_mito_infect_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Cardamone_MitoInfec_randomly_selected_1_times.438_induced_genes.least_distance'))
Demeyer_NUP214_ABL1 <- density_plot(read_tss_least_distant_genes('observed/Demeyer_NUP214-ABL1_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Demeyer_NUP214-ABL1_randomly_selected_1_times.5323_induced_genes.least_distance'))
Ebisuya_FGF_stimulated <- density_plot(read_tss_least_distant_genes('observed/Ebisuya_FGF_stimulated_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Ebisuya_FGFstimulated_randomly_selected_1_times.380_induced_genes.least_distance'))
Esnault_Serum <- density_plot(read_tss_least_distant_genes('observed/Esnault_Serum_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Ensault_Serum_randomly_selected_1_times.1439_induced_genes.least_distance'))
Ferrari_Serum_starvation <- density_plot(read_tss_least_distant_genes('observed/Ferrari_Serum-starvation_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Ferrari_SerumStarvation_randomly_selected_1_times.1274_induced_genes.least_distance'))
Franco_E2 <- density_plot(read_tss_least_distant_genes('observed/Franco_E2_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Franco_E2_randomly_selected_1_times.591_induced_genes.least_distance'))
Franco_E2_TNFa <- density_plot(read_tss_least_distant_genes('observed/Franco_E2_TNFa_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Franco_E2_TNFa_randomly_selected_1_times.777_induced_genes.least_distance'))
Franco_TNFa <- density_plot(read_tss_least_distant_genes('observed/Franco_TNFa_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Franco_TNFa_randomly_selected_1_times.397_induced_genes.least_distance'))
Gualdrini_SRF <- density_plot(read_tss_least_distant_genes('observed/Gualdrini_SRF_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Gualdrini_Starvation_randomly_selected_1_times.1597_induced_genes.least_distance'))
Hancock_Serum_starvation <- density_plot(read_tss_least_distant_genes('observed/Hancock_Serum-starvation_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Hancok_Serum_randomly_selected_1_times.1217_induced_genes.least_distance'))
Hogan_IL1B <- density_plot(read_tss_least_distant_genes('observed/Hogan_IL1B_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Hogan_IL1b_randomly_selected_1_times.820_induced_genes.least_distance'))
Hogan_TNFa <- density_plot(read_tss_least_distant_genes('observed/Hogan_TNFa_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Hogan_TNFa_randomly_selected_1_times.674_induced_genes.least_distance'))
Jin_TNFa <- density_plot(read_tss_least_distant_genes('observed/Jin_TNFa_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Jin_TNFa_randomly_selected_1_times.613_induced_genes.least_distance'))
Jubb_GR_hMDM <- density_plot(read_tss_least_distant_genes('observed/Jubb_GR_hMDM_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Jubb_GR.hMDM_randomly_selected_1_times.183_induced_genes.least_distance'))
Jubb_GR_mBMDM <- density_plot(read_tss_least_distant_genes('observed/Jubb_GR_mBMDM_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Jubb_GR.mBMDM_randomly_selected_1_times.126_induced_genes.least_distance'))
Santiago_IFN <- density_plot(read_tss_least_distant_genes('observed/Santiago_IFN_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Santiago_IFN_randomly_selected_1_times.482_induced_genes.least_distance'))
Kusnadi_TNF_SREBP <- density_plot(read_tss_least_distant_genes('observed/Kusnadi_TNF_SREBP_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Kusnadi_TNF_randomly_selected_1_times.1044_induced_genes.least_distance'))
Langlais_IFNg <- density_plot(read_tss_least_distant_genes('observed/Langlais_IFNg_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Langlais_IFNg_randomly_selected_1_times.561_induced_genes.least_distance'))
Lyu_HS <- density_plot(read_tss_least_distant_genes('observed/Lyu_HS_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Lyu_HS_randomly_selected_1_times.1155_induced_genes.least_distance'))
Mahat_HS <- density_plot(read_tss_least_distant_genes('observed/Mahat_HS_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Mahat_HS_randomly_selected_1_times.421_induced_genes.least_distance'))
Mancino_LPS <- density_plot(read_tss_least_distant_genes('observed/Mancino_LPS_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Mancino_LPS_randomly_selected_1_times.252_induced_genes.least_distance'))
Park_IFN_TNF <- density_plot(read_tss_least_distant_genes('observed/Park_IFN_TNF_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Park_TNF_IFN_randomly_selected_1_times.409_induced_genes.least_distance'))
Park_TNF <- density_plot(read_tss_least_distant_genes('observed/Park_TNF_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Park_TNF_randomly_selected_1_times.997_induced_genes.least_distance'))
Phanstiel_PMA <- density_plot(read_tss_least_distant_genes('observed/Phanstiel_PMA_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Phanstiel_PMA_randomly_selected_1_times.4738_induced_genes.least_distance'))
Piccolo_IFNg <- density_plot(read_tss_least_distant_genes('observed/Piccolo_IFNg_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Piccolo_IFNg_randomly_selected_1_times.896_induced_genes.least_distance'))
Piccolo_IL4 <- density_plot(read_tss_least_distant_genes('observed/Piccolo_IL4_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Piccolo_IL4_randomly_selected_1_times.370_induced_genes.least_distance'))
Porter_DSB <- density_plot(read_tss_least_distant_genes('observed/Porter_DSB_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Porter_HS_randomly_selected_1_times.1155_induced_genes.least_distance'))
Schmidt_TNF <- density_plot(read_tss_least_distant_genes('observed/Schmidt_TNF_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Schmidt_TNF_randomly_selected_1_times.2899_induced_genes.least_distance'))
Vierbuchen_Serum <- density_plot(read_tss_least_distant_genes('observed/Vierbuchen_Serum_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Vierbuchen_Serum_randomly_selected_1_times.1011_induced_genes.least_distance'))
Viervaara_HS <- density_plot(read_tss_least_distant_genes('observed/Viervaara_HS_least_distant_induced_genes.distance'), read_tss_least_distant_genes('randomly_selected/Viervaahra_HS_randomly_selected_1_times.778_induced_genes.least_distance'))
#blank plot
p00 <- ggplot() + theme_void()

p <- ggarrange (Santiago_IFN, Viervaara_HS, Camps_Hypoxia, Franco_E2, Franco_E2_TNFa,Franco_TNFa,
                Porter_DSB, Ferrari_Serum_starvation, Hancock_Serum_starvation, Hogan_IL1B, Hogan_TNFa, Jin_TNFa, 
                Kusnadi_TNF_SREBP, Lyu_HS, Park_IFN_TNF, Park_TNF, Phanstiel_PMA, Brown_TNFa, 
                Schmidt_TNF, Jubb_GR_hMDM, Jubb_GR_mBMDM, Biddie_Dex,Biddie_Dex_Tetracycline, Cardamone_mito_infect,
                Demeyer_NUP214_ABL1,Ebisuya_FGF_stimulated,Esnault_Serum, Gualdrini_Starvation, Vierbuchen_Serum, Mahat_HS, 
                p00,Mancino_LPS, Langlais_IFNg, Piccolo_IFNg, Piccolo_IL4,
                nrow=6,ncol=6,
                labels = c('Santiago_IFN', 'Viervaara_HS', 'Camps_Hypoxia', 'Franco_E2', 'Franco_E2_TNFa','Franco_TNFa',
                           'Porter_DSB', 'Ferrari_Serum-starvation', 'Hancock_Serum-starvation', 'Hogan_IL1B', 'Hogan_TNFa', 'Jin_TNFa', 
                           'Kusnadi_TNF_SREBP', 'Lyu_HS', 'Park_IFN_TNF', 'Park_TNF', 'Phanstiel_PMA', 'Brown_TNFa', 
                           'Schmidt_TNF', 'Jubb_GR_hMDM', 'Jubb_GR_mBMDM', 'Biddie_Dex','Biddie_Dex_Tetracycline', 'Cardamone_mito_infect',
                           'Demeyer_NUP214-ABL1','Ebisuya_FGF_stimulated','Esnault_Serum', 'Gualdrini_Starvation', 'Vierbuchen_Serum', 'Mahat_HS', 
                           '','Mancino_LPS', 'Langlais_IFNg', 'Piccolo_IFNg', 'Piccolo_IL4'),
                #hjust = -2.5, vjust =2,
                font.label = list(size=11),common.legend = T)#,legend="top")
p
ggsave(p, file="tinku_testing_observed_vs_expected.png", width = 28, height = 28, dpi = 100)
