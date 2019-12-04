'
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('DESeq2')
'

library(ggplot2)
library(gridExtra)
library(ggpubr)
bubble_plot_with_size<-function(filename){
  df <- read.csv(filename,sep="\t",header = T)
  df$Resources <- with(df,paste0(df$TF,"_",df$Condition,"_",df$Publication))
  names(df)[3:5] <- c("Total.No.Of.Clusters","ClustersTFBinding1P","Epromoter.likeness.score")
  p1 <- ggplot(df,aes(x=ClustersTFBinding1P, y=Resources,color=Epromoter.likeness.score,size=Total.No.Of.Clusters))
  p1<-p1+geom_point()+scale_colour_viridis_c()+labs(x='Clusters TF Bining on One Promoter',y='')
  
  p2 <- ggplot(df,aes(x=Total.No.Of.Clusters, y=Resources,color=Epromoter.likeness.score,size=ClustersTFBinding1P))
  p2<-p2+geom_point()+scale_colour_viridis_c()+labs(x='Total No of Clusters',y='')
  
  p <- grid.arrange(p2,p1,nrow=1,ncol=2)
  ggsave('bubble_plot.png',p,width = 20,height=8)
  print('bubble_plot.png is generated...')
}
bubble_plot_with_size('test.csv')



bubble_plot_without_size<-function(filename){
  df <- read.csv(filename,sep="\t",header = T)
  df$Resources_TF <- with(df,paste0(df$TF,"_",df$Condition,"_",df$Publication))
  df$Resources_Condtion <- with(df,paste0(df$Condition,"_",df$TF,"_",df$Publication))
  
  names(df)[3:5] <- c("Total.No.Of.Clusters","ClustersTFBinding1P","Epromoter.likeness.score")
  
  p1 <- ggplot(df,aes(x=Total.No.Of.Clusters, y=Resources_Condtion,color=Epromoter.likeness.score))
  p1<-p1+geom_point()+scale_colour_viridis_c()+labs(x='Total No of Clusters',y='')
  
  p2 <- ggplot(df,aes(x=Total.No.Of.Clusters, y=Resources_TF,color=Epromoter.likeness.score))
  p2<-p2+geom_point()+scale_colour_viridis_c()+labs(x='Total No of Clusters',y='')
  
  p <- grid.arrange(p1,p2,nrow=1,ncol=2)
  ggsave('bubble_plot.png',p,width = 20,height=8)
  print('bubble_plot.png is generated...')
}

bubble_plot_without_size('database_epromoter_like.csv')

