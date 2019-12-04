library(circlize)
library(readxl)

my_data <- read_excel("Him_Circos data file_GWAS_Nov 2019.xlsx", sheet = "Complete SNP information")
#snp = data.frame(my_data$CHR,my_data$Location,my_data$SNP,my_data$pValue)
snp = data.frame(my_data$CHR,my_data$Location,my_data$Location,my_data$SNP)
snp$my_data.CHR <- gsub("^","chr",snp$my_data.CHR)
colnames(snp)<-c('chr','start','end','snp')
head(snp,2)

nrow(snp)
bed=generateRandomBed(nr=62005,n=2)
nrow(bed)

data = data.frame(name=snp$snp,start=snp$start,end=snp$start,value1=bed$value1,value2=bed$value2)

circos.clear()

circos.par("gap.degree" = rep(c(2, 4), 11))#gap of 2 in 22 chromosomes
circos.initializeWithIdeogram(chromosome.index = paste0("chr", c(1:22)))#by default it use hg19, to change species='hg38'
'
circos.genomicTrackPlotRegion(data, ylim = c(0, 1),#add rectangle. To add extra layer repeat same command
                              panel.fun = function(region, value, ...) {
                                circos.genomicPoints(region, value, ...)
                              })
'
bed <- data.frame(snp$chr,data$start,data$end,data$value1)
colnames(bed)<-c('chr','start','end','value1')
head(bed,2)
bgcol = rep(c("#EFEFEF", "#CCCCCC"), 11)

circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, pch = 21, cex = 0.01,bg = bgcol, ...)
})


circos.genomicTrackPlotRegion(data, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, ...)
})

