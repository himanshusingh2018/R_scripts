#https://www.bioconductor.org/packages/release/bioc/vignettes/rGREAT/inst/doc/rGREAT.html
#GREAT Analysis
BiocManager::install("rGREAT")
library('rGREAT')
set.seed(123)


promoter=read.table('promoter.prox.bed',sep="\t")
promoter=promoter[1:3]
colnames(promoter) <- c('chr','start','end')

#background data
bg = read.table('overlap.promoter.bed',sep="\t")
bg = bg[1:3]
colnames(bg) <- c('chr','start','end')

job1 = submitGreatJob(promoter,bg,species = 'mm9')


