#!/usr/bin/R
library(stringr) 
library(triplex)
library(BSgenome.Hsapiens.UCSC.hg38)


triplex_search <- function(chrNo,start,end){
  t <- triplex.search(Hsapiens[[chrNo]][start:end],min_score=17,min_len=8)
  writeXStringSet(as(t, "DNAStringSet"), file="triplex.seq.fa", format="fasta")
  
  for(i in 1:(length(readLines(file("triplex.seq.fa")))/2)){
    print(paste0('triplex.',i,'.fa'))
    writeXStringSet(as(triplex.diagram(ts[1]), "DNAStringSet"), 
                    file=paste0('triplex.',i,'.fa'), 
                    format="fasta")
  }  
}
triplex_search('chr9',69035251,69079076)

#create dataframe to generate all triplex table
for(num in 1:(length(readLines(file("triplex.seq.fa")))/2)){
  #print(paste0("triplex.",num,".fa"))
  x <- readLines(paste0("triplex.",num,".fa")) 
  nCC=0; nTC=0;nCT=0;nTT=0;nC=0;nT=0

  for(i in seq(from=2, to=8, by=2)){
    #print(i)
    nCC = str_count(x[[i]], pattern = "CC")+nCC
    nTC = str_count(x[[i]], pattern = "TC")+nTC
    nCT = str_count(x[[i]], pattern = "CT")+nCT
    nTT = str_count(x[[i]], pattern = "TT")+nTT
    nC = str_count(x[[i]], pattern = "C")+nC
    nT = str_count(x[[i]], pattern = "T")+nT
  }
  
  pH=7
  dH = -4.9*(nCC)-8.9*(nTC + nCT)-7.4*(nTT)
  dG= -3.00*(nC)-0.65*(nT)+1.65*(nCC)+6.0+(nC)*(pH-5.0)*(1.26-0.08*(nCC))
  Tm=(310*dH)/dH-dG-310*0.00198*1.386
  df1 = rbind(df,data.frame(num,dH,dG,Tm))
}
print(df1)

