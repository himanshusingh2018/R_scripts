dna_triplex_humangenome <- function(Input,Output){
  ### ------ load library ----------------------------------
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager",verbose = 3.12)
    BiocManager::install("org.Hs.eg.db")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  
  if(!require(triplex, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("triplex")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  
  if(!require(BSgenome.Hsapiens.UCSC.hg38,quietly = TRUE)){
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("Bsgenome.Hsapiens.UCSC.hg38")
    require(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)
  }
  
  if(!require(tidyverse, quietly = TRUE)) {
    if(!requireNamespace("tidyverse", quietly = TRUE))
      install.packages("tidyverse")
    install.packages("tidyverse")
    require(tidyverse, quietly = TRUE)
  }
  
  ### ------ THEORY OF PRIDICTION OF TRIPLEXES AND ITS MELTING TEMPERATURE -------------
  
  # A rough prediction of the melting temperature (Tm) of DNA-triplex in the genome by 
  # Roberts and Crothers empirical equations (Roberts and Crothers, 1996; Go?i et al., 2004):
  #   Tm=(310???H???)/??H????????G??????310?R?ln?(4/Tcp)
  #     where CTP is the concentration of DNA- triplex and R is gas constant.
  #  The enthalpy (??H?) is evaluated (in kcal/mol) using equation (ii):
  #   ??H??? = ???4.9 (CC) ??? 8.9 (TC +CT)???7.4 (TT)
  # The Gibbs free energy (??G?) at 37?C was computed using equation (iii):
  #   ??G???37=???3.00(C)???0.65(T)+1.65(CC)+ 6.0 + (C)(pH???5.0)[1.26???0.08(CC)]
  #     where (X) and (Y) represents cytosine and thymine
  
  ### ------- USER INPUT -----------
  user_input <- function(){
    data <- read_tsv(Input)
    #gas_constant_R <- as.double(readline("GAS Constant [Kcal K-1 mol-1]: "))
    gas_constant_R = 0.002
    #pH <- as.numeric(readline("pH: "))
    pH = 7.4
    #Tcp <- as.numeric(readline("Concertration proportion of third strand: "))
    Tcp = 1
    return(list(gas_constant_R = gas_constant_R,pH = pH,Tcp = Tcp, data = data))
  }
  ip = user_input()
  ip$data <- head(ip$data,200)
  data <- tibble(chr = as.character(),
                 start = as.numeric(),
                 end = as.numeric(),
                 description = as.character(),
                 width = as.numeric(),
                 score = as.numeric(),
                 pvalue = as.numeric(),
                 dH = as.numeric(),
                 dG37 = as.numeric(),
                 Tm = as.numeric(),
                 seq = as.character()
  )
  
  ### ----- CALCULATE TRIPLEX SEQUENCE & DNT TRIPLEX MELTHING TEMPERATURE ----------------
  #create 'output' folder, if it does not exist
  ifelse(!dir.exists(file.path('output')),dir.create(file.path('output')),FALSE)
  for(i in 1:dim(ip$data)[1]){
    #SEARCH TRIPLEX MOTIF
    triplex_motif <- triplex.search(Hsapiens[[ip$data[[1]][i]]][ip$data[[2]][i]:ip$data[[3]][i]],min_score=17,min_len=8)
    if(length(triplex_motif)>0){ #if triplex motif is present
      for(j in 1:length(triplex_motif)){ #iterate all motifs one by one
        triplex <- triplex.alignment(triplex_motif[j]) #triplex strand and loop sequence
        #save the triplex schematic diagram; name: chr_start_end_triplexMotifNumber
        png(paste0(Output,'/',
                   ip$data[['#chrom']][i],'_',
                   ip$data[['chromStart']][i],'_',
                   ip$data[['chromEnd']][i],'_',j,'.png')) # diagram file name
        
        triplex.diagram(triplex_motif[j]) #schematic diagram
        dev.off()
        
        #count nucleotides frequencies
        CC = str_count(paste(triplex[[1]],triplex[[2]],triplex[[3]],triplex[[4]]), "CC")
        TC = str_count(paste(triplex[[1]],triplex[[2]],triplex[[3]],triplex[[4]]), "TC")
        CT = str_count(paste(triplex[[1]],triplex[[2]],triplex[[3]],triplex[[4]]), "CT")
        TT = str_count(paste(triplex[[1]],triplex[[2]],triplex[[3]],triplex[[4]]), "TT")
        C = str_count(paste(triplex[[1]],triplex[[2]],triplex[[3]],triplex[[4]]), "C")
        
        #thermodynamics parameters calculations
        dH = -4.9*(CC) - 8.9*(TC+CT) - 7.4*(TT)
        dG37 = -3.00*(C)-0.65*(T)+1.65*(CC)+ 6.0 + (C)*(ip$pH-5.0)*(1.26-0.08*(CC))
        Tm = (310*dH)/(dH-dG37-310*ip$gas_constant_R*log((4/ip$Tcp),base = exp(1))) 
        seq = paste(triplex[[1]],triplex[[2]],triplex[[3]],triplex[[4]])
        
        #save all data in the dataframe
        data <- data%>%
          add_row(chr = ip$data[['#chrom']][i],
                  start = ip$data[['chromStart']][i],
                  end = ip$data[['chromEnd']][i],
                  description = ip$data[['description']][i],
                  width = width(triplex_motif[j]),
                  score = score(triplex_motif[j]),
                  pvalue = pvalue(triplex_motif[j]),
                  dH = dH,
                  dG37 = dG37,
                  Tm = Tm,
                  seq = seq)
      }
    }
    write_tsv(data, paste0(Output,"/","dna.triplex.GRCh38.ENCODE.cCREs.11.05.2021.tsv"))#write output
  }
  #data %>%View()
  print('Task is completed successfully !!!')
  return(data)
}

data <- dna_triplex_humangenome('11.05.2021.GRCh38.ENCODE.cCREs.txt',
                                'output')

#data <- dna_triplex_humangenome('C:\\Users\\himanshu\\Dropbox\\short_articles\\11.05.2021.GRCh38.ENCODE.cCREs.txt',
#                        'C:\\Users\\himanshu\\Dropbox\\short_articles\\output')


