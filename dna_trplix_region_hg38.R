if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Load necessary libraries and genomes.
BiocManager::install("triplex")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(triplex)
library(BSgenome.Hsapiens.UCSC.hg38)

triplex_search_genome_region <- function(chrNo,chrStart,chrEnd){
  #read genomic region
  t <- triplex.search(Hsapiens[[chrNo]][chrStart:chrEnd],min_score=17,min_len=8)
  ts <- t[order(score(t),decreasing=TRUE)]
  #Save triple 2D diagram
  jpeg('first_triplex_2D_diagram.jpg',res = 100)
  triplex.diagram(ts[1])
  dev.off()
  print('first_triplex_2D_diagram.jpg is generated...')  
  #triplex.3D(ts[1])
  
  # Export all triplexes into a GFF3 format file.
  library(rtracklayer)
  export(as(t, "GRanges"),"all_triplex_forming_regions.gff", version="3")
  print('all_triplex_forming_regions.gff is generated...')
  
  # Export all triplexes into a FASTA format file.
  writeXStringSet(as(t, "DNAStringSet"), file="all_triplex_forming_regions.fa", format="fasta")
  print('all_triplex_forming_regions.gff is generated...')
  
  # Show histogram for score distribution of detected triplexes.
  jpeg('triplex_distribution.histogram.jpg',res=100)
  hist(score(t), breaks=20)
  dev.off()
  print('triplex_distribution.histogram.jpg is generated...')
  
  # Show triplex distribution along the chromosome or other analysed sequence.
  jpeg('triplex_coverage.jpg',res=100)
  plot(coverage(ts[0:length(t)]), type="s", col="grey75")
  dev.off()
  print('triplex_coverage.jpg is generated...')
  return(ts)
}

ts = triplex_search_genome_region('chr1',1,10000)
ts

calculate_triplex_stability<-function(){
  
}































t <- triplex.search(Hsapiens[["chrX"]][1:100000],min_score=17,min_len=8)
ts <- t[order(score(t),decreasing=TRUE)]
triplex.diagram(ts[1])
triplex.3D(ts[1])

# Export all triplexes into a GFF3 format file.
library(rtracklayer)
export(as(t, "GRanges"),"test.gff", version="3")

# Export all triplexes into a FASTA format file.
writeXStringSet(as(t, "DNAStringSet"), file="test.fa", format="fasta")
# Show histogram for score distribution of detected triplexes.
hist(score(t), breaks=20)
# Show triplex distribution along the chromosome or other analysed sequence.
plot(coverage(ts[0:length(t)]), type="s", col="grey75")
