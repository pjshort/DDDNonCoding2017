# process de novos and various annotations using GRanges

require(GenomicRanges)

genomicus_target_gene <- function(de_novos, genomicus){
  dn = GRanges(seqnames=Rle(paste0("chr", de_novos$chr)), ranges = IRanges(start = de_novos$pos, end = de_novos$pos + 1))
  gen = GRanges(seqnames=Rle(paste0("chr", genomicus$CREs_chr.1)), ranges = IRanges(start = genomicus$Start.hg19.6., end = genomicus$End.hg19.7.))
  mcols(gen) = genomicus[ , c("Predicted.target.20.")]

  # find overlap between denovos and genomicus annotation
  hits = findOverlaps(dn, gen)
  dn_hits_idx = queryHits(hits) # get index of de novos with genomicus hit
  gen_hits_idx = subjectHits(hits) # get index of genomicus ranges to retrieve target gene

  target_genes = rep("NONE", nrow(de_novos))

  # all de novos with no genom. target will remain NONE
  target_genes[dn_hits_idx] = gsub(" $", "", as.character(genomicus$Predicted.target.20.[gen_hits_idx]))

  return(target_genes)
}

count_genomicus_hits <- function(de_novos, genomicus){
  dn = GRanges(seqnames=Rle(paste0("chr", de_novos$chr)), ranges = IRanges(start = de_novos$pos, end = de_novos$pos + 1))
  gen = GRanges(seqnames=Rle(paste0("chr", genomicus$CREs_chr.1.)), ranges = IRanges(start = genomicus$Start.hg19.6., end = genomicus$End.hg19.7.))

  # find overlap between denovos and genomicus annotation
  hits = findOverlaps(dn, gen)
  dn_hits_idx = unique(queryHits(hits)) # get index of de novos with genomicus hit

  count = length(dn_hits_idx)

  return(count)
}

genomicus_conservation_scores <- function(de_novos, genomicus){
  dn = GRanges(seqnames=Rle(paste0("chr", de_novos$chr)), ranges = IRanges(start = de_novos$pos, end = de_novos$pos + 1))
  gen = GRanges(seqnames=Rle(paste0("chr", genomicus$CREs_chr.1.)), ranges = IRanges(start = genomicus$Start.hg19.6., end = genomicus$End.hg19.7.))

  # find overlap between denovos and genomicus annotation
  hits = findOverlaps(dn, gen)
  genomicus_CNEs_idx = subjectHits(hits) # get index of de novos with genomicus hit
  conservation_scores = genomicus$Conservation.score..12.[genomicus_CNEs_idx]

  return(conservation_scores)

}

count_bed_hits <- function(de_novos, bed){
  # assumes that the bed input has at least three columns with chr, start, end as the first three (all others ignored)
  dn = GRanges(seqnames=Rle(paste0("chr", de_novos$chr)), ranges = IRanges(start = de_novos$pos, end = de_novos$pos + 1))
  gen = GRanges(seqnames=Rle(bed[,1]), ranges = IRanges(start = bed[,2], end = bed[,3]))

  # find overlap between denovos and genomicus annotation
  hits = findOverlaps(dn, gen)
  dn_hits_idx = queryHits(hits) # get index of de novos with genomicus hit

  count = length(dn_hits_idx)

  return(count)
}

filter_with_bed <- function(de_novos, bed){
  # assumes that the bed input has at least three columns with chr, start, end as the first three (all others ignored)

  if (any(!grepl("^chr", de_novos$chr))) {
    de_novos$chr = paste0("chr", de_novos$chr)
  }
  if (any(!grepl("^chr", bed[,1]))) {
    bed[ ,1] = paste0("chr", bed[ ,1])
  }

  dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$pos, end = de_novos$pos))
  gen = GRanges(seqnames=Rle(bed[,1]), ranges = IRanges(start = bed[,2], end = bed[,3]))

  # find overlap between denovos and genomicus annotation
  hits = findOverlaps(dn, gen)
  dn_hits_idx = unique(queryHits(hits)) # get index of de novos with genomicus hit

  de_novos = de_novos[dn_hits_idx, ]

  return(de_novos)
}

overlap_with_bed <- function(regions, bed){
  # assumes that the bed input has at least three columns with chr, start, end as the first three (all others ignored)
  reg = GRanges(seqnames=Rle(regions$chr), ranges = IRanges(start = regions$start, end = regions$stop))
  bed = GRanges(seqnames=Rle(bed[,1]), ranges = IRanges(start = bed[,2], end = bed[,3]))

  # intersect regions with bed range (i.e. genomicus)
  reg_bed_overlap = intersect(reg, bed)
  reg = data.frame(chr = as.character(reg_bed_overlap@seqnames),
                   start = as.integer(reg_bed_overlap@ranges@start),
                   stop = as.integer(reg_bed_overlap@ranges@start + reg_bed_overlap@ranges@width - 1))
  reg$region_id <- paste(reg$chr, reg$start, reg$stop, sep = ".")

  return(reg)
}


split_to_json <- function(de_novos, by){
  probands_by_factor = split(de_novos$person_stable_id, de_novos[, by])
  probands_by_factor = Filter(function(x) length(x) > 0, probands_by_factor)
  region_json = toJSON(probands_by_factor, pretty = TRUE)
  save_name = sprintf("../data/probands_by_%s.json", by)
  sink(save_name)
  cat(region_json)
  sink()
}

# get pLI based on exac data set
exac_pli = read.table("../data/exac_pLI.txt", sep = "\t", header = TRUE)

gene_pli = function(gene_name){
  if (gene_name %in% exac_pli$gene){
    pli = as.numeric(subset(exac_pli, gene == gene_name)$pLI)
  } else {
    pli = 0
  }
  return(pli)
}

get_closest_gene <- function(de_novos, CNEs){
  # assumes that the CNE input has at least three columns with chr, start, end as the first three (all others ignored)

  if (any(!grepl("^chr", de_novos$chr))) {
    vars$chr = paste0("chr", vars$chr)
  }
  if (any(!grepl("^chr", CNEs$chr))) {
    CNEs$chr = paste0("chr", CNEs$chr)
  }

  if ("end" %in% colnames(vars)){ # region instead of de novos - use the first position of the region to get closest gene!
    dn = GRanges(seqnames=Rle(vars$chr), ranges = IRanges(start = vars$start, end = vars$start + 1))
  } else {
    dn = GRanges(seqnames=Rle(vars$chr), ranges = IRanges(start = vars$pos, end = vars$pos + 1))
  }

  cne = GRanges(seqnames=Rle(CNEs$chr), ranges = IRanges(start = CNEs$start, end = CNEs$stop))

  # find overlap between denovos and annotated CNEs
  hits = findOverlaps(dn, cne)
  dn_hits_idx = queryHits(hits) # get index of de novos
  CNE_hits_idx = subjectHits(hits) # get index of CNEs

  closest_gene = CNEs$closest_gene[CNE_hits_idx]

  return(closest_gene)
}

get_CNE_annotation <- function(de_novos, CNEs){
  # assumes that the CNE input has at least three columns with chr, start, end as the first three (all others ignored)

  if (any(!grepl("^chr", de_novos$chr))) {
    de_novos$chr = paste0("chr", de_novos$chr)
  }
  if (any(!grepl("^chr", CNEs$chr))) {
    CNEs$chr = paste0("chr", CNEs$chr)
  }

  if ("end" %in% colnames(de_novos)){ # region instead of de novos - use the first position of the region to get closest gene!
    dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$start, end = de_novos$start + 1))
  } else {
    dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$pos, end = de_novos$pos + 1))
  }

  cne = GRanges(seqnames=Rle(CNEs$chr), ranges = IRanges(start = CNEs$start, end = CNEs$stop))

  # find overlap between denovos and annotated CNEs
  hits = findOverlaps(dn, cne)
  dn_hits_idx = queryHits(hits) # get index of de novos
  CNE_hits_idx = subjectHits(hits) # get index of CNEs

  enhancer = CNEs$enhancer[CNE_hits_idx]
  conserved = CNEs$conserved[CNE_hits_idx]
  heart = CNEs$heart[CNE_hits_idx]

  # assumption in these cases is that information quality is enhancer > heart > conserved
  case1 = enhancer == 1 # call this enhancers
  case2 = conserved == 1 & enhancer == 0 | conserved == 1 & heart == 0 # these are only high phastcons score
  case3 = heart == 1 & enhancer == 0 # heart set - implied that conserved could be 1 or 0

  annotation = rep("NONE", length(CNE_hits_idx))
  annotation[case1] = "Enhancer"
  annotation[case2] = "Conserved"
  annotation[case3] = "Heart"

  return(annotation)
}

get_region_id <- function(de_novos, CNEs){
  # assumes that the CNE input has at least three columns with chr, start, end as the first three (all others ignored)

  if (any(!grepl("^chr", de_novos$chr))) {
    de_novos$chr = paste0("chr", de_novos$chr)
  }
  if (any(!grepl("^chr", CNEs$chr))) {
    CNEs$chr = paste0("chr", CNEs$chr)
  }

  if ("end" %in% colnames(de_novos)){ # region instead of de novos - use the first position of the region to get closest gene!
    dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$start, end = de_novos$start))
  } else {
    dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$pos, end = de_novos$pos))
  }

  cne = GRanges(seqnames=Rle(CNEs$chr), ranges = IRanges(start = CNEs$start, end = CNEs$stop))

  # find overlap between denovos and annotated CNEs
  hits = findOverlaps(dn, cne)
  dn_hits_idx = queryHits(hits) # get index of de novos
  CNE_hits_idx = subjectHits(hits) # get index of CNEs


  return(CNEs$region_id[CNE_hits_idx])
}

get_gene <- function(de_novos, genes){
  # assumes that the CNE input has at least three columns with chr, start, end as the first three (all others ignored)
  
  if (any(!grepl("^chr", de_novos$chr))) {
    de_novos$chr = paste0("chr", de_novos$chr)
  }
  if (any(!grepl("^chr", genes$chr))) {
    genes$chr = paste0("chr", genes$chr)
  }
  
  if ("end" %in% colnames(de_novos)){ # region instead of de novos - use the first position of the region to get closest gene!
    dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$start, end = de_novos$start + 1))
  } else {
    dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$pos, end = de_novos$pos + 1))
  }
  
  cne = GRanges(seqnames=Rle(genes$chr), ranges = IRanges(start = genes$start, end = genes$stop))
  
  # find overlap between denovos and annotated genes
  hits = findOverlaps(dn, cne)
  dn_hits_idx = queryHits(hits) # get index of de novos
  hits_idx = subjectHits(hits) # get index of genes
  
  
  return(genes$gene[hits_idx])
}


get_chromHMM <- function(vars, bed){
  # assumes that the CNE input has at least three columns with chr, start, end as the first three (all others ignored)
  
  if (any(!grepl("^chr", vars$chr))) {
    vars$chr = paste0("chr", vars$chr)
  }
  if (any(!grepl("^chr", bed$chr))) {
    bed$chr = paste0("chr", bed$chr)
  }
  
  if ("stop" %in% colnames(vars) | "end" %in% colnames(vars)){ # region instead of de novos - use the first position of the region to get closest gene!
    v = GRanges(seqnames=Rle(vars$chr), ranges = IRanges(start = vars$start, end = vars$start + 1))
  } else {
    v = GRanges(seqnames=Rle(vars$chr), ranges = IRanges(start = vars$pos, end = vars$pos))
  }
  
  b = GRanges(seqnames=Rle(bed$chr), ranges = IRanges(start = bed$start, end = bed$stop-1))
  
  # find overlap between denovos and annotated CNEs
  hits = findOverlaps(v, b)
  var_hits_idx = queryHits(hits) # get index of de novos
  bed_hits_idx = subjectHits(hits) # get index of CNEs
  
  states = bed$chromHMM[bed_hits_idx]
  
  return(states)
}

get_sequence <- function(chr, start, stop, version = "hg19") {

  # input: (multiple) chr, start, stop, hg version (defaults to hg19)
  # output: list of sequences as DNAStrings object for each input

  if (version == "hg19"){
    library(BSgenome.Hsapiens.UCSC.hg19)
  } else if (version == "hg18"){
    library(BSgenome.Hsapiens.UCSC.hg18) # TODO need to add download of hg18 to build.R
  }

  if (!all(grepl(pattern = "^chr", chr))){  # assert that chromosome column have chr in front
    warning("Not all entries in the chromosome column start with \"chr\" - try reformatting this column e.g. \"chrX\" instead of \"X\" with paste0(\"chr\",chr_number")
    chr = paste0("chr", chr)
  }

  seqs = getSeq(Hsapiens, chr, start, stop)
  return(seqs)
}


sim_observed_slice = function(de_novos, sim_de_novos, split_by_diagnosis = TRUE, factor1 = NULL, factor2 = NULL){
 u = subset(de_novos, diagnosed == FALSE)
 d = subset(de_novos, diagnosed == TRUE)

 if (!is.null(factor1)){
   u_factor1 = subset(u, u[,factor] == TRUE)
   u_not_factor1 = subset(u, u[,factor] == FALSE)
   d_factor1 = subset(u, u[,factor] == TRUE)
   d_not_factor1 = subset(u, u[,factor] == FALSE)
   if (!is.null(factor2)){
     u_factor12 = subset(u_factor1, u[,factor] == TRUE)
     u_factor1only = subset(u, u[,factor] == TRUE)
     u_factor2only = subset(u, u[,factor] == TRUE)
     u_no_factors = subset(u, u[,factor] == TRUE)
     u_factor12 = subset(u, u[,factor] == TRUE)
     u_factor1only = subset(u, u[,factor] == TRUE)
     u_factor2only = subset(u, u[,factor] == TRUE)
     u_no_factors = subset(u, u[,factor] == TRUE)


   }
 }



}



get_distance_to_closest_gene <- function(variants, gencode = gencode){
  # variant input must have 'closest_gene' annotation
  genes = variants$closest_gene
  
  gencode_start = gencode$start[match(genes, gencode$gene)]
  gencode_stop = gencode$stop[match(genes, gencode$gene)]
  gencode_strand = gencode$strand[match(genes, gencode$gene)]
  
  # if strand is '+', then take gencode START to find distance. if strand is '-' then take gencode STOP to find distance.
  gencode_TSS = sapply(seq(1,length(gencode_strand)), function(i) ifelse(gencode_strand[i] == "+", gencode_start[i], gencode_stop[i]))
  
  if ("pos" %in% colnames(variants)) {
    distance = abs(variants$pos - gencode_TSS)
  } else {
    distance = abs(variants$start - gencode_TSS)
  }
  
  return(distance)
}

get_element_closest_gene <- function(elements, gencode){
  # assumes that the CNE input has at least three columns with chr, start, end as the first three (all others ignored)
  
  # load gencode annotation with strand
  #gencode = read.table("../data/gencode_protein_coding_genes_v19_+strand.txt", header = TRUE, sep = "\t")
  
  if (any(!grepl("^chr", elements$chr))) {
    elements$chr = paste0("chr", elements$chr)
  }
  if (any(!grepl("^chr", gencode$chr))) {
    gencode$chr = paste0("chr", gencode$chr)
  }
  
  if ('pos' %in% colnames(elements)){
    e = GRanges(seqnames=Rle(elements$chr), ranges = IRanges(start = elements$pos, end = elements$pos))
  } else {
    e = GRanges(seqnames=Rle(elements$chr), ranges = IRanges(start = elements$start, end = elements$start))
  }
  genes = GRanges(seqnames=Rle(gencode$chr), ranges = IRanges(start = gencode$start, end = gencode$stop), strand = gencode$strand)
  
  closest_gene_idx = nearest(e, genes)
  closest_gene = as.character(gencode[closest_gene_idx, "gene"])
    
  return(closest_gene)
}

get_distance_to_nearest_exon <- function(elements, exons){
  # assumes that the CNE input has at least three columns with chr, start, end as the first three (all others ignored)
  
  # load gencode annotation with strand
  #gencode = read.table("../data/gencode_protein_coding_genes_v19_+strand.txt", header = TRUE, sep = "\t")
  
  if (any(!grepl("^chr", elements$chr))) {
    elements$chr = paste0("chr", elements$chr)
  }
  if (any(!grepl("^chr", exons$chr))) {
    exons$chr = paste0("chr", exons$chr)
  }
  
  if ('pos' %in% colnames(elements)){
    e = GRanges(seqnames=Rle(elements$chr), ranges = IRanges(start = elements$pos, end = elements$pos))
  } else {
    e = GRanges(seqnames=Rle(elements$chr), ranges = IRanges(start = elements$start, end = elements$start))
  }
  genes = GRanges(seqnames=Rle(exons$chr), ranges = IRanges(start = exons$start, end = exons$stop), strand = exons$strand)
  
  closest_gene_idx = nearest(e, genes)
  closest_gene = as.character(exons[closest_gene_idx, "gene"])
  distance = exons[closest_gene_idx, "start"] - elements$start
  
  return(distance)
}


