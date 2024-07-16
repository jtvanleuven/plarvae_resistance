library(Rsamtools)
library(readxl)

#read excel file for sample names
get_sample_names <- function(){
  samples <- read_excel("rawdata/3388_JamesVanLeuven_06092023.xlsx", skip = 1)
  samples <- samples %>%
    rename("ID"="SeqCoast Tube ID", "NAME"="Sample Name") %>%
    filter(!NAME=="Y3650") #reference 
  
  samples$NAME <- str_replace(samples$NAME, '.r.Y3650 ', "_")
  samples$NAME <- str_replace(samples$NAME, ".r.Vx ", "_")
  samples$NAME <- str_replace(samples$NAME, "Soct", "Scot")
  samples$NAME <- str_replace(samples$NAME, " \\(1/13\\)", "")
  samples$NAME <- str_replace(samples$NAME, " ", "")
  return(data.frame(samples))
}

#import data from excel
get_data <- function(){
  sheets <- 3:20 #resistants
  grow <- data.frame(matrix(ncol=4))
  names(grow) <- c("Time", "isolate", "rep", "OD600")
  for(i in sheets){
    dat <- read_excel('rawdata/Resistant Host Growth Rate Curve.xlsx', sheet=i)
    dat.l <- dat %>%
      select(c(1:4)) %>%
      pivot_longer(!Time, names_to = "isolate", values_to = "OD600")
    dat.l$rep <- str_split(dat.l$isolate, " ", simplify = T)[,3]
    dat.l$isolate <- substr(dat.l$isolate,1,nchar(dat.l$isolate)-2)
    dat.l$isolate <- str_replace_all(dat.l$isolate, "\\s+", "_")
    dat.l <- dat.l[,names(grow)]
    grow <- rbind(grow, dat.l)
  }
  grow <- na.omit(grow)
  rm(dat, dat.l)
  
  anc <- read_excel('rawdata/Resistant Host Growth Rate Curve.xlsx', sheet=2, skip = 2)
  anc <- anc[,1:13]
  names(anc) <- str_split(names(anc), "\\.", simplify=T)[,1]
  names(anc)[1] <- "Time"
  anc.l <- anc %>%
    pivot_longer(!Time, names_to = "isolate", values_to = "OD600")
  anc.l$rep <- substr(anc.l$isolate,nchar(anc.l$isolate)-1, nchar(anc.l$isolate))
  anc.l$isolate <- substr(anc.l$isolate,1,nchar(anc.l$isolate)-2)
  anc.l <- anc.l %>%
    group_by(Time, isolate) %>%
    mutate(rep = 1:6)
  
  grow.all <- rbind(grow, anc.l)
}


bamcoverage <- function (bamfile) {
  # read in the bam file
  bam <- scanBam(bamfile)[[1]] # the result comes in nested lists
  # filter reads without match position
  ind <- ! is.na(bam$pos)
  ## remove non-matches, they are not relevant to us
  bam <- lapply(bam, function(x) x[ind])
  ranges <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))
  ## names of the bam data frame:
  ## "qname"  "flag"   "rname"  "strand" "pos"    "qwidth"
  ## "mapq"   "cigar"  "mrnm"   "mpos"   "isize"  "seq"    "qual"
  ## construc: genomic ranges object containing all reads
  ranges <- GRanges(seqnames=Rle(bam$rname), ranges=ranges, strand=Rle(bam$strand), flag=bam$flag, readid=bam$rname )
  ## returns a coverage for each reference sequence (aka. chromosome) in the bam file
  return (mean(coverage(ranges)))      
}
