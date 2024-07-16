##for summarizing breseq data
library(stringr)
library(tidyverse)
source('R_functions.R')
library(readxl)
library(seqinr)


paths<-list.dirs("rawdata/breseq_try2/", recursive = F)
vcf_all <- data.frame(matrix(ncol = 11))
names(vcf_all) <- c("SAMP", "POS",  "REF",  "ALT", "AApos", "REFcod", "ALTcod", "FREQ", "GENE", "GENEproduct", "CAT")
for(p in paths){
  fannot <- paste(p, "/output/evidence/annotated.gd", sep="")
  annot <- read.table(file = fannot, sep = "\t",comment.char = "#", fill = T, col.names = 1:100) ##annoying. new format has tons of columns
  annot.short <- annot %>%
    filter(X1 %in% c("DEL", "SNP", "INS", "MOB"))  #descriptions here: https://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/output.html#mutation-display
  tmp <- data.frame(matrix(nrow=nrow(annot.short), ncol=11))
  collapsed <-  data.frame(col1 = do.call(paste, c(annot.short, sep=";")))
  names(tmp) <- c("SAMP", "POS",  "REF",  "ALT", "AApos", "REFcod", "ALTcod", "FREQ", "GENE", "GENEproduct", "CAT")
  tmp$SAMP <- tail(as.vector(str_split(p, "\\/", simplify = T)), n=1)
  tmp$POS <- str_extract(collapsed$col1, "hybrid\\;(\\d+)", group=1)
  tmp$REF <- str_extract(collapsed$col1, "aa_ref_seq\\=(\\w)", group=1)
  tmp$ALT <- str_extract(collapsed$col1, "aa_new_seq\\=(\\w)", group=1)
  tmp$AApos <- str_extract(collapsed$col1, "aa_position\\=(\\d+)", group=1)
  tmp$ALTcod <- str_extract(collapsed$col1, "codon_new_seq\\=(\\w+)", group=1)
  tmp$REFcod <- str_extract(collapsed$col1, "codon_ref_seq\\=(\\w+)", group=1)
  tmp$FREQ <- str_extract(collapsed$col1, "frequency\\=([\\d\\.\\w\\-]+)\\;", group=1)
  tmp$GENE <- str_extract(collapsed$col1, "gene_name\\=([\\w\\s]+)", group=1)
  tmp$GENEproduct <- str_extract(collapsed$col1, "gene_product\\=([\\w\\s\\d]+)", group=1)
  tmp$CAT <- str_extract(collapsed$col1, "mutation_category\\=([\\w\\s\\d\\_]+)", group=1)
  
  vcf_all <- rbind(vcf_all, tmp)
}
vcf.filt <- vcf_all %>%
  filter(!is.na(SAMP)) %>%
  mutate(FREQ = as.numeric(FREQ), POS = as.integer(POS), AApos = as.integer(AApos)) %>%
  filter(FREQ > 0.9) %>%
  filter(!POS %in% c(364936, 2269414)) %>% #remove variants present in ancestor
  mutate(SAMPnum = as.numeric(str_split(SAMP, "_", simplify = T)[,2])) %>%
  arrange(SAMPnum, POS)

#sample submission sheet for SeqCoast
samps <- read_excel("rawdata/3388_JamesVanLeuven_06092023.xlsx", sheet = 1, skip = 1) %>%
  select("SeqCoast Tube ID", "Sample Name")

names(samps) <- c("ID", "NAME")
  
##we had some extra smaples in this run. use just these:
samps <- samps %>%
  filter(ID %in% c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34))

row.names(samps) <- samps$ID

vcf.filt <- merge(x = vcf.filt, y=samps, by.x="SAMPnum", by.y="ID")
vcf.filt[is.na(vcf.filt$GENEproduct),]$GENEproduct <- vcf.filt[is.na(vcf.filt$GENEproduct),]$CAT
write.csv(vcf.filt, file="results/filtered_variants.csv", quote = F, row.names = F)


##which samples have phage stuff going on?
unique(vcf.filt[str_detect(vcf.filt$GENEproduct, '(?i)phage'),]$NAME)
##these mutations in all the Fern and Vegas resistant variants
##make heatmap to see what mutations are in what samples
##might consider consolitaing individual mutations by regions
##these samples 
#3388_28_S21_R - V
#3388_26_S19_R - V
#3388_2_S12_R - F
#3388_4_S29_R - F
#3388_27_S20_R - V
#3388_3_S23_R - F
phage_muts <- vcf.filt[str_detect(vcf.filt$GENEproduct, '(?i)phage'),]
phage_muts_wide <- phage_muts %>%
  select(NAME, POS, GENEproduct) %>%
  pivot_wider(names_from = NAME, values_from = GENEproduct) %>%
  arrange(as.numeric(POS))

p <- ggplot(phage_muts, aes(x=as.factor(POS), y=NAME, fill=GENEproduct, color=GENEproduct)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_jama() +
  scale_fill_jama() +
  labs(x="Genome position", y="")

ggsave("figs/phage_muts.png", p, units = "in", width = 10, height=4)
ggsave("figs/phage_muts.pdf", p, units = "in", width = 10, height=4)

##genome regions to extract for comparison. these coordinates are from
###breseq
#1327332 (1327000) to 1332482 (1332500) = 5500 bp
#and
#2421901 to 2423651 = 1750 bp
#and
#2925106 (2925000) to 2927691 (2927800) = 2800 bp

#blast hits elsewhere in genome. There are several smaller (~2.5kb) hits 
##across the genome, but the three below are the largest.
##plan is to extract these regions and run stuff through clinker
#1300000-1350000  - 50kb
#2655000-2675000 - 20kb
#2910000-2950000 - 40kb


##this appears to be quite a mess. spent a lot of time trying to figure this out
##compared running breseq on 
    ## 1 - plarvae genome
    ## 2 - plarvae genome and challenge phage genome (multifasta file)
    ## 3 - challenge phage genome appended plarvae genome
##for Fern yb (3388_2_S12_R) this pretty much fixed all issue. In the appended run, there were no polymorphic sites
##and the fern genome has reads mapping to it with some changes (~300 snps)
##coverage for phage genome appears to be slightly higher (~200) than for plarvae (~150)

##for other challenges, appending the phage genome to plarvae did not fix the issue
##tested for all the V phages

#for example, the F-resistant isolate x still has lots of polymorphic sites in phage-annotated regions
#also happens to have no apparent fixed mutations
#looked into this. Seemed to find reads that support the presence of another phage or bacterial isolate
#see SNP 26. checked to see if reads that support polymorphism match phage F (they do not) or 3650 anc (they do not)
#blasting a few of these reads shows that they match 100% to other plarvae phages or plarva isolates (100% identity for 150 bp)

##there are far fewer differences in the V phage genome (from reference to what we are using)
##only around 20 SNPs
##V coverage is very close to plarvae coverage

##scanned alignments by eye and identified one place in Fern-resistant and Vegas-
##resistant isolates that there is some sort of recombination.
##it looks like Fern and Vegas integrated. 
##the Fern looks quite clear. Both breseq and my alignmnets support this.
##for Vegas is is a little weirder because the coverage of Vegas is not even.
##I copied the phage genomes into the bacterial genomes and re-aligned. 
##coverage is perfectly even across the junctions. 




ref <- read.fasta("ref_seq/Plarvae_3650.fasta")
seq1 <- paste(ref$`Paenibacillus_larvae_strain_3650_(hybrid)`[1300000:1350000], collapse = "")
seq2 <- paste(ref$`Paenibacillus_larvae_strain_3650_(hybrid)`[2655000:2675000], collapse = "")
seq3 <- paste(ref$`Paenibacillus_larvae_strain_3650_(hybrid)`[2910000:2950000], collapse = "")

ref_regions <- paste(seq1, seq2, seq3)
ref_regions <- str_replace_all(ref_regions, " ", "")
#ref_regions <- paste(seq1, seq2, seq3, collapse = "")

write.fasta(ref_regions, "plarvae_1300000-1350000_2655000-2675000_2910000-2950000","ref_seq/clinker/plarvae_3650_selected.fasta", as.string = T)


vcf.nophage <- vcf.filt[!str_detect(vcf.filt$GENEproduct, '(?i)phage'),]

write.csv(vcf.nophage, file="results/filtered_variants_nophage.csv", quote = F, row.names = F)

isolates <- vcf.nophage %>%
  group_by(NAME) %>%
  summarize(n())
View(isolates)

length(table(vcf.nophage$POS)) ###only 23 unique mutations observed

vcf.nophage %>%
  group_by(POS) %>%
  summarize(n=n()) %>%
  arrange(desc(n))
