library(stringr)
library(tidyverse)
library(readxl)
library(gcplyr)
library(lme4)
library(moments)   #for calculating skewness
library(survival)
library(survminer)
source('R_functions.R')

##compare growth rates of different isolates
####resistant strain growth data
grow_all <- get_data()

#transformations do not like zeros
grow_all$OD600 <- grow_all$OD600 + 0.001
##have decided to leave out Hal because I cannot confirm its identity
##also, Scot y and Scot a for some reason have identical growth
##measurements. Dropping Scot a
grow_all <- grow_all %>%
  filter(!isolate %in% c("Scot_a", "Hal_ab", "Hal_Bb", "Hal_yb"))


ex_dat_mrg <- mutate(group_by(grow_all, isolate, rep),
                     deriv = calc_deriv(x = Time, y = OD600))


# Now let's plot the derivative
ggplot(data = dplyr::filter(ex_dat_mrg, isolate=="Scot_B"),
       aes(x = Time, y = deriv)) +
  geom_line() +
  facet_wrap(~rep, scales = "free")


ex_dat_mrg <- mutate(group_by(grow_all, isolate, rep),
                     deriv_percap5 = calc_deriv(x = Time, y = OD600, 
                                                percapita = TRUE, blank = 0,
                                                window_width_n = 5, trans_y = "log"))


ex_dat_mrg <- mutate(group_by(grow_all, isolate, rep),
                     deriv_percap5 = calc_deriv(x = Time, y = OD600, 
                                                percapita = TRUE, blank = 0,
                                                window_width_n = 5, trans_y = "log"),
                     doub_time = doubling_time(y = deriv_percap5))

ggplot(data = dplyr::filter(ex_dat_mrg, isolate=="Scot_B"),
       aes(x = Time, y = deriv_percap5)) +
  geom_line() +
  facet_wrap(~rep, scales = "free")

#look at data in plot format, compare replicates using these values
ggplot(data = dplyr::filter(ex_dat_mrg),
       aes(x = Time, y = deriv_percap5, color=rep)) +
  geom_line() +
  facet_wrap(~isolate, scales = "free")

#look at data in plot format, compare replicates using these values
ggplot(data = dplyr::filter(ex_dat_mrg),
       aes(x = Time, y = doub_time, color=rep)) +
  geom_line() +
  facet_wrap(~isolate, scales = "free")


ex_dat_mrg_sum <- 
  summarize(group_by(ex_dat_mrg, isolate, rep),
            n=n(),
            max_dens = max_gc(OD600, na.rm = TRUE),
            max_time = extr_val(Time, which_max_gc(OD600)),
            lag_time = lag_time(y = OD600, x = Time, 
                                deriv = deriv_percap5),
            max_percap = max_gc(deriv_percap5),
            max_percap_time = Time[which_max_gc(deriv_percap5)],
            max_percap_dens = OD600[which_max_gc(deriv_percap5)],
            doub_time = doubling_time(y = max_percap),                      ###max growth rate
            auc = auc(x = Time, y = OD600))

ex_dat_mrg_sum$phage <- str_split(ex_dat_mrg_sum$isolate, "_", simplify = T)[,1]
ex_dat_mrg_sum[ex_dat_mrg_sum$phage=="Y3650 Anc",]$phage <- "Y3650"

##next steps
###summarize replicates/phage level
###perform stats: https://m-clark.github.io/mixed-models-with-R/random_intercepts.html

ggplot(ex_dat_mrg_sum, aes(x=isolate, y=max_dens, color=isolate)) +
  geom_jitter() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  facet_wrap(~phage, nrow=1, scales = "free_x")

ggplot(ex_dat_mrg_sum, aes(x=isolate, y=max_percap_dens, color=isolate)) +
  geom_jitter() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  facet_wrap(~phage, nrow=1, scales = "free_x")

ggplot(ex_dat_mrg_sum, aes(x=isolate, y=doub_time, color=isolate)) +
  geom_jitter() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  facet_wrap(~phage, nrow=1, scales = "free_x")

ggplot(ex_dat_mrg_sum, aes(x=isolate, y=auc, color=isolate)) +
  geom_jitter() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  facet_wrap(~phage, nrow=1, scales = "free_x")

ggplot(ex_dat_mrg_sum, aes(x=isolate, y=lag_time, color=isolate)) +
  geom_jitter() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  facet_wrap(~phage, nrow=1, scales = "free_x")



###do some stats

#think we should use these metrics:
  #max density
  #doubling time
  #lag time


#try this: https://statsandr.com/blog/anova-in-r/
mod_dat <- ex_dat_mrg_sum %>%
  select(isolate, phage, rep, n, max_dens, lag_time, doub_time) %>%
  mutate(isolate=str_remove_all(isolate, " Anc")) %>%
  mutate(phage=as.factor(phage)) %>%
  mutate(isolate=as.factor(isolate))


mod_dat$phage <- relevel(mod_dat$phage, ref="Y3650")
mod_dat$isolate <- relevel(mod_dat$isolate, ref="Y3650")

#code to extract sig p-values from stats below
pvals <- data.frame(matrix(nrow=length(unique(mod_dat$isolate))-1, ncol=7)) #subtract one because we don't want to count the reference
names(pvals) <- c("isolate", 
                  "coeff_dens", 
                  "pval_dens", 
                  "coeff_doub", 
                  "pval_doub",
                  "coeff_lag",
                  "pval_lag")

##data frame to save predictor variable significance

p <- 0.05 ###change pvalue here
#max density
#anova to see if phage and isolate are important
res_aov <- aov(max_dens ~ phage+isolate, data=mod_dat)
summary(res_aov) #yes, both are significant
#linear model to get coeffciencts for each group
res_lm <- lm(max_dens ~ phage, data=mod_dat)
summary(res_lm)
res_lm2 <- lm(max_dens ~ isolate, data=mod_dat)
summary(res_lm2)
coeffs <- data.frame(summary(res_lm2)$coefficients)
coeffs <- coeffs %>% filter(!row.names(coeffs)=="(Intercept)") #drop intercept 
pvals[,1] <- row.names(coeffs)
sig <- which(coeffs$Pr...t.. < p)
pvals$coeff_dens[sig] <- coeffs$Estimate[sig]
pvals$pval_dens[sig] <- coeffs$Pr...t..[sig]

#lag
#anova to see if phage and isolate are important
res_aov <- aov(lag_time ~ phage+isolate, data=mod_dat)
summary(res_aov) 
#linear model to get coeffciencts for each group
res_lm <- lm(lag_time ~ phage, data=mod_dat)
summary(res_lm)
res_lm2 <- lm(lag_time ~ isolate, data=mod_dat)
summary(res_lm2)
coeffs <- data.frame(summary(res_lm2)$coefficients)
coeffs <- coeffs %>% filter(!row.names(coeffs)=="(Intercept)") #drop intercept 
sig <- which(coeffs$Pr...t.. < p) 
pvals$coeff_lag[sig] <- coeffs$Estimate[sig]
pvals$pval_lag[sig] <- coeffs$Pr...t..[sig]


#double time
#anova to see if phage and isolate are important
res_aov <- aov(doub_time ~ phage+isolate, data=mod_dat)
summary(res_aov) #yes, both are significant
#linear model to get coeffciencts for each group
res_lm <- lm(doub_time ~ phage, data=mod_dat)
summary(res_lm)
res_lm2 <- lm(doub_time ~ isolate, data=mod_dat)
summary(res_lm2)
coeffs <- data.frame(summary(res_lm2)$coefficients)
coeffs <- coeffs %>% filter(!row.names(coeffs)=="(Intercept)") #drop intercept 
sig <- which(coeffs$Pr...t.. < p) 
pvals$coeff_doub[sig] <- coeffs$Estimate[sig]
pvals$pval_doub[sig] <- coeffs$Pr...t..[sig]

pvals$isolate <- str_remove(pvals$isolate, "isolate")


#want to compare to model with gene identity
#need to import more data
vcf <- read_excel("rawdata/vcf_all.xlsx")
table(vcf$GENE)
#one of these two mutations show up in all the samples
ancs <- c(841945, 3045798)
vcf_short <- vcf %>%
  filter(!POS %in% ancs)
length(unique(vcf_short$SAMP)) #three samples have only these mutations. hmm. seems like maybe we are missing some regions in the reference.

samples <- get_sample_names()
samples <- samples %>%
  select("ID", "NAME")
vcf_short <- merge(vcf_short, samples, by.x="SAMP", by.y="ID")

mod_dat$inter <- NA #intergenic
mod_dat$prsa <- NA #prsA
mod_dat$dnaj <- NA #dnaJ
mod_dat$trans <- NA #transcriptional regulator-like protein
mod_dat$pin <- NA #putative PIN and TRAM-domain containing protein precursor
for(i in 1:nrow(mod_dat)){
  muts <- vcf_short %>%
    filter(NAME == mod_dat$isolate[i]) %>%
    select(GENE)
  if("intergenic" %in% muts$GENE){
    mod_dat[i,]$inter <- 1
  }
  if("prsA" %in% muts$GENE){
    mod_dat[i,]$prsa <- 1
  }
  if("dnaJ" %in% muts$GENE){
    mod_dat[i,]$dnaj <- 1
  }
  if("transcriptional regulator-like protein" %in% muts$GENE){
    mod_dat[i,]$trans <- 1
  }
  if("putative PIN and TRAM-domain containing protein precursor" %in% muts$GENE){
    mod_dat[i,]$pin <- 1
  }
}
mod_dat[is.na(mod_dat)] <- 0

#Unity_y - have growth data for, but no sequence
#Unity_x - Emma says something is wrong with this sample. I don't see any mutations in it. growth looks reduced.

#compare models
#linear model to get coeffciencts for each group
#comparing AIC values suggests that "Gene" is an important predictor variable,
#but I really don't see this in plotting the data. 
{
phage_lm <- lm(max_dens ~ phage, data=mod_dat)
summary(phage_lm)
isolate_lm <- lm(max_dens ~ isolate, data=mod_dat)
summary(isolate_lm)
genes_lm <- lm(max_dens ~ inter+prsa+dnaj+trans+pin, data=mod_dat)
summary(genes_lm)

AIC(phage_lm)
AIC(isolate_lm)
AIC(genes_lm)

models <- data.frame(matrix(nrow=3, ncol = 3))
names(models) <- c("model", "AdjR2", "AIC")
models$model <- c("Challenge phage", "Isolate", "Gene mutated")
models$AIC <- c(AIC(phage_lm),
                AIC(isolate_lm),
                AIC(genes_lm))
models$AdjR2 <- c(summary(phage_lm)$adj.r.squared,
               summary(isolate_lm)$adj.r.squared,
               summary(genes_lm)$adj.r.squared)

plot_dat <- mod_dat %>%
  mutate(none=str_detect(phage, "Y3650")) %>%
  mutate(none=as.integer(as.logical(none))) %>%
  pivot_longer(cols=c("inter","prsa", "dnaj", "trans", "pin", "none")) %>%
  filter(value == 1)

ggplot(plot_dat, aes(x=isolate, y=max_dens, color=isolate)) +
  geom_jitter() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  facet_wrap(~name, nrow=1, scales = "free_x")
}


##look at Emma's plating effciency stuff. 
##I think there is something interesting here, but really hard to analyze in the Excel format
##First thing is to simply look at EOP of original phage stock. We expect to see a very low EOP. 
##Probably plaques are phage variants that can grow on resistant hosts. 
##for this, need to extract certain data from Emma's excel file
titers <- read_excel("rawdata/Picked plaques Titers Comparison.xlsx", sheet = 3) ##must copy to working excel file to Project directory.
titers[titers==0] <- 1 #pseudocount to avoid infinity

titers$eop <- titers$resistantHost / titers$titer3650

eop.l <- titers %>%
  drop_na() %>%
  pivot_longer(cols=c("titer3650", "resistantHost"), values_to = "titer")

ggplot(eop.l, aes(x=name, log(titer, base = 10), group=host)) +
  geom_point() +
  geom_line() +
  facet_wrap(~evolvedPlaque) +
  scale_x_discrete(limits=c("titer3650", "resistantHost"))


mean(na.omit(titers[titers$evolvedPlaque=="n",]$eop))
sd(titers[titers$evolvedPlaque=="n",]$eop, na.rm = T)
mean(na.omit(titers[titers$evolvedPlaque=="y",]$eop))
sd(titers[titers$evolvedPlaque=="y",]$eop, na.rm = T)
min(na.omit(titers[titers$evolvedPlaque=="y",]$eop))



titers_cross <- read_excel("rawdata/Picked plaques Titers Comparison.xlsx", sheet = 4) ##must copy to working excel file to Project directory.
#drop challenged hosts
get <- str_split(titers_cross$host, " ", simplify = T)[,1] == titers_cross$phage
titers_cross_short <- titers_cross %>%
  filter(!get)   #have pseudocount of 1 in there
lookup <- titers_cross %>%
  filter(get)   #have pseudocount of 1 in there
row.names(lookup) <- lookup$host
hist(-1 * titers_cross_short$diffFromAnc, breaks=20) #a change in the positive direction
#indicates that the phage plates better on the resistant host than the phage
#used to challenge that host.
shapiro.test(titers_cross_short$diffFromAnc)
skewness(titers_cross$diffFromAnc)
#try t-test
titers_cross_short$compEOP <- NA
for(i in 1:nrow(titers_cross_short)){
  titers_cross_short[i,]$compEOP <- lookup[titers_cross_short[i,]$host,]$logEOP
}

  
t.test(x=titers_cross_short$relativeEOP, y = titers_cross_short$compEOP, paired = T, p.adjust="BH")
t.test(x=titers_cross_short$logEOP, y = titers_cross_short$compEOP, paired = T, p.adjust="BH")
p.adjust(7.617E-5, method="BH", n = 40)
hist(titers_cross_short$relativeEOP, breaks = 20)  #definate not normally distributed
hist(titers_cross_short$logEOP, breaks = 20)   #looks normal
#going to report the paired t-test with BH correction. I think that this is the most
#appropriate statistical test, however, it does change the result
#from not significant (Shapiro) to significant (t-test)

#survival data from July 15
#survival data from July 15
surv <- read.csv("rawdata/phage_resistant_survival.csv")
names(surv) <- c("treat", "date","0", "1", "2", "3")
surv$dose <- str_split(surv$treat, "\\(", simplify = T)[,2]
surv$dose <- str_replace(surv$dose, "\\)", "")
surv$dose[1] <- "r1"
surv[str_detect(surv$treat, "larvaeAndPhage"),]$dose <- "r1"

surv.long <- surv %>% 
  pivot_longer(cols=c("0","1","2","3"), names_to = "day", values_to = "count") %>% #make long
  mutate(treat = str_split(treat, "\\ \\(", simplify = T)[,1]) %>%  ##fix the treatment names
  group_by(treat,date,dose)  %>% #calculate proportions
  mutate(max = max(count), freq = count / max)     #calculate proportions

big.tab <- surv.long %>% 
  mutate(dead=max-count, alive=count) %>%
  select(-count,-freq) %>%
  gather(status, count, c(alive,dead)) %>%
  slice(rep(1:n(), count)) %>%
  transform(surv=ifelse(status=="alive", 0, 1), status=NULL) %>%
  select(-count) %>%
  arrange(treat, dose, day, surv) %>%
  mutate(day=as.numeric(day), treat=as.factor(treat))

sfit <- survfit(Surv(day, surv) ~ treat, data=big.tab)
summary(sfit)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE)

coxph(Surv(day, surv) ~ treat, data=big.tab)

pairwise_survdiff(Surv(day, surv) ~ treat, data=big.tab)



##was interested to see if I could see the presence of phages in resistant isolates
##aligned reads to a multifasta file containing all the phage genomes and 3650
##processed sam files with "run_variant_calling.sh" script and then open
##sorted bams with this bit of code. It calculates average coverage of each
##genome

samples <- get_sample_names()
samples <- samples %>%
  select("ID", "NAME")
bams <- list.files(path="rawdata/phage_alignments/bam")
bams <- bams[str_detect(bams,"sorted.bam")]
bams <- bams[!str_detect(bams,"bai")]

fastas <- c("KT361649", "KT361650", "KT361652", "KT361654", "MH460824", "MH460825", "MH460826", "MH460827", "Paenibacillus_larvae_strain_3650_(hybrid)")
bamTab <- data.frame(matrix(nrow=length(bams), ncol=length(fastas)+1))
names(bamTab) <- c("sample", fastas)

for(i in 1:nrow(bamTab)){
  bamName <- paste("rawdata/phage_alignments/bam/", bams[i], sep="")
  covs <- bamcoverage(bamName)
  bamTab[i,] <- covs[names(bamTab)]
  bamTab[i,1] <- bams[i]
  print(bamName)
  print(covs)
}
names(bamTab) <- c("sample", "Fern", "Willow", "Xenia", "Vegas", "Unity", "Scottie", "Heath", "Halcyone", "B3650")
bamTab$sampleNum <- str_split(bamTab$sample, "_", simplify = T)[,2]
bamTab2 <- merge(bamTab, samples, by.x="sampleNum", by.y="ID", all=T)  
bamTab <- bamTab2 %>%
  select("sample", "NAME", "Fern", "Willow", "Xenia", "Vegas", "Unity", "Scottie", "Heath", "Halcyone", "B3650")
bamTab[1,2] <- "anc"

write.csv(bamTab, file="results/allPhageW3650_coverages.csv", quote = F, col.names = T, row.names = F)
