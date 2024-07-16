#make figures
library(stringr)
library(tidyverse)
library(readxl)
source('R_functions.R')
library(ggsci)
library(ggpubr)


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

surv.sum <- surv.long %>%
  group_by(treat,day) %>%
  summarise(mean=mean(freq),
            sd=sd(freq),
            se=sd/sqrt(n()))

p <- ggplot(surv.sum, aes(x=day, y=mean, color=treat, group=treat, shape=treat)) +
  geom_point(size=2) +
  geom_line() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  theme_classic() +
  scale_color_npg() +
  scale_shape_manual(values = c(16, 15, 17, 18)) +
  labs(y="Survival (%)", x="Days post-exposure")
p

#surv.long %>%
#  filter(!treat=="Control") %>%
#  mutate(treat=as.factor(treat)) %>%
#  mutate(day=as.factor(day)) %>%
#  compare_means(count ~ treat, data=., group.by = "day")
ggsave("figs/survive.png", plot = p, width=4, height=2.5, units="in")
ggsave("figs/survive.pdf", plot = p, width=4, height=2.5, units="in")





#function to import data from excel
grow.all <- get_data()

grow.sum <- grow.all %>%
  group_by(Time, isolate) %>%
  summarize(mean=mean(OD600), sd=sd(OD600), n=n()) %>%
  filter(!isolate=="Y3650 Anc")

grow.sum$ref <- str_detect(grow.sum$isolate, "Y3650")
grow.sum$phage <- str_split(grow.sum$isolate, "_", simplify = T)[,1]
tmp <- grow.sum[grow.sum$isolate=="Y3650",]
grow.sum <- grow.sum %>%
  filter(!isolate=="Y3650")
tmp <- rbind(
  cbind(tmp[,-7],"F"),
  cbind(tmp[,-7],"Hal"),
  cbind(tmp[,-7],"Scot"),
  cbind(tmp[,-7],"Unity"),
  cbind(tmp[,-7],"V"),
  cbind(tmp[,-7],"XIII")
)
names(tmp) <- names(grow.sum)
grow.sum <- rbind(grow.sum, tmp)
grow.sum$isolate <- as.factor(grow.sum$isolate)
levels(grow.sum$isolate)
#####drop Unity_x - no mutations in it
#grow.sum <- grow.sum %>%
#  filter(!isolate=="Unity_x")
cols <- c("#1b9e77",
          "#1b9e77",
          "#1b9e77",
          "#d95f02",
          "#d95f02",
          "#d95f02",
          "#7570b3",
          "#7570b3",
          "#7570b3",
          "#e7298a",
          "#e7298a",
          "#e7298a",
          "#e7298a",
          "#66a61e",
          "#66a61e",
          "#66a61e",
          "#e6ab02",
          "#e6ab02",
          "black"
)
p <- grow.sum %>%
  filter(Time <= 48) %>%
ggplot(aes(x=Time, y=mean, ymin=mean-sd, ymax=mean+sd, fill=isolate)) +
  geom_line(aes(linewidth=ref)) +
  geom_ribbon(alpha=0.25) +
  scale_linewidth_manual(values=c(0.5,2)) +
  facet_wrap(~phage, ncol = 2, scales="free_y") +
  scale_fill_manual(values=cols) +
  theme_light() +
  theme(legend.position="none") +
  labs(x="Time (hr)", y="Cell growth (OD600)")
p

ggsave("figs/cellgrowth.png", p, units = "in", width = 4, height=5.5)


p2 <- grow.sum %>%
  filter(Time <= 40) %>%
  ggplot(aes(x=Time, y=mean, ymin=mean-sd, ymax=mean+sd, fill=isolate)) +
  geom_line(aes(linewidth=ref)) +
  geom_ribbon(alpha=0.25) +
  scale_linewidth_manual(values=c(0.5,2)) +
  facet_wrap(~phage, ncol = 3, scales="free_y") +
  scale_fill_manual(values=cols) +
  theme_light() +
  theme(legend.position="none") +
  labs(x="Time (hr)", y="Cell growth (OD600)")
p2

ggsave("figs/cellgrowth_horz.png", p2, units = "in", width = 7.5, height=4.75)
ggsave("figs/cellgrowth_horz.pdf", p2, units = "in", width = 7.5, height=4.75)


###EOP plot
titers <- read_excel("rawdata/Picked plaques Titers Comparison.xlsx", sheet = 3) ##must copy to working excel file to Project directory.
titers[titers==0] <- 1 #pseudocount to avoid infinity

titers$eop <- titers$resistantHost / titers$titer3650

eop.l <- titers %>%
  drop_na() %>%
  pivot_longer(cols=c("titer3650", "resistantHost"), values_to = "titer")

p <- ggplot(eop.l, aes(x=name, log(titer, base = 10), group=host)) +
  geom_point() +
  geom_line() +
  facet_wrap(~evolvedPlaque) +
  scale_x_discrete(limits=c("titer3650", "resistantHost")) +
  labs(x="", y="Titer (log10 pfu/mL)") +
  theme_light()
p
ggsave("figs/eop.png", plot = p, width=3.5, height=3.5, units="in")
ggsave("figs/eop.pdf", plot = p, width=3.5, height=3.5, units="in")


##############################################
#second EOP plot
#tries to show relative plating efficiency of survivor phages
#on alternative hosts
#In addition to testing plaques from the host that was challenged with a phage,
#Emma measured the titer of these survivor phages on other permissive hosts.
#example: Phage F used to evolve resistance in 3650 (F-r-3650). A survivor phage picked.
#survivor phage stock titer tested on 3650 and F resistant 3650. Change in EOP from
#initial phage stock to survivor phage stock measured (see above figure)
#F is also able to infect Wa resistant hosts, so change in EOP measured
#for initial F phage stock and F survivor phages, tested on Wa resistant 3650
#cannot calculate this change if there is missing data (no plaques). 
#this is a hard number to figure out what it is. It is a ratio or ratios that 
#compares the change in plating efficiency of challenge and alt phages. It attempts
#to answer the question: Do other phages that can grow on a host experience similar 
#improvements after a round of selection as the challenge phage?
##############################################

relEOP <- read_excel("rawdata/Picked plaques Titers Comparison.xlsx", sheet = 6) ##must copy to working excel file to Project directory.
reEOP_long <- relEOP %>%
  pivot_longer(!c("challenge_phage", "alt_phage"), names_to = "time", values_to = "relEOP")

pd <- position_dodge(0.1)
p <- ggplot(reEOP_long, aes(x=time, y=log(relEOP, base = 10), group=interaction(challenge_phage, alt_phage))) +
  geom_point(position=pd, alpha=0.75) +
  geom_line(position=pd, alpha=0.75) +
  scale_x_discrete(limits=c("start", "end")) +
  labs(x="", y=expression(Delta~EOP)) +
  theme_light()
p
ggsave("figs/DeltaEOP.png", plot = p, width=2, height=3.5, units="in")
ggsave("figs/DeltaEOP.pdf", plot = p, width=2, height=3.5, units="in")


##############################################
#Matrix plot showing cross resistance - use up-to-date excel file
##############################################
{
cross <- read_excel("rawdata/P. larvae Resistant Hosts - Cross Resistance Tables.xlsx", sheet=3, skip=3)  
cross <- cross %>%
  select("...2","F", "Wa", "XIII", "Heath", "Halcyone", "Scottie", "Unity", "V...10")

names(cross) <- c("strain","F", "Wa", "XIII", "Heath", "Halcyone", "Scottie", "Unity", "V")
cross.l <- cross %>%
  pivot_longer(!strain, names_to = "phage", values_to = "titer")

lookup <- data.frame(ref=c("Complete Clearing","TMTC","Cleared", "Almost completely", "Almost cleared", "Almost completely cleared"),
                     data=c("1E7", "1E7", "1E7", "1E7", "1E7", "1E7"))
row.names(lookup) <- lookup$ref

for(i in 1:nrow(cross.l)){
  if(cross.l$titer[i] %in% lookup$ref){
    cross.l$titer[i] <- lookup[cross.l$titer[i],]$data
  }
  if(!is.na(cross.l$titer[i])){
    if(str_detect(cross.l$titer[i], "[F-Zf-z]")){
      print("Found a comment that didn't get fixed")
      exit()
    }
  }
}
cross.l$titer <- as.numeric(cross.l$titer)
cross.l$phage <- factor(cross.l$phage)
cross.l$strain <- factor(cross.l$strain)

cols <- c("#1b9e77",   
          "#d95f02",   
          "#7570b3",   
          "#e7298a",   
          "#66a61e",   
          "#e6ab02",   
          "black",
          "grey"
)


p <- ggplot(cross.l, aes(x=phage, y=rev(strain), alpha=titer, color=phage)) +
  geom_point(size=3.25) +
  scale_color_manual(values=cols) +
  theme_classic() +
  labs(x="", y="") +
  theme(legend.position = "none")+
  scale_x_discrete(limits=c("F", "Wa", "XIII", "Heath",  "Scottie", "Unity", "V", "Halcyone") ,position = "top") +
  scale_y_discrete(limits = rev(levels(cross.l$strain)))
p
ggsave(p, file="figs/cross_resistance.png", width = 4, height=6, units="in")
ggsave(p, file="figs/cross_resistance.pdf", width = 4, height=6, units="in")
}
