#title: "Nevertheless sex persisted: facultative sex is common but costly in a parthenogenetic wasp"
#author: "Rebecca Boulton"
#date: "2023-08-03"
#subtitle: Asexual females get the worst of both worlds when exposed to male conspecifics
#from a sexual population in *Lysiphlebus fabarum*
#  ---

library(survival)
library(survminer)
library(dplyr)
library(formattable)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(gtable)
library(glmmTMB)
library(Rmisc)
library(tidyr)
library(Hmisc)
library(formattable)
library(knitr)
library(lme4)
library(ggplot2)
library(ResourceSelection)
library(binom)
library(DHARMa)
library(MASS)
library(Rmisc)
library(ggbeeswarm)


G1.raw <- read.delim("~/G1-raw-check.txt")
G2.raw <- read.delim("~/G2-raw-check.txt")
G1.mating.rates <- read.delim("~/G1.mating.rates.txt")

# put asexual line names as full names
G1.raw <- G1.raw |> 
  mutate(Line = case_when(Line == "84" ~ "CV17-84", 
                          Line == "658" ~ "IL06-658",
                          Line == "348" ~ "IL09-348",
                          Line == "64" ~ "IL07-64",
                          Line == "66" ~ "CHC17-66",
                          Line == "554" ~ "IL09-554",
                          Line == "402" ~ "IL09-402",
                          Line == "Sexual" ~ "Sexual"))

G2.raw <- G2.raw |> 
  mutate(Line = case_when(Line == "84" ~ "CV17-84", 
                          Line == "658" ~ "IL06-658",
                          Line == "348" ~ "IL09-348",
                          Line == "64" ~ "IL07-64",
                          Line == "66" ~ "CHC17-66",
                          Line == "554" ~ "IL09-554",
                          Line == "402" ~ "IL09-402",
                          Line == "Sexual" ~ "Sexual"))




G1 <- subset(G1.raw, Exclude..lost.or.killed.mother. == 0)
G1$gen <- as.factor("G1")
G2 <- subset(G2.raw, Exclude..lost.or.killed.mother. == 0)
G2 <- subset(G2.raw, Exclude..no.G2. == 0)
G2$gen <- as.factor("G2")

G1 <- G1 %>% mutate(Sexual.Asexual = ifelse(Line == "Sexual", "Sexual", "Asexual"))
G2 <- G2 %>% mutate(Sexual.Asexual = ifelse(Line == "Sexual", "Sexual", "Asexual"))

G1 <- G1 %>% mutate(survival.cens = ifelse(survival..hours. >=72, 0, 1))
G2 <- G2 %>% mutate(survival.cens = ifelse(Survival..hours. >=72, 0, 1))

G1$notwasps <- G1$Mummies-G1$Wasps
G2$notwasps <- G2$Mummies-G2$Wasps

G1.Sex <- subset(G1, Line == "Sexual")
G1.Sex <- G1.Sex %>% 
  filter(!(Mated.Virgin == "Virgin" & Females >= 1))
G1.ASex <- subset(G1, Line != "Sexual")

G2.Sex <- subset(G2, Line == "Sexual")
G2.Sex <- G2.Sex %>% 
  filter(!(Mated.Virgin == "Virgin" & Females >= 1))
G2.ASex <- subset(G2, Line != "Sexual")


#subset mated asexual females

G1.ASex.Mated <- subset(G1.ASex, Mated.Virgin == "Mated")
G2.ASex.Mated <- subset(G2.ASex, Mated.Virgin == "Mated")

G1.suc <- subset(G1, Failed.to.produce.offspring == 0)
G1.suc <- G1.suc %>% 
  filter(!(Mated.Virgin == "Virgin" & Sexual.Asexual == "Sexual" & Females >= 1))
G2.suc <- subset(G2, Failed.to.produce.offspring == 0)

G1.suc.ASex <- subset(G1.suc, Line != "Sexual")
G2.suc.ASex <- subset(G2.suc, Line != "Sexual")

G1.mating.rates <- G1.mating.rates[!is.na(G1.mating.rates$No.attempt),]


New_nums <- c(138:147)
G1[119:128, "Number"] <- New_nums
G1$ID <- as.numeric(G1$Number)
G1$Block.fac <- as.factor(G1$Block)

G2[4, "Number.G1"] <- 138
G2[5, "Number.G1"] <- 142
New_nums <- c(138,139,141,143,148,145,146,147,148,149)
G2[96:105, "Number.G1"] <- New_nums
G2$ID <- as.numeric(G2$Number.G1)
G2 <- G2[-c(46), ]
G2$Block.fac.1 <- as.factor(G2$Block.1)
G2$Block.fac.2 <- as.factor(G2$Block.2)

#remove virgin sexual females that produced daughters
G1 <- G1 %>% filter(!(Mated.Virgin == "Virgin" & Sexual.Asexual == "Sexual" & Females >= 1))
G2 <- G2 %>% filter(!(Mated.Virgin == "Virgin" & Sexual.Asexual == "Sexual" & Females >= 1))

#remove all virgin sexual females
G1_no_VS <- G1 %>% filter(!(Mated.Virgin == "Virgin" & Sexual.Asexual == "Sexual"))
G2_no_VS <- G2 %>% filter(!(Mated.Virgin == "Virgin" & Sexual.Asexual == "Sexual"))
G1_sex_no_VS <- G1.Sex %>% filter(!(Mated.Virgin == "Virgin" & Sexual.Asexual == "Sexual"))


# get sample sizes for each generation
# total number of females for G1 by line and mated status

total.G1  <- G1 %>%
  dplyr::summarise( 
    n())

# total that produced offspring in G1

total.G1.suc  <- G1.suc %>%
  dplyr::summarise( 
    n())

#### Mating rates for G1 wasps 

Mating.summaries <- G1.mating.rates %>%
  dplyr::summarise(  #need to state dplyr is plyr is attached or this doesn't work
    Success = sum(Mated),
    Virgin = sum(Virgin),
    Failed.to.find = sum(No.attempt),
    Refused = sum(Refuse),
    Total = sum(Mated + No.attempt + Refuse),
    proportion = Success / Total)

G1.mating.rates <- G1.mating.rates %>%
  mutate(Line = case_when(Line == "84" ~ "CV17-84", 
                          Line == "658" ~ "IL06-658",
                          Line == "348" ~ "IL09-348",
                          Line == "64" ~ "IL07-64",
                          Line == "66" ~ "CHC17-66",
                          Line == "554" ~ "IL09-554",
                          Line == "402" ~ "IL09-402",
                          Line == "Sexual" ~ "Sexual"))


Mating.summaries.line.mating <- G1.mating.rates %>%
  group_by(Line) %>%
  dplyr::summarise(  #need to state dplyr is plyr is attached or this doesn't work
    "Remained virgin" = sum(Virgin),
    "Attempted matings" =sum(Mated + No.attempt + Refuse),
    "Successful matings" = sum(Mated),
    "No attempt" = sum(No.attempt),
    "Rejected" = sum(Refuse),
    "proportion attempts successful" = round(sum(Mated)/sum(Mated+No.attempt+Refuse), digits = 2))

table1 <- formattable(Mating.summaries.line.mating, align=c("l", "c", "c", "c", "c", "c", "r"),
                      list(`Line` = formatter(
                        "span", style = ~ style(color = "grey",font.weight = "bold")), 
                        `proportion attempts successful` = color_bar(color = "lightgreen")) 
)


#mated females that used sperm

G1.ASex.Mated <- G1.ASex.Mated[!is.na(G1.ASex.Mated$SA),]

Sperm.summaries <- G1.ASex.Mated %>%
  dplyr::summarise(  #need to state dplyr is plyr is attached or this doesn't work
    "Used sperm" = sum(SA),
    "Total" = n(),
    "Did not use sperm" = sum(n()-sum(SA)))

Sperm.summaries.Line <- G1.ASex.Mated %>%
  group_by(Line) %>%
  dplyr::summarise(  #need to state dplyr is plyr is attached or this doesn't work
    "Sexual allele present" = sum(SA),
    "Sexual allele absent" = sum(n()-sum(SA)),
    "Proportion used sperm" = round(sum(sum(SA)/n()), digits = 2))



table3 <- formattable(Sperm.summaries.Line, align=c("l", "c", "c", "r"),
                      list(`Line` = formatter(
                        "span", style = ~ style(color = "grey",font.weight = "bold")), 
                        `Proportion used sperm` = color_bar(color = "lightgreen")) 
)

# total that produced offspring in G1

G2.suc  <- G2.suc %>%
  dplyr::summarise( 
    n())


#parasitised but failed to produce adult offspring

G1.offs.fail <- subset(G1, Failed.to.produce.offspring == 1)

G1.para.but.offs.fail <- subset(G1.offs.fail, Failed.to.parasitise == 0)

G1.para.but.offs.fail.n <- G1.para.but.offs.fail %>%
  dplyr::summarise( 
    n())

#survival analyses G1

G1 <- G1 %>% mutate(survival.cens = ifelse(survival..hours. >=72, 0, 1))

survived <- subset(G1, survival.cens == 0)
ded <- subset(G1, survival.cens == 1)

n.survivors <- survived %>% 
  dplyr::summarise(
    n())

n.ded <- ded %>%
  dplyr::summarise(
    n())

surv_G1 <- Surv(time = G1$survival..hours., event = G1$survival.cens)
G1.surv <- survfit(surv_G1 ~ Sexual.Asexual + Mated.Virgin, data = G1)
leg_labs.sex.mated <- list("Asexual mated", "Asexual virgin", "Sexual mated", "Sexual virgin")
G1.surv.plot <- ggsurvplot(G1.surv, data = G1, conf.int = TRUE, 
                           pval = TRUE, legend.title="", legend.labs = leg_labs.sex.mated, title="", surv.plot.height = 1) + theme_survminer(base_size = 0.5) + xlab(" ") + ylab(" ")

G1.surv.plot <- G1.surv.plot$plot + facet_grid(~Mated.Virgin)

#log rank
survdiff(Surv(survival..hours.,, survival.cens) ~ Sexual.Asexual + Mated.Virgin, data=G1)

```{r, include=FALSE, message=FALSE, warning = FALSE}


#test overall effect of mating or not on both sexual and asexual wasps

G1.Para.fail <- glmer(Failed.to.parasitise ~ Sexual.Asexual*Mated.Virgin + N.nymphs + survival..hours. + (1|Sexual.Asexual/Line), family = binomial(link = "logit"), data = G1)
simulationOutputG1parafail <- simulateResiduals(fittedModel = G1.Para.fail, plot = T, use.u = T)
testZeroInflation(simulationOutputG1parafail)
para.fail.F <- car::Anova(G1.Para.fail, test = "Chisq")
summary(G1.Para.fail)


result1 <- G1 %>%
  group_by(Sexual.Asexual, Mated.Virgin) %>%
  dplyr::summarise(  #need to state dplyr is plyr is attached or this doesn't work
    Count_fail = sum(Failed.to.parasitise),
    Total = length(Failed.to.parasitise),
    Count_suc = Total-Count_fail,
    Percent_fail = Count_fail/Total*100,
    Percent_suc = Count_suc/Total*100
  )

result1.long <- result1 %>%
  pivot_longer(cols = starts_with("Percent_"), names_to = "Suc.fail", values_to = "Percent")

long.leg.labs <- list("Failed", "Succeeded")

G1.fail.plot <- ggplot(result1.long, aes(x = Sexual.Asexual, y = Percent, fill=factor(Suc.fail))) +
  geom_bar(stat = "identity") + labs(x = "Mated status", y = "% successful or failed parasitism attempts", fill = "") +
  scale_fill_manual(values=c('#fc8d59', '#ffffbf'), labels=c('Failed', 'Succeeded')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_grid(cols = vars(Mated.Virgin)) +
  theme(strip.text = element_text(face="bold.italic", size = "14"), axis.title=element_text(face="bold",  size="14", color="black"),
 axis.text = element_text(size = "11")) + theme(legend.text = element_text(size="11"))


# subset successful mothers out
G1.suc <- subset(G1, Failed.to.produce.offspring == 0)
G1.ASex.suc <- subset(G1.ASex, Failed.to.produce.offspring == 0)

G1.suc.tot <- G1.suc %>% 
  dplyr::summarise(
    n())

G1.ASex.suc.tot <- G1.ASex.suc %>% 
  dplyr::summarise(
    n())


#test overall effect of mating or not on both sexual and asexual wasps

G1.Mummies <- glmmTMB(Mummies~ Sexual.Asexual*Mated.Virgin + N.nymphs + survival..hours. + (1|Sexual.Asexual/Line), family = nbinom1, data=G1.suc)
simulationOutputG1Mummies <- simulateResiduals(fittedModel = G1.Mummies, plot = T, use.u = T)
G1.mummies.F <- car::Anova(G1.Mummies)
summary(G1.Mummies)

dodge <- position_dodge(width = 1)


G1.mummies.plot <- ggplot(G1.suc, aes(x = Mated.Virgin, y=Mummies, fill=Mated.Virgin)) + labs(x = "", y = "# Mummies", fill = "")
+ geom_beeswarm(aes(colour = Mated.Virgin, fill = Mated.Virgin), shape = 21, size = 7, colour = "black", alpha = 0.7, cex = 5) + 
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange", colour = "black", size = 1) + scale_fill_manual(values=c('#af8dc3', '#7fbf7b'),
 labels=c('Mated', 'Virgin')) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(strip.text = element_text(face="bold.italic", size = "18"), 
axis.title=element_text(face="bold",  size="18", color="black"), axis.text = element_text(size = "14")) + theme(legend.position="none") + facet_grid(cols = vars(Sexual.Asexual)) 


G1.Wasps.bin <- glmmTMB(cbind(Wasps,Unemerged)~ Sexual.Asexual*Mated.Virgin + N.nymphs + survival..hours. + (1|Line) + (1|Number), family = binomial(link = "logit"), data=G1.suc)
simulationOutputG1wasps <- simulateResiduals(fittedModel = G1.Wasps.bin, plot = T, use.u = T)
G1.wasps.F <- car::Anova(G1.Wasps.bin)
summary(G1.Wasps.bin)

G1.wasps.plot <- ggplot(G1.suc, aes(x = Mated.Virgin, y=(Wasps/Mummies * 100), fill=Mated.Virgin)) + labs(x = "", y = " % Adult offspring that emerged from mummies",fill = "") 
+ geom_beeswarm(aes(colour = Mated.Virgin, fill = Mated.Virgin), shape = 21, size = 7, colour = "black", alpha = 0.7, cex = 5) + 
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange", colour = "black", size = 1) + scale_fill_manual(values=c('#af8dc3', '#7fbf7b'), labels=c('Mated', 'Virgin')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(strip.text = element_text(face="bold.italic", size = "18"), axis.title=element_text(face="bold",  size="18", color="black"),
  axis.text = element_text(size = "16")) + theme(legend.position="none") + stat_summary(fun.data = "mean_cl_boot", geom = "pointrange", colour = "black") + facet_grid(cols = vars(Sexual.Asexual))

G1_no_VS$Sex_mating <- paste(G1_no_VS$Sexual.Asexual, G1_no_VS$Mated.Virgin)

G1.daughters <- glmmTMB(Females ~ Sexual.Asexual*Males + N.nymphs + survival..hours. + (1|Sexual.Asexual/Line), family = nbinom1, data=G1_no_VS)
simulationOutputG1Females <- simulateResiduals(fittedModel = G1.daughters, plot = T, use.u = T)
G1.daughters.F <- car::Anova(G1.daughters)
summary(G1.daughters)


G1.daughters.ASex <- glmmTMB(Females ~ Mated.Virgin*Males + N.nymphs + survival..hours. + (1|Line), family = nbinom1, data=G1.ASex.suc)
simulationOutputG1Females <- simulateResiduals(fittedModel = G1.daughters, plot = T, use.u = T)
G1.daughters.ASex.F <- car::Anova(G1.daughters.ASex)
summary(G1.daughters.ASex)

##G2

total.G2  <- G2.ASex %>%
  dplyr::summarise( 
    n())

Total.mated.virgin.G2 <-  G2.ASex %>%
  group_by(Mated.Virgin) %>%
  dplyr::summarise( 
    n())

total.G2.suc  <- G2.suc.ASex %>%
  dplyr::summarise( 
    n())


Total.mated.virgin.G2.suc <-  G2.suc.ASex %>%
  group_by(Mated.Virgin) %>%
  dplyr::summarise( 
    n())


G2.ASex$"Mated or Virgin mother*" <- G2.ASex$Mated.Virgin

total.Mated.line.G2 <- G2.ASex %>%
  group_by(Line, `Mated or Virgin mother*`) %>%
  dplyr::summarise( 
    "N" = n())


table2 <- formattable(total.Mated.line.G2, align=c("l", "c", "r"),
                      list(`Line` = formatter(
                        "span", style = ~ style(color = "grey",font.weight = "bold"))
                      ))


G2.offs.fail <- subset(G2.ASex, Failed.to.produce.offspring == 1)

G2.para.but.offs.fail <- subset(G2.offs.fail, Failed.to.parasitise == 0)

G2.para.but.offs.fail.n <- G2.para.but.offs.fail %>%
  dplyr::summarise( 
    n())

G2.survived <- subset(G2.ASex, survival.cens == 0)
G2.ded <- subset(G2.ASex, survival.cens == 1)

G2.n.survivors <- G2.survived %>% 
  dplyr::summarise(
    n())

G2.n.ded <- G2.ded %>%
  dplyr::summarise(
    n())

surv_G2 <- Surv(time = G2.ASex$Survival..hours., event = G2.ASex$survival.cens)
G2.surv <- survfit(surv_G2 ~ Mated.Virgin, data = G2.ASex)
leg_labs.G2 <- list("Mated", "Virgin")
G2.surv.plot <- ggsurvplot(G2.surv, data = G2.ASex, conf.int = TRUE, 
                           pval = TRUE, legend.title="", legend.labs = leg_labs.G2, title="", surv.plot.height = 1) + theme_survminer(base_size = 0.5) + xlab(" ") + ylab(" ")

#logrank test
survdiff(Surv(Survival..hours.,, survival.cens) ~ Mated.Virgin, data=G2.ASex)

G2.ASex <- G2.ASex[!is.na(G2.ASex$Failed.to.parasitise),]

G2.suc.tot <- G2.suc.ASex %>% 
  dplyr::summarise(
    n())


G2.Para.fail <- glmer(Failed.to.parasitise ~ Mated.Virgin + N.nymphs + Survival..hours. + as.character(SA) + (1|Line), family = binomial(link = "logit"), data = G2.ASex)
simulationOutputG2parafail <- simulateResiduals(fittedModel = G2.Para.fail, plot = T, use.u = T)
testZeroInflation(simulationOutputG2parafail)
G2.para.fail.F <- car::Anova(G2.Para.fail)
summary(G2.Para.fail)


result2 <- G2.ASex %>%
  group_by(Mated.Virgin) %>%
  dplyr::summarise(  #need to state dplyr is plyr is attached or this doesn't work
    Count_fail = sum(Failed.to.parasitise),
    Total = length(Failed.to.parasitise),
    Count_suc = Total-Count_fail,
    Percent_fail = Count_fail/Total*100,
    Percent_suc = Count_suc/Total*100
  )

result2.long <- result2 %>%
  pivot_longer(cols = starts_with("Percent_"), names_to = "Suc.fail", values_to = "Percent")

long.leg.labs <- list("Failed", "Succeeded")


G2.fail.plot <- ggplot(result2.long, aes(x = Mated.Virgin, y = Percent, fill=factor(Suc.fail))) + 
  geom_bar(stat = "identity") + labs(x = "Mated status", y = "% successful or failed parasitism attempts", fill = "") + 
  scale_fill_manual(values=c('#fc8d59', '#ffffbf'), labels=c('Failed', 'Succeeded')) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(strip.text = element_text(face="bold.italic", size = "14"), axis.title=element_text(face="bold",  size="14", color="black"),
     axis.text = element_text(size = "11")) + theme(legend.text = element_text(size="11"))

G2.ASex.Mated <- G2.ASex.Mated[!is.na(G2.ASex.Mated$SA),]


result3 <- G2.ASex.Mated %>%
  group_by(SA) %>%
  dplyr::summarise(  #need to state dplyr is plyr is attached or this doesn't work
    Count_fail = sum(Failed.to.parasitise),
    Total = length(Failed.to.parasitise),
    Count_suc = Total-Count_fail,
    Percent_fail = Count_fail/Total*100,
    Percent_suc = Count_suc/Total*100
  )

result3.long <- result3 %>%
  pivot_longer(cols = starts_with("Percent_"), names_to = "Suc.fail", values_to = "Percent")

long.leg.labs <- list("Failed", "Succeeded")


G2.fail.plot.mated <- ggplot(result3.long, aes(x = SA, y = Percent, fill=factor(Suc.fail))) + 
  geom_bar(stat = "identity") + labs(x = "Sexual allele present/absent", y = "% successful or failed parasitism attempts", fill = "") + scale_fill_discrete(labels=c('Failed', 'Succeeded'))  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


SA.present <- subset(G2.ASex.Mated, SA == 1)
SA.absent <- subset(G2.ASex.Mated, SA == 0)

SAP <- as.vector(SA.present$Failed.to.parasitise)
SAA <- as.vector(SA.absent$Failed.to.parasitise)

m <- cbind(table(SAP), table(SAA))
chisq.test(m)

G2.Mummies <- glmmTMB(Mummies~ Mated.Virgin + N.nymphs + Survival..hours. + (1|Line), family = nbinom1, data=G2.suc.ASex)
simulationOutputG2Mummies <- simulateResiduals(fittedModel = G2.Mummies, plot = T, use.u = T)
G2.mummies.F <- car::Anova(G2.Mummies)
summary(G2.Mummies)

G2.mummies.plot <- ggplot(G2.suc.ASex, aes(x = Mated.Virgin, y=Mummies, fill=Mated.Virgin)) + 
  labs(x = "", y = "# Mummies", fill = "") + geom_beeswarm(aes(colour = Mated.Virgin, fill = Mated.Virgin), 
shape = 21, size = 7, colour = "black", alpha = 0.7, cex = 5) + stat_summary(fun.data = "mean_cl_boot", geom = "pointrange", colour = "black", size = 1) + 
  scale_fill_manual(values=c('#af8dc3', '#7fbf7b'), labels=c('Mated', 'Virgin')) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(strip.text = element_text(face="bold.italic", size = "18"), axis.title=element_text(face="bold",  size="18", color="black"), axis.text = element_text(size = "14")) + theme(legend.position="none") 



G2.suc.ASex$Wasps<- as.numeric(G2.suc.ASex$Wasps)
G2.suc.ASex$Unemerged <- as.numeric(G2.suc.ASex$Unemerged)

G2.Wasps.bin <- glmmTMB(cbind(Wasps,Unemerged)~ Mated.Virgin + N.nymphs + Survival..hours. + (1|Line) + (1|Number.G1), family = binomial(link = "logit"), data=G2.suc.ASex)
simulationOutputG2wasps <- simulateResiduals(fittedModel = G2.Wasps.bin, plot = T, use.u = T)
G2.wasps.F <- car::Anova(G2.Wasps.bin)
summary(G2.Wasps.bin)


G2.wasps.plot <- ggplot(G2.suc.ASex, aes(x = Mated.Virgin, y=(Wasps/Mummies * 100), fill=Mated.Virgin)) + 
  labs(x = "", y = "% Adult offspring that emerged from mummies", fill = "") + 
  geom_beeswarm(aes(colour = Mated.Virgin, fill = Mated.Virgin), shape = 21, size = 7, colour = "black", alpha = 0.7, cex = 5) +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange", colour = "black", size = 1) + scale_fill_manual(values=c('#af8dc3', '#7fbf7b'), labels=c('Mated', 'Virgin')) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(strip.text = element_text(face="bold.italic", size = "18"), 
axis.title=element_text(face="bold",  size="18", color="black"), axis.text = element_text(size = "14")) + theme(legend.position="none")


G2.daughters <- glmmTMB(Females~ Mated.Virgin*Males + N.nymphs + Survival..hours. + (1|Line), family = nbinom1, data=G2.suc.ASex)
simulationOutputG2daughters <- simulateResiduals(fittedModel = G2.daughters, plot = T, use.u = T)
G2.daughters.F <- car::Anova(G2.daughters)
summary(G2.daughters)

G2.sperm <- glmmTMB(SA~ Mated.Virgin, family = binomial, data=G2.suc.ASex)
simulationOutputG2sperm <- simulateResiduals(fittedModel = G2.sperm, plot = T, use.u = T)
car::Anova(G2.sperm)
summary(G2.sperm)
