####"State Project 2021: *A. cervicornis* & *A. palmata* Long-term vs. CBASS"####
#libraries
library(car)
library(corrplot)
library(lmerTest)
library(ggcorrplot)
library(ggproto)
library(emmeans)
library(ggplot2)
library(ggpubr) 
library(Hmisc)
library(dplyr)
library(FSA)
library(FactoMineR)
library(RColorBrewer)
library(PerformanceAnalytics)
library(factoextra)
library(MASS)
library(reshape2)
library(cowplot)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis",force=TRUE)
library(pairwiseAdonis)
BiocManager::install('PCAtools')
library(PCAtools)
library(vegan)
library(ggdendro)
library(tidyr)

###Import data and set factors and levels###
# data set
trait <- read.csv("20221114_2021State_MasterSheet.csv", header = T, na.strings=c("","NA"))
str(trait)
summary(trait)

#set factor levels
trait<-transform(trait,
                Genotype=as.factor(Genotype),
                Treatment=as.factor(Treatment),
                Tank=as.factor(Tank),
                Experiment=as.factor(Experiment),
                Exp=as.factor(Exp),
                Region=as.factor(Region),
                Rep=as.factor(Rep))

head(trait)
levels(trait$Genotype)<-c("LK31","LK41","LK50","LK62","LK7", "UK12", "UK19", "UK70","UK76","UK80")
trait$Exp=factor(trait$Exp,levels(trait$Exp)[c(9,1,2,3,4,5,6,7,8,10)])
trait$Treatment=factor(trait$Treatment,levels(trait$Treatment)[c(9,1,2,3,4,5,6,7,8,10)])

#subset new df based on experiment
trait<-droplevels(trait[!trait$Exp == 'CBASS-37' ,])
trait$Treatment <- factor(trait$Treatment) 
LT<-trait[trait$Experiment=='LT',]
CBASS<-trait[trait$Experiment=='CBASS',]

#for working with heat samples only
trait2 <- droplevels(trait[!trait$Treatment == 'Control' ,])
trait2$Treatment <- factor(trait2$Treatment) 
trait3 <- droplevels(trait2[!trait2$Treatment == '27.5' ,])
trait3$Treatment <- factor(trait3$Treatment) 

#####NMDS of LT, CBASS#####
#select columns of interest
traitdata<-dplyr::select(trait, Exp, Genotype, FvFm, HP, TotalChl,ED) 
traitdata<-na.omit(traitdata) #rm NAs

## Create a community data frame and a predictor variable data frame 
trait_com<-dplyr::select(traitdata,FvFm, HP, TotalChl)
#standardize by centering and scaling
trait_com<-trait_com %>% mutate_if(is.numeric, scale)
trait_meta<-dplyr::select(traitdata,Exp, Genotype)

#This produce the coordinates for each point
# Calculate distance matrix non-rarefied
distmat <- vegdist(trait_com, method = "euclidean",na.rm=TRUE)
#Run PERMANOVA
adonis2(distmat~Exp*Genotype, data=trait_meta)
#getOption('max.print')
pairwise.adonis(trait_com,factors=trait_meta$Exp,sim.function='euclidean',p.adjust.m='bonferroni')

#can't get a pairwise 35.6 and LT-OW, trying a workaround
trait4<-trait[trait$Exp=='LT-OW'| trait$Exp=='CBASS-35.6',]
trait_com2<-dplyr::select(trait4,FvFm, HP, TotalChl)
#standardize by centering and scaling
trait_com2<-trait_com2 %>% mutate_if(is.numeric, scale)
trait_meta2<-dplyr::select(trait4,Exp, Genotype)
distmat <- vegdist(trait_com2, method = "euclidean",na.rm=TRUE)
#Run PERMANOVA
adonis2(distmat~Exp, data=trait_meta2)
#read in pvals to bonferroni adjust
pvals<-read.delim("2021_PERMANOVApvals.txt", header=T)
p.adjust(pvals$p.value, method="bonferroni")

#This produce the coordinates for each point
#trait_com<-as.matrix(trait_com)
ord_IR<-metaMDS(distmat)
ord_IR$stress
stressplot(ord_IR)

data.scores = as.data.frame(scores(ord_IR, display="site"))
data.scores <- cbind(data.scores,trait_meta)

fit<-envfit(ord_IR, trait_com, perm=999,na.rm=TRUE)
fit
scores(fit, "vectors")
plot(ord_IR,axes=FALSE, ylim=c(-.25,0.35), xlim=c(-0.2,0.2))
plot(fit, p.max = 0.95, col = "black") + abline(h = 0, v = 0, col = "black", lty = 2)

vec.df <- as.data.frame(fit$vectors$arrows*sqrt(fit$vectors$r))
vec.df$variables <- rownames(vec.df)

pNMDS<-ggplot(data.scores, aes(x = NMDS1, y = NMDS2, col=Exp)) + 
       scale_colour_brewer(palette = 'RdYlBu',direction=-1) +
       stat_ellipse(geom = "polygon", alpha = 0.25, size= 0.5, level = 0.95, aes(fill=Exp)) + 
       scale_fill_brewer(palette = 'RdYlBu',direction=-1) +
       geom_point(size = 2, aes(shape=Exp)) + 
       scale_shape_manual(values=c(17,19,19,19,19,19,19,19,17)) +
 geom_segment(data = vec.df,
              aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
              arrow = arrow(length = unit(.3, "cm")),
              colour="black", size=1,
              inherit.aes = FALSE) + 
       labs("Experiment")+
       coord_fixed() +                                              # same axis scaling
       theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
       geom_vline(xintercept=0, linetype="dashed") + geom_hline(yintercept=0, linetype="dashed") +
       theme(text = element_text(size=15, face="bold")) + 
       theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
       theme(axis.ticks.length=unit(.2,"cm")) 

ggsave(pNMDS, height = 6 , width = 8, filename = "./20221118_NMDS_byExperiment-sites.pdf", useDingbats=FALSE)
          
####Correlation matrix of traits across similar treatments####
#CBASS 34.3
CBASS34.3<-subset(traitdata, Exp=="CBASS-34.3"  ) 
CBASS34.3<-dplyr::select(CBASS34.3, FvFm, HP, TotalChl,ED) 
C<-cor(CBASS34.3,use = "complete.obs")
p.mat<-cor_pmat(CBASS34.3)
CBASS34.3plot<-ggcorrplot(C,hc.order=TRUE,outline.col="white",p.mat=p.mat,type='upper',lab=TRUE)
ggsave("./20230321_corrmatrix-34.3.pdf", width = 6, height = 6)

# CBASS 35.6
CBASS35.6<-subset(traitdata, Exp=="CBASS-35.6"  ) 
CBASS35.6<-dplyr::select(CBASS35.6, FvFm, HP, TotalChl,ED) 
C<-cor(CBASS35.6,use = "complete.obs")
p.mat<-cor_pmat(CBASS35.6)
CBASS35.6plot<-ggcorrplot(C,hc.order=TRUE,outline.col="white",p.mat=p.mat,type='upper',lab=TRUE)
ggsave("./20230321_corrmatrix-35.6.pdf", width = 6, height = 6)

####Correlation of genotype correlations by trait####
#since there are so many NAs, first averaging replicates then working with mean values
#LT2<-LT %>%
#  group_by(Experiment,Exp,Genotype,Treatment) %>% 
#  summarise_at(vars("FvFm", "HP", "TotalChl","ED","initialFvFm"), list(mean,sd), na.rm=TRUE)
#CBASS2<-CBASS %>%
#  group_by(Experiment,Exp,Genotype,Treatment) %>% 
#  summarise_at(vars("FvFm", "HP", "TotalChl","ED"), mean, na.rm=TRUE)
#
#
HotCBASS2<-ddply(HotCBASS, .(Genotype), transform, deltFvFmHot=(FvFm[Treatment=="34.3"]-FvFm[Treatment=='27.5'])/FvFm[Treatment=='27.5'],
                deltChlHot=(TotalChl[Treatment=="34.3"]-TotalChl[Treatment=='27.5'])/TotalChl[Treatment=='27.5'], 
                deltHPHot=(HP[Treatment=="34.3"]-HP[Treatment=='27.5'])/HP[Treatment=='27.5'])
##HHotCBASS<-ddply(HHotCBASS, .(Genotype), transform, deltFvFmHHot=FvFm[Treatment=="35.6"]-FvFm[Treatment=='27.5'],deltChlHHot=TotalChl[Treatment=="35.6"]-TotalChl[Treatment=='27.5'], deltHPHHot=HP[Treatment=="35.6"]-HP[Treatment=='27.5'])
#LT2<-ddply(LT2, .(Genotype), transform, deltFvFmLT=(FvFm[Treatment=="OW"]-initialFvFm[Treatment=="OW"])-(FvFm[Treatment=="Control"]-initialFvFm[Treatment=="Control"]),
#           deltChlLT=(TotalChl[Treatment=="OW"]-TotalChl[Treatment=='Control'])/TotalChl[Treatment=='Control'], deltHPLT=(HP[Treatment=="OW"]-HP[Treatment=='Control'])/HP[Treatment=='Control'])
HotCBASS<-CBASS[CBASS$Treatment=='27.5' | CBASS$Treatment=='34.3',]

#not averaging first to get confidence intervals
HotCBASS<-ddply(HotCBASS, .(Rep), transform, deltFvFmHot=FvFm[Treatment=="34.3"]-FvFm[Treatment=='27.5'], deltChlHot=TotalChl[Treatment=="34.3"]-TotalChl[Treatment=='27.5'], deltHPHot=HP[Treatment=="34.3"]-HP[Treatment=='27.5'])
LT2<-ddply(LT, .(Rep), transform, deltFvFmLT=(FvFm[Treatment=="OW"]-initialFvFm[Treatment=="OW"])-(FvFm[Treatment=="Control"]-initialFvFm[Treatment=="Control"]),
           deltChlLT=(TotalChl[Treatment=="OW"]-TotalChl[Treatment=='Control'])/TotalChl[Treatment=='Control'], deltHPLT=HP[Treatment=="OW"]/HP[Treatment=='Control'])

HotCBASS<-CBASS2[CBASS2$Treatment=='34.3',]
LT3<-LT2[LT2$Treatment=='OW',]

#combine two df
geno_corr<-merge(HotCBASS2,LT3,  by=c('Genotype'))
geno_corr<-dplyr::rename(geno_corr, ED50=ED.x)
geno_corr<-dplyr::rename(geno_corr, ED25=ED.y)
geno_corr2<-dplyr::select(geno_corr, Genotype, ED50, ED25, deltFvFmHot, deltChlHot, deltHPHot, deltFvFmLT, deltHPLT, deltChlLT) 

# correlation - ED values 
CBASS.vs.LT_hightemp_lm <- lm(geno_corr2$ED25 ~ geno_corr2$ED50)
summary(CBASS.vs.LT_hightemp_lm) #R-squared:0.4194 #pval 0.3056559 

# correlation - ED50 and LTFvFm
CBASS.vs.LT_hightemp_lm <- lm(geno_corr2$deltFvFmLT ~ geno_corr2$ED50)
summary(CBASS.vs.LT_hightemp_lm) #R-squared:0.4162 # pval 0.04290505 

# correlation - ED50 and LT TotalChl 
CBASS.vs.LT_hightemp_lm <- lm(geno_corr2$deltChlLT~ geno_corr2$ED50)
summary(CBASS.vs.LT_hightemp_lm) #R-squared:0.504 #pval 0.352

# correlation -ED50 and LT HP 
CBASS.vs.LT_hightemp_lm <- lm(geno_corr2$deltHPLT ~ geno_corr2$ED50)
summary(CBASS.vs.LT_hightemp_lm) #R-squared: 0.2581, #pval 0.134

# correlation - ED50 and CBASS FvFm 
CBASS.vs.LT_hightemp_lm <- lm(geno_corr2$deltFvFmHot~ geno_corr2$ED50)
summary(CBASS.vs.LT_hightemp_lm) #R-squared:0.9723 #pval 1.62e-07

# correlation - ED50 and CBASS TotalChl 
CBASS.vs.LT_hightemp_lm <- lm(geno_corr2$deltChlHot~ geno_corr2$ED50)
summary(CBASS.vs.LT_hightemp_lm) #R-squared:0.108 #pval 0.352

# correlation - ED50 and CBASS HP
CBASS.vs.LT_hightemp_lm <- lm(geno_corr2$deltHPHot~ geno_corr2$ED50)
summary(CBASS.vs.LT_hightemp_lm) #R-squared:0.108 #pval 0.1015333

# correlation - LT FvFM and Chl
CBASS.vs.LT_hightemp_lm <- lm(geno_corr2$deltChlLT ~ geno_corr2$deltFvFmLT)
summary(CBASS.vs.LT_hightemp_lm) #R-squared:0.6427 #pval 0.00528 

# correlation - CBASS FvFM and Chl
CBASS.vs.LT_hightemp_lm <- lm(geno_corr2$deltChlHot ~ geno_corr2$deltFvFmHot)
summary(CBASS.vs.LT_hightemp_lm) #R-squared:0.1164 #pval 0.3346

####corr bn CBASS 35.6 and LT####
# correlation model - FvFm 35.6 
CBASS.vs.LT_corr <- cor.test(x = geno_corr$deltFvFmHHot, y = geno_corr$deltFvFmLT, method = "pearson")
CBASS.vs.LT_corr$p.value #0.6387923
CBASS.vs.LT_hightemp_lm <- lm(geno_corr$deltFvFmLT ~ geno_corr$deltFvFmHHot)
summary(CBASS.vs.LT_hightemp_lm) #R-squared:0.02888

# correlation model - TotalChl 35.6 
CBASS.vs.LT_corr <- cor.test(x = geno_corr$deltTotalChlHHot, y = geno_corr$deltTotalChlLT, method = "pearson")
CBASS.vs.LT_corr$p.value #0.1603066
CBASS.vs.LT_hightemp_lm <- lm(geno_corr$deltTotalChlLT ~ geno_corr$deltTotalChlHHot)
summary(CBASS.vs.LT_hightemp_lm) #R-squared:0.2304

# correlation model - HP 35.6 
CBASS.vs.LT_corr <- cor.test(x = geno_corr$deltHPHHot, y = geno_corr$deltHPLT, method = "pearson")
CBASS.vs.LT_corr$p.value #0.9099472
CBASS.vs.LT_hightemp_lm <- lm(geno_corr$deltHPLT ~ geno_corr$deltHPHHot)
summary(CBASS.vs.LT_hightemp_lm) #R-squared:0.0017 


####corrplot with genotype points####
#for averaging reps, probably not needed
geno_corr2<-dplyr::select(geno_corr, Genotype, ED50, ED25, deltFvFmLT, deltHPLT, deltChlLT) 
geno_corr3 <- geno_corr2 %>% group_by(Genotype,ED50) %>% mean_table(deltChlLT)
geno_corr3<-dplyr::rename(geno_corr3, ED50=group_2_cat)
geno_corr3<-dplyr::rename(geno_corr3, Genotype=group_1_cat)

#  xlab("\u0394Fv/Fm CBASS 34.3°C") + ylab("\u0394Fv/Fm LT") +
#  xlab("\u0394Host Protein CBASS 35.6°C") + ylab("\u0394Host Protein LT") +
#scale_x_continuous(n.breaks=5, labels = scales::label_number(accuracy = 1)) + scale_y_continuous(n.breaks=5) +

Pgeno_corr <- ggplot(geno_corr2, aes(x=deltFvFmHot,y=deltChlHot)) +
  stat_smooth(method = "lm",formula = y ~ x, se = TRUE, level = 0.95, na.rm = TRUE, colour="grey40")+
  geom_point(aes(color=Genotype), size = 8, alpha=0.8)+
  xlab("CBASS \u0394Fv/Fm") + ylab("CBASS \u0394Total Chl") +
  geom_label(aes(x = -0.37, y = -.79), 
             label = paste("Pearson R =",signif(summary(CBASS.vs.LT_hightemp_lm)$r.squared, 4),
                           "\np =",signif(summary(CBASS.vs.LT_hightemp_lm)$coef[2,4], 4)), label.size=NA, hjust='left', fill=alpha(c("white"),0)) +
  geom_text(aes(label=Genotype), vjust=0.5, size=3, fontface=c("bold")) +
  scale_x_continuous(n.breaks=5) + scale_y_continuous(n.breaks=5) +
  theme_classic() +
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.2 , "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position="none"
)
ggsave("./20230321_genocorr-CBASSFvFm_CBASSChl.pdf", width = 6, height = 6,device=cairo_pdf)

####Correlation of Fv/Fm with Total Chl across heat####
# correlation model - FvFm vs TotalChl 
CBASS34.3_corr <- cor.test(x = CBASS34.3$TotalChl, y = CBASS34.3$FvFm, method = "pearson")
CBASS34.3_corr$p.value #0.07701758
CBASS34.3_corr_lm <- lm(CBASS34.3$FvFm ~ CBASS34.3$TotalChl)
summary(CBASS34.3_corr_lm) # Adjusted R-squared: -0.02082  

CBASS35.6_corr <- cor.test(x = CBASS35.6$TotalChl, y = CBASS35.6$FvFm, method = "pearson")
CBASS35.6_corr$p.value #0.4440274
CBASS35.6_corr_lm <- lm(CBASS35.6$FvFm ~ CBASS35.6$TotalChl)
summary(CBASS35.6_corr_lm) # Adjusted R-squared: 0.1238 

LTOW<-subset(traitdata, Exp=="LT-OW" ) 
LTOW_corr <- cor.test(x = LTOW$TotalChl, y = LTOW$FvFm, method = "pearson")
LTOW_corr$p.value #0.0002271994
LTOW_corr_lm <- lm(LTOW$FvFm ~ LTOW$TotalChl)
summary(LTOW_corr_lm) # Adjusted R-squared: 0.3902 

# correlation model - ED50 vs TotalChl 
geno_corr34.3 <- cor.test(x = geno_corr$TotalChl, y = geno_corr$ED50, method = "pearson")
geno_corr34.3$p.value #0.4836629
geno_corr34.3_lm <- lm(geno_corr$ED50 ~ geno_corr$TotalChl)
summary(geno_corr34.3_lm) # Adjusted R-squared: 0.02763 

#CBASS35.6_corr <- cor.test(x = CBASS35.6$TotalChl, y = CBASS35.6$FvFm, method = "pearson")
#CBASS35.6_corr$p.value #0.4440274
#CBASS35.6_corr_lm <- lm(CBASS35.6$FvFm ~ CBASS35.6$TotalChl)
#summary(CBASS35.6_corr_lm) # Adjusted R-squared: 0.1238 
LT3<-LT3%>%drop_na(TotalChl)
LTOW_corr <- cor.test(x = LT3$TotalChl, y = LT3$ED, method = "pearson")
LTOW_corr$p.value #0.6366283
LTOW_corr_lm <- lm(LT3$ED ~ LT3$TotalChl)
summary(LTOW_corr_lm) # Adjusted R-squared: 0.007059


# plot FvFm vs Total chl for LT, CBASS 34.3, CBASS 35.6
P1 <- ggplot(CBASS34.3,aes(x=TotalChl,y=FvFm)) +
  theme_classic() +
  stat_smooth(method = "lm",formula = y ~ x, se = TRUE, level = 0.95, na.rm = TRUE, colour="grey40")+
  geom_point(size = 2, alpha=0.8)+
  theme(legend.position = 'none')+
  xlab("TotalChl") + ylab(label = expression(paste(italic(F[V]/F[M])))) + ggtitle("CBASS 34.3°C") +
  geom_label(aes(x = 1.5, y = 0.425), hjust="left",
             label =paste("Adj R2 = ",signif(summary(CBASS34.3_corr_lm)$adj.r.squared, 5),
                           "\nP =",signif(summary(CBASS34.3_corr_lm)$coef[2,4], 5)), label.size=NA) +
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.2 , "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)
  )
P2 <- ggplot(CBASS35.6,aes(x=TotalChl,y=FvFm)) +
  theme_classic() +
  stat_smooth(method = "lm",formula = y ~ x, se = TRUE, level = 0.95, na.rm = TRUE, colour="grey40")+
  geom_point(size = 2, alpha=0.8)+
  theme(legend.position = 'none')+
  xlab("TotalChl ") + ylab(label = expression(paste(italic(F[V]/F[M])))) + ggtitle("CBASS 35.6°C") +
  geom_label(aes(x = 0.2, y = 0.15), hjust="left",
             label =paste("Adj R2 = ",signif(summary(CBASS35.6_corr_lm)$adj.r.squared, 5),
                          "\nP =",signif(summary(CBASS35.6_corr_lm)$coef[2,4], 5)), label.size=NA) +
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.2 , "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)
  )
P3 <- ggplot(LTOW,aes(x=TotalChl,y=FvFm)) +
  theme_classic() +
  stat_smooth(method = "lm",formula = y ~ x, se = TRUE, level = 0.95, na.rm = TRUE, colour="grey40")+
  geom_point(size = 2, alpha=0.8)+
  theme(legend.position = 'none')+
  xlab("TotalChl") + ylab(label = expression(paste(italic(F[V]/F[M])))) + ggtitle("LT-OW") +
  geom_label(aes(x = 0.1, y = 0.585), hjust="left",
             label =paste("Adj R2 = ",signif(summary(LTOW_corr_lm)$adj.r.squared, 5),
                          "\nP =",signif(summary(LTOW_corr_lm)$coef[2,4], 5)), label.size=NA) +
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.2 , "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)
  )
ggarrange(P3,P1,P2, labels=c("A","B","C"), font.label=list(size=16), ncol=3, nrow=1)
ggsave("./20221118_FvFmvsTotalChl.pdf", width = 10, height = 6,device=cairo_pdf)

#####LTvsCBASS FvFm######
#Both experiments 
bartlett.test(FvFm~Exp, data=trait) #Genotype equal variance, exp unequal
#welch one-way with unequal variance
oneway.test(FvFm~Exp, data=trait, var.equal=FALSE)
pairwise.wilcox.test(trait$FvFm, trait$Exp, p.adjust.method = "none", exact=FALSE)   
FvFm_lmer<-lmer(FvFm ~ Exp+ (1|Tank),data=trait)
anova(FvFm_lmer)
rand(FvFm_lmer)

#CBASS 34.3 and LT-OW
FvFmheat<-trait[trait$Exp=='LT-OW'| trait$Exp=='CBASS-34.3'  ,]
#CBASSheat<-trait[trait$Exp=='CBASS-34.3'  ,]
bartlett.test(deltFvFm~Exp, data=trait2) #Genotype and exp equal
#welch one-way with unequal variance
#FvFmcont_welch<-oneway.test(FvFm~Exp, data=FvFmcont, var.equal=FALSE)
#pairwise.wilcox.test(FvFmheat$FvFm, FvFmheat$Exp, p.adjust.method = "bonf", exact=FALSE)    
FvFm_lmer<-lmer(deltFvFm ~ Exp*Genotype  + (1|Tank),data=trait2)
sjPlot::plot_model(FvFm_lmer, type="diag")
step(FvFm_lmer)
FvFm_lmer<-lmer(deltFvFm ~ Exp*Region + (1|Tank),data=trait2)
anova(FvFm_lmer)
rand(FvFm_lmer)
print(emmeans(FvFm_lmer, list(pairwise ~Region|Exp)), adjust = c("tukey"))

#all CBASS trts
#welch one-way with unequal variance
bartlett.test(FvFm~Exp, data=CBASS) # exp unequal
FvFm_welch<-oneway.test(FvFm~Exp, data=CBASS)
pairwise.wilcox.test(CBASS$FvFm, CBASS$Exp, p.adjust.method = "bonf", exact=FALSE)   
FvFm_lmer<-lmer(FvFm ~ Treatment*Genotype + (1|Tank),data=CBASS)
#FvFm_lmer<-lmer(FvFm ~Genotype + (1|Tank),data=CBASSheat)

anova(FvFm_lmer)
ranova(FvFm_lmer)
print(emmeans(FvFm_lmer, list(pairwise ~Exp)), adjust = c("tukey"))

#LT Stats overall model testing Exp and genotype with tank as random
FvFm_lmer<-lmer(FvFm ~ Treatment+Genotype+Treatment*Genotype + (1|Tank) + (1|Rep),data=CBASS)
anova(FvFm_lmer)
rand(FvFm_lmer)
step(FvFm_lmer, reduce.random=TRUE)
#region
FvFm_lmer<-lmer(FvFm ~ Treatment*Region + (1|Tank),data=CBASS)
anova(FvFm_lmer)
rand(FvFm_lmer)
print(emmeans(FvFm_lmer, list(pairwise ~ Region|Treatment)), adjust = c("tukey"))

#photoacclimation
LTtime<-transform(LTtime,
                  Genotype=as.factor(Genotype),
                  Treatment=as.factor(Treatment),
                  Tank=as.factor(Tank))
LTtimecont<-LTtime[LTtime$Treatment=='Control' ,]
LTtimeheat<-LTtime[LTtime$Treatment=='OW' ,]

bartlett.test(FvFm.I~Treatment, data=LTtimecont) # time unequal
oneway.test(FvFm.I~Time, data=LTtimecont)

FvFm_lmer<-lmer(deltFvFm ~ Treatment*Genotype + (1|Tank),data=LTtime)
anova(FvFm_lmer)
rand(FvFm_lmer)
sjPlot::plot_model(FvFm_lmer, type="diag")
#prints nonsensical comparisons ie cu80 37 vs ac31 34.3
print(emmeans(FvFm_lmer, list(pairwise ~ Genotype|Treatment)), adjust = c("tukey"))
print(emmeans(FvFm_lmer, list(pairwise ~ Exp)), adjust = c("tukey"))

#or, to get letters of diff pairwise
leastsquare=lsmeans(FvFm_lmer,pairwise~Exp,adjust='tukey')
multcomp::cld(leastsquare,alpha=0.05,Letters=letters,adjust='tukey')

# Model fitting and assumptions diagnostic 
x = residuals(FvFm_lmer)
plot(fitted(FvFm_lmer), x) 
leveneTest(FvFm ~ Genotype * Exp, data=trait3, center=mean) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
bartlett.test(FvFm~Exp, data=trait3) #Genotype equal variance, exp unequal
qqnorm(x); qqline(x) # qq plot to check for normal distribution of residuals
hist(x) # histogram of residuals to check for normal distribution of residuals
shapiro.test(x) # formal statistical test
#boxcox doesn't like zeros, add .015 to all data values
bc<-boxcox(FvFm.I~Genotype*Time*Treatment, data=LTtime)
lambda <- bc$x[which.max(bc$y)]
#lambda=2
LTtime$FvFmsq<-LTtime$FvFm.I^2

#summary for delta point plot, FvFm Time
LTtime$deltFvFm<- (LTtime$FvFm.E-LTtime$FvFm.I)
LTtime$Grtime<-paste(LTtime$Genotype, LTtime$Time, sep=".")
sum<-Summarize(FvFm~Treatment+Grtime, data=LTtime, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$Genotype<-gl(10,4,labels=c("LK31","LK41","LK50","LK62","LK7", "UK12", "UK19", "UK70","UK76","UK80"))
sum$Time<-gl(2,2,40,labels=c("Final","Initial"))
sum$Time<-factor(sum$Time,levels(sum$Time)[c(2,1)])

geom_size <- 1
pd=position_dodge(.75)
plot.fvfm<-ggplot(sum,aes(x=Time,y=mean,color=Genotype)) + 
  geom_point(size=4,position=pd) +
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=Genotype), width=.2, size=geom_size,position=pd) + 
  geom_line(aes(group=Genotype),position=pd,size=.75) +
  facet_wrap(~Treatment) +
  ylab(label = expression(italic(F[V]/F[M]))) + 
  xlab("Time") +
  #theme_classic() + 
  theme(line= element_line(size = 0.6),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.2 , "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = "gray80",color="black",linewidth=1),
        strip.text= element_text(size= 12),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "bottom",
        legend.key = element_rect(fill = NA)) 

ggsave(plot.fvfm, height = 6 , width = 8, filename = "./2023_acer_LTfvfmtime.pdf",device=cairo_pdf)
#legend.position = c(1, .6),
#legend.text = element_text(size=10),
#legend.title = element_blank()
#summary for point plot
sum<-Summarize(deltFvFm~Exp, data=trait2, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
Region<-factor(c("Lower Keys","Upper Keys"))
sum$Region<-rep(Region,times=2,  each=5)
SubRegion<-factor(c("CBASS Lower Keys","CBASS Upper Keys","LT Lower Keys","LT Upper Keys"))
sum$SubRegion<-rep(SubRegion,  each=5)
Regionmean <- sum %>% group_by(SubRegion) %>% summarise(count = n(),mean_val=mean(mean),    
                                                        LB_FvFm = mean_val - qnorm(0.975) * sd(mean) / sqrt(count),
                                                        UB_FvFm = mean_val + qnorm(0.975) * sd(mean) / sqrt(count) )
sum$Regionmean<-c(0.2710,0.2710,0.2710,0.2710,0.2710,0.3184,0.3184,0.3184,0.3184,0.3184,0.3690,0.3690,0.3690,0.3690,0.3690,0.3242,0.3242,0.3242,0.3242,0.3242)
sum$RegionLB<-c(0.2164,0.2164,0.2164,0.2164,0.2164,0.2830,0.2830,0.2830,0.2830,0.2830,0.3329,0.3329,0.3329,0.3329,0.3329,0.2828,0.2828,0.2828,0.2828,0.2828)
sum$RegionUB<-c(0.3255,0.3255,0.3255,0.3255,0.3255,0.3537,0.3537,0.3537,0.3537,0.3537,0.4050,0.4050,0.4050,0.4050,0.4050,0.3655,0.3655,0.3655,0.3655,0.3655)

plot.pam<-ggplot(sum,aes(x=Genotype,y=mean,color=Exp,shape=Region)) + 
  facet_grid(Exp~Region,scales="free_x") + 
  geom_ribbon(aes(group=SubRegion,ymin=RegionLB,ymax=RegionUB,fill=Exp),alpha=.2,colour=NA) +
  geom_line(aes(group=SubRegion,y=Regionmean),size=2,alpha=.4) +
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=Exp), width=.3, size=1.5) + 
  geom_point(size=5,stroke=2) +
  scale_color_manual(values=c("#FDAE61","#A50026")) +
  scale_fill_manual(values=c("#FDAE61","#A50026"), guide='none') +
  scale_shape_manual(values=c(19,17,19,17)) +
  labs(color= "Experiment") +
  #ylim(-0.05,0.6) +
  ylab(label = expression(paste(italic(F[V]/F[M])))) + 
  xlab("Genotype") + 
  theme_classic() + 
  theme(line= element_line(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.2 , "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.spacing.x=unit(0, "lines"),
        strip.background.x = element_blank(), 
        strip.text.y=element_blank(),
        strip.text.x=element_text(size=14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = 'bottom',
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) 
ggsave(plot.pam, height = 6 , width = 7, bg='transparent',filename = "./20230310_CBASS34.3vsLTpam-regionXExp.pdf", useDingbats=FALSE)

#boxplot of CBASS and LT
trait$Treatment=factor(trait$Treatment,levels(trait$Treatment)[c(2,3,4,5,6,7,8,1,9)])

brewer.pal(n = 10, name = "RdYlBu")
#       
cols<- c("CBASS-27.5" = "#4575B4", "CBASS-28.9" = "#74ADD1","CBASS-30.2" = "#ABD9E9", "CBASS-31.6" = "#E0F3F8", 
         "CBASS-32.9" = "#FEE090", "CBASS-34.3" = "#FDAE61", "CBASS-35.6" = "#F46D43" , "LT-Control" = "#313695", "LT-OW" = "#D73027") 
FvFmbox <- ggplot(data=trait, aes(x=Treatment, y=FvFm, label=Treatment, fill=Exp)) +
  scale_fill_manual(values = cols) +
  stat_boxplot(geom ='errorbar', width = 0.7, lwd=1) +
  geom_boxplot(width=0.7, lwd=1, fatten=1, outlier.shape=NA) +
  geom_jitter(color="black", width=0.1, alpha=0.9) +
  expand_limits(y = 0)+
  theme_bw() +
  ylab(label = expression(paste(italic(F[V]/F[M])))) +
  #annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, colour = "black", size=1) +
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.2 , "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "none")

ggsave(FvFmbox, height = 6 , width = 10, filename = "./20221116_FvFm_experiment-boxplot.pdf", useDingbats=FALSE)

#####LTvsCBASS heat Host Proteins#####
###stats###
# Statistical tests for maximum quantum yield retained at end of assay 
# Analysis of variance (2-way) ANOVA dat (mixed linear model)
# Mixed Effects model origin_dest- random effect of colony nested in tank
HP.model <- lmer(HPsqrt ~ Genotype*Exp+ (1|Tank),data=trait)
HP.model <- lmer(HPsqrt ~ Treatment*Region+ (1|Tank),data=CBASS)

anova(HP.model)
rand(HP.model)
step(HP.model, reduce.random=FALSE)
# statistical pairwise comparisons
print(emmeans(HP.model, list(pairwise ~ Exp)), adjust = c("tukey"))
print(emmeans(HP.model, list(pairwise ~ Treatment|Region)), adjust = c("tukey"))

# Model fitting and assumptions diagnostic 
sjPlot::plot_model(HP.model, type="diag")
x = residuals(HP.model)
shapiro.test(x) # formal statistical test
leveneTest(HPsqrt ~ Genotype * Exp, data=LT, center=mean) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
#find optimal lambda for Box-Cox transformation 
#boxcox doesn't like zeros, add .015 to all data values
trait$HPpos <- trait$HP+0.015
bc<-boxcox(HPpos~Genotype*Exp, data=trait)
lambda <- bc$x[which.max(bc$y)]
#lambda=0.3
FvFmheat$HPsqrt <- sqrt(FvFmheat$HP)
LT$HPsqrt <- sqrt(LT$HP)
CBASS$HPsqrt <- sqrt(CBASS$HP)

#welchs one-way since Exp is unequal var
bartlett.test(HPsqrt~Exp, data=trait)
oneway.test(HPsqrt~Exp, data=trait)
pairwise.wilcox.test(trait$HPsqrt, trait$Exp, p.adjust.method = "none", exact=FALSE)            

#CBASS 34.3 and 35.6 and LT-OW
bartlett.test(HPsqrt~Exp, data=FvFmheat) #Equal variance with HPsqrt
HP_lmer<-lmer(HPsqrt ~ Exp*Region + (1|Tank) ,data=FvFmheat)
sjPlot::plot_model(HP_lmer, type="diag")
step(HP_lmer)
anova(HP_lmer)
rand(HP_lmer)
# statistical pairwise comparisons
print(emmeans(HP_lmer, list(pairwise ~ Exp|Region)), adjust = c("tukey"))

###plots###
###boxplot by experiment###
Protein <- ggplot(data=trait, aes(x=Treatment, y=HP, label=Treatment, fill=Exp)) +
  scale_fill_manual(values = cols) +
  stat_boxplot(geom ='errorbar', width = 0.7, lwd=1) +
  geom_boxplot(width=0.7, lwd=1, fatten=1, outlier.shape=NA) +
  geom_jitter(color="black", width=0.1, alpha=0.9) +
  expand_limits(y = 0)+
  theme_bw() +
  ylab(bquote('Protein (ug/cm'^'2'*')')) +
  #annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, colour = "black", size=1) +
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.2 , "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "none")
ggsave(Protein, height = 6 , width = 10, filename = "./20221213_HP_byexperiment-boxplot.pdf", useDingbats=FALSE)

###points###
sum<-Summarize(HP~Genotype+Exp, data=FvFmheat, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
Region<-factor(c("Lower Keys","Upper Keys"))
sum$Region<-rep(Region,times=2,  each=5)
SubRegion<-factor(c("CBASS Lower Keys","CBASS Upper Keys","LT Lower Keys","LT Upper Keys"))
sum$SubRegion<-rep(SubRegion,  each=5)
Regionmean <- sum %>% group_by(SubRegion) %>% summarise(count = n(),mean_val=mean(mean),    
                                                        LB_HP = mean_val - qnorm(0.975) * sd(mean) / sqrt(count),
                                                        UB_HP = mean_val + qnorm(0.975) * sd(mean) / sqrt(count) )
mean_val<-c(44.3306,39.2504,22.2306,14.8160)
sum$Regionmean<-rep(mean_val, each=5)
RegionLB<-c(31.866581,23.205951,15.660585,7.031679)
sum$RegionLB<-rep(RegionLB, each=5)
RegionUB<-c(56.79462,55.29485,28.80061,22.60032)
sum$RegionUB<-rep(RegionUB, each=5)

#points with sd or ci 95% bars, comment out which you want
# #FDAE61",'#F46D43',"#A50026"
plot.HP<-ggplot(sum,aes(x=Genotype,y=mean,color=Exp,shape=Region)) + 
  facet_grid(Exp~Region,scales="free_x") + 
  geom_ribbon(aes(group=SubRegion,ymin=RegionLB,ymax=RegionUB,fill=Exp),alpha=.2,colour=NA) +
  geom_line(aes(group=SubRegion,y=Regionmean),size=2,alpha=.4) +
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=Exp), width=.3, size=1.5) + 
  geom_point(size=5,stroke=2) +
  scale_color_manual(values=c("#FDAE61","#A50026")) +
  scale_fill_manual(values=c("#FDAE61","#A50026"), guide='none') + 
  scale_shape_manual(values=c(19,17,19,17)) +
  labs(color= "Experiment") +
  #ylim(-0.05,0.6) +
  ylab(bquote('Protein (ug/cm'^'2'*')')) + 
  xlab("Genotype") + 
  theme_classic() + 
  theme(line= element_line(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.2 , "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.spacing.x=unit(0, "lines"),
        strip.background.x = element_blank(), 
        strip.text.y=element_blank(),
        strip.text.x=element_text(size=14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = 'bottom',
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) 
ggsave(plot.HP, height = 6 , width = 7, bg='transparent', filename = "./20230310_CBASS34.3vsLTHP-regionXExp.pdf", useDingbats=FALSE)

#####LTvsCBASS heat Chlorophylls#####
# Mixed Effects model origin_dest- random effect of colony nested in tank
TotalChl.model <- lmer(TotalChlsqrt ~ Genotype*Exp+ (1|Tank),data=trait)
TotalChl.model <- lmer(TotalChlsqrt ~ Genotype*Treatment + (1|Tank),data=CBASS)

anova(TotalChl.model)
rand(TotalChl.model)
step(TotalChl.model, reduce.random=FALSE)
# statistical pairwise comparisons
print(emmeans(TotalChl.model, list(pairwise ~ Exp)), adjust = c("tukey"))
print(emmeans(TotalChl.model, list(pairwise ~ Treatment)), adjust = c("tukey"))

# Model fitting and assumptions diagnostic 
sjPlot::plot_model(TotalChl.model, type="diag")
x = residuals(TotalChl.model)
shapiro.test(x) # formal statistical test
leveneTest(TotalChl ~ Genotype * Exp, data=LT, center=mean) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
#find optimal lambda for Box-Cox transformation 
#boxcox doesn't like zeros, add .015 to all data values
trait$TotalChlpos <- trait$TotalChl+0.015
bc<-boxcox(TotalChlpos~Genotype*Exp, data=trait)
lambda <- bc$x[which.max(bc$y)]
#lambda=0.3
FvFmheat$TotalChlsqrt <- sqrt(FvFmheat$TotalChl)
CBASS$TotalChlsqrt <- sqrt(CBASS$TotalChl)
FvFmheat$TotalChlsqrt <- sqrt(FvFmheat$TotalChl)

#welchs one-way since Exp is unequal var
bartlett.test(TotalChlsqrt~Exp, data=CBASS) #Exp unequal geno equal
oneway.test(TotalChlsqrt~Exp, data=CBASS)
pairwise.wilcox.test(CBASS$TotalChlsqrt, CBASS$Exp, p.adjust.method = "bonf", exact=FALSE)            

#CBASS 34.3 and 35.6 and LT-OW
bartlett.test(TotalChlsqrt~Exp, data=FvFmheat) #Equal variance with TotalChlsqrt
TotalChl_lmer<-lmer(TotalChlsqrt ~ Exp*Genotype+ (1|Tank),data=FvFmheat)
sjPlot::plot_model(TotalChl_lmer, type="diag")
step(TotalChl_lmer)
anova(TotalChl_lmer)
rand(TotalChl_lmer)
# statistical pairwise comparisons
print(emmeans(TotalChl_lmer, list(pairwise ~ Exp|Region)), adjust = c("tukey"))

###plots###
###boxplot by experiment###
TotalChl <- ggplot(data=trait, aes(x=Treatment, y=TotalChl, label=Treatment, fill=Exp)) +
  scale_fill_manual(values = cols) +
  stat_boxplot(geom ='errorbar', width = 0.7, lwd=1) +
  geom_boxplot(width=0.7, lwd=1, fatten=1, outlier.shape=NA) +
  geom_jitter(color="black", width=0.1, alpha=0.9) +
  expand_limits(y = 0)+
  theme_bw() +
  ylab(bquote('Total Chlorophyll (ug/cm'^'2'*')')) +
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.2 , "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "none")
ggsave(TotalChl, height = 6 , width = 10, filename = "./20221213_TotalChl_byexperiment-boxplot.pdf", useDingbats=FALSE)

###points###
sum<-Summarize(TotalChl~Genotype+Exp, data=FvFmheat, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
Region<-factor(c("Lower Keys","Upper Keys"))
sum$Region<-rep(Region,times=2,  each=5)
SubRegion<-factor(c("CBASS Lower Keys","CBASS Upper Keys","LT Lower Keys","LT Upper Keys"))
sum$SubRegion<-rep(SubRegion,  each=5)
Regionmean <- sum %>% group_by(SubRegion) %>% summarise(count = n(),mean_val=mean(mean),    
                                                        LB_TotalChl = mean_val - qnorm(0.975) * sd(mean) / sqrt(count),
                                                        UB_TotalChl = mean_val + qnorm(0.975) * sd(mean) / sqrt(count) )
mean_val<-c(5.8672,6.1482,1.3842,0.4868)
sum$Regionmean<-rep(mean_val, each=5)
RegionLB<-c(4.3398407,4.0046507,0.9307978,0.2688353)
sum$RegionLB<-rep(RegionLB, each=5)
RegionUB<-c(7.3945593,8.2917493,1.8376022,0.7047647)
sum$RegionUB<-rep(RegionUB, each=5)

#points with sd or ci 95% bars, comment out which you want
# #FDAE61",'#F46D43',"#A50026"
plot.TotalChl<-ggplot(sum,aes(x=Genotype,y=mean,color=Exp,shape=Region)) + 
  facet_grid(Exp~Region,scales="free_x") + 
  geom_ribbon(aes(group=SubRegion,ymin=RegionLB,ymax=RegionUB,fill=Exp),alpha=.2,colour=NA) +
  geom_line(aes(group=SubRegion,y=Regionmean),size=2,alpha=.4) +
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=Exp), width=.3, size=1.5) + 
  geom_point(size=5,stroke=2) +
  scale_color_manual(values=c("#FDAE61","#A50026")) +
  scale_fill_manual(values=c("#FDAE61","#A50026"), guide='none') + 
  scale_shape_manual(values=c(19,17,19,17)) +
  labs(color= "Experiment") +
  #ylim(-0.05,0.6) +
  ylab(bquote('Total Chlorophyll (ug/cm'^'2'*')')) +   
  xlab("Genotype") + 
  theme_classic() + 
  theme(line= element_line(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.2 , "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.spacing.x=unit(0, "lines"),
        strip.background.x = element_blank(), 
        strip.text.y=element_blank(),
        strip.text.x=element_text(size=14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = 'bottom',
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) 
ggsave(plot.TotalChl, height = 6 , width = 7, bg='transparent', filename = "./20230310_CBASS34.3vsLTTotalChl-regionXExp.pdf", useDingbats=FALSE)

#### TotalChlCBASS2.0####
# Statistical tests for TotalChl at of assay 
# Analysis of variance (2-way) ANOVA dat (mixed linear model)
# Mixed Effects model origin_dest- random effect of colony nested in tank
# Since the model is completely balanced the SS type I, II or III will provide the same result
TotalChl.model <- lmer(TotalChl ~ Treatment  * Genotype + (1|Tank) , data = cbass)
#on transformed data
TotalChl.model <- lmer(TotalChlsqrt ~ Treatment  * Genotype + (1|Tank) , data = cbass)

#summary(TotalChl.model)
anova(TotalChl.model)
rand(TotalChl.model)
#confint(TotalChl.model,oldNames=F)

x = residuals(TotalChl.model)
plot(fitted(TotalChl.model), x) 
leveneTest(x ~ Genotype * Treatment, data=cbass, center=mean) # formal statistical test for homogeinity of variance (not recommed due the small sample size)
qqnorm(x); qqline(x) # qq plot to check for normal distribution of residuals
hist(x) # histogram of residuals to check for normal distribution of residuals
shapiro.test(x) # formal statistical test
outlier(cbass$TotalChl) #OF125
cbass<- tibble::rownames_to_column(cbass, "VALUE") #rownames retained from original df, rowname to column to find out true row, here it's 125
cbass<-cbass[-c(125), ] 

#find optimal lambda for Box-Cox transformation 
bc<-boxcox(TotalChl~Genotype*Treatment, data=cbass)
lambda <- bc$x[which.max(bc$y)]
#lambda=.5 sqrt
cbass <- cbind(cbass, sqrt(cbass$TotalChl))
dim(cbass)
names(cbass)[46] <- "TotalChlsqrt"
##### statistical pairwise comparisons
# Genotype*treatment
TotalChl.emms <- emmeans(TotalChl.model, pairwise ~ Genotype|Treatment, weights = "proportional", adjust="none")
# P.value adjustment of the Bonferroni
rbind(TotalChl.emms$contrasts, adjust="tukey")
#Genotype
TotalChl.emms <- emmeans(TotalChl.model, pairwise ~ Genotype, weights = "proportional", adjust="none")
rbind(TotalChl.emms$contrasts, adjust="tukey")
#treatment
TotalChl.emms <- emmeans(TotalChl.model, pairwise ~ Treatment, weights = "proportional", adjust="none")
rbind(TotalChl.emms$contrasts, adjust="tukey")

### Plots
##### points over time
sum<-Summarize(TotalChl~Genotype+Treatment, data=cbass, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$Treatment=factor(sum$Treatment,levels(sum$Treatment)[c(1,2,4,3)])
sum$Genotype=factor(sum$Genotype,levels(sum$Genotype)[c(1,6,8,10,12,2,3,7,9,11,4,5)])

#points with sd or ci 95% bars, comment out which you want
pd=position_dodge(.5)
plot.TotalChl<-ggplot(sum,aes(x=Treatment,y=mean,group=Genotype)) + 
  geom_line(aes(color=Genotype)) +
  geom_point() +
  scale_color_brewer(palette = "Paired") +
  ylab("Total Chlorophyll (ug/cm^2)") + 
  xlab("Treatment") +
  theme_classic() + 
  theme(axis.text.x=element_text(size=rel(0.7),angle = -45, hjust = 0)) +
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) 
ggsave(plot.TotalChl, height = 6 , width = 10, filename = "./202201_cbass2_TotalChl.pdf", useDingbats=FALSE)

####ACER LT weight####
#pivot table to long format, need to subset by trait of interest
acertopivot<-dplyr::select(acer,Genotype,Treatment,Tank,IW,MW,EW)
acerlong<-pivot_longer(data=acertopivot,cols = -c(Genotype,Treatment,Tank), names_to = "time", values_to = "weight")

# Statistical tests for maximum quantum yield retained at end of assay 
# Analysis of variance (2-way) ANOVA dat (mixed linear model)
# Mixed Effects model origin_dest- random effect of colony nested in tank
# Since the model is completely balanced the SS type I, II or III will provide the same result
weight.model <- lmer(I.Erate ~ Treatment  * Genotype + (1|Tank) , data = acer)
weighttime.model <- lmer(weight ~ Treatment  * Genotype * time + (1|Tank) , data = acerlong)

#summary(weight.model)
anova(weighttime.model)
rand(weighttime.model)
#confint(weight.model,oldNames=F)

# Model fitting and assumptions diagnostic 
x = residuals(weight.model)
plot(fitted(weight.model), x) 
leveneTest(x ~ Genotype * Treatment, data=acer, center=mean) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
qqnorm(x); qqline(x) # qq plot to check for normal distribution of residuals
hist(x) # histogram of residuals to check for normal distribution of residuals
shapiro.test(x) # formal statistical test
#acer outliers > ac31 oaow row 11, cu12 OA  row 106, remove
acer<-acer[-c(11,108), ]

# statistical pairwise comparisons
# genotype
weight.emms <- emmeans(weighttime.model, pairwise ~ Genotype|time|Treatment, weights = "proportional", adjust="none")
summary(weight.emms$emmeans)
# P.value adjustment of the Bonferroni
rbind(weight.emms$contrasts, adjust="tukey")
# treatment
weight.emms <- emmeans(weight.model, pairwise ~ Genotype, weights = "proportional", adjust="none")
summary(weight.emms$emmeans)
rbind(weight.emms$contrasts, adjust="tukey")
#or
print(emmeans(weighttime.model, list(pairwise ~ time)), adjust = c("tukey"))

###plot###
###points over time
acerlong$grtime<-interaction(acerlong$Treatment,acerlong$time) 

sum<-Summarize(weight~Genotype+grtime, data=acerlong, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$Treatment<-gl(4,10,labels=c("Control","OW","OA","OAOW"))
sum$Time<-gl(3,40,labels=c("End","Initial","Mid"))
sum$Time=factor(sum$Time,levels(sum$Time)[c(2,3,1)])
summary(sum)

#points with sd or ci 95% bars, comment out which you want
pd=position_dodge(.75)
plot.weight<-ggplot(sum,aes(x=Time,y=mean,color=Genotype,group=Genotype)) + 
  geom_point(size=3,position=pd) +
  geom_line(position=pd) +
  scale_color_brewer(palette = "Paired") +
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=Genotype), width=.2, position=pd) + 
  ylab("Dry Weight (g)") + 
  xlab("Time point") +
  theme_bw() + 
  facet_wrap(~Treatment) +
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) 
ggsave(plot.weight, height = 6 , width = 10, filename = "./figures/20210701_acer_weightbytrt.pdf", useDingbats=FALSE)

#points as end growth rate
sum<-Summarize(I.Erate~Genotype+Treatment, data=acer, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$Treatment=factor(sum$Treatment,levels(sum$Treatment)[c(1,4,2,3)])

summary(sum)

plot.EIgrate<-ggplot(sum,aes(x=Genotype,y=mean,color=Genotype)) + 
  geom_point(size=3) +
  scale_color_brewer(palette = "Paired") +
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=Genotype), width=.2) + 
  ylab("Growth Rate (g dry weigh/day)") + 
  xlab("Genotype") +
  theme_bw() + 
  facet_wrap(~Treatment) +
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) 
ggsave(plot.EIgrate, height = 6 , width = 10, filename = "./figures/20210706_acer_endgrowthrate.pdf", useDingbats=FALSE)


####ACER LT FvFm####
#pivot table to long format, need to subset by trait of interest
acertopivot<-dplyr::select(acer,Genotype,Treatment,Tank,IFvFm,MFvFm,EFvFm)
acerlong<-pivot_longer(data=acertopivot,cols = -c(Genotype,Treatment,Tank), names_to = "time", values_to = "fvfm")

# Statistical tests for maximum quantum yield retained at end of assay 
# Analysis of variance (2-way) ANOVA dat (mixed linear model)
# Mixed Effects model origin_dest- random effect of colony nested in tank
# Since the model is completely balanced the SS type I, II or III will provide the same result
fvfm.model <- lmer(EFvFm ~ Treatment  * Genotype + (1|Tank) , data = acer)
fvfmtime.model <- lmer(fvfm ~ Treatment  * Genotype * time + (1|Tank) , data = acerlong)

#summary(fvfm.model)
anova(fvfmtime.model)
rand(fvfmtime.model)
#confint(fvfm.model,oldNames=F)

# Model fitting and assumptions diagnostic 
x = residuals(fvfm.model)
plot(fitted(fvfm.model), x) 
leveneTest(x ~ Genotype * Treatment, data=acer, center=mean) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
qqnorm(x); qqline(x) # qq plot to check for normal distribution of residuals
hist(x) # histogram of residuals to check for normal distribution of residuals
shapiro.test(x) # formal statistical test

# statistical pairwise comparisons
# genotype
fvfm.emms <- emmeans(fvfmtime.model, pairwise ~ Genotype|time|Treatment, weights = "proportional", adjust="none")
summary(fvfm.emms$emmeans)
# P.value adjustment of the Bonferroni
rbind(fvfm.emms$contrasts, adjust="tukey")
# treatment
fvfm.emms <- emmeans(fvfm.model, pairwise ~ Genotype, weights = "proportional", adjust="none")
summary(fvfm.emms$emmeans)
rbind(fvfm.emms$contrasts, adjust="tukey")
#or
print(emmeans(fvfmtime.model, list(pairwise ~ time)), adjust = c("tukey"))

###plot###
###points over time
acerlong$grtime<-interaction(acerlong$Treatment,acerlong$time) 

sum<-Summarize(fvfm~Genotype+grtime, data=acerlong, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$Treatment<-gl(4,10,labels=c("Control","OW","OA","OAOW"))
sum$Time<-gl(3,40,labels=c("End","Initial","Mid"))
sum$Time=factor(sum$Time,levels(sum$Time)[c(2,3,1)])
summary(sum)

#points with sd or ci 95% bars, comment out which you want
pd=position_dodge(.5)
plot.fvfm<-ggplot(sum,aes(x=Time,y=mean,color=Genotype,group=Genotype)) + 
  geom_point(size=3,position=pd) +
  scale_color_brewer(palette = "Paired") +
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=Genotype), width=.2, position=pd) + 
  ylab("Fv/Fm") + 
  xlab("Time point") +
  theme_bw() + 
  facet_wrap(~Treatment) +
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) 
ggsave(plot.fvfm, height = 6 , width = 10, filename = "./figures/20210706_acer_fvfmbytrt.pdf", useDingbats=FALSE)

#points as end Fv/Fm normalized to start
acer$mqy<-(1-(acer$IFvFm-acer$EFvFm)/acer$IFvFm)
sum<-Summarize(IFvFm~Genotype+Treatment, data=acer, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$Treatment=factor(sum$Treatment,levels(sum$Treatment)[c(1,4,2,3)])

summary(sum)

plot.endfvfm<-ggplot(sum,aes(x=Genotype,y=mean,color=Genotype)) + 
  geom_point(size=3) +
  scale_color_brewer(palette = "Paired") +
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=Genotype), width=.2) + 
  ylab("Fv/Fm Retained") + 
  xlab("Genotype") +
  theme_bw() + 
  facet_wrap(~Treatment) +
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) 
ggsave(plot.endfvfm, height = 6 , width = 10, filename = "./figures/20210706_acer_fvfmretention.pdf", useDingbats=FALSE)


####APAL LT weight####
#pivot table to long format
apaltopivot<-dplyr::select(apal,Genotype,Treatment,Tank,IW,MW,EW)
apallong<-pivot_longer(data=apaltopivot,cols = -c(Genotype,Treatment,Tank), names_to = "time", values_to = "weight")

# Statistical tests for maximum quantum yield retained at end of assay 
# Analysis of variance (2-way) ANOVA dat (mixed linear model)
# Mixed Effects model origin_dest- random effect of colony nested in tank
# Since the model is completely balanced the SS type I, II or III will provide the same result
weight.model <- lmer(I.Erate ~ Treatment  * Genotype + (1|Tank) , data = apal)
weighttime.model <- lmer(weight ~ Treatment  * Genotype * time + (1|Tank) , data = apallong)
#summary(weight.model)
anova(weight.model)
rand(weight.model)
anova(weighttime.model)
rand(weighttime.model)

#confint(weight.model,oldNames=F)

# Model fitting and assumptions diagnostic 
x = residuals(weight.model)
plot(fitted(weight.model), x) 
leveneTest(x ~ Genotype * Treatment, data=apal, center=mean) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
plot(weight.model) # Residual vs Fitted values
qqnorm(x); qqline(x) # qq plot to check for normal distribution of residuals
hist(x) # histogram of residuals to check for normal distribution of residuals
shapiro.test(x) # formal statistical test (not recommended due the small sample size)
#apal outliers > row ,   row , remove
#apal<-apal[-c(11,106), ]

# statistical pairwise comparisons
# genotype
weight.emms <- emmeans(weight.model, pairwise ~ Genotype|Treatment, weights = "proportional", adjust="none")
summary(weight.emms$emmeans)
# P.value adjustment of the Bonferroni
rbind(weight.emms$contrasts, adjust="tukey")
# treatment
weight.emms <- emmeans(weight.model, pairwise ~ Treatment, weights = "proportional", adjust="none")
summary(weight.emms$emmeans)
rbind(weight.emms$contrasts, adjust="tukey")
#or
print(emmeans(weight.model, list(pairwise ~ Treatment)), adjust = c("tukey"))

###plot###
###points over time
apallong$grtime<-interaction(apallong$Treatment,apallong$time) 

sum<-Summarize(weight~Genotype+grtime, data=apallong, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$Treatment<-gl(4,5,labels=c("Control","OW","OA","OAOW"))
sum$Time<-gl(3,20,labels=c("End","Initial","Mid"))
sum$Time=factor(sum$Time,levels(sum$Time)[c(2,3,1)])
summary(sum)

#points with sd or ci 95% bars, comment out which you want
pd=position_dodge(.75)
plot.weight<-ggplot(sum,aes(x=Time,y=mean,color=Genotype,group=Genotype)) + 
  geom_point(size=3,position=pd,shape=17) +
  geom_line(position=pd) +
  scale_color_brewer(palette = "Paired") +
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=Genotype), width=.2, position=pd) + 
  ylab("Dry Weight (g)") + 
  xlab("Time point") +
  theme_bw() + 
  facet_wrap(~Treatment) +
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) 
ggsave(plot.weight, height = 6 , width = 10, filename = "./figures/20210701_apal_weightbytrt.pdf", useDingbats=FALSE)

#points as end growth rate
sum<-Summarize(I.Erate~Genotype+Treatment, data=apal, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$Treatment=factor(sum$Treatment,levels(sum$Treatment)[c(1,4,2,3)])

summary(sum)

plot.EIgrate<-ggplot(sum,aes(x=Genotype,y=mean,color=Genotype)) + 
  geom_point(size=3) +
  scale_color_brewer(palette = "Paired") +
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=Genotype), width=.2) + 
  ylab("Growth Rate (g dry weigh/day)") + 
  xlab("Genotype") +
  theme_bw() + 
  facet_wrap(~Treatment) +
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) 
ggsave(plot.EIgrate, height = 6 , width = 10, filename = "./figures/20210706_apal_endgrowthrate.pdf", useDingbats=FALSE)

####APAL LT FvFm####
#pivot table to long format, need to subset by trait of interest
apaltopivot<-dplyr::select(apal,Genotype,Treatment,Tank,IFvFm,MFvFm,EFvFm)
apallong<-pivot_longer(data=apaltopivot,cols = -c(Genotype,Treatment,Tank), names_to = "time", values_to = "fvfm")

# Statistical tests for maximum quantum yield retained at end of assay 
# Analysis of variance (2-way) ANOVA dat (mixed linear model)
# Mixed Effects model origin_dest- random effect of colony nested in tank
# Since the model is completely balanced the SS type I, II or III will provide the same result
fvfm.model <- lmer(EFvFm ~ Treatment  * Genotype + (1|Tank) , data = apal)
fvfmtime.model <- lmer(fvfm ~ Treatment  * Genotype * time + (1|Tank) , data = apallong)

#summary(fvfm.model)
anova(fvfm.model)
rand(fvfm.model)
#confint(fvfm.model,oldNames=F)

# Model fitting and assumptions diagnostic 
x = residuals(fvfm.model)
plot(fitted(fvfm.model), x) 
leveneTest(x ~ Genotype * Treatment, data=apal, center=mean) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
qqnorm(x); qqline(x) # qq plot to check for normal distribution of residuals
hist(x) # histogram of residuals to check for normal distribution of residuals
shapiro.test(x) # formal statistical test

# statistical pairwise comparisons
# genotype
fvfm.emms <- emmeans(fvfm.model, pairwise ~ Genotype|Treatment, weights = "proportional", adjust="none")
summary(fvfm.emms$emmeans)
# P.value adjustment of the Bonferroni
rbind(fvfm.emms$contrasts, adjust="tukey")
# treatment
fvfm.emms <- emmeans(fvfm.model, pairwise ~ Treatment, weights = "proportional", adjust="none")
summary(fvfm.emms$emmeans)
rbind(fvfm.emms$contrasts, adjust="tukey")
#or
print(emmeans(fvfm.model, list(pairwise ~ Treatment)), adjust = c("tukey"))

###plot###
###points over time
apallong$grtime<-interaction(apallong$Treatment,apallong$time) 

sum<-Summarize(fvfm~Genotype+grtime, data=apallong, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$Treatment<-gl(4,5,labels=c("Control","OW","OA","OAOW"))
sum$Time<-gl(3,20,labels=c("End","Initial","Mid"))
sum$Time=factor(sum$Time,levels(sum$Time)[c(2,3,1)])
summary(sum)

#points with sd or ci 95% bars, comment out which you want
pd=position_dodge(.5)
plot.fvfm<-ggplot(sum,aes(x=Time,y=mean,color=Genotype,group=Genotype)) + 
  geom_point(size=3,position=pd) +
  geom_line(position=pd) +
  scale_color_brewer(palette = "Paired") +
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=Genotype), width=.2, position=pd) + 
  ylab("Fv/Fm") + 
  xlab("Time point") +
  theme_bw() + 
  facet_wrap(~Treatment) +
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) 
ggsave(plot.fvfm, height = 6 , width = 10, filename = "./figures/20210706_apal_fvfmbytrt.pdf", useDingbats=FALSE)

#points as end Fv/Fm normalized to start
apal$mqy<-(1-(apal$IFvFm-apal$EFvFm)/apal$IFvFm)
sum<-Summarize(EFvFm~Genotype+Treatment, data=apal, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$Treatment=factor(sum$Treatment,levels(sum$Treatment)[c(1,4,2,3)])

summary(sum)

plot.endfvfm<-ggplot(sum,aes(x=Genotype,y=mean,color=Genotype)) + 
  geom_point(size=3) +
  scale_color_brewer(palette = "Paired") +
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=Genotype), width=.2) + 
  ylab("Fv/Fm Retained") + 
  xlab("Genotype") +
  theme_bw() + 
  facet_wrap(~Treatment) +
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) 
ggsave(plot.endfvfm, height = 6 , width = 10, filename = "./figures/20210706_apal_endfvfm.pdf", useDingbats=FALSE)

####Alluvial plot to show change in ranking among traits####
dat<-read.csv("202306_alluvialdata.csv", header=T)
dat<-transform(dat,
                 Genotype=as.factor(Genotype),
                 Treatment=as.factor(Treatment),
                 Tank=as.factor(Tank),
                 Experiment=as.factor(Experiment),
                 Exp=as.factor(Exp),
                 Region=as.factor(Region))
dat$deltFvFm<-as.numeric(dat$deltFvFm)
dat2 <-dplyr::select(dat,Genotype,Experiment,Tank, deltFvFm,HP,TotalChl)
LTdat<-dat2[dat2$Experiment=='LT',]
CBASSdat<-dat2[dat2$Experiment=='CBASS',]
##LT-OW##
#pivot table to long format
sum<-LTdat%>%
  group_by(Genotype) %>% 
  summarise_at(vars("deltFvFm", "TotalChl", "HP"), mean, na.rm=TRUE)
LTdatlong<-pivot_longer(data=sum,cols = -c(Genotype), names_to = "trait", values_to = "values")
#make sure traits and Genotypes stay in the specified order
col_order <- c("deltFvFm","TotalChl","HP")
LTdatlong$trait <- factor(LTdatlong$trait, levels = col_order)
#scale values by each trait level independently
LTdatlong <- transform(LTdatlong, ScaledTrait = ave(values, trait, FUN = scale))

##CBASS 34.3##
CBASSsum<-CBASSdat%>%
  group_by(Genotype) %>% 
  summarise_at(vars("deltFvFm", "TotalChl", "HP"), mean, na.rm=TRUE)
CBASSdatlong<-pivot_longer(data=CBASSsum,cols = -c(Genotype), names_to = "trait", values_to = "values")
#make sure traits and Genotypes stay in the specified order
col_order <- c("deltFvFm","TotalChl","HP")
CBASSdatlong$trait <- factor(CBASSdatlong$trait, levels = col_order)
#scale values by each trait level independently
CBASSdatlong <- transform(CBASSdatlong, ScaledTrait = ave(values, trait, FUN = scale))

LTalluvial<-ggplot(LTdatlong, aes(x = trait, stratum = Genotype, alluvium = Genotype, y = ScaledTrait, 
                      fill = Genotype, label = Genotype)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow(decreasing = TRUE) +
  geom_stratum(decreasing = TRUE,alpha = .5) +
  geom_text(stat = "stratum", size = 3,decreasing = TRUE) +
  theme_classic() +
  theme(legend.position = "right") +
  ylab("Scaled Trait Value") + 
  xlab("Trait")

CBASSalluvial<-ggplot(CBASSdatlong, aes(x = trait, stratum = Genotype, alluvium = Genotype, y = ScaledTrait, 
                                  fill = Genotype, label = Genotype)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow(decreasing = TRUE) +
  geom_stratum(decreasing = TRUE,alpha = .5) +
  geom_text(stat = "stratum", size = 3,decreasing = TRUE) +
  theme_classic() +
  theme(legend.position = "right") +
  ylab("Scaled Trait Value") + 
  xlab("Trait")

alluvial<-ggarrange(LTalluvial,CBASSalluvial, common.legend=TRUE, legend="right",labels= c("A. LT-OW", "B. CBASS 34.3°C"))
ggsave(alluvial,height = 6 , width = 10, filename ="20230610_alluvial-genetrankingbyexperiment.pdf", useDingbats=FALSE)
