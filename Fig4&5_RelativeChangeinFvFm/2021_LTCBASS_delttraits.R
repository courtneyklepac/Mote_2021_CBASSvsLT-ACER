library(tidyverse)
library(glmmTMB)    # for glmmTMB
library(sjPlot)     # for outputs and partial plots
library(car)        # for regression diagnostics
library(DHARMa)     # for residual diagnostics
library(effects)    # for partial effects plots
library(ggeffects)  # for effects plots in ggplot
library(multcomp)
library(dplyr)
library(MASS)
library(MuMIn)
library(rstatix)
library(emmeans)
library(lmtest)
library(multcompView)
library(patchwork)  # for arranging multiple plots
library(ggplot2)    # for graphing
library(ggfortify)  # for model diagnostics
library(gridExtra)
library(extrafont)  # for setting fonts in ggplot
loadfonts(device = "win") # needed to apply the fonts

trait$ID<-paste(trait$Genotype,trait$Tank,sep=".")
#####LTvsCBASS deltFvFm######
#Both experiments 
bartlett.test(deltFvFm~Genotype, data=trait) #Genotype equal variance, exp unequal
#welch one-way with unequal variance
oneway.test(deltFvFm~Exp, data=trait, var.equal=FALSE)
pairwise.wilcox.test(trait$deltFvFm, trait$Exp, p.adjust.method = "none", exact=FALSE)   
deltFvFm_lmer<-lmer(deltFvFm ~ Exp+ (1|Tank),data=trait)
anova(deltFvFm_lmer)
rand(deltFvFm_lmer)

#CBASS 34.3 and LT-OW
trait<-trait[trait$Exp=='LT-OW'| trait$Exp=='CBASS-34.3'  ,]
#CBASSheat<-trait[trait$Exp=='CBASS-34.3'  ,]
bartlett.test(deltFvFm~Exp, data=trait) #Genotype and exp equal
#welch one-way with unequal variance
#deltFvFmcont_welch<-oneway.test(deltFvFm~Exp, data=deltFvFmcont, var.equal=FALSE)
#pairwise.wilcox.test(trait$deltFvFm, trait$Exp, p.adjust.method = "bonf", exact=FALSE)    
deltFvFm_lmer<-lmer(deltFvFm ~ Exp*Region + (1|Rep) + (1|Tank),data=trait)
sjPlot::plot_model(deltFvFm_lmer, type="diag")
step(deltFvFm_lmer)
deltFvFm_lmer<-lmer(deltFvFm ~ Exp*Region + (1|Tank),data=trait)
anova(deltFvFm_lmer)
rand(deltFvFm_lmer)
print(emmeans(deltFvFm_lmer, list(pairwise ~Region|Exp)), adjust = c("tukey"))

#all CBASS trts
#welch one-way with unequal variance
bartlett.test(deltFvFm~Exp, data=CBASS) # exp unequal
deltFvFm_welch<-oneway.test(deltFvFm~Exp, data=CBASS)
pairwise.wilcox.test(CBASS$deltFvFm, CBASS$Exp, p.adjust.method = "bonf", exact=FALSE)   
deltFvFm_lmer<-lmer(deltFvFm ~ Treatment*Genotype + (1|Tank),data=CBASS)
#deltFvFm_lmer<-lmer(deltFvFm ~Genotype + (1|Tank),data=CBASSheat)

anova(deltFvFm_lmer)
ranova(deltFvFm_lmer)
print(emmeans(deltFvFm_lmer, list(pairwise ~Exp)), adjust = c("tukey"))

#LT Stats overall model testing Exp and genotype with tank as random
deltFvFm_lmer<-lmer(deltFvFm ~ Treatment+Genotype+Treatment*Genotype + (1|Tank) + (1|Rep),data=CBASS)
anova(deltFvFm_lmer)
rand(deltFvFm_lmer)
step(deltFvFm_lmer, reduce.random=TRUE)
#region
deltFvFm_lmer<-lmer(deltFvFm ~ Treatment*Region + (1|Tank),data=CBASS)
anova(deltFvFm_lmer)
rand(deltFvFm_lmer)
print(emmeans(deltFvFm_lmer, list(pairwise ~ Region|Treatment)), adjust = c("tukey"))

#or, to get letters of diff pairwise
leastsquare=lsmeans(deltFvFm_lmer,pairwise~Exp,adjust='tukey')
multcomp::cld(leastsquare,alpha=0.05,Letters=letters,adjust='tukey')

# Model fitting and assumptions diagnostic 
x = residuals(deltFvFm_lmer)
plot(fitted(deltFvFm_lmer), x) 
leveneTest(deltFvFm ~ Genotype * Exp, data=trait3, center=mean) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
bartlett.test(deltFvFm~Exp, data=trait3) #Genotype equal variance, exp unequal
qqnorm(x); qqline(x) # qq plot to check for normal distribution of residuals
hist(x) # histogram of residuals to check for normal distribution of residuals
shapiro.test(x) # formal statistical test
#boxcox doesn't like zeros, add .015 to all data values
bc<-boxcox(deltFvFm.I~Genotype*Time*Treatment, data=LTtime)
lambda <- bc$x[which.max(bc$y)]
#lambda=2
LTtime$deltFvFmsq<-LTtime$deltFvFm.I^2

#legend.position = c(1, .6),
#legend.text = element_text(size=10),
#legend.title = element_blank()
#summary for point plot
sum<-Summarize(deltFvFm~Genotype+Exp, data=trait, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$Region<-gl(2,5,label=c("Lower Keys","Upper Keys"))
sum$SubRegion<-gl(4,5,label=c("CBASS Lower Keys","CBASS Upper Keys","LT Lower Keys","LT Upper Keys"))
Regionmean <- sum %>% group_by(SubRegion) %>% dplyr::summarise(count = n(),mean_val=mean(mean),    
                                                        LB_deltFvFm = mean_val - qnorm(0.975) * sd(mean) / sqrt(count),
                                                        UB_deltFvFm = mean_val + qnorm(0.975) * sd(mean) / sqrt(count) )
mean_val<-c(-0.304,-0.256,-0.049,-0.115)
sum$Regionmean<-rep(mean_val, each=5)
RegionLB<-c(-0.356 ,-0.290 ,-0.0895 ,-0.156)
sum$RegionLB<-rep(RegionLB, each=5)
RegionUB<-c(-0.252,-0.222 ,-0.00853,-0.0732 )
sum$RegionUB<-rep(RegionUB, each=5)

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
        strip.background = element_rect(fill = "gray80",color="black",linewidth=1),
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

#####LTvsCBASS heat Host Proteins#####
trait1 <- trait %>% 
  group_by(Exp) %>% 
  identify_outliers(deltHP) %>%
  filter(is.outlier == "TRUE") #  4 outliers identified

trait2 <- trait %>% filter(!(ID=="UK76.C9")) # outliers removed
modelProt = glmmTMB(deltHP ~ Exp*Genotype + (1|Tank), data = trait2) # gaussian family unless specified
ModelProt.resid = simulateResiduals(modelProt, plot = TRUE) # something there isn't great
# so test for zero inflation
testZeroInflation(ModelProt.resid) # p value =1, not zero inflated.
# check for dispersion
testDispersion(ModelProt.resid) # not over/underdispersed either. 
# So issues likely due to complete lack of protein in last 2 timepoints. 
# partial plot
ggemmeans(modelProt, ~Genotype) %>% plot(add.data=TRUE) 
# R2 performance of model
performance::r2(modelProt) # 0.48
# Model outputs
summary(modelProt)
# post-hoc comparisons with emmeans 
emmeans(modelProt, ~ Genotype|Exp)%>%
  pairs()
###plots###

###points###
sum<-Summarize(deltHP~Genotype+Exp, data=trait2, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$Region<-gl(2,5,label=c("Lower Keys","Upper Keys"))
sum$SubRegion<-gl(4,5,label=c("CBASS Lower Keys","CBASS Upper Keys","LT Lower Keys","LT Upper Keys"))
Regionmean <- sum %>% group_by(SubRegion) %>% dplyr::summarise(count = n(),mean_val=mean(mean),    
                                                        LB_deltHP = mean_val - qnorm(0.975) * sd(mean) / sqrt(count),
                                                        UB_deltHP = mean_val + qnorm(0.975) * sd(mean) / sqrt(count) )
mean_val<-c(-0.0722,-0.209  ,-0.498,-0.581)
sum$Regionmean<-rep(mean_val, each=5)
RegionLB<-c(-0.439,-0.487,-0.676,-0.825 )
sum$RegionLB<-rep(RegionLB, each=5)
RegionUB<-c( 0.294,0.0694,-0.319,-0.337)
sum$RegionUB<-rep(RegionUB, each=5)

#points with sd or ci 95% bars, comment out which you want
# #FDAE61",'#F46D43',"#A50026"
plot.deltHP<-ggplot(sum,aes(x=Genotype,y=mean,color=Exp,shape=Region)) + 
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
ggsave(plot.deltHP, height = 6 , width = 7, bg='transparent', filename = "./20230324_CBASS34.3vsLTdeltHP-regionXExp.pdf", useDingbats=FALSE)

#####LTvsCBASS heat Chlorophylls#####
# Mixed Effects model origin_dest- random effect of colony nested in tank
deltChl.model <- lmer(deltChlsqrt ~ Genotype*Exp+ (1|Tank),data=trait)
deltChl.model <- lmer(deltChlsqrt ~ Genotype*Treatment + (1|Tank),data=CBASS)

anova(deltChl.model)
rand(deltChl.model)
step(deltChl.model, reduce.random=FALSE)
# statistical pairwise comparisons
print(emmeans(deltChl.model, list(pairwise ~ Exp)), adjust = c("tukey"))
print(emmeans(deltChl.model, list(pairwise ~ Treatment)), adjust = c("tukey"))

# Model fitting and assumptions diagnostic 
sjPlot::plot_model(deltChl.model, type="diag")
x = residuals(deltChl.model)
shapiro.test(x) # formal statistical test
leveneTest(deltChl ~ Genotype * Exp, data=LT, center=mean) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
#find optimal lambda for Box-Cox transformation 
#boxcox doesn't like zeros, add .015 to all data values
trait$deltChlpos <- trait$deltChl+0.015
bc<-boxcox(deltChlpos~Genotype*Exp, data=trait)
lambda <- bc$x[which.max(bc$y)]
#lambda=0.3
trait$deltChlsqrt <- sqrt(trait$deltChl)
CBASS$deltChlsqrt <- sqrt(CBASS$deltChl)
trait$deltChlsqrt <- sqrt(trait$deltChl)

#welchs one-way since Exp is unequal var
bartlett.test(deltChlsqrt~Exp, data=CBASS) #Exp unequal geno equal
oneway.test(deltChlsqrt~Exp, data=CBASS)
pairwise.wilcox.test(CBASS$deltChlsqrt, CBASS$Exp, p.adjust.method = "bonf", exact=FALSE)            

#CBASS 34.3 and 35.6 and LT-OW
bartlett.test(deltChlsqrt~Exp, data=trait) #Equal variance with deltChlsqrt
deltChl_lmer<-lmer(deltChlsqrt ~ Exp*Genotype+ (1|Tank),data=trait)
sjPlot::plot_model(deltChl_lmer, type="diag")
step(deltChl_lmer)
anova(deltChl_lmer)
rand(deltChl_lmer)
# statistical pairwise comparisons
print(emmeans(deltChl_lmer, list(pairwise ~ Exp|Region)), adjust = c("tukey"))

###plots###
###points###
sum<-Summarize(deltChl~Genotype+Exp, data=trait, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$Region<-gl(2,5,label=c("Lower Keys","Upper Keys"))
sum$SubRegion<-gl(4,5,label=c("CBASS Lower Keys","CBASS Upper Keys","LT Lower Keys","LT Upper Keys"))
Regionmean <- sum %>% group_by(SubRegion) %>% dplyr::summarise(count = n(),mean_val=mean(mean),    
                                                        LB_deltChl = mean_val - qnorm(0.975) * sd(mean) / sqrt(count),
                                                        UB_deltChl = mean_val + qnorm(0.975) * sd(mean) / sqrt(count) )
mean_val<-c(-0.447,-0.430,-0.885,-0.96)
sum$Regionmean<-rep(mean_val, each=5)
RegionLB<-c(-0.526,-0.630,-0.908,-0.972 )
sum$RegionLB<-rep(RegionLB, each=5)
RegionUB<-c(-0.367,-0.229,-0.861,-0.948)
sum$RegionUB<-rep(RegionUB, each=5)

#points with sd or ci 95% bars, comment out which you want
# #FDAE61",'#F46D43',"#A50026"
plot.deltChl<-ggplot(sum,aes(x=Genotype,y=mean,color=Exp,shape=Region)) + 
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
ggsave(plot.deltChl, height = 6 , width = 7, bg='transparent', filename = "./20230310_CBASS34.3vsLTdeltChl-regionXExp.pdf", useDingbats=FALSE)

