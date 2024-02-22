##Statistical analyses for PAM at T1
library(lmerTest)
library(emmeans)
library(sjPlot)
library(drc)
library(ggplot2)
library(Rmisc)
library(brms)
library(tidyr)

#need to set wd to source file location
#setwd("~/Dropbox/2021 State funded Project (Apal and Acer)/CBASS/PAM")
#Load file
pam <- read.csv("20221114_2021State_MasterSheet.csv", header = T, na.strings=c("","NA"))
pam<-read.csv("FvFmheatloading_CBASS.csv")
pamOW<-read.csv("FvFmheatloading_LT.csv")
#set factor levels
pam<-transform(pam,
                 Exp=as.factor(Exp),
                 Treatment=as.factor(Treatment),
                 Genotype=as.factor(Genotype),
                 Tank=as.factor(Tank),
                 Heatload30.5=as.factor(Heatload30.5))

levels(pam$Genotype)<-c("LK31","LK41","LK50","LK62","LK7", "UK12", "UK19", "UK70","UK76","UK80")
pam$Exp=factor(pam$Exp,levels(pam$Exp)[c(9,1,2,3,4,5,6,7,8,10)])
pam$Treatment=factor(pam$Treatment,levels(pam$Treatment)[c(9,1,2,3,4,5,6,7,8,10)])
pam<-subset(pam,Experiment=="CBASS")
pam$Heatload30.5<-as.character(pam$Heatload30.5)
pam$Heatload30.5<-  as.numeric(pam$Heatload30.5)
pam<-na.omit(pam)

pamOW<-transform(pamOW,
               Genotype=as.factor(Genotype),
               Tank=as.factor(Tank),
               Exp=as.factor(Exp),
               RealTemp=as.factor(RealTemp),
               Region=as.factor(Region),
               Heatload30.5=as.factor(Heatload30.5))
pamOW$Genotype=factor(pamOW$Genotype,levels(pamOW$Genotype)[c(4,5,6,7,1,2,3,8,9,10)])
levels(pamOW$Genotype)<-c("LK31","LK41","LK50","LK62","LK7", "UK12", "UK19", "UK70","UK76","UK80")
pamOW$Heatload30.5<-as.character(pamOW$Heatload30.5)
pamOW$Heatload30.5<-  as.numeric(pamOW$Heatload30.5)
pamOW<-subset(pamOW,Exp=="OW")
pamOW<-na.omit(pamOW)
#Subset data
#pam1<-subset(pam, Rep=="1")
#pam2<-subset(pam, Rep=="2")
#pam3<-subset(pam, Rep=="3")

############################################################################################################
###### Compute and plot Fv/Fm ED50 for CBASS Genotypes
### Treatment as numeric for drc fitting
pam$Temp<-as.character(pam$Temp)
pam$Temp<-as.numeric(pam$Temp)

pam
####CBASS-Generate curves and ED50s ####
DRC <- drm(FvFm ~ Treatment, data = pam, curveid = Genotype, fct = LL.3(names = c('hill', 'max', 'ed50')))
mselect(DRC, list(LL.2(), LL.4(), LL.5(), LL2.2(), LL2.3(), LL2.4(), LL2.5(), AR.2(), AR.3(), EXD.2(), EXD.3()), icfct = AIC)
modelFit(DRC)
summary(DRC)
plot(DRC)
DRC_coeff<-data.frame(ED(DRC, c(50)))

#extract ED50s for each genet
geno_coeff<-DRC$coefficients[21:30]
coeffs<-data.frame(geno_coeff)
coeffs[,'geno_coeff']=round(coeffs[,'geno_coeff'],2) #round to two decimals
coeffs<-rownames_to_column(coeffs)

pam<-pamCBASS

#Run full model
CBASS <- drm(FvFm ~ Temp, data = pam, fct = LL.3()) #34.37
summary(CBASS)
plot(CBASS)

#overall DRC predictors
CBASS_preddata = data.frame(Temp = seq(27,38, length.out = 100))
CBASS_pred = as.data.frame(predict(CBASS, newdata = CBASS_preddata, interval = 'confidence'))
CBASS_preddata = data.frame(CBASS_preddata, FvFm = CBASS_pred$Prediction, Lower = CBASS_pred$Lower, Upper = CBASS_pred$Upper)

#### PLOT  DRC for experiment
DRC_plot<- ggplot() +
  geom_jitter(data = pam, aes(x = Temp, y = FvFm), size = 1, width = 0.25) +
  scale_x_continuous(limits=c(26,38), breaks=c(27,29,31,33,35,37)) +
  scale_y_continuous(limits=c(-0.1, 0.7), breaks=c(0, 0.2, 0.4, 0.6)) +
  
  geom_line(data = CBASS_preddata, aes(x = Temp, y = FvFm), color = 'blue3', show.legend = FALSE) +
  geom_ribbon(data = CBASS_preddata, aes(x = Temp, ymin=Lower, ymax=Upper), color = 'blue3', linetype=2, alpha = 0.2) +
  geom_vline(aes(xintercept = 34.45), color = 'blue3', show.legend = FALSE) +
  #geom_text(aes(label = '34.37', x = 35.2, y = 0.6), color='blue') +
  
  ylab(label = expression(paste(italic(F[V]/F[M])))) +
  xlab("Temperature (°C)") +
  theme_classic()

ggsave(DRC_plot, height = 6 , width = 10, filename = "20221206_CBASS_nurseryDRC.pdf", useDingbats=FALSE)


####CBASS-Run model for each geno####
LK31drc <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="LK31",], fct = LL.3())
summary(LK31drc)
plot(LK31drc)

LK41drc <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="LK41",], fct = LL.3())
summary(LK41drc)
plot(LK41drc)

LK50drc <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="LK50",], fct = LL.3())
summary(LK50drc)
plot(LK50drc)

LK62drc <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="LK62",], fct = LL.3())
summary(LK62drc)
plot(LK62drc)

LK7drc <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="LK7",], fct = LL.3())
summary(LK7drc)
plot(LK7drc)

UK12drc <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="UK12",], fct = LL.3())
summary(UK12drc)
plot(UK12drc)

UK19drc <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="UK19",], fct = LL.3())
summary(UK19drc)
plot(UK19drc)

UK70drc <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="UK70",], fct = LL.3())
summary(UK70drc)
plot(UK70drc)

UK76drc <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="UK76",], fct = LL.3())
summary(UK76drc)
plot(UK76drc)

UK80drc <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="UK80",], fct = LL.3())
summary(UK80drc)
plot(UK80drc)

#### Combine ED50 data plus predict curves from models for plotting ####
#Genotype predictors
geno_preddata = data.frame(Temp = seq(27,38, length.out = 100))
LK31_pred = as.data.frame(predict(LK31drc, newdata = geno_preddata, interval = 'confidence'))
LK31_preddata = data.frame(geno_preddata, Y = LK31_pred$Prediction, Lower = LK31_pred$Lower, Upper = LK31_pred$Upper)

LK41_pred = as.data.frame(predict(LK41drc, newdata = geno_preddata, interval = 'confidence'))
LK41_preddata = data.frame(geno_preddata, Y = LK41_pred$Prediction, Lower = LK41_pred$Lower, Upper = LK41_pred$Upper)

LK50_pred = as.data.frame(predict(LK50drc, newdata = geno_preddata, interval = 'confidence'))
LK50_preddata = data.frame(geno_preddata, Y = LK50_pred$Prediction, Lower = LK50_pred$Lower, Upper = LK50_pred$Upper)

LK62_pred = as.data.frame(predict(LK62drc, newdata = geno_preddata, interval = 'confidence'))
LK62_preddata = data.frame(geno_preddata, Y = LK62_pred$Prediction, Lower = LK62_pred$Lower, Upper = LK62_pred$Upper)

LK7_pred = as.data.frame(predict(LK7drc, newdata = geno_preddata, interval = 'confidence'))
LK7_preddata = data.frame(geno_preddata, Y = LK7_pred$Prediction, Lower = LK7_pred$Lower, Upper = LK7_pred$Upper)

UK12_pred = as.data.frame(predict(UK12drc, newdata = geno_preddata, interval = 'confidence'))
UK12_preddata = data.frame(geno_preddata, Y = UK12_pred$Prediction, Lower = UK12_pred$Lower, Upper = UK12_pred$Upper)

UK19_pred = as.data.frame(predict(UK19drc, newdata = geno_preddata, interval = 'confidence'))
UK19_preddata = data.frame(geno_preddata, Y = UK19_pred$Prediction, Lower = UK19_pred$Lower, Upper = UK19_pred$Upper)

UK70_pred = as.data.frame(predict(UK70drc, newdata = geno_preddata, interval = 'confidence'))
UK70_preddata = data.frame(geno_preddata, Y = UK70_pred$Prediction, Lower = UK70_pred$Lower, Upper = UK70_pred$Upper)

UK76_pred = as.data.frame(predict(UK76drc, newdata = geno_preddata, interval = 'confidence'))
UK76_preddata = data.frame(geno_preddata, Y = UK76_pred$Prediction, Lower = UK76_pred$Lower, Upper = UK76_pred$Upper)

UK80_pred = as.data.frame(predict(UK80drc, newdata = geno_preddata, interval = 'confidence'))
UK80_preddata = data.frame(geno_preddata, Y = UK80_pred$Prediction, Lower = UK80_pred$Lower, Upper = UK80_pred$Upper)

####CBASS-Plot DRC for each geno####
library(scales)
hex<-hue_pal()(10)
hex #"#F8766D" "#D89000" "#A3A500" "#39B600" "#00BF7D" "#00BFC4" "#00B0F6" "#9590FF" "#E76BF3" "#FF62BC"
DRC_plot<- ggplot() +
  geom_jitter(data = pam, aes(x = Treatment, y = FvFm, color=Genotype), size = 1, width = 0.25) +
  scale_x_continuous(limits=c(26,38), breaks=c(27,29,31,33,35,37)) +
  scale_y_continuous(limits=c(-0.1, 0.7), breaks=c(0, 0.2, 0.4, 0.6)) +
  
  geom_line(data = LK31_preddata, aes(x = Treatment, y = FvFm), color = '#F8766D', show.legend = FALSE) +
  geom_ribbon(data = LK31_preddata, aes(x = Treatment, ymin=Lower, ymax=Upper), color = '#F8766D', linetype=2, alpha = 0.2) +
  geom_vline(data = coeffs, aes(xintercept = geno_coeff[1]), color = '#F8766D', show.legend = FALSE) +
  geom_text(data = coeffs, aes(label = "LK31: 34.00"), x = 36, y = 0.33, show.legend = FALSE, color = '#F8766D', hjust='left') +

  geom_line(data = LK41_preddata, aes(x = Treatment, y = FvFm), color = '#D89000', show.legend = FALSE) +
  geom_vline(data = coeffs, aes(xintercept = geno_coeff[2]), color = '#D89000', show.legend = FALSE) +
  geom_text(data = coeffs, aes(label = "LK41:34.44"), x = 36, y = 0.49, show.legend = FALSE, color = '#D89000', hjust='left') +
  
  geom_line(data = LK50_preddata, aes(x = Treatment, y = FvFm), color = '#A3A500', show.legend = FALSE) +
  geom_vline(data = coeffs, aes(xintercept = geno_coeff[3]), color = '#A3A500', show.legend = FALSE) +
  geom_text(data = coeffs, aes(label = "LK50: 34.04"), x = 36, y = 0.37, show.legend = FALSE, color = '#A3A500', hjust='left') +
 
  geom_line(data = LK62_preddata, aes(x = Treatment, y = FvFm), color = '#39B600', show.legend = FALSE) +
  geom_vline(data = coeffs, aes(xintercept = geno_coeff[4]), color = '#39B600', show.legend = FALSE) +
  geom_text(data = coeffs, aes(label = "LK62: 34.15"), x = 36, y = 0.41, show.legend = FALSE, color = '#39B600', hjust='left') +

  geom_line(data = LK7_preddata, aes(x = Treatment, y = FvFm), color = '#00BF7D', show.legend = FALSE) +
  geom_vline(data = coeffs, aes(xintercept = geno_coeff[5]), color = '#00BF7D', show.legend = FALSE) +
  geom_text(data = coeffs, aes(label = "LK7: 34.58"), x = 36, y = 0.61, show.legend = FALSE, color = '#00BF7D', hjust='left') +
  
  geom_line(data = UK12_preddata, aes(x = Treatment, y = FvFm), color = '#00BFC4', show.legend = FALSE) +
  geom_vline(data = coeffs, aes(xintercept = geno_coeff[6]), color = '#00BFC4', show.legend = FALSE) +
  geom_text(data = coeffs, aes(label = "UK12: 34.72"), x = 36, y = 0.69, show.legend = FALSE, color = '#00BFC4', hjust='left') +
  
  geom_line(data = UK19_preddata, aes(x = Treatment, y = FvFm), color = '#00B0F6', show.legend = FALSE) +
  geom_vline(data = coeffs, aes(xintercept = geno_coeff[7]), color = '#00B0F6', show.legend = FALSE) +
  geom_text(data = coeffs, aes(label = "UK19: 34.62"), x = 36, y = 0.65, show.legend = FALSE, color = '#00B0F6', hjust='left') +
  
  geom_line(data = UK70_preddata, aes(x = Treatment, y = FvFm), color = '#9590FF', show.legend = FALSE) +
  geom_vline(data = coeffs, aes(xintercept = geno_coeff[8]), color = '#9590FF', show.legend = FALSE) +
  geom_text(data = coeffs, aes(label = "UK70: 34.54"), x = 36, y = 0.57, show.legend = FALSE, color = '#9590FF', hjust='left') +
  
  geom_line(data = UK76_preddata, aes(x = Treatment, y = FvFm), color = '#E76BF3', show.legend = FALSE) +
  geom_vline(data = coeffs, aes(xintercept = geno_coeff[9]), color = '#E76BF3', show.legend = FALSE) +
  geom_text(data = coeffs, aes(label = "UK76: 34.47"), x = 36, y = 0.53, show.legend = FALSE, color = '#E76BF3', hjust='left') +
  
  geom_line(data = UK80_preddata, aes(x = Treatment, y = FvFm), color = '#FF62BC', show.legend = FALSE) +
  geom_vline(data = coeffs, aes(xintercept = geno_coeff[10]), color = '#FF62BC', show.legend = FALSE) +
  geom_text(data = coeffs, aes(label = "UK80: 34.23"), x = 36, y = 0.45, show.legend = FALSE, color = '#FF62BC', hjust='left') +
  
  ylab(label = expression(paste(italic(F[V]/F[M])))) +
  xlab("Temperature (°C)") +
  theme_classic()+ theme(legend.position = "none")

DRC_plot

ggsave(DRC_plot, height = 6 , width = 10, filename = "20221206_CBASS_genoDRCs.pdf", useDingbats=FALSE)

####model testing####
#compare ED50 across region
coeffs$Region<-factor(c("LowerKeys","LowerKeys","LowerKeys","LowerKeys","LowerKeys","UpperKeys","UpperKeys","UpperKeys","UpperKeys","UpperKeys"))
bartlett.test(geno_coeff~Region, data=coeffs) 
t.test(geno_coeff~Region, data=coeffs, var.equal = TRUE)
#ED25
geno_drcLT_coeff_long$Region<-factor(c("LowerKeys","LowerKeys","LowerKeys","LowerKeys","LowerKeys","UpperKeys","UpperKeys","UpperKeys","UpperKeys","UpperKeys"))
geno_drcLT_coeff_long$Genotype<-factor(c("LK31","LK41","LK50","LK62","LK7","UK12","UK19","UK70","UK76","UK80"))
t.test(geno_coeff~Region, data=geno_drc_coeff_long, var.equal = TRUE)

#Stats model
pam$Treatment<-as.factor(pam$Treatment)

pamCBASS<-lmer(FvFm ~ Treatment*Genotype + (1|Tank) , data=pam)
sjPlot::plot_model(pamCBASS, type="diag")
step(pamCBASS, reduce.random=FALSE)
# Model fitting and assumptions diagnostic 
leveneTest(FvFm ~ Genotype*Treatment, data=pam, center=mean) # Genotype=.9488 Treatment=2.2e-16 int=1.17e-16
shapiro.test(x) # formal statistical test
anova(pamCBASS)
rand(pamCBASS)
print(emmeans(pamCBASS, list(pairwise ~ Genotype)), adjust = c("tukey"))

####compare pam CBASS vs. long term####
pam34<-pam[pam$Treatment=='34.3',]
pam34<-subset(pam34, select=c(2:5))
#names(pam34)[2]<-"EFvFm"
pamOW<-pamOW[pamOW$Exp=='OW',]
pamOW<-pamOW[pamOW$Day=='51',]
pamOW<-subset(pamOW, select=c(3:6))
spoing<-rbind(pam34, pamOW)
levels(spoing$Exp)[levels(spoing$Exp)=='OW'] <- 'LT'

####LT-ED25 DRC-RealTemp####
DRC <- drm(FvFm ~ RealTemp, data = pamOW, curveid = Genotype, fct = LL.4())
mselect(DRC, list(LL.2(), LL.4(), LL.5(), LL2.2(), LL2.3(), LL2.4(), LL2.5(), AR.2(), AR.3(), EXD.2(), EXD.3()), icfct = AIC)
modelFit(DRC)
summary(DRC)
plot(DRC)

#extract coeffs
DRC25_coeff<-data.frame(ED(DRC, c(25)))

DRC25_coeff_mean<-DRC25_coeff[,1]
DRC25_coeff_lower<-DRC25_coeff[,1] - 1.96*DRC25_coeff[,2]
DRC25_coeff_upper<-DRC25_coeff[,1] + 1.96*DRC25_coeff[,2]

#Run full model
LTdrc <- drm(FvFm ~ RealTemp, data = pamOW, fct = LL.4()) 
summary(LTdrc)
plot(LTdrc)
LT_coeff<-ED(LTdrc, c(25))#31.43

#overall drcLT predictors
LT_preddata = data.frame(RealTemp = seq(27,38, length.out = 100))
LT_pred = as.data.frame(predict(LTdrc, newdata = LT_preddata, interval = 'confidence'))
LT_preddata = data.frame(LT_preddata, FvFm = LT_pred$Prediction, Lower = LT_pred$Lower, Upper = LT_pred$Upper)

#### PLOT ED25 DRC for LT
DRC_plot<- ggplot() +
  geom_jitter(data = pamOW, aes(x = RealTemp, y = FvFm), size = 1, width = 0.25) +
  scale_x_continuous(limits=c(26,33), breaks=c(27,29,31,33)) +
  scale_y_continuous(limits=c(0, 0.7), breaks=c(0, 0.2, 0.4, 0.6)) +
  
  geom_line(data = LT_preddata, aes(x = RealTemp, y = FvFm), color = 'blue3', show.legend = FALSE) +
  geom_ribbon(data =LT_preddata, aes(x = RealTemp, ymin=Lower, ymax=Upper), color = 'blue3', linetype=2, alpha = 0.2) +
  geom_vline(aes(xintercept = 31.43), color = 'blue3', show.legend = FALSE) +
  geom_text(aes(label = 'LT ED25:\n31.43', x = 32.2, y = 0.6), color='blue') +
  
  ylab(label = expression(paste(italic(F[V]/F[M])))) +
  xlab("Temperature (°C)") +
  theme_classic()

ggsave(DRC_plot, height = 6 , width = 10, filename = "20230320_LT_ED25DRC.pdf", useDingbats=FALSE)

####LT-Run model for each geno####
#### LK31 LT ED25 ####
#model
LK31drcLT  <- drm(FvFm ~ RealTemp, data = pamOW[pamOW$Genotype=="LK31",], fct = LL.4())
summary(LK31drcLT)
plot(LK31drcLT)
#extract coeffs 31.04732
LK31drcLT_coeff<-data.frame(ED(LK31drcLT, c(25)))

LK31drcLT_coeff_mean<-LK31drcLT_coeff[,1]
LK31drcLT_coeff_lower<-LK31drcLT_coeff[,1] - 1.96*LK31drcLT_coeff[,2]
LK31drcLT_coeff_upper<-LK31drcLT_coeff[,1] + 1.96*LK31drcLT_coeff[,2]


#### LK41 LT ED25 ####
#model
LK41drcLT  <- drm(FvFm ~ RealTemp, data = pamOW[pamOW$Genotype=="LK41",], fct = LL.4())
summary(LK41drcLT)
plot(LK41drcLT)
#extract coeffs 31.10161
LK41drcLT_coeff<-data.frame(ED(LK41drcLT, c(25)))

LK41drcLT_coeff_mean<-LK41drcLT_coeff[,1]
LK41drcLT_coeff_lower<-LK41drcLT_coeff[,1] - 1.96*LK41drcLT_coeff[,2]
LK41drcLT_coeff_upper<-LK41drcLT_coeff[,1] + 1.96*LK41drcLT_coeff[,2]


#### LK50 LT ED25 ####
#model
LK50drcLT  <- drm(FvFm ~ RealTemp, data = pamOW[pamOW$Genotype=="LK50",], fct = LL.4())
summary(LK50drcLT)
plot(LK50drcLT)
#extract coeffs 31.0905
LK50drcLT_coeff<-data.frame(ED(LK50drcLT, c(25)))

LK50drcLT_coeff_mean<-LK50drcLT_coeff[,1]
LK50drcLT_coeff_lower<-LK50drcLT_coeff[,1] - 1.96*LK50drcLT_coeff[,2]
LK50drcLT_coeff_upper<-LK50drcLT_coeff[,1] + 1.96*LK50drcLT_coeff[,2]


#### LK62 LT ED25 ####
#model
LK62drcLT  <- drm(FvFm ~ RealTemp, data = pamOW[pamOW$Genotype=="LK62",], fct = LL.4())
summary(LK62drcLT)
plot(LK62drcLT)
#extract coeffs 31.298
LK62drcLT_coeff<-data.frame(ED(LK62drcLT, c(25)))

LK62drcLT_coeff_mean<-LK62drcLT_coeff[,1]
LK62drcLT_coeff_lower<-LK62drcLT_coeff[,1] - 1.96*LK62drcLT_coeff[,2]
LK62drcLT_coeff_upper<-LK62drcLT_coeff[,1] + 1.96*LK62drcLT_coeff[,2]


#### LK7 LT ED25 ####
#model
LK7drcLT  <- drm(FvFm ~ RealTemp, data = pamOW[pamOW$Genotype=="LK7",], fct = LL.4())
summary(LK7drcLT)
plot(LK7drcLT)
#extract coeffs 30.97742
LK7drcLT_coeff<-data.frame(ED(LK7drcLT, c(25)))

LK7drcLT_coeff_mean<-LK7drcLT_coeff[,1]
LK7drcLT_coeff_lower<-LK7drcLT_coeff[,1] - 1.96*LK7drcLT_coeff[,2]
LK7drcLT_coeff_upper<-LK7drcLT_coeff[,1] + 1.96*LK7drcLT_coeff[,2]


#### UK12 LT ED25 ####
#model
UK12drcLT  <- drm(FvFm ~ RealTemp, data = pamOW[pamOW$Genotype=="UK12",], fct = LL.4())
summary(UK12drcLT)
plot(UK12drcLT)
#extract coeffs 31.3669
UK12drcLT_coeff<-data.frame(ED(UK12drcLT, c(25)))

UK12drcLT_coeff_mean<-UK12drcLT_coeff[,1]
UK12drcLT_coeff_lower<-UK12drcLT_coeff[,1] - 1.96*UK12drcLT_coeff[,2]
UK12drcLT_coeff_upper<-UK12drcLT_coeff[,1] + 1.96*UK12drcLT_coeff[,2]


#### UK19 LT ED25 ####
#model
UK19drcLT  <- drm(FvFm ~ RealTemp, data = pamOW[pamOW$Genotype=="UK19",], fct = LL.4())
summary(UK19drcLT)
plot(UK19drcLT)
#extract coeffs 31.01445 
UK19drcLT_coeff<-data.frame(ED(UK19drcLT, c(25)))

UK19drcLT_coeff_mean<-UK19drcLT_coeff[,1]
UK19drcLT_coeff_lower<-UK19drcLT_coeff[,1] - 1.96*UK19drcLT_coeff[,2]
UK19drcLT_coeff_upper<-UK19drcLT_coeff[,1] + 1.96*UK19drcLT_coeff[,2]


#### UK70 LT ED25 ####
#model
UK70drcLT  <- drm(FvFm ~ RealTemp, data = pamOW[pamOW$Genotype=="UK70",], fct = LL.4())
summary(UK70drcLT)
plot(UK70drcLT)
#extract coeffs 31.02729
UK70drcLT_coeff<-data.frame(ED(UK70drcLT, c(25)))

UK70drcLT_coeff_mean<-UK70drcLT_coeff[,1]
UK70drcLT_coeff_lower<-UK70drcLT_coeff[,1] - 1.96*UK70drcLT_coeff[,2]
UK70drcLT_coeff_upper<-UK70drcLT_coeff[,1] + 1.96*UK70drcLT_coeff[,2]


#### UK76 LT ED25 ####
#model
UK76drcLT  <- drm(FvFm ~ RealTemp, data = pamOW[pamOW$Genotype=="UK76",], fct = LL.4())
summary(UK76drcLT)
plot(UK76drcLT)
#extract coeffs 31.149
UK76drcLT_coeff<-data.frame(ED(UK76drcLT, c(25)))

UK76drcLT_coeff_mean<-UK76drcLT_coeff[,1]
UK76drcLT_coeff_lower<-UK76drcLT_coeff[,1] - 1.96*UK76drcLT_coeff[,2]
UK76drcLT_coeff_upper<-UK76drcLT_coeff[,1] + 1.96*UK76drcLT_coeff[,2]


#### UK80 LT ED25 ####
#model
UK80drcLT  <- drm(FvFm ~ RealTemp, data = pamOW[pamOW$Genotype=="UK80",], fct = LL.4())
summary(UK80drcLT)
plot(UK80drcLT)
#extract coeffs 30.86556
UK80drcLT_coeff<-data.frame(ED(UK80drcLT, c(25)))

UK80drcLT_coeff_mean<-UK80drcLT_coeff[,1]
UK80drcLT_coeff_lower<-UK80drcLT_coeff[,1] - 1.96*UK80drcLT_coeff[,2]
UK80drcLT_coeff_upper<-UK80drcLT_coeff[,1] + 1.96*UK80drcLT_coeff[,2]

####LT- geno ED25 predictors and confidence intervals####
geno_drcLT_coeff_means<- data.frame(LK31drcLT_coeff_mean, LK41drcLT_coeff_mean, LK50drcLT_coeff_mean, LK62drcLT_coeff_mean, LK7drcLT_coeff_mean, UK12drcLT_coeff_mean, UK19drcLT_coeff_mean, UK70drcLT_coeff_mean, UK76drcLT_coeff_mean, UK80drcLT_coeff_mean)
geno_drcLT_coeff_long<-pivot_longer(cols=everything(),geno_drcLT_coeff_means, names_to="name",values_to="geno_coeff")
geno_drcLT_coeff_long$name<-as.factor(geno_drcLT_coeff_long$name) 
geno_drcLT_coeff_long<- subset(geno_drcLT_coeff_long, select = -c(name))
geno_drcLT_coeff_long[,'geno_coeff']=round(geno_drcLT_coeff_long[,'geno_coeff'],2) #round to two decimals

geno_drcLT_coeff_lowers<-data.frame(LK31drcLT_coeff_lower,LK41drcLT_coeff_lower,LK50drcLT_coeff_lower,LK62drcLT_coeff_lower,LK7drcLT_coeff_lower,UK12drcLT_coeff_lower,UK19drcLT_coeff_lower,UK70drcLT_coeff_lower,UK76drcLT_coeff_lower,UK80drcLT_coeff_lower )
geno_drcLT_coeff_uppers<-data.frame(LK31drcLT_coeff_upper,LK41drcLT_coeff_upper,LK50drcLT_coeff_upper,LK62drcLT_coeff_upper,LK7drcLT_coeff_upper,UK12drcLT_coeff_upper,UK19drcLT_coeff_upper,UK70drcLT_coeff_upper,UK76drcLT_coeff_upper,UK80drcLT_coeff_upper )

LK31drcLT_preddata = data.frame(temp = seq(26,38, length.out = 100))
LK31drcLT_pred = as.data.frame(predict(LK31drcLT, newdata = LK31drcLT_preddata, interval = 'confidence'))
LK31drcLT_preddata = data.frame(LK31drcLT_preddata, fvfm = LK31drcLT_pred$Prediction, Lower = LK31drcLT_pred$Lower, Upper = LK31drcLT_pred$Upper)

LK41drcLT_preddata = data.frame(temp = seq(26,38, length.out = 100))
LK41drcLT_pred = as.data.frame(predict(LK41drcLT, newdata = LK41drcLT_preddata, interval = 'confidence'))
LK41drcLT_preddata = data.frame(LK41drcLT_preddata, fvfm = LK41drcLT_pred$Prediction, Lower = LK41drcLT_pred$Lower, Upper = LK41drcLT_pred$Upper)

LK50drcLT_preddata = data.frame(temp = seq(26,38, length.out = 100))
LK50drcLT_pred = as.data.frame(predict(LK50drcLT, newdata = LK50drcLT_preddata, interval = 'confidence'))
LK50drcLT_preddata = data.frame(LK50drcLT_preddata, fvfm = LK50drcLT_pred$Prediction, Lower = LK50drcLT_pred$Lower, Upper = LK50drcLT_pred$Upper)

LK62drcLT_preddata = data.frame(temp = seq(26,38, length.out = 100))
LK62drcLT_pred = as.data.frame(predict(LK62drcLT, newdata = LK62drcLT_preddata, interval = 'confidence'))
LK62drcLT_preddata = data.frame(LK62drcLT_preddata, fvfm = LK62drcLT_pred$Prediction, Lower = LK62drcLT_pred$Lower, Upper = LK62drcLT_pred$Upper)

LK7drcLT_preddata = data.frame(temp = seq(26,38, length.out = 100))
LK7drcLT_pred = as.data.frame(predict(LK7drcLT, newdata = LK7drcLT_preddata, interval = 'confidence'))
LK7drcLT_preddata = data.frame(LK7drcLT_preddata, fvfm = LK7drcLT_pred$Prediction, Lower = LK7drcLT_pred$Lower, Upper = LK7drcLT_pred$Upper)

UK12drcLT_preddata = data.frame(temp = seq(26,38, length.out = 100))
UK12drcLT_pred = as.data.frame(predict(UK12drcLT, newdata = UK12drcLT_preddata, interval = 'confidence'))
UK12drcLT_preddata = data.frame(UK12drcLT_preddata, fvfm = UK12drcLT_pred$Prediction, Lower = UK12drcLT_pred$Lower, Upper = UK12drcLT_pred$Upper)

UK19drcLT_preddata = data.frame(temp = seq(26,38, length.out = 100))
UK19drcLT_pred = as.data.frame(predict(UK19drcLT, newdata = UK19drcLT_preddata, interval = 'confidence'))
UK19drcLT_preddata = data.frame(UK19drcLT_preddata, fvfm = UK19drcLT_pred$Prediction, Lower = UK19drcLT_pred$Lower, Upper = UK19drcLT_pred$Upper)

UK70drcLT_preddata = data.frame(temp = seq(26,38, length.out = 100))
UK70drcLT_pred = as.data.frame(predict(UK70drcLT, newdata = UK70drcLT_preddata, interval = 'confidence'))
UK70drcLT_preddata = data.frame(UK70drcLT_preddata, fvfm = UK70drcLT_pred$Prediction, Lower = UK70drcLT_pred$Lower, Upper = UK70drcLT_pred$Upper)

UK76drcLT_preddata = data.frame(temp = seq(26,38, length.out = 100))
UK76drcLT_pred = as.data.frame(predict(UK76drcLT, newdata = UK76drcLT_preddata, interval = 'confidence'))
UK76drcLT_preddata = data.frame(UK76drcLT_preddata, fvfm = UK76drcLT_pred$Prediction, Lower = UK76drcLT_pred$Lower, Upper = UK76drcLT_pred$Upper)

UK80drcLT_preddata = data.frame(temp = seq(26,38, length.out = 100))
UK80drcLT_pred = as.data.frame(predict(UK80drcLT, newdata = UK80drcLT_preddata, interval = 'confidence'))
UK80drcLT_preddata = data.frame(UK80drcLT_preddata, fvfm = UK80drcLT_pred$Prediction, Lower = UK80drcLT_pred$Lower, Upper = UK80drcLT_pred$Upper)

####LT-plot ED25 DRC for genos####
hex<-hue_pal()(10)
hex #"#F8766D" "#D89000" "#A3A500" "#39B600" "#00BF7D" "#00BFC4" "#00B0F6" "#9590FF" "#E76BF3" "#FF62BC"
DRC_plot<- ggplot() +
  geom_jitter(data = pamOW, aes(x = RealTemp, y = FvFm, color=Genotype), size = 1, width = 0.25) +
  scale_x_continuous(limits=c(26,33), breaks=c(27,29,31,33)) +
  scale_y_continuous(limits=c(0.0, 0.7), breaks=c(0, 0.2, 0.4, 0.6)) +
  
  geom_line(data = LK31drcLT_preddata, aes(x = temp, y = fvfm), color = '#F8766D', show.legend = FALSE) +
  #geom_ribbon(data = LK31drcLT_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = '#F8766D', linetype=2, alpha = 0.2) +
  geom_vline(data = geno_drcLT_coeff_long, aes(xintercept = geno_coeff[1]), color = '#F8766D', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drcLT_coeff_long, aes(label = paste("LK31:",geno_coeff[1])), x = 32, y = 0.55, show.legend = FALSE, color = '#F8766D', hjust='left') +
  geom_line(data = LK41drcLT_preddata, aes(x = temp, y = fvfm), color = '#D89000', show.legend = FALSE) +
  geom_vline(data = geno_drcLT_coeff_long, aes(xintercept = geno_coeff[2]), color = '#D89000', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drcLT_coeff_long, aes(label = paste("LK41:",geno_coeff[2])), x = 32, y = 0.61, show.legend = FALSE, color = '#D89000', hjust='left') +
  geom_line(data = LK50drcLT_preddata, aes(x = temp, y = fvfm), color = '#A3A500', show.legend = FALSE) +
  geom_vline(data = geno_drcLT_coeff_long, aes(xintercept = geno_coeff[3]), color = '#A3A500', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drcLT_coeff_long, aes(label = paste("LK50:",geno_coeff[3])), x = 32, y = 0.58, show.legend = FALSE, color = '#A3A500', hjust='left') +
  geom_line(data = LK62drcLT_preddata, aes(x = temp, y = fvfm), color = '#39B600', show.legend = FALSE) +
  geom_vline(data = geno_drcLT_coeff_long, aes(xintercept = geno_coeff[4]), color = '#39B600', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drcLT_coeff_long, aes(label = paste("LK62:",geno_coeff[4])), x = 32, y = 0.67, show.legend = FALSE, color = '#39B600', hjust='left') +
  geom_line(data = LK7drcLT_preddata, aes(x = temp, y = fvfm), color = '#00BF7D', show.legend = FALSE) +
  geom_vline(data = geno_drcLT_coeff_long, aes(xintercept = geno_coeff[5]), color = '#00BF7D', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drcLT_coeff_long, aes(label = paste("LK7:",geno_coeff[5])), x = 32, y = 0.46, show.legend = FALSE, color = '#00BF7D', hjust='left') +
  geom_line(data = UK12drcLT_preddata, aes(x = temp, y = fvfm), color = '#00BFC4', show.legend = FALSE) +
  geom_vline(data = geno_drcLT_coeff_long, aes(xintercept = geno_coeff[6]), color = '#00BFC4', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drcLT_coeff_long, aes(label = paste("UK12:",geno_coeff[6])), x = 32, y = 0.7, show.legend = FALSE, color = '#00BFC4', hjust='left') +
  geom_line(data = UK19drcLT_preddata, aes(x = temp, y = fvfm), color = '#00B0F6', show.legend = FALSE) +
  geom_vline(data = geno_drcLT_coeff_long, aes(xintercept = geno_coeff[7]), color = '#00B0F6', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drcLT_coeff_long, aes(label = paste("UK19:",geno_coeff[7])), x = 32, y = 0.49, show.legend = FALSE, color = '#00B0F6', hjust='left') +
  geom_line(data = UK70drcLT_preddata, aes(x = temp, y = fvfm), color = '#9590FF', show.legend = FALSE) +
  geom_vline(data = geno_drcLT_coeff_long, aes(xintercept = geno_coeff[8]), color = '#9590FF', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drcLT_coeff_long, aes(label = paste("UK70:",geno_coeff[8])), x = 32, y = 0.52, show.legend = FALSE, color = '#9590FF', hjust='left') +
  geom_line(data = UK76drcLT_preddata, aes(x = temp, y = fvfm), color = '#E76BF3', show.legend = FALSE) +
  geom_vline(data = geno_drcLT_coeff_long, aes(xintercept = geno_coeff[9]), color = '#E76BF3', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drcLT_coeff_long, aes(label = paste("UK76:",geno_coeff[9])), x = 32, y = 0.64, show.legend = FALSE, color = '#E76BF3', hjust='left') +
  geom_line(data = UK80drcLT_preddata, aes(x = temp, y = fvfm), color = '#FF62BC', show.legend = FALSE) +
  geom_vline(data = geno_drcLT_coeff_long, aes(xintercept = geno_coeff[10]), color = '#FF62BC', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drcLT_coeff_long, aes(label = paste("UK80:",geno_coeff[10])), x = 32, y = 0.43, show.legend = FALSE, color = '#FF62BC', hjust='left') +
  
  ylab(label = expression(paste(italic(F[V]/F[M])))) +
  xlab("Temperature (°C)") +
  theme_classic()+ theme(legend.position = "none")

DRC_plot

ggsave(DRC_plot, height = 6 , width = 10, filename = "20230321_LT_ED25genoDRCs.pdf", useDingbats=FALSE)
####CBASS-ED25 DRC-Temp####
#### LK31 CBASS ED25 
#model
LK31drc  <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="LK31",], fct = LL.3())
#extract coeffs 33.29234 
LK31drc_coeff<-data.frame(ED(LK31drc, c(25)))
LK31drc_coeff_mean<-LK31drc_coeff[,1]
LK31drc_coeff_lower<-LK31drc_coeff[,1] - 1.96*LK31drc_coeff[,2]
LK31drc_coeff_upper<-LK31drc_coeff[,1] + 1.96*LK31drc_coeff[,2]

#### LK41 CBASS ED25
#model
LK41drc  <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="LK41",], fct = LL.3())
#extract coeffs 33.74078
LK41drc_coeff<-data.frame(ED(LK41drc, c(25)))
LK41drc_coeff_mean<-LK41drc_coeff[,1]
LK41drc_coeff_lower<-LK41drc_coeff[,1] - 1.96*LK41drc_coeff[,2]
LK41drc_coeff_upper<-LK41drc_coeff[,1] + 1.96*LK41drc_coeff[,2]

#### LK50 CBASS ED25
#model
LK50drc  <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="LK50",], fct = LL.3())
#extract coeffs 33.28866
LK50drc_coeff<-data.frame(ED(LK50drc, c(25)))
LK50drc_coeff_mean<-LK50drc_coeff[,1]
LK50drc_coeff_lower<-LK50drc_coeff[,1] - 1.96*LK50drc_coeff[,2]
LK50drc_coeff_upper<-LK50drc_coeff[,1] + 1.96*LK50drc_coeff[,2]

#### LK62 CBASS ED25
#model
LK62drc  <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="LK62",], fct = LL.3())
#extract coeffs 33.381987
LK62drc_coeff<-data.frame(ED(LK62drc, c(25)))
LK62drc_coeff_mean<-LK62drc_coeff[,1]
LK62drc_coeff_lower<-LK62drc_coeff[,1] - 1.96*LK62drc_coeff[,2]
LK62drc_coeff_upper<-LK62drc_coeff[,1] + 1.96*LK62drc_coeff[,2]

#### LK7 CBASS ED25 
#model
LK7drc  <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="LK7",], fct = LL.3())
#extract coeffs 33.80948
LK7drc_coeff<-data.frame(ED(LK7drc, c(25)))
LK7drc_coeff_mean<-LK7drc_coeff[,1]
LK7drc_coeff_lower<-LK7drc_coeff[,1] - 1.96*LK7drc_coeff[,2]
LK7drc_coeff_upper<-LK7drc_coeff[,1] + 1.96*LK7drc_coeff[,2]

#### UK12 CBASS ED25 
#model
UK12drc  <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="UK12",], fct = LL.3())
#extract coeffs 33.96815
UK12drc_coeff<-data.frame(ED(UK12drc, c(25)))
UK12drc_coeff_mean<-UK12drc_coeff[,1]
UK12drc_coeff_lower<-UK12drc_coeff[,1] - 1.96*UK12drc_coeff[,2]
UK12drc_coeff_upper<-UK12drc_coeff[,1] + 1.96*UK12drc_coeff[,2]

#### UK19 CBASS ED25 
#model
UK19drc  <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="UK19",], fct = LL.3())
#extract coeffs 33.89115 
UK19drc_coeff<-data.frame(ED(UK19drc, c(25)))
UK19drc_coeff_mean<-UK19drc_coeff[,1]
UK19drc_coeff_lower<-UK19drc_coeff[,1] - 1.96*UK19drc_coeff[,2]
UK19drc_coeff_upper<-UK19drc_coeff[,1] + 1.96*UK19drc_coeff[,2]

#### UK70 CBASS ED25
#model
UK70drc  <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="UK70",], fct = LL.3())
#extract coeffs 33.68033
UK70drc_coeff<-data.frame(ED(UK70drc, c(25)))
UK70drc_coeff_mean<-UK70drc_coeff[,1]
UK70drc_coeff_lower<-UK70drc_coeff[,1] - 1.96*UK70drc_coeff[,2]
UK70drc_coeff_upper<-UK70drc_coeff[,1] + 1.96*UK70drc_coeff[,2]

#### UK76 CBASS ED25 
#model
UK76drc  <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="UK76",], fct = LL.3())
#extract coeffs 33.7162 
UK76drc_coeff<-data.frame(ED(UK76drc, c(25)))
UK76drc_coeff_mean<-UK76drc_coeff[,1]
UK76drc_coeff_lower<-UK76drc_coeff[,1] - 1.96*UK76drc_coeff[,2]
UK76drc_coeff_upper<-UK76drc_coeff[,1] + 1.96*UK76drc_coeff[,2]

#### UK80 CBASS ED25 
#model
UK80drc  <- drm(FvFm ~ Treatment, data = pam[pam$Genotype=="UK80",], fct = LL.3())
#extract coeffs 33.34052
UK80drc_coeff<-data.frame(ED(UK80drc, c(25)))
UK80drc_coeff_mean<-UK80drc_coeff[,1]
UK80drc_coeff_lower<-UK80drc_coeff[,1] - 1.96*UK80drc_coeff[,2]
UK80drc_coeff_upper<-UK80drc_coeff[,1] + 1.96*UK80drc_coeff[,2]

####pulling geno ED25 predictors and confidence intervals####
geno_drc_coeff_means<- data.frame(LK31drc_coeff_mean, LK41drc_coeff_mean, LK50drc_coeff_mean, LK62drc_coeff_mean, LK7drc_coeff_mean, UK12drc_coeff_mean, UK19drc_coeff_mean, UK70drc_coeff_mean, UK76drc_coeff_mean, UK80drc_coeff_mean)
geno_drc_coeff_long<-pivot_longer(cols=everything(),geno_drc_coeff_means, names_to="name",values_to="geno_coeff")
geno_drc_coeff_long$name<-as.factor(geno_drc_coeff_long$name) 
geno_drc_coeff_long<- subset(geno_drc_coeff_long, select = -c(name))
geno_drc_coeff_long[,'geno_coeff']=round(geno_drc_coeff_long[,'geno_coeff'],2) #round to two decimals
geno_drc_coeff_lowers<-data.frame(LK31drc_coeff_lower,LK41drc_coeff_lower,LK50drc_coeff_lower,LK62drc_coeff_lower,LK7drc_coeff_lower,UK12drc_coeff_lower,UK19drc_coeff_lower,UK70drc_coeff_lower,UK76drc_coeff_lower,UK80drc_coeff_lower )
geno_drc_coeff_uppers<-data.frame(LK31drc_coeff_upper,LK41drc_coeff_upper,LK50drc_coeff_upper,LK62drc_coeff_upper,LK7drc_coeff_upper,UK12drc_coeff_upper,UK19drc_coeff_upper,UK70drc_coeff_upper,UK76drc_coeff_upper,UK80drc_coeff_upper )

LK31drc_preddata = data.frame(temp = seq(26,38, length.out = 100))
LK31drc_pred = as.data.frame(predict(LK31drc, newdata = LK31drc_preddata, interval = 'confidence'))
LK31drc_preddata = data.frame(LK31drc_preddata, fvfm = LK31drc_pred$Prediction, Lower = LK31drc_pred$Lower, Upper = LK31drc_pred$Upper)

LK41drc_preddata = data.frame(temp = seq(26,38, length.out = 100))
LK41drc_pred = as.data.frame(predict(LK41drc, newdata = LK41drc_preddata, interval = 'confidence'))
LK41drc_preddata = data.frame(LK41drc_preddata, fvfm = LK41drc_pred$Prediction, Lower = LK41drc_pred$Lower, Upper = LK41drc_pred$Upper)

LK50drc_preddata = data.frame(temp = seq(26,38, length.out = 100))
LK50drc_pred = as.data.frame(predict(LK50drc, newdata = LK50drc_preddata, interval = 'confidence'))
LK50drc_preddata = data.frame(LK50drc_preddata, fvfm = LK50drc_pred$Prediction, Lower = LK50drc_pred$Lower, Upper = LK50drc_pred$Upper)

LK62drc_preddata = data.frame(temp = seq(26,38, length.out = 100))
LK62drc_pred = as.data.frame(predict(LK62drc, newdata = LK62drc_preddata, interval = 'confidence'))
LK62drc_preddata = data.frame(LK62drc_preddata, fvfm = LK62drc_pred$Prediction, Lower = LK62drc_pred$Lower, Upper = LK62drc_pred$Upper)

LK7drc_preddata = data.frame(temp = seq(26,38, length.out = 100))
LK7drc_pred = as.data.frame(predict(LK7drc, newdata = LK7drc_preddata, interval = 'confidence'))
LK7drc_preddata = data.frame(LK7drc_preddata, fvfm = LK7drc_pred$Prediction, Lower = LK7drc_pred$Lower, Upper = LK7drc_pred$Upper)

UK12drc_preddata = data.frame(temp = seq(26,38, length.out = 100))
UK12drc_pred = as.data.frame(predict(UK12drc, newdata = UK12drc_preddata, interval = 'confidence'))
UK12drc_preddata = data.frame(UK12drc_preddata, fvfm = UK12drc_pred$Prediction, Lower = UK12drc_pred$Lower, Upper = UK12drc_pred$Upper)

UK19drc_preddata = data.frame(temp = seq(26,38, length.out = 100))
UK19drc_pred = as.data.frame(predict(UK19drc, newdata = UK19drc_preddata, interval = 'confidence'))
UK19drc_preddata = data.frame(UK19drc_preddata, fvfm = UK19drc_pred$Prediction, Lower = UK19drc_pred$Lower, Upper = UK19drc_pred$Upper)

UK70drc_preddata = data.frame(temp = seq(26,38, length.out = 100))
UK70drc_pred = as.data.frame(predict(UK70drc, newdata = UK70drc_preddata, interval = 'confidence'))
UK70drc_preddata = data.frame(UK70drc_preddata, fvfm = UK70drc_pred$Prediction, Lower = UK70drc_pred$Lower, Upper = UK70drc_pred$Upper)

UK76drc_preddata = data.frame(temp = seq(26,38, length.out = 100))
UK76drc_pred = as.data.frame(predict(UK76drc, newdata = UK76drc_preddata, interval = 'confidence'))
UK76drc_preddata = data.frame(UK76drc_preddata, fvfm = UK76drc_pred$Prediction, Lower = UK76drc_pred$Lower, Upper = UK76drc_pred$Upper)

UK80drc_preddata = data.frame(temp = seq(26,38, length.out = 100))
UK80drc_pred = as.data.frame(predict(UK80drc, newdata = UK80drc_preddata, interval = 'confidence'))
UK80drc_preddata = data.frame(UK80drc_preddata, fvfm = UK80drc_pred$Prediction, Lower = UK80drc_pred$Lower, Upper = UK80drc_pred$Upper)

####plot ED25 DRC for CBASS genos####
hex #"#F8766D" "#D89000" "#A3A500" "#39B600" "#00BF7D" "#00BFC4" "#00B0F6" "#9590FF" "#E76BF3" "#FF62BC"
DRC_plot<- ggplot() +
  geom_jitter(data = pam, aes(x = Treatment, y = FvFm, color=Genotype), size = 1, width = 0.25) +
  scale_x_continuous(limits=c(26,37), breaks=c(27,29,31,33,35,37)) +
  scale_y_continuous(limits=c(0.0, 0.7), breaks=c(0, 0.2, 0.4, 0.6)) +
  
  geom_line(data = LK31drc_preddata, aes(x = temp, y = fvfm), color = '#F8766D', show.legend = FALSE) +
  #geom_ribbon(data = LK31drc_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = '#F8766D', linetype=2, alpha = 0.2) +
  geom_vline(data = geno_drc_coeff_long, aes(xintercept = geno_coeff[1]), color = '#F8766D', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drc_coeff_long, aes(label = paste("LK31:",geno_coeff[1])), x = 36, y = 0.43, show.legend = FALSE, color = '#F8766D', hjust='left') +
  geom_line(data = LK41drc_preddata, aes(x = temp, y = fvfm), color = '#D89000', show.legend = FALSE) +
  geom_vline(data = geno_drc_coeff_long, aes(xintercept = geno_coeff[2]), color = '#D89000', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drc_coeff_long, aes(label = paste("LK41:",geno_coeff[2])), x = 36, y = 0.61, show.legend = FALSE, color = '#D89000', hjust='left') +
  geom_line(data = LK50drc_preddata, aes(x = temp, y = fvfm), color = '#A3A500', show.legend = FALSE) +
  geom_vline(data = geno_drc_coeff_long, aes(xintercept = geno_coeff[3]), color = '#A3A500', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drc_coeff_long, aes(label = paste("LK50:",geno_coeff[3])), x = 36, y = 0.46, show.legend = FALSE, color = '#A3A500', hjust='left') +
  geom_line(data = LK62drc_preddata, aes(x = temp, y = fvfm), color = '#39B600', show.legend = FALSE) +
  geom_vline(data = geno_drc_coeff_long, aes(xintercept = geno_coeff[4]), color = '#39B600', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drc_coeff_long, aes(label = paste("LK62:",geno_coeff[4])), x = 36, y = 0.52, show.legend = FALSE, color = '#39B600', hjust='left') +
  geom_line(data = LK7drc_preddata, aes(x = temp, y = fvfm), color = '#00BF7D', show.legend = FALSE) +
  geom_vline(data = geno_drc_coeff_long, aes(xintercept = geno_coeff[5]), color = '#00BF7D', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drc_coeff_long, aes(label = paste("LK7:",geno_coeff[5])), x = 36, y = 0.64, show.legend = FALSE, color = '#00BF7D', hjust='left') +
  geom_line(data = UK12drc_preddata, aes(x = temp, y = fvfm), color = '#00BFC4', show.legend = FALSE) +
  geom_vline(data = geno_drc_coeff_long, aes(xintercept = geno_coeff[6]), color = '#00BFC4', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drc_coeff_long, aes(label = paste("UK12:",geno_coeff[6])), x = 36, y = 0.7, show.legend = FALSE, color = '#00BFC4', hjust='left') +
  geom_line(data = UK19drc_preddata, aes(x = temp, y = fvfm), color = '#00B0F6', show.legend = FALSE) +
  geom_vline(data = geno_drc_coeff_long, aes(xintercept = geno_coeff[7]), color = '#00B0F6', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drc_coeff_long, aes(label = paste("UK19:",geno_coeff[7])), x = 36, y = 0.67, show.legend = FALSE, color = '#00B0F6', hjust='left') +
  geom_line(data = UK70drc_preddata, aes(x = temp, y = fvfm), color = '#9590FF', show.legend = FALSE) +
  geom_vline(data = geno_drc_coeff_long, aes(xintercept = geno_coeff[8]), color = '#9590FF', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drc_coeff_long, aes(label = paste("UK70:",geno_coeff[8])), x = 36, y = 0.55, show.legend = FALSE, color = '#9590FF', hjust='left') +
  geom_line(data = UK76drc_preddata, aes(x = temp, y = fvfm), color = '#E76BF3', show.legend = FALSE) +
  geom_vline(data = geno_drc_coeff_long, aes(xintercept = geno_coeff[9]), color = '#E76BF3', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drc_coeff_long, aes(label = paste("UK76:",geno_coeff[9])), x = 36, y = 0.58, show.legend = FALSE, color = '#E76BF3', hjust='left') +
  geom_line(data = UK80drc_preddata, aes(x = temp, y = fvfm), color = '#FF62BC', show.legend = FALSE) +
  geom_vline(data = geno_drc_coeff_long, aes(xintercept = geno_coeff[10]), color = '#FF62BC', show.legend = FALSE, lwd=1) +
  geom_text(data = geno_drc_coeff_long, aes(label = paste("UK80:",geno_coeff[10])), x = 36, y = 0.49, show.legend = FALSE, color = '#FF62BC', hjust='left') +
  
  ylab(label = expression(paste(italic(F[V]/F[M])))) +
  xlab("Temperature (°C)") +
  theme_classic()+ theme(legend.position = "none")

DRC_plot

ggsave(DRC_plot, height = 6 , width = 10, filename = "20230321_CBASS_ED25genoDRCs.pdf", useDingbats=FALSE)

####ED25 DRC-HL####
pam_clean<-subset(pam, select=c(Genotype, FvFm, Heatload30.5))
pamOW_clean<-subset(pamOW, select=c(Genotype, FvFm, Heatload30.5))
#Full_data<-rbind(pam_clean,pamOW_clean)

#model for CBASS HL
CBASS_HL <- drm(FvFm ~ Heatload30.5, curveid = Genotype, data = pam_clean, fct = LL.3())
summary(CBASS_HL)
plot(CBASS_HL)

#extract coeffs
CBASS_HL_coeff<-data.frame(ED(CBASS_HL, c(25)))

CBASS_HL_coeff_mean<-CBASS_HL_coeff[,1]
CBASS_HL_coeff_lower<-CBASS_HL_coeff[,1] - 1.96*CBASS_HL_coeff[,2]
CBASS_HL_coeff_upper<-CBASS_HL_coeff[,1] + 1.96*CBASS_HL_coeff[,2]

#model for LT HL
LT_HL <- drm(FvFm ~ Heatload30.5, curveid = Genotype, data = pamOW_clean, fct = LL.3())
mselect(LT_HL, list(LL.2(), LL.4(), LL.5(), LL2.2(), LL2.3(), LL2.4(), LL2.5(), AR.2(), AR.3(), EXD.2(), EXD.3()), icfct = AIC)
modelFit(LT_HL)
summary(LT_HL)
plot(LT_HL)

#extract coeffs
LT_HL_coeff<-data.frame(ED(LT_HL, c(25)))

LT_HL_coeff_mean<-LT_HL_coeff[,1]
LT_HL_coeff_lower<-LT_HL_coeff[,1] - 1.96*LT_HL_coeff[,2]
LT_HL_coeff_upper<-LT_HL_coeff[,1] + 1.96*LT_HL_coeff[,2]

#combine ED25 data plus predict curves from models
Acr_HL_coeff_means<-data.frame(CBASS_HL_coeff_mean, LT_HL_coeff_mean)
Acr_HL_coeff_lowers<-data.frame(CBASS_HL_coeff_lower,LT_HL_coeff_lower)
Acr_HL_coeff_uppers<-data.frame(CBASS_HL_coeff_upper,LT_HL_coeff_upper)

Acr_exp_HL_LT_preddata = data.frame(Treatment = seq(0,365, length.out = 100))
Acr_exp_HL_LT_pred = as.data.frame(predict(LT_HL, newdata = Acr_exp_HL_LT_preddata, interval = 'confidence'))
Acr_exp_HL_LT_preddata = data.frame(Acr_exp_HL_LT_preddata, fvfm = Acr_exp_HL_LT_pred$Prediction, Lower = Acr_exp_HL_LT_pred$Lower, Upper = Acr_exp_HL_LT_pred$Upper)

Acr_exp_HL_CBASS_preddata = data.frame(Treatment = seq(0,10, length.out = 100))
Acr_exp_HL_CBASS_pred = as.data.frame(predict(CBASS_HL, newdata = Acr_exp_HL_CBASS_preddata, interval = 'confidence'))
Acr_exp_HL_CBASS_preddata = data.frame(Acr_exp_HL_CBASS_preddata, fvfm = Acr_exp_HL_CBASS_pred$Prediction, Lower = Acr_exp_HL_CBASS_pred$Lower, Upper = Acr_exp_HL_CBASS_pred$Upper)

##CBASS-plot ED25 DRC-HL####
HL_Acr_CBASS<- ggplot() +
  geom_jitter(data = pam_clean, aes(x = Heatload30.5, y = FvFm, color = Genotype), size = 2) + geom_jitter() +
  scale_x_continuous(limits=c(0, 6), breaks=c(0,1,2,3, 4, 5, 6)) +
  scale_y_continuous(limits=c(-0.02, 0.75), breaks=c(0, 0.2, 0.4, 0.6)) +
  geom_line(data = Acr_exp_HL_CBASS_preddata, aes(x = Treatment, y = fvfm), color = 'tomato', show.legend = FALSE) +
  geom_ribbon(data = Acr_exp_HL_CBASS_preddata, aes(x = Treatment, ymin=Lower, ymax=Upper), color = 'tomato', linetype=2, alpha = 0.2) +
  geom_vline(data = Acr_HL_coeff_means, aes(xintercept = CBASS_HL_coeff_mean), color = 'tomato', show.legend = FALSE) +
  annotate("rect", xmin=Acr_HL_coeff_lowers$CBASS_HL_coeff_lower, xmax=Acr_HL_coeff_uppers$CBASS_HL_coeff_upper, ymin=-0.02, ymax=0.75, fill= 'tomato',  alpha = 0.2) +
  geom_text(data = Acr_HL_coeff_means, aes(label=round(CBASS_HL_coeff_mean, digits = 2)), x = 5.65, y = 0.6, show.legend = FALSE, color = 'tomato',size=6) +
  scale_color_brewer(palette = "Paired") +
  ggtitle("CBASS") +
  ylab(label = expression(paste(italic(F[V]/F[M])))) +
  xlab("Heatload (hrs) above 30.5°C") +
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
        legend.position = "none",
        legend.text = element_text(size=14),
        legend.title = element_text(size=16)
  )
ggsave(HL_Acr_CBASS, height = 6 , width = 6, filename = "./202206_CBASS_ED25-HL.pdf", useDingbats=FALSE)

##LT-plot ED25 DRC-HL####
HL_Acr_LT<- ggplot() +
  geom_jitter(data = pamOW_clean, aes(x = Heatload30.5, y = FvFm, color = Genotype), size = 2) + geom_jitter() +
  scale_x_continuous(limits=c(0, 365), breaks=c(0,100,200,300,400)) +
  scale_y_continuous(limits=c(-0.02, 0.75), breaks=c(0, 0.2, 0.4, 0.6)) +
  geom_line(data = Acr_exp_HL_LT_preddata, aes(x = Treatment, y = fvfm), color = 'blue', show.legend = FALSE) +
  geom_ribbon(data = Acr_exp_HL_LT_preddata, aes(x = Treatment, ymin=Lower, ymax=Upper), color = 'blue', linetype=2, alpha = 0.2) +
  geom_vline(data = Acr_HL_coeff_means, aes(xintercept = LT_HL_coeff_mean), color = 'blue', show.legend = FALSE) +
  annotate("rect", xmin=Acr_HL_coeff_lowers$LT_HL_coeff_lower, xmax=Acr_HL_coeff_uppers$LT_HL_coeff_upper, ymin=-Inf, ymax=Inf, fill= 'blue',  alpha = 0.2) +
  geom_text(data = Acr_HL_coeff_means, aes(label=round(LT_HL_coeff_mean, digits = 2)), x = 250, y = 0.6, show.legend = FALSE, color = 'blue',size=5) +
  scale_color_brewer(palette = "Paired") +
  ggtitle("LT") +
  ylab(label = expression(paste(italic(F[V]/F[M])))) +
  xlab("Heatload (hrs) above 30.5°C") +
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
        legend.position = "none",
        legend.text = element_text(size=14),
        legend.title = element_text(size=16)
  )
ggsave(HL_Acr_LT, height = 6 , width = 6, filename = "./202206_LT_ED25-HL.pdf", useDingbats=FALSE)
