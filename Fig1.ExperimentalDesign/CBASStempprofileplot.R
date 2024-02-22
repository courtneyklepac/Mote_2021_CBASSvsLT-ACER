#### Inkbird plotting ####
library(ggplot2)
library(reshape2)
library(dplyr)
library(openair)
library(FSA)
library(plyr)

#SD card files
#merge and clean temp files, temp1 one system/SD card, temp2 second system/SD
#merge and clean temp files, temp1 one system, temp2 second
temp1<-read.csv("20210621_CBASS2_RW5N_1-4.txt")
temp1$time<-strptime(paste(temp1$Date,paste(temp1$Th,temp1$Tm,temp1$Ts,sep=":"),sep=" "),format="%Y_%B_%d %H:%M:%S")
temp1<-temp1[order(temp1$time),]
names(temp1)
temp1<-subset(temp1, select=c(24,8,12,16,20))

temp2<-read.csv("20210621_CBASS2_RW6N_5-8.txt")
temp2$time<-strptime(paste(temp2$Date,paste(temp2$Th,temp2$Tm,temp2$Ts,sep=":"),sep=" "),format="%Y_%B_%d %H:%M:%S")
temp2<-temp2[order(temp2$time),]
temp2<-subset(temp2, select=c(24,8,12,16,20))

library(plyr)
#6/19/21 35.6, 32.9, 27.5, 30.2, 34.3, 28.9, 37, 31.6
temp1<-rename(temp1, c("T1inT"="35.6", "T2inT"="32.9", "T3inT"="27.5", "T4inT"="30.2"))
temp2<-rename(temp2, c("T1inT"="34.3", "T2inT"="28.9", "T3inT"="37", "T4inT"="31.6"))

#6/20/21 34.3, 30.2, 35.6, 28.9, 31.6, 27.5, 32.9, 37
temp1<-rename(temp1, c("T1inT"="34.3", "T2inT"="30.2", "T3inT"="35.6", "T4inT"="28.9"))
temp2<-rename(temp2, c("T1inT"="31.6", "T2inT"="27.5", "T3inT"="32.9", "T4inT"="37"))

#6/21/21 37, 28.9, 34.3, 30.2, 35.6, 31.6, 27.5, 32.9
temp1<-rename(temp1, c("T1inT"="37", "T2inT"="28.9", "T3inT"="34.3", "T4inT"="30.2"))
temp2<-rename(temp2, c("T1inT"="35.6", "T2inT"="31.6", "T3inT"="27.5", "T4inT"="32.9"))

library(reshape2)
temp1$time<-as.POSIXct(temp1$time)
temp1.1<-melt(temp1, id.vars=c("time"))

temp2$time<-as.POSIXct(temp2$time)
temp2.1<-melt(temp2, id.vars=c("time"))

spoing<-merge(temp1, temp2, all=T)
CBASS_full <- melt(spoing, id="time")

names(CBASS_full)[names(CBASS_full) == "time"] <- "date"
names(CBASS_full)[names(CBASS_full) == "variable"] <- "Tank"
names(CBASS_full)[names(CBASS_full) == "value"] <- "Temp"
#6/19
CBASS_full$Tank=factor(CBASS_full$Tank,levels(CBASS_full$Tank)[c(3,6,4,8,2,5,1,7)])
#6/20
CBASS_full$Tank=factor(CBASS_full$Tank,levels(CBASS_full$Tank)[c(6,4,2,5,7,1,3,8)])
#6/21
CBASS_full$Tank=factor(CBASS_full$Tank,levels(CBASS_full$Tank)[c(7,2,4,6,8,3,5,1)])

library(openair)
CBASS_full$date<-as.POSIXct(CBASS_full$date)
CBASS_full <- na.omit(CBASS_full)
CBASS_full$Temp<-as.numeric(CBASS_full$Temp)

#Select for experimental period within larger file, average data to 15-min intervals for plotting
CBASS_av<-timeAverage(CBASS_full, avg.time = "15 min", type = c("Tank","Temp"))

#calculate mean temp during the hold
CBASS_hold<-selectByDate(CBASS_full, start = "2021-06-02", end = "2021-06-02",hour = 16:19)
CBASS_hold$Temp<-as.numeric(CBASS_hold$Temp)
CBASS_summary_data <- ddply(CBASS_hold, c("Tank"), summarise, N= length(Temp), mean = mean(Temp),
                            sd   = sd(Temp),
                            se   = sd / sqrt(N))
#plotting Fig. 1 8temps
cols <- c("27.5" = "navyblue", "28.9" = "blue", "30.2" = "yellow2", "31.6" = "goldenrod", 
          "32.9" = "darkorange2", "34.3" = "tomato", "35.6" = "red3", "37" = "firebrick")
CBASS_plot<-ggplot(CBASS_full, aes(x=date, y=Temp, col=Tank)) + scale_colour_manual(values=cols) +
  geom_rect(aes(xmin = as.POSIXct("2021-06-21 19:00"), xmax = as.POSIXct("2021-06-22 07:00"), 
                ymin = -Inf, ymax = Inf), color=NA, fill = "lightgray") +
  geom_line(aes(linetype=Tank),lwd=1.5) +
  scale_linetype_manual(values=c("solid","solid","solid","solid","solid","solid","solid","solid")) +
  scale_y_continuous(limits=c(27.5,37), breaks = c(27,29,31,33,35,37,38)) + 
  scale_x_datetime(name = "Hours", date_labels = "%Hh", date_breaks = "2 hour", limits=c(as.POSIXct("2021-06-21 12:00:01"), as.POSIXct("2021-06-22 08:00:00"))) +
  geom_vline(xintercept=as.POSIXct("2021-06-21 20:00"),linetype="solid",lwd=0.7,color="gray50") +
  geom_vline(xintercept=as.POSIXct("2021-06-21 22:00"),linetype="solid",lwd=0.7,color="gray50") +
  annotate("text", x = as.POSIXct("2021-06-21 21:00"), y = 27.7, angle = 90, adj = 0, size = 5,
           label = expression(paste(italic(F[V]/F[M]))))+
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none') + 
  theme(text = element_text(size=14, face="bold")) + 
  labs(x = "Time (hours)", y = "Temperature (°C)") +
  theme(axis.ticks.length=unit(.3,"cm")) 
ggsave(CBASS_plot, height = 6 , width = 10, filename = "20221111_CBASS2_tempprofileplot.pdf", useDingbats=FALSE)

#########Calculate hours over threshold########
#average by min
spoimg_av<-timeAverage(CBASS_full, avg.time = "1 min",type='Tank')
#subset into Tanks
spoing37<-spoimg_av[spoimg_av$Tank=='37',]
spoing35.6<-spoimg_av[spoimg_av$Tank=='35.6',]
spoing34.3<-spoimg_av[spoimg_av$Tank=='34.3',]
spoing32.9<-spoimg_av[spoimg_av$Tank=='32.9',]
spoing31.6<-spoimg_av[spoimg_av$Tank=='31.6',]
spoing30.2<-spoimg_av[spoimg_av$Tank=='30.2',]

#total number of hours logged over threshold
spoing37_30.5<- (sum(spoing37$Temp > 30.5,na.rm=T)/60) #6.08
spoing35.6_30.5<- (sum(spoing35.6$Temp > 30.5,na.rm=T)/60) #5.73
spoing34.3_30.5<- (sum(spoing34.3$Temp > 30.5,na.rm=T)/60) #5.4
spoing32.9_30.5<- (sum(spoing32.9$Temp > 30.5,na.rm=T)/60) #4.92
spoing31.6_30.5<- (sum(spoing31.6$Temp > 30.5,na.rm=T)/60) #4.28
spoing30.2_30.5<- (sum(spoing30.2$Temp > 30.5,na.rm=T)/60) #0

#### LT plotting ####
tempLT<-read.csv("APAL_LTtemps2.csv")
tempLT$Date<-as.POSIXct(paste(tempLT$Date), format = "%m/%d/%y") 
tempLT$Temp<-as.numeric(tempLT$Temp)
tempLT$date<-as.Date(tempLT$Date) 
tempLT$date<-as.factor(tempLT$date) 
tempLT$Tank<-as.factor(tempLT$Tank)
tempLT$Treatment<-as.factor(tempLT$Treatment)
tempLT$Days<-as.factor(tempLT$Days)

AP2<- tempLT[tempLT$Treatment=="Control" | tempLT$Treatment== "OW",]

tempLT_summary<- Summarize(Temp~Treatment+date, data=AP2, digits=3)
tempLT_summary$date<-as.Date(tempLT_summary$date)
#plotting Treatments
cols <- c("Control" = "navyblue", "OW" = "firebrick")
tempLT_plot<-ggplot(tempLT, aes(x=date, y=Temp, col=Treatment)) + scale_colour_manual(values=cols) +
  geom_line(aes(linetype=Treatment),size=1.5) +
  scale_linetype_manual(values=c("solid","solid")) +
  scale_y_continuous(limits=c(27.5,31.5), breaks = c(27,28,29,30,31)) + 
  scale_x_date(date_breaks = "1 week",date_labels='%d')+
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none', panel.background = element_rect(fill='transparent'), #transparent panel bg
                          plot.background = element_rect(fill='transparent', color=NA)) + 
  theme(text = element_text(size=14, face="bold")) + 
  theme(axis.ticks.length=unit(.3,"cm")) +
  ylab("Temperature (°C)") + xlab("Days") +
  geom_vline(xintercept=as.numeric(tempLT$date[1]),linetype="solid",lwd=0.7,color="gray50") +
  geom_vline(xintercept=as.numeric(tempLT$date[8]),linetype="solid",lwd=0.7,color="gray50") +
  geom_vline(xintercept=as.numeric(tempLT$date[59]),linetype="solid",lwd=0.7,color="gray50") +
  geom_vline(xintercept=as.numeric(tempLT$date[65]),linetype="solid",lwd=0.7,color="gray50") +
  annotate("text", x = as.Date("2021-03-18"), y = 27.7, angle = 90, adj = 0, size = 5,
           label = expression(paste(italic(F[V]/F[M])))) +
  annotate("text", x = as.Date("2021-05-15"), y = 27.7, angle = 90, adj = 0, size = 5,
           label = expression(paste(italic(F[V]/F[M])))) +
  geom_hline(yintercept=c(27.5,30.5,31.5), linetype="dashed", color = c("navyblue","red","firebrick"), size=0.5) 

ggsave(tempLT_plot, bg='transparent', height = 6 , width = 10, filename = "20221111_LT-tempprofileplot.pdf",useDingbats=FALSE)

#plotting each tank
tempLT_plot<-ggplot(tempLT, aes(x=Date, y=Temp, col=Tank)) +
  geom_line(aes(linetype=Tank)) +
  scale_linetype_manual(values=c("dashed","dashed","dashed","dashed","dashed","solid","solid","solid","solid","solid")) +
  scale_y_continuous(limits=c(26.5,32), breaks = c(27,29,31)) + 
  #scale_x_datetime(name = "03/23-05/13/2021", date_breaks = "2 week") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.3,"cm")) +
  ylab("Temperature (°C)") + 
  geom_hline(yintercept=c(27,30.5,31.5), linetype="dashed", size=0.5) 

ggsave(tempLT_plot, height = 6 , width = 10, filename = "20210621_tempLT_tempprofileplot-tank.pdf", useDingbats=FALSE)

#########Calculate hours over threshold########
#subset into Treatment
tempC <-tempLT_summary[tempLT_summary$Treatment=='Control',]
tempOW<-tempLT_summary[tempLT_summary$Treatment=='OW',]
#subset into timepoints
tempC_mid<-  tempC[tempC$Date >= "2021-03-23" & tempC$Date <= "2021-04-23", ]
tempOW_mid<-tempOW[tempOW$Date >= "2021-03-23" &tempOW$Date<= "2021-04-23", ]

#total number of hours logged over threshold
tempC_end30.5<- (sum(tempC$mean > 30.5,na.rm=T)) 
tempOW_end30.5<- (sum(tempOW$mean > 30.5,na.rm=T))  #13 days = 312h
tempC_mid30.5<- (sum(tempC_mid$mean > 30.5,na.rm=T)) 
tempOW_mid30.5<- (sum(tempOW_mid$mean > 30.5,na.rm=T))  #4 days = 96h

#####PAR levels#####
PARAP<-read.csv("2021_apal_PARlevels.csv")
PARAP$PAR<-(as.numeric(PARAP$PAR))
PARAP<-na.omit(PARAP)
PARAP$Date<-as.POSIXct(paste(PARAP$Date), format = "%m/%d/%y") 
PARAP$date<-as.Date(PARAP$Date) 
PARAP$date<-as.factor(PARAP$date) 
PARAP$Tank<-as.factor(PARAP$Tank)
PARAP$Treatment<-as.factor(PARAP$Treatment)

#acer Control - 2, 4, 6, 7, 10; OW - 11, 14, 15, 17, 19
#apal Control - 12, 13, 15, 18, 20; OW - 1, 2, 5, 7, 9
parAC <- parAC[parAC$Tank=='2'| parAC$Tank=='4'| parAC$Tank=='6'| parAC$Tank=='7'| parAC$Tank=='10'| parAC$Tank=='11'| parAC$Tank=='14'|parAC$Tank=='15'| parAC$Tank=='17'| parAC$Tank=='19',]
parAP <- parAP[parAP$Tank=='12'| parAP$Tank=='13'| parAP$Tank=='15'| parAP$Tank=='18'| parAP$Tank=='20'| parAP$Tank=='1'| parAP$Tank=='2'|parAP$Tank=='5'| parAP$Tank=='7'| parAP$Tank=='9',]

PARAC2<- PARAC[PARAC$Treatment=="Control" | PARAC$Treatment== "OW", ]
PARAP2<- PARAP[PARAP$Treatment=="Control" | PARAP$Treatment== "OW",]

#plotting Treatments
PARAP_summary<- Summarize(PAR~Treatment+date, data=PARAP2, digits=3)

cols <- c("Control" = "blue", "OW" = "red")
PARLT_plot<-ggplot(PARAP_summary, aes(x=date, y=mean, col=Treatment, group=Treatment)) + scale_colour_manual(values=cols) +
  geom_line(aes(linetype=Treatment),size=1.5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color=Treatment, width=.4), alpha=.4)+
  scale_linetype_manual(values=c("solid","solid")) +
  scale_y_continuous(limits=c(50,650), breaks = c(100,200,300,400,500,600)) + 
  #scale_x_datetime(name = "03/23-05/13/2021", date_breaks = "2 week") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.3,"cm")) +
  labs(y = expression ("PAR (µmol photon"~m^2*")")) 

ggsave(PARLT_plot, height = 6 , width = 10, filename = "20210621_apal_PARprofileplot.pdf", useDingbats=FALSE)

#plotting each tank
PARLT_plot<-ggplot(parAC, aes(x=Date, y=PAR, col=Tank)) +
  geom_line(aes(linetype=Tank)) +
  scale_linetype_manual(values=c("dashed","dashed","dashed","dashed","dashed","solid","solid","solid","solid","solid")) +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.3,"cm")) +
  labs(y = expression ("PAR (µmol photon"~m^2*")")) 

ggsave(PARLT_plot, height = 6 , width = 10, filename = "20210621_acer_PARprofileplot-tank.pdf", useDingbats=FALSE)

#APAL
PARLT_plot<-ggplot(parAP, aes(x=Date, y=PAR, col=Tank)) +
  geom_line(aes(linetype=Tank)) +
  scale_linetype_manual(values=c("dashed","dashed","dashed","dashed","dashed","solid","solid","solid","solid","solid")) +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.3,"cm")) +
  labs(y = expression ("PAR (µmol photon"~m^2*")")) 

ggsave(PARLT_plot, height = 6 , width = 10, filename = "20210621_apal_PARprofileplot-tank.pdf", useDingbats=FALSE)

