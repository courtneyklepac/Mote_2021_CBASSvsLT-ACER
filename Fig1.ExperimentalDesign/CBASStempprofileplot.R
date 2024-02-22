#### CBASS plotting ####
library(ggplot2)
#SD card files
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
#plotting Fig. 1c
library(ggplot2)
cols <- c("27.5" = "navyblue", "28.9" = "blue", "30.2" = "yellow2", "31.6" = "goldenrod", 
          "32.9" = "darkorange2", "34.3" = "tomato", "35.6" = "red3", "37" = "firebrick")
CBASS_plot<-ggplot(CBASS_full, aes(x=date, y=Temp, col=Tank)) + scale_colour_manual(values=cols) +
  geom_line(aes(linetype=Tank)) +
  scale_linetype_manual(values=c("solid","solid","solid","solid","solid","solid","solid","solid")) +
  scale_y_continuous(limits=c(26.5,38), breaks = c(27,29,31,33,35,37,39)) + 
  scale_x_datetime(name = "06/21/2021", date_labels = "%Hh", date_breaks = "2 hour", limits=c(as.POSIXct("2021-06-21 11:40:01"), as.POSIXct("2021-06-22 09:00:01"))) +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.3,"cm")) +
  geom_hline(yintercept=c(27.5, 28.9,30.2,31.6,32.9,34.3,35.6,37), linetype="dashed", color = c("navyblue","blue","yellow2","goldenrod","darkorange2","tomato","red3","firebrick"), size=0.5) 

ggsave(CBASS_plot, height = 6 , width = 10, filename = "20210621_CBASS2_tempprofileplot.pdf", useDingbats=FALSE)
