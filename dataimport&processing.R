library(zoo)
library(dplyr)
library(plyr)

###############################################################################
##### DATA IMPORT AND PROCESSING BEFORE APPLYING THE IC-IS METHODOLOGY ########
###############################################################################

############################# Participant VS SIP Time stamps ##################

nbboPart<-read.csv('./inputs/cq20161003nbboPart.csv', sep = ',', header = TRUE)
nbboSIP<-read.csv('./inputs/cq20161003nbboSIP.csv', sep = ',', header = TRUE)

############################### IBM Natural Time 1-sec resolution #################
# in each second interval I take the median value 

IBM.part<-nbboPart[which(nbboPart[,2]=="IBM"),] 
IBM.sip<-nbboSIP[which(nbboSIP[,2]=="IBM"),] 

# create data set from 9.30 am to 16 pm 
h_am <- 3600*9.5 
h_pm <- 3600*16
Nr <- h_pm - h_am + 1 # number of rows of the array
Nc <- 3 # number of columns of the array

# participant data
array.IBMpart<-matrix(replicate(Nr*Nc,0), nrow=Nr, ncol = Nc)
colnames(array.IBMpart)<-c("timeMid","bBid","bOfr")
array.IBMpart[,1]<-c(h_am:h_pm)                                               

bBid.IBM<-IBM.part[which(!is.na(IBM.part[,"bBid"])),c("timeMid","bBid")]
bOfr.IBM<-IBM.part[which(!is.na(IBM.part[,"bOfr"])),c("timeMid","bOfr")]

# take median values in each second interval (in alternative, take the last available obs. in a given second interval)
for (i in 1:(Nr-1)) {
  array.IBMpart[i+1,"bBid"]<-median(bBid.IBM[which(bBid.IBM[,"timeMid"]>= array.IBMpart[i,"timeMid"] & bBid.IBM[,"timeMid"]<= array.IBMpart[i+1,1]),"bBid"])
  array.IBMpart[i+1,"bOfr"]<-median(bOfr.IBM[which(bOfr.IBM[,"timeMid"]>= array.IBMpart[i,"timeMid"] & bOfr.IBM[,"timeMid"]<= array.IBMpart[i+1,1]),"bOfr"]) 
  
}
array.IBMpart<-array.IBMpart[-1,]

# sip data
array.IBMsip<-matrix(replicate(Nr*Nc,0), nrow=Nr, ncol = Nc)
colnames(array.IBMsip)<-c("timeMid","bBid","bAsk")
array.IBMsip[,1]<-c(h_am:h_pm)                                               

bBid.IBM<-IBM.sip[which(!is.na(IBM.sip[,"bBid"])),c("timeMid","bBid")]
bAsk.IBM<-IBM.sip[which(!is.na(IBM.sip[,"bAsk"])),c("timeMid","bAsk")]

# take median values in each second interval (in alternative, take the last available obs. in a given second interval)
for (i in 1:(Nr-1)) {
  array.IBMsip[i+1,"bBid"]<-median(bBid.IBM[which(bBid.IBM[,"timeMid"]>= array.IBMsip[i,"timeMid"] & bBid.IBM[,"timeMid"]<= array.IBMsip[i+1,"timeMid"]),"bBid"])
  array.IBMsip[i+1,"bAsk"]<-median(bAsk.IBM[which(bAsk.IBM[,"timeMid"]>= array.IBMsip[i,"timeMid"] & bAsk.IBM[,"timeMid"]<= array.IBMsip[i+1,"timeMid"]),"bAsk"]) 
  
}
array.IBMsip<-array.IBMsip[-1,]

part.sip.IBMdata <- cbind.data.frame(array.IBMpart,array.IBMsip[,c("bBid","bAsk")])
df.partsip.1sec.IBM <- na.locf(part.sip.IBMdata) # data set for the analysis at the 1-sec level
colnames(df.partsip.1sec.IBM) <- c("timeMid","NBBpart","NBOpart","NBBsip","NBOsip")
write.csv(df.partsip.1sec.IBM, file = './data processed/df.partsip.1sec.IBM.csv', row.names = FALSE)

############################################ IBM Event-time ##########################################
#the time-counter t is incremented whenever there is an update to any variable 
#involved in the system.I increase the time counter whenever there is an update to ANY variable,
#and I propagate simultaneously the prices of the other variables which have no update. 
#Thus missing values are not excluded but replaced by previous prices. 

events.IBM <- full_join(IBM.part, IBM.sip, by = "timeMid")
events.IBM <- events.IBM[order(events.IBM$timeMid),]
events.IBM <- events.IBM[,c(3,4,5,8,9)]
events.IBM <- na.locf(events.IBM)
df.partsip.events.IBM <- events.IBM[which(events.IBM[,"timeMid"] > h_am & events.IBM[,"timeMid"]< h_pm),]
rownames(df.partsip.events.IBM) <- 1:nrow(df.events.IBM)
colnames(df.partsip.events.IBM) <- c("timeMid","NBBpart","NBOpart","NBBsip","NBOsip")

write.csv(df.partsip.events.IBM, file='./data processed/df.partsip.events.IBM.csv', row.names = FALSE)
###########################################################################################################

#################################### Primary Listing VS Others ################################

IBMraw.exc<-read.csv('./inputs/IBM20161003CQLexPart.csv', sep = ',', header = TRUE)

#################################  IBM Event-Time ##########################################
df.exc.event.IBM <- na.locf(IBMraw.exc)
write.csv(df.exc.event.IBM, file = './data processed/df.exc.events.IBM.csv', row.names = FALSE)

#################################### IBM Natural Time 1-sec ########################################
# create data set from 9.30 am to 16 pm 
h_am <- 3600*9.5 
h_pm <- 3600*16
Nr <- h_pm - h_am + 1 # number of rows of the array
Nc <- 5 # number of columns of the array

array.IBMexc<-matrix(replicate(Nr*Nc,0), nrow=Nr, ncol = Nc)
array.IBMexc[,1]<-c(h_am:h_pm)
colnames(array.IBMexc)<-c("timeMid","NBBother","NBOother","NBBlist","NBOlist")

NBBother<-IBMraw.exc[which(!is.na(IBMraw.exc[,4])),c(3,4)]
NBOother<-IBMraw.exc[which(!is.na(IBMraw.exc[,5])),c(3,5)]
NBBlist<-IBMraw.exc[which(!is.na(IBMraw.exc[,6])),c(3,6)]
NBOlist<-IBMraw.exc[which(!is.na(IBMraw.exc[,7])),c(3,7)]

for (i in 1:(Nr-1)) {
  array.IBMexc[i+1,"NBBother"]<-median(NBBother[which(NBBother[,1]>= array.IBMexc[i,1] & NBBother[,1]<= array.IBMexc[i+1,1]),2])
  array.IBMexc[i+1,"NBOother"]<-median(NBOother[which(NBOother[,1]>= array.IBMexc[i,1] & NBOother[,1]<= array.IBMexc[i+1,1]),2]) 
  array.IBMexc[i+1,"NBBlist"]<-median(NBBlist[which(NBBlist[,1]>= array.IBMexc[i,1] & NBBlist[,1]<= array.IBMexc[i+1,1]),2])
  array.IBMexc[i+1,"NBOlist"]<-median(NBOlist[which(NBOlist[,1]>= array.IBMexc[i,1] & NBOlist[,1]<= array.IBMexc[i+1,1]),2]) 
}

array.IBMexc<-as.data.frame(array.IBMexc[-1,])
df.exc.1sec.IBM <- na.locf(array.IBMexc)
write.csv(df.exc.1sec.IBM, file = './data processed/df.exc.1sec.IBM.csv', row.names = FALSE)

###############################################################################################

#################################### Quotes VS Lit/Dark Trades ################################

raw.qt<-read.csv('./inputs/ct20161003Part.csv', sep = ',', header = TRUE)

######################################## Event-Time ###########################################
IBM.part <- nbboPart[which(nbboPart[,2]=="IBM"),] ### participant time stamps for IBM
events.IBMpart <- IBM.part[,c(3,4,5)]
part.data <- na.locf(events.IBMpart)
part.data <- part.data[-1,]

darktrades <- raw.qt[which(raw.qt[,2]=="D" & raw.qt[,3]=="IBM"),c(4,5)]
colnames(darktrades) <- c("D","timeMid")

littrades <- raw.qt[which(raw.qt[,2]!="D" & raw.qt[,3]=="IBM"),c(4,5)]
colnames(littrades) <- c("Lit","timeMid")

trades <- full_join(littrades, darktrades, by = "timeMid")
qt <- full_join(trades, part.data, by = "timeMid")
qt <- qt[order(qt$timeMid),]
qt <- qt[which(qt[,"timeMid"] > h_am & qt[,"timeMid"]< h_pm),c(2,1,3,4,5)]
rownames(qt) <- c(1:nrow(qt))
df.qt.events.IBM <- na.locf(qt)

write.csv(df.qt.events.IBM, file = './data processed/df.qt.events.IBM.csv', row.names = FALSE)


######################################## Natural time ###########################################
array.IBMqt <- matrix(replicate(Nr*Nc,0), nrow=Nr, ncol = Nc)
colnames(array.IBMqt)<-c("timeMid","Lit","D","bBid","bOfr")
array.IBMqt[,1]<-c(h_am:h_pm)

for (i in 1:(Nr-1)) {
  array.IBMqt[i+1,"bBid"] <- median(bBid.IBM[which(bBid.IBM[,1]>= array.IBMqt[i,1] & bBid.IBM[,1]<= array.IBMqt[i+1,1]),"bBid"])
  array.IBMqt[i+1,"bOfr"] <- median(bOfr.IBM[which(bOfr.IBM[,1]>= array.IBMqt[i,1] & bOfr.IBM[,1]<= array.IBMqt[i+1,1]),"bOfr"])  
  array.IBMqt[i+1,"Lit"] <- median(littrades[which(littrades[,2]>= array.IBMqt[i,1] & littrades[,2]<= array.IBMqt[i+1,1]),"Lit"])
  array.IBMqt[i+1,"D"] <- median(darktrades[which(darktrades[,2]>= array.IBMqt[i,1] & darktrades[,2]<= array.IBMqt[i+1,1]),"D"]) 
}
array.IBMqt<-array.IBMqt[-1,]
df.qt.1sec.IBM<-as.data.frame(na.locf(array.IBMqt))

write.csv(df.qt.1sec.IBM, file = './data processed/df.qt.1sec.IBM.csv', row.names = FALSE)

rm(array.IBMexc, IBM.part, IBM.sip, NBBlist, NBBother, NBOlist, NBOother,
   array.IBMpart,array.IBMsip,bAsk.IBM,bBid.IBM,bOfr.IBM,IBM.part,IBM.sip,part.sip.IBMdata,
   littrades,darktrades,trades,part.data,raw.qt,array.IBMqt,IBMraw.exc,nbboPart,nbboSIP,events.IBMpart)
