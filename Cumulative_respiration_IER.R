#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#####################################Cumulative respiration#############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Define exposition periods
resp_all$Year<-as.numeric(substr(resp_all$time, 1,4))
resp_all$Month<-as.numeric(substr(resp_all$time, 7, 8))
resp_all$Wsl_cat<-ifelse(resp_all$Wsl<0, "Decreasing", "Increasing")
env$Year<-as.numeric(substr(env$time, 1,4))
env$Month<-as.numeric(substr(env$time, 6,7))
env$Day<-as.numeric(substr(env$time, 9,10))

#PLO_1
plo1<-subset(resp_all, id=="PLO_1")
plo1[(plo1$Year==2018 & plo1$Tair<15 & plo1$resp_cor>0.3), "outliers"]<-"YES"

##statistical model
plo1$Temp<-1/(plo1$Tair+273.15)/8.314
lm_plo1<-lm(log(resp_corr)~Temp+Year, plo1, subset = outliers=="NO")
summary(lm_plo1)

##Figure
plo1$pred<-exp(predict(lm_plo1, newdata = plo1))
ggplot(subset(plo1, outliers=="NO"), aes(Tair, resp_corr))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

##Environmental data
plo1_env<-env[, c("time", "pl1_t3", "Year", "Month", "Day")]

#IER exposition period
plo1_env$Period<-vector("character", length=nrow(plo1_env))
for(i in 1:nrow(plo1_env)){
  if(plo1_env$Year[i]==2017 & plo1_env$Month[i]>=5 & plo1_env$Day[i]>18){
    plo1_env[i, "Period"]<-c("Summer")
  }else{
    if(plo1_env$Year[i]==2017 & plo1_env$Month[i]<=11 & plo1_env$Day[i]<9){
      plo1_env[i, "Period"]<-c("Summer")
    }else{
      if(plo1_env$Year[i]==2017 & plo1_env$Month[i]>=11 & plo1_env$Day[i]>=9){
        plo1_env[i, "Period"]<-c("Winter")
      }
      if(plo1_env$Year[i]==2018 & plo1_env$Month[i]<=5 & plo1_env$Day[i]<22){
        plo1_env[i, "Period"]<-c("Summer")
      }else{
        
      }
    }
  }
}

summary(as.factor(plo1_env$Period))
nrow(plo1_env[(plo1_env$Year[i]==2018 & plo1_env$Month[i]>5),])
