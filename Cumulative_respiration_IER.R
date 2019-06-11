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
env$PL_period<-vector(length=nrow(env))
env$CT_period<-vector(length=nrow(env))

env[(plo1_env$Year==2017 & plo1_env$Month==5 & plo1_env$Day>18), "PL_period"]<-c("Summer")
env[(plo1_env$Year==2017 & plo1_env$Month>5 & plo1_env$Month<11), "PL_period"]<-c("Summer")
env[(plo1_env$Year==2017 & plo1_env$Month==11 & plo1_env$Day<9), "PL_period"]<-c("Summer")
env[(plo1_env$Year==2017 & plo1_env$Month==11 & plo1_env$Day>=9), "PL_period"]<-c("Winter")
env[(plo1_env$Year==2017 & plo1_env$Month==12), "PL_period"]<-c("Winter")
env[(plo1_env$Year==2018 & plo1_env$Month<5), "PL_period"]<-c("Winter")
env[(plo1_env$Year==2018 & plo1_env$Month==5 & plo1_env$Day<22), "PL_period"]<-c("Winter")
env[(plo1_env$Year==2018 & plo1_env$Month==5 & plo1_env$Day>=22), "PL_period"]<-c("Summer")
env[(plo1_env$Year==2018 & plo1_env$Month>5 & plo1_env$Month<11), "PL_period"]<-c("Summer")

summary(as.factor(env$PL_period))

env[(plo1_env$Year==2017 & plo1_env$Month==5 & plo1_env$Day>19), "CT_period"]<-c("Summer")
env[(plo1_env$Year==2017 & plo1_env$Month>5 & plo1_env$Month<11), "CT_period"]<-c("Summer")
env[(plo1_env$Year==2017 & plo1_env$Month>=11), "CT_period"]<-c("Winter")
env[(plo1_env$Year==2018 & plo1_env$Month<5), "CT_period"]<-c("Winter")
env[(plo1_env$Year==2018 & plo1_env$Month==5 & plo1_env$Day<17), "CT_period"]<-c("Winter")
env[(plo1_env$Year==2018 & plo1_env$Month==5 & plo1_env$Day>=17), "CT_period"]<-c("Summer")
env[(plo1_env$Year==2018 & plo1_env$Month>5 & plo1_env$Month<10), "CT_period"]<-c("Summer")
env[(plo1_env$Year==2018 & plo1_env$Month==10 & plo1_env$Day<30), "CT_period"]<-c("Summer")

summary(as.factor(env$CT_period))

#PLO_1
plo1<-subset(resp_all, id=="PLO_1")
plo1[(plo1$Year==2018 & plo1$Tair<15 & plo1$resp_corr>0.3), "outliers"]<-"YES"

##statistical model
plo1$Temp<-1/(plo1$Tair+273.15)/8.314
lm_plo1<-lm(log(resp)~Temp+Year, plo1, subset = outliers=="NO")
summary(lm_plo1)

##Figure
plo1$pred<-exp(predict(lm_plo1, newdata = plo1))
ggplot(subset(plo1, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

##Calculating cumulative respiration
###g(CO2)/m2/15 minutes (which is the measurement interval)
plo1_env<-env[, c("pl1_t3", "Year", "PL_period")]
plo1_env$Temp<-1/(plo1_env$pl1_t3+273.15)/8.314

plo1_env$resp_pred<-exp(predict(lm_plo1, newdata = plo1_env))*60*15#per 15 minutes of the measurement
plo1_env$resp_pred.se<-exp(predict(lm_plo1, newdata = plo1_env)+
                             predict(lm_plo1, newdata = plo1_env, se.fit = T)$se.fit)*60*15-
  exp(predict(lm_plo1, newdata = plo1_env))*60*15#standard error
###Summer 2017
plo1_envS2017<-subset(plo1_env, Year==2017 & PL_period=="Summer")
####sum the produced CO2
plo1_S2017<-sum(plo1_envS2017$resp_pred)#574.7831 µmol CO2/m2
plo1_S2017.se<-sum(plo1_envS2017$resp_pred.se)#22.99 µmol CO2/m2
###recalculate to mmols of C per transplant 
plo1_S2017<-plo1_S2017*(67/10000)
plo1_S2017.se<-plo1_S2017.se*(67/10000)


###################################Generating all predictions models###################################
#PLO_2
plo2<-subset(resp_all, id=="PLO_2")
#plo1[(plo1$Year==2018 & plo1$Tair<15 & plo1$resp_cor>0.3), "outliers"]<-"YES"

##statistical model
plo2$Temp<-1/(plo2$Tair+273.15)/8.314
lm_plo2<-lm(log(resp)~Temp+Year, plo2, subset = outliers=="NO")
summary(lm_plo2)

##Figure
plo2$pred<-exp(predict(lm_plo2, newdata = plo2))
ggplot(subset(plo2), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLO_3
plo3<-subset(resp_all, id=="PLO_3")
#plo1[(plo1$Year==2018 & plo1$Tair<15 & plo1$resp_cor>0.3), "outliers"]<-"YES"

##statistical model
plo3$Temp<-1/(plo3$Tair+273.15)/8.314
lm_plo3<-lm(log(resp)~Temp+Year, plo3, subset = outliers=="NO")
summary(lm_plo3)

##Figure
plo3$pred<-exp(predict(lm_plo3, newdata = plo3))
ggplot(subset(plo3), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLO_4
plo4<-subset(resp_all, id=="PLO_4")
plo4[(plo4$Year==2017 & plo4$Tair<20 & plo4$Tair>10 & plo4$resp_corr>0.5), "outliers"]<-"YES"

##statistical model
plo4$Temp<-1/(plo4$Tair+273.15)/8.314
lm_plo4<-lm(log(resp)~Temp+Year, plo4, subset = outliers=="NO")
summary(lm_plo4)

##Figure
plo4$pred<-exp(predict(lm_plo4, newdata = plo4))
ggplot(subset(plo4, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLO_5
plo5<-subset(resp_all, id=="PLO_5")
#plo4[(plo4$Year==2017 & plo4$Tair<20 & plo4$Tair>10 & plo4$resp_cor>0.5), "outliers"]<-"YES"

##statistical model
plo5$Temp<-1/(plo5$Tair+273.15)/8.314
lm_plo5<-lm(log(resp)~Temp, plo5, subset = outliers=="NO")
summary(lm_plo5)

##Figure
plo5$pred<-exp(predict(lm_plo5, newdata = plo5))
ggplot(subset(plo5, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLO_6
plo6<-subset(resp_all, id=="PLO_6")
#plo4[(plo4$Year==2017 & plo4$Tair<20 & plo4$Tair>10 & plo4$resp_cor>0.5), "outliers"]<-"YES"

##statistical model
plo6$Temp<-1/(plo6$Tair+273.15)/8.314
lm_plo6<-lm(log(resp)~Temp+Year, plo6, subset = outliers=="NO")
summary(lm_plo6)

##Figure
plo6$pred<-exp(predict(lm_plo6, newdata = plo6))
ggplot(subset(plo6, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLO_7
plo7<-subset(resp_all, id=="PLO_7")
#plo4[(plo4$Year==2017 & plo4$Tair<20 & plo4$Tair>10 & plo4$resp_cor>0.5), "outliers"]<-"YES"

##statistical model
plo7$Temp<-1/(plo7$Tair+273.15)/8.314
lm_plo7<-lm(log(resp)~Temp+Year, plo7, subset = outliers=="NO")
summary(lm_plo7)

##Figure
plo7$pred<-exp(predict(lm_plo7, newdata = plo7))
ggplot(subset(plo7, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLO_8
plo8<-subset(resp_all, id=="PLO_8" & !is.na(Tsoil))
plo8[(plo8$Year==2018 & plo8$Tair<20 & plo8$resp_corr>0.2), "outliers"]<-"YES"

##statistical model
plo8$Temp<-1/(plo8$Tair+273.15)/8.314
lm_plo8<-lm(log(resp)~Temp+Year, plo8, subset = outliers=="NO")
summary(lm_plo8)

##Figure
plo8$pred<-exp(predict(lm_plo8, newdata = plo8))
ggplot(subset(plo8, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLO_9
plo9<-subset(resp_all, id=="PLO_9" & !is.na(Tsoil))
#plo9[(plo8$Year==2018 & plo8$Tair<20 & plo8$resp>0.2), "outliers"]<-"YES"

##statistical model
plo9$Temp<-1/(plo9$Tair+273.15)/8.314
lm_plo9<-lm(log(resp)~Temp+Year, plo9, subset = outliers=="NO")
summary(lm_plo9)

##Figure
plo9$pred<-exp(predict(lm_plo9, newdata = plo9))
ggplot(subset(plo9, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLO_10
plo10<-subset(resp_all, id=="PLO_10")
plo10[(plo10$Year==2018 & plo10$Tair<20 & plo10$Tair>10 & plo10$resp_corr>0.275), "outliers"]<-"YES"

##statistical model
plo10$Temp<-1/(plo10$Tair+273.15)/8.314
lm_plo10<-lm(log(resp)~Temp, plo10, subset = outliers=="NO")
summary(lm_plo10)

##Figure
plo10$pred<-exp(predict(lm_plo10, newdata = plo10))
ggplot(subset(plo10, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLA_1
pla1<-subset(resp_all, id=="PLA_1")
pla1[(pla1$Year==2018 & pla1$resp_corr>1), "outliers"]<-"YES"

##statistical model
pla1$Temp<-1/(pla1$Tair+273.15)/8.314
lm_pla1<-lm(log(resp)~Temp, pla1, subset = outliers=="NO")
summary(lm_pla1)

##Figure
pla1$pred<-exp(predict(lm_pla1, newdata = pla1))
ggplot(subset(pla1, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLA_2
pla2<-subset(resp_all, id=="PLA_2")
#pla1[(pla1$Year==2018 & pla1$resp>1), "outliers"]<-"YES"

##statistical model
pla2$Temp<-1/(pla2$Tair+273.15)/8.314
lm_pla2<-lm(log(resp)~Temp+Year, pla2, subset = outliers=="NO")
summary(lm_pla2)

##Figure
pla2$pred<-exp(predict(lm_pla2, newdata = pla2))
ggplot(subset(pla2, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLA_3
pla3<-subset(resp_all, id=="PLA_3")
pla3[(pla3$Tair<20 & pla3$resp_corr>0.4), "outliers"]<-"YES"

##statistical model
pla3$Temp<-1/(pla3$Tair+273.15)/8.314
lm_pla3<-lm(log(resp)~Temp, pla3, subset = outliers=="NO")
summary(lm_pla3)

##Figure
pla3$pred<-exp(predict(lm_pla3, newdata = pla3))
ggplot(subset(pla3, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLA_4
pla4<-subset(resp_all, id=="PLA_4")
pla4[(pla4$Tair<20 & pla4$resp_corr>0.9), "outliers"]<-"YES"

##statistical model
pla4$Temp<-1/(pla4$Tair+273.15)/8.314
lm_pla4<-lm(log(resp)~Temp+Year, pla4, subset = outliers=="NO")
summary(lm_pla4)

##Figure
pla4$pred<-exp(predict(lm_pla4, newdata = pla4))
ggplot(subset(pla4, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLA_5
pla5<-subset(resp_all, id=="PLA_5")
#pla4[(pla4$Tair<20 & pla4$resp>0.9), "outliers"]<-"YES"

##statistical model
pla5$Temp<-1/(pla5$Tair+273.15)/8.314
lm_pla5<-lm(log(resp)~Temp, pla5, subset = outliers=="NO")
summary(lm_pla5)

##Figure
pla5$pred<-exp(predict(lm_pla5, newdata = pla5))
ggplot(subset(pla5, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLA_6
pla6<-subset(resp_all, id=="PLA_6")
#pla4[(pla4$Tair<20 & pla4$resp>0.9), "outliers"]<-"YES"

##statistical model
pla6$Temp<-1/(pla6$Tair+273.15)/8.314
lm_pla6<-lm(log(resp)~Temp, pla6, subset = outliers=="NO")
summary(lm_pla6)

##Figure
pla6$pred<-exp(predict(lm_pla6, newdata = pla6))
ggplot(subset(pla6, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLA_7
pla7<-subset(resp_all, id=="PLA_7")
#pla4[(pla4$Tair<20 & pla4$resp>0.9), "outliers"]<-"YES"

##statistical model
pla7$Temp<-1/(pla7$Tair+273.15)/8.314
lm_pla7<-lm(log(resp)~Temp, pla7, subset = outliers=="NO")
summary(lm_pla7)

##Figure
pla7$pred<-exp(predict(lm_pla7, newdata = pla7))
ggplot(subset(pla7, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLA_8
pla8<-subset(resp_all, id=="PLA_8")
#pla4[(pla4$Tair<20 & pla4$resp>0.9), "outliers"]<-"YES"

##statistical model
pla8$Temp<-1/(pla8$Tair+273.15)/8.314
lm_pla8<-lm(log(resp)~Temp+Year, pla8, subset = outliers=="NO")
summary(lm_pla8)

##Figure
pla8$pred<-exp(predict(lm_pla8, newdata = pla8))
ggplot(subset(pla8, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLA_9
pla9<-subset(resp_all, id=="PLA_9")
#pla4[(pla4$Tair<20 & pla4$resp>0.9), "outliers"]<-"YES"

##statistical model
pla9$Temp<-1/(pla9$Tair+273.15)/8.314
lm_pla9<-lm(log(resp)~Temp, pla9, subset = outliers=="NO")
summary(lm_pla9)

##Figure
pla9$pred<-exp(predict(lm_pla9, newdata = pla9))
ggplot(subset(pla9, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLA_10
pla10<-subset(resp_all, id=="PLA_10")
#pla4[(pla4$Tair<20 & pla4$resp>0.9), "outliers"]<-"YES"

##statistical model
pla10$Temp<-1/(pla10$Tair+273.15)/8.314
lm_pla10<-lm(log(resp)~Temp+Year, pla10, subset = outliers=="NO")
summary(lm_pla10)

##Figure
pla10$pred<-exp(predict(lm_pla10, newdata = pla10))
ggplot(subset(pla10, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLC_1
plc1<-subset(resp_all, id=="PLC_1")
#pla4[(pla4$Tair<20 & pla4$resp>0.9), "outliers"]<-"YES"

##statistical model
plc1$Temp<-1/(plc1$Tair+273.15)/8.314
lm_plc1<-lm(log(resp)~Temp, plc1, subset = outliers=="NO")
summary(lm_plc1)

##Figure
plc1$pred<-exp(predict(lm_plc1, newdata = plc1))
ggplot(subset(plc1, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLC_2
plc2<-subset(resp_all, id=="PLC_2")
#pla4[(pla4$Tair<20 & pla4$resp>0.9), "outliers"]<-"YES"

##statistical model
plc2$Temp<-1/(plc2$Tair+273.15)/8.314
lm_plc2<-lm(log(resp)~Temp+Year, plc2, subset = outliers=="NO")
summary(lm_plc2)

##Figure
plc2$pred<-exp(predict(lm_plc2, newdata = plc2))
ggplot(subset(plc2, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLC_3
plc3<-subset(resp_all, id=="PLC_3")
#pla4[(pla4$Tair<20 & pla4$resp>0.9), "outliers"]<-"YES"

##statistical model
plc3$Temp<-1/(plc3$Tair+273.15)/8.314
lm_plc3<-lm(log(resp)~Temp, plc3, subset = outliers=="NO")
summary(lm_plc3)

##Figure
plc3$pred<-exp(predict(lm_plc3, newdata = plc3))
ggplot(subset(plc3, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLC_4
plc4<-subset(resp_all, id=="PLC_4")
plc4[(plc4$Tair>20 & plc4$Year==2017 & plc4$resp_corr<0.9), "outliers"]<-"YES"

##statistical model
plc4$Temp<-1/(plc4$Tair+273.15)/8.314
lm_plc4<-lm(log(resp)~Temp, plc4, subset = outliers=="NO")
summary(lm_plc4)

##Figure
plc4$pred<-exp(predict(lm_plc4, newdata = plc4))
ggplot(subset(plc4, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#PLC_5
plc5<-subset(resp_all, id=="PLC_5")
plc5[(plc5$Tair>25 & plc5$Year==2018 & plc5$resp_corr<1.25), "outliers"]<-"YES"
plc5[(plc5$Tair<15 & plc5$Year==2018 & plc5$resp_corr<0.5), "outliers"]<-"YES"

##statistical model
plc5$Temp<-1/(plc5$Tair+273.15)/8.314
lm_plc5<-lm(log(resp)~Temp+Year, plc5, subset = outliers=="NO")
summary(lm_plc5)

##Figure
plc5$pred<-exp(predict(lm_plc5, newdata = plc5))
ggplot(subset(plc5, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTO_1
cto1<-subset(resp_all, id=="CTO_1")
cto1[(cto1$Tair>25 & cto1$resp_corr<0.1), "outliers"]<-"YES"

##statistical model
cto1$Temp<-1/(cto1$Tair+273.15)/8.314
lm_cto1<-lm(log(resp)~Temp+Year, cto1, subset = outliers=="NO")
summary(lm_cto1)

##Figure
cto1$pred<-exp(predict(lm_cto1, newdata = cto1))
ggplot(subset(cto1, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTO_2
cto2<-subset(resp_all, id=="CTO_2")
cto2[(cto2$Tair>25 & cto2$resp<0.2 & cto2$Year==2018), "outliers"]<-"YES"

##statistical model
cto2$Temp<-1/(cto2$Tair+273.15)/8.314
lm_cto2<-lm(log(resp)~Temp+Year, cto2, subset = outliers=="NO")
summary(lm_cto2)

##Figure
cto2$pred<-exp(predict(lm_cto2, newdata = cto2))
ggplot(subset(cto2, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTO_3
cto3<-subset(resp_all, id=="CTO_3")
#cto2[(cto2$Tair>25 & cto2$resp<0.2 & cto2$Year==2018), "outliers"]<-"YES"

##statistical model
cto3$Temp<-1/(cto3$Tair+273.15)/8.314
lm_cto3<-lm(log(resp)~Temp+Year, cto3, subset = outliers=="NO")
summary(lm_cto3)

##Figure
cto3$pred<-exp(predict(lm_cto3, newdata = cto3))
ggplot(subset(cto3, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTO_4
cto4<-subset(resp_all, id=="CTO_4")
#cto2[(cto2$Tair>25 & cto2$resp<0.2 & cto2$Year==2018), "outliers"]<-"YES"

##statistical model
cto4$Temp<-1/(cto4$Tair+273.15)/8.314
lm_cto4<-lm(log(resp)~Temp+Year, cto4, subset = outliers=="NO")
summary(lm_cto4)

##Figure
cto4$pred<-exp(predict(lm_cto4, newdata = cto4))
ggplot(subset(cto4, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTO_5
cto5<-subset(resp_all, id=="CTO_5")
#cto2[(cto2$Tair>25 & cto2$resp<0.2 & cto2$Year==2018), "outliers"]<-"YES"

##statistical model
cto5$Temp<-1/(cto5$Tair+273.15)/8.314
lm_cto5<-lm(log(resp)~Temp, cto5, subset = outliers=="NO")
summary(lm_cto5)

##Figure
cto5$pred<-exp(predict(lm_cto5, newdata = cto5))
ggplot(subset(cto5, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTO_6
cto6<-subset(resp_all, id=="CTO_6")
#cto2[(cto2$Tair>25 & cto2$resp<0.2 & cto2$Year==2018), "outliers"]<-"YES"

##statistical model
cto6$Temp<-1/(cto6$Tair+273.15)/8.314
lm_cto6<-lm(log(resp)~Temp+Year, cto6, subset = outliers=="NO")
summary(lm_cto6)

##Figure
cto6$pred<-exp(predict(lm_cto6, newdata = cto6))
ggplot(subset(cto6, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTO_7
cto7<-subset(resp_all, id=="CTO_7")
#cto2[(cto2$Tair>25 & cto2$resp<0.2 & cto2$Year==2018), "outliers"]<-"YES"

##statistical model
cto7$Temp<-1/(cto7$Tair+273.15)/8.314
lm_cto7<-lm(log(resp)~Temp+Year, cto7, subset = outliers=="NO")
summary(lm_cto7)

##Figure
cto7$pred<-exp(predict(lm_cto7, newdata = cto7))
ggplot(subset(cto7, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTO_8
cto8<-subset(resp_all, id=="CTO_8")
cto8[(cto8$Tair<20 & cto8$resp_corr>0.4 & cto8$Year==2017), "outliers"]<-"YES"

##statistical model
cto8$Temp<-1/(cto8$Tair+273.15)/8.314
lm_cto8<-lm(log(resp)~Temp+Year, cto8, subset = outliers=="NO")
summary(lm_cto8)

##Figure
cto8$pred<-exp(predict(lm_cto8, newdata = cto8))
ggplot(subset(cto8, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTO_9
cto9<-subset(resp_all, id=="CTO_9")
#cto8[(cto8$Tair<20 & cto8$resp>0.4 & cto8$Year==2017), "outliers"]<-"YES"

##statistical model
cto9$Temp<-1/(cto9$Tair+273.15)/8.314
lm_cto9<-lm(log(resp)~Temp+Year, cto9, subset = outliers=="NO")
summary(lm_cto9)

##Figure
cto9$pred<-exp(predict(lm_cto9, newdata = cto9))
ggplot(subset(cto9, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTO_10
cto10<-subset(resp_all, id=="CTO_10")
#cto8[(cto8$Tair<20 & cto8$resp>0.4 & cto8$Year==2017), "outliers"]<-"YES"

##statistical model
cto10$Temp<-1/(cto10$Tair+273.15)/8.314
lm_cto10<-lm(log(resp)~Temp, cto10, subset = outliers=="NO")
summary(lm_cto10)

##Figure
cto10$pred<-exp(predict(lm_cto10, newdata = cto10))
ggplot(subset(cto10, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTA_1
cta1<-subset(resp_all, id=="CTA_1")
cta1[(cta1$resp_corr>0.55), "outliers"]<-"YES"

##statistical model
cta1$Temp<-1/(cta1$Tair+273.15)/8.314
lm_cta1<-lm(log(resp)~Temp+Year, cta1, subset = outliers=="NO")
summary(lm_cta1)

##Figure
cta1$pred<-exp(predict(lm_cta1, newdata = cta1))
ggplot(subset(cta1, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTA_2
cta2<-subset(resp_all, id=="CTA_2")
#cta1[(cta1$resp>0.55), "outliers"]<-"YES"

##statistical model
cta2$Temp<-1/(cta2$Tair+273.15)/8.314
lm_cta2<-lm(log(resp)~Temp, cta2, subset = outliers=="NO")
summary(lm_cta2)

##Figure
cta2$pred<-exp(predict(lm_cta2, newdata = cta2))
ggplot(subset(cta2, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTA_3
cta3<-subset(resp_all, id=="CTA_3")
cta3[(cta3$resp_corr>0.5 & cta3$Year==2018), "outliers"]<-"YES"

##statistical model
cta3$Temp<-1/(cta3$Tair+273.15)/8.314
lm_cta3<-lm(log(resp)~Temp+Year, cta3, subset = outliers=="NO")
summary(lm_cta3)

##Figure
cta3$pred<-exp(predict(lm_cta3, newdata = cta3))
ggplot(subset(cta3, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTA_4
cta4<-subset(resp_all, id=="CTA_4")
cta4[(cta4$resp_corr>0.75 & cta4$Year==2017), "outliers"]<-"YES"

##statistical model
cta4$Temp<-1/(cta4$Tair+273.15)/8.314
lm_cta4<-lm(log(resp)~Temp+Year, cta4, subset = outliers=="NO")
summary(lm_cta4)

##Figure
cta4$pred<-exp(predict(lm_cta4, newdata = cta4))
ggplot(subset(cta4, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTA_5
cta5<-subset(resp_all, id=="CTA_5")
cta5[(cta5$resp_corr<0.4 & cta5$Year==2017 & cta5$Tair>20), "outliers"]<-"YES"

##statistical model
cta5$Temp<-1/(cta5$Tair+273.15)/8.314
lm_cta5<-lm(log(resp)~Temp+Year, cta5, subset = outliers=="NO")
summary(lm_cta5)

##Figure
cta5$pred<-exp(predict(lm_cta5, newdata = cta5))
ggplot(subset(cta5, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTA_6
cta6<-subset(resp_all, id=="CTA_6")
#cta5[(cta5$resp<0.4 & cta5$Year==2017 & cta5$Tair>20), "outliers"]<-"YES"

##statistical model
cta6$Temp<-1/(cta6$Tair+273.15)/8.314
lm_cta6<-lm(log(resp)~Temp+Year, cta6, subset = outliers=="NO")
summary(lm_cta6)

##Figure
cta6$pred<-exp(predict(lm_cta6, newdata = cta6))
ggplot(subset(cta6, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTA_7
cta7<-subset(resp_all, id=="CTA_7")
#cta5[(cta5$resp<0.4 & cta5$Year==2017 & cta5$Tair>20), "outliers"]<-"YES"

##statistical model
cta7$Temp<-1/(cta7$Tair+273.15)/8.314
lm_cta7<-lm(log(resp)~Temp+Year, cta7, subset = outliers=="NO")
summary(lm_cta7)

##Figure
cta7$pred<-exp(predict(lm_cta7, newdata = cta7))
ggplot(subset(cta7, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTA_8
cta8<-subset(resp_all, id=="CTA_8")
#cta5[(cta5$resp<0.4 & cta5$Year==2017 & cta5$Tair>20), "outliers"]<-"YES"

##statistical model
cta8$Temp<-1/(cta8$Tair+273.15)/8.314
lm_cta8<-lm(log(resp)~Temp+Year, cta8, subset = outliers=="NO")
summary(lm_cta8)

##Figure
cta8$pred<-exp(predict(lm_cta8, newdata = cta8))
ggplot(subset(cta8, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTA_9
cta9<-subset(resp_all, id=="CTA_9")
#cta5[(cta5$resp<0.4 & cta5$Year==2017 & cta5$Tair>20), "outliers"]<-"YES"

##statistical model
cta9$Temp<-1/(cta9$Tair+273.15)/8.314
lm_cta9<-lm(log(resp)~Temp+Year, cta9, subset = outliers=="NO")
summary(lm_cta9)

##Figure
cta9$pred<-exp(predict(lm_cta9, newdata = cta9))
ggplot(subset(cta9, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTA_10
cta10<-subset(resp_all, id=="CTA_10")
#cta5[(cta5$resp<0.4 & cta5$Year==2017 & cta5$Tair>20), "outliers"]<-"YES"

##statistical model
cta10$Temp<-1/(cta10$Tair+273.15)/8.314
lm_cta10<-lm(log(resp)~Temp+Year, cta10, subset = outliers=="NO")
summary(lm_cta10)

##Figure
cta10$pred<-exp(predict(lm_cta10, newdata = cta10))
ggplot(subset(cta10, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTC_1
ctc1<-subset(resp_all, id=="CTC_1")
#cta5[(cta5$resp<0.4 & cta5$Year==2017 & cta5$Tair>20), "outliers"]<-"YES"

##statistical model
ctc1$Temp<-1/(ctc1$Tair+273.15)/8.314
lm_ctc1<-lm(log(resp)~Temp, ctc1, subset = outliers=="NO")
summary(lm_ctc1)

##Figure
ctc1$pred<-exp(predict(lm_ctc1, newdata = ctc1))
ggplot(subset(ctc1, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTC_2
ctc2<-subset(resp_all, id=="CTC_2")
#cta5[(cta5$resp<0.4 & cta5$Year==2017 & cta5$Tair>20), "outliers"]<-"YES"

##statistical model
ctc2$Temp<-1/(ctc2$Tair+273.15)/8.314
lm_ctc2<-lm(log(resp)~Temp, ctc2, subset = outliers=="NO")
summary(lm_ctc2)

##Figure
ctc2$pred<-exp(predict(lm_ctc2, newdata = ctc2))
ggplot(subset(ctc2, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTC_3
ctc3<-subset(resp_all, id=="CTC_3")
#cta5[(cta5$resp<0.4 & cta5$Year==2017 & cta5$Tair>20), "outliers"]<-"YES"

##statistical model
ctc3$Temp<-1/(ctc3$Tair+273.15)/8.314
lm_ctc3<-lm(log(resp)~Temp+Year, ctc3, subset = outliers=="NO")
summary(lm_ctc3)

##Figure
ctc3$pred<-exp(predict(lm_ctc3, newdata = ctc3))
ggplot(subset(ctc3, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTC_4
ctc4<-subset(resp_all, id=="CTC_4")
#cta5[(cta5$resp<0.4 & cta5$Year==2017 & cta5$Tair>20), "outliers"]<-"YES"

##statistical model
ctc4$Temp<-1/(ctc4$Tair+273.15)/8.314
lm_ctc4<-lm(log(resp)~Temp+Year, ctc4, subset = outliers=="NO")
summary(lm_ctc4)

##Figure
ctc4$pred<-exp(predict(lm_ctc4, newdata = ctc4))
ggplot(subset(ctc4, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#CTC_5
ctc5<-subset(resp_all, id=="CTC_5")
#cta5[(cta5$resp<0.4 & cta5$Year==2017 & cta5$Tair>20), "outliers"]<-"YES"

##statistical model
ctc5$Temp<-1/(ctc5$Tair+273.15)/8.314
lm_ctc5<-lm(log(resp)~Temp, ctc5, subset = outliers=="NO")
summary(lm_ctc5)

##Figure
ctc5$pred<-exp(predict(lm_ctc5, newdata = ctc5))
ggplot(subset(ctc5, outliers=="NO"), aes(Tair, resp))+geom_point(cex=4, aes(colour=as.factor(Year)))+
  geom_line(aes(Tair, pred, colour=as.factor(Year)))

#######################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Defining the list of models
Mlist<-data.frame(Catchment = c(rep("Plesne", 25), rep("Certovo", 25)),
                  Horizon = c(rep("Litter", 10), rep("Organic soil", 10), rep("Control", 5),
                              rep("Litter", 10), rep("Organic soil", 10), rep("Control", 5)),
                  Soil = c(rep("Plesne", 5), rep("Certovo", 5), rep("Plesne", 5), rep("Certovo", 5), rep("Plesne", 5),
                           rep("Plesne", 5), rep("Certovo", 5), rep("Plesne", 5), rep("Certovo", 5), rep("Certovo", 5)),
                  Senzor = c(rep(c("pl1_t3", "pl2_t3", "pl3_t3", "pl4_t3", "pl5_t3"), times = 5),
                             rep(c("ct1_t3", "ct2_t3", "ct3_t3", "ct4_t3", "ct5_t3"), times = 5)),
                  Area = c(rep(67, 10), rep(111, 10), rep(78, 5),
                           rep(67, 10), rep(111, 10), rep(78, 5)),
                  ID = c("PLO_1", "PLO_2", "PLO_3", "PLO_4", "PLO_5",
                         "CTO_1", "CTO_2", "CTO_3", "CTO_4", "CTO_5",
                         "PLA_1", "PLA_2", "PLA_3", "PLA_4", "PLA_5",
                         "CTA_1", "CTA_2", "CTA_3", "CTA_4", "CTA_5",
                         "PLC_1", "PLC_2", "PLC_3", "PLC_4", "PLC_5",
                         "PLO_6", "PLO_7", "PLO_8", "PLO_9", "PLO_10",
                         "CTO_6", "CTO_7", "CTO_8", "CTO_9", "CTO_10",
                         "PLA_6", "PLA_7", "PLA_8", "PLA_9", "PLA_10",
                         "CTA_6", "CTA_7", "CTA_8", "CTA_9", "CTA_10",
                         "CTC_1", "CTC_2", "CTC_3", "CTC_4", "CTC_5"))
Mlist$Origin <- ifelse(Mlist$Catchment==Mlist$Soil, "Native", "Transplanted")

Model_list<-list(lm_plo1, lm_plo2, lm_plo3, lm_plo4, lm_plo5,
                 lm_cto1, lm_cto2, lm_cto3, lm_cto4, lm_cto5,
                 lm_pla1, lm_pla2, lm_pla3, lm_pla4, lm_pla5,
                 lm_cta1, lm_cta2, lm_cta3, lm_cta4, lm_cta5,
                 lm_plc1, lm_plc2, lm_plc3, lm_plc4, lm_plc5,
                 lm_plo6, lm_plo7, lm_plo8, lm_plo9, lm_plo1,
                 lm_cto6, lm_cto7, lm_cto8, lm_cto9, lm_cto1,
                 lm_pla6, lm_pla7, lm_pla8, lm_pla9, lm_pla1,
                 lm_cta6, lm_cta7, lm_cta8, lm_cta9, lm_cta10,
                 lm_ctc1, lm_ctc2, lm_ctc3, lm_ctc4, lm_ctc5)
#Checking env data frame
##Removing all undefined rows
env_calc<-env[!(env$PL_period=="FALSE" | env$CT_period=="FALSE"), ]
summary(env_calc)


#Do the calculation 
##Defining the data frame storing the results
Cmin_field<-data.frame(Cmin = vector("numeric"),
                       Cmic.se = vector("numeric"),
                       Period = vector("character"),
                       Year = vector("numeric"),
                       Soil = vector("character"),
                       Catchment = vector("character"),
                       Horizon = vector("character"),
                       Origin = vector("character"),
                       ID = vector("character"))

for(i in 1:nrow(Mlist)){
  if(Mlist$Catchment[i]=="Plesne"){
    ##Calculating cumulative respiration
    ###g(CO2)/m2/15 minutes (which is the measurement interval)
    e<-data.frame(Tair = env_calc[, as.character(Mlist$Senzor[i])],
                  Year = env_calc[, "Year"],
                  PL_period = env_calc[, "PL_period"])
    e$Temp<-1/(e$Tair+273.15)/8.314
    
    e$resp_pred<-exp(predict(Model_list[[i]], newdata = e))/4#per 15 minutes of the measurement
    e$resp_pred.se<-exp(predict(Model_list[[i]], newdata = e)+
                          predict(Model_list[[i]], newdata = e, se.fit = T)$se.fit)/4-
      exp(predict(Model_list[[i]], newdata = e))/4#standard error
    ###Summer 2017
    S2017d<-sum(e[(e$Year==2017 & e$PL_period=="Summer"), "resp_pred"], na.rm=T)/44*1e3
    S2017se<-sum(e[(e$Year==2017 & e$PL_period=="Summer"), "resp_pred.se"], na.rm=T)/44*1e3
    
    ###Winter 2017 - 2018
    W20178d<-(sum(e[(e$Year==2017 & e$PL_period=="Winter"), "resp_pred"], na.rm=T)+
                sum(e[(e$Year==2018 & e$PL_period=="Winter"), "resp_pred"], na.rm=T))/44*1e3
    W20178se<-(sum(e[(e$Year==2017 & e$PL_period=="Winter"), "resp_pred.se"], na.rm=T)+
                 sum(e[(e$Year==2018 & e$PL_period=="Winter"), "resp_pred.se"], na.rm=T))/44*1e3
    
    ###Summer 2018
    S2018d<-sum(e[(e$Year==2018 & e$PL_period=="Summer"), "resp_pred"], na.rm=T)/44*1e3
    S2018se<-sum(e[(e$Year==2018 & e$PL_period=="Summer"), "resp_pred.se"], na.rm=T)/44*1e3
    
    #Store the results
    Cmin_field<-rbind(Cmin_field, data.frame(Cmin = c(S2017d, W20178d, S2018d),
                                             Cmic.se = c(S2017se, W20178se, S2018se),
                                             Period = c("Summer", "Winter", "Summer"),
                                             Year = c(2017, 2018, 2018),
                                             Soil = rep(Mlist$Soil[i], times=3),
                                             Catchment = rep(Mlist$Catchment[i], times=3),
                                             Horizon = rep(Mlist$Horizon[i], times=3),
                                             Origin = rep(Mlist$Origin[i], times=3),
                                             ID = rep(Mlist$ID[i], times=3)))
  }else{
    ##Calculating cumulative respiration
    ###g(CO2)/m2/15 minutes (which is the measurement interval)
    e<-data.frame(Tair = env_calc[, as.character(Mlist$Senzor[i])],
                  Year = env_calc[, "Year"],
                  CT_period = env_calc[, "CT_period"])
    e$Temp<-1/(e$Tair+273.15)/8.314
    
    e$resp_pred<-exp(predict(Model_list[[i]], newdata = e))/4#per 15 minutes of the measurement
    e$resp_pred.se<-exp(predict(Model_list[[i]], newdata = e)+
                          predict(Model_list[[i]], newdata = e, se.fit = T)$se.fit)/4-
      exp(predict(Model_list[[i]], newdata = e))/4#standard error
    ###Summer 2017
    S2017d<-sum(e[(e$Year==2017 & e$CT_period=="Summer"), "resp_pred"], na.rm=T)/44*1e3
    S2017se<-sum(e[(e$Year==2017 & e$CT_period=="Summer"), "resp_pred.se"], na.rm=T)/44*1e3
    
    ###Winter 2017 - 2018
    W20178d<-(sum(e[(e$Year==2017 & e$CT_period=="Winter"), "resp_pred"], na.rm=T)+
                sum(e[(e$Year==2018 & e$CT_period=="Winter"), "resp_pred"], na.rm=T))/44*1e3
    W20178se<-(sum(e[(e$Year==2017 & e$CT_period=="Winter"), "resp_pred.se"], na.rm=T)+
                 sum(e[(e$Year==2018 & e$CT_period=="Winter"), "resp_pred.se"], na.rm=T))/44*1e3
    
    ###Summer 2018
    S2018d<-sum(e[(e$Year==2018 & e$CT_period=="Summer"), "resp_pred"], na.rm=T)/44*1e3
    S2018se<-sum(e[(e$Year==2018 & e$CT_period=="Summer"), "resp_pred.se"], na.rm=T)/44*1e3
    
    #Store the results
    Cmin_field<-rbind(Cmin_field, data.frame(Cmin = c(S2017d, W20178d, S2018d),
                                             Cmic.se = c(S2017se, W20178se, S2018se),
                                             Period = c("Summer", "Winter", "Summer"),
                                             Year = c(2017, 2018, 2018),
                                             Soil = rep(Mlist$Soil[i], times=3),
                                             Catchment = rep(Mlist$Catchment[i], times=3),
                                             Horizon = rep(Mlist$Horizon[i], times=3),
                                             Origin = rep(Mlist$Origin[i], times=3),
                                             ID = rep(Mlist$ID[i], times=3)))
  }
}


############################################IER data#############################################
IER<-read.csv("IER_field.csv")

All_field<-merge(Cmin_field, IER, by.x = c("Year", "Period", "ID"),
                 by.y = c("Year", "Period", "ID"))
All_field$outliers<-"NO"
All_field[(All_field$Origin=="Native" & All_field$Period=="Summer" &
             All_field$Horizon=="Litter" & All_field$Pmin<10), "outliers"]<-"YES"
All_field[(All_field$Origin=="Native" & All_field$Period=="Summer" &
             All_field$Horizon=="Organic soil" & All_field$Pmin<10), "outliers"]<-"YES"

ggplot(All_field[All_field$outliers=="NO", ], aes(Cmin/1000, Pmin))+geom_point(cex=6, aes(colour = Horizon))+
  geom_errorbarh(aes(xmin = Cmin/1000-Cmic.se/1000, xmax=Cmin/1000+Cmic.se/1000))+
  facet_grid(Soil~Origin+Period)+scale_y_log10()+scale_x_log10()+
  stat_smooth(method = lm, se = F)

Ccontrols<-All_field[(All_field$Soil=="Certovo" & All_field$Horizon=="Control"), ]
Ccontrols$Origin<-"Transplanted"
Ccontrols$Soil<-"Plesne"

Pcontrols<-All_field[(All_field$Soil=="Plesne" & All_field$Horizon=="Control"), ]
Pcontrols$Origin<-"Transplanted"
Pcontrols$Soil<-"Certovo"

All_field2<-rbind(All_field, Ccontrols, Pcontrols)
All_field2$Legend<-All_field2$Horizon

ggplot(All_field2, aes(Cmin/1000, Pmin))+geom_point(cex=6, aes(colour = Legend))+
  geom_errorbarh(aes(xmin = Cmin/1000-Cmic.se/1000, xmax=Cmin/1000+Cmic.se/1000))+
  facet_grid(Soil~Origin+Period, scales = "free")+scale_y_log10()+scale_x_log10()+
  stat_smooth(method = lm, se = F)+
  geom_point(data = All_field2[(All_field2$Origin=="Transplanted" & All_field2$Horizon=="Control"), ],
             aes(Cmin/1000, Pmin), alpha=0.5, cex=6, colour="grey")+
  theme_min+
  ylab(expression(paste(PO[4], " (", mu,"mol ", box^{-1}, ")")))+
  xlab(expression(paste(CO[2], " (mmol ", box^{-1}, ")")))+
  geom_text(data=lbplNSa, label = paste("R^2 == ", round(summary(plNS)$adj.r.squared, 2)), parse = T, cex=6)+
  geom_text(data=lbplNSb, label = "p < 0.001", cex=6)+
  geom_text(data=lbplTSa, label = paste("R^2 == ", round(summary(plTS)$adj.r.squared, 2)), parse = T, cex=6)+
  geom_text(data=lbplTSb, label = "p < 0.001", cex=6)+
  geom_text(data=lbplNWa, label = paste("R^2 == ", round(summary(plNW)$adj.r.squared, 2)), parse = T, cex=6)+
  geom_text(data=lbplNWb, label = "p < 0.001", cex=6)+
  geom_text(data=lbplTWa, label = paste("R^2 == ", round(summary(plTW)$adj.r.squared, 2)), parse = T, cex=6)+
  geom_text(data=lbplTWb, label = "p < 0.001", cex=6)+
  geom_text(data=lbctNSa, label = paste("R^2 == ", round(summary(ctNS)$adj.r.squared, 2)), parse = T, cex=6)+
  geom_text(data=lbctNSb, label = "p < 0.001", cex=6)+
  geom_text(data=lbctTSa, label = paste("R^2 == ", round(summary(ctTS)$adj.r.squared, 2)), parse = T, cex=6)+
  geom_text(data=lbctTSb, label = "p < 0.001", cex=6)+
  geom_text(data=lbctNWa, label = paste("R^2 == ", round(summary(ctNW)$adj.r.squared, 2)), parse = T, cex=6)+
  geom_text(data=lbctNWb, label = "p < 0.001", cex=6)+
  geom_text(data=lbctTWa, label = paste("R^2 == ", round(summary(ctTW)$adj.r.squared, 2)), parse = T, cex=6)+
  geom_text(data=lbctTWb, label = "p < 0.001", cex=6)

#Statistics and labels
##Plesne
plNS<-lm(log(Pmin)~log(Cmin), data=subset(All_field2, Soil=="Plesne" & Origin=="Native" & Period=="Summer"))
summary(plNS)
lbplNSa<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 2, Origin = factor("Native", levels = c("Native", "Transplanted")),
                    Soil = factor("Plesne", levels = c("Plesne", "Certovo")),
                    Period = factor("Summer", levels = c("Summer", "Winter")))
lbplNSb<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 1.2, Origin = factor("Native", levels = c("Native", "Transplanted")),
                    Soil = factor("Plesne", levels = c("Plesne", "Certovo")),
                    Period = factor("Summer", levels = c("Summer", "Winter")))

plNW<-lm(log(Pmin)~log(Cmin), data=subset(All_field2, Soil=="Plesne" & Origin=="Native" & Period=="Winter"))
summary(plNW)

lbplNWa<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 2, Origin = factor("Native", levels = c("Native", "Transplanted")),
                    Soil = factor("Plesne", levels = c("Plesne", "Certovo")),
                    Period = factor("Winter", levels = c("Summer", "Winter")))
lbplNWb<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 1.2, Origin = factor("Native", levels = c("Native", "Transplanted")),
                    Soil = factor("Plesne", levels = c("Plesne", "Certovo")),
                    Period = factor("Winter", levels = c("Summer", "Winter")))

plTS<-lm(log(Pmin)~log(Cmin), data=subset(All_field2, Soil=="Plesne" & Origin=="Transplanted" & Period=="Summer"))
summary(plTS)

lbplTSa<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 2, Origin = factor("Transplanted", levels = c("Native", "Transplanted")),
                    Soil = factor("Plesne", levels = c("Plesne", "Certovo")),
                    Period = factor("Summer", levels = c("Summer", "Winter")))
lbplTSb<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 1.2, Origin = factor("Transplanted", levels = c("Native", "Transplanted")),
                    Soil = factor("Plesne", levels = c("Plesne", "Certovo")),
                    Period = factor("Summer", levels = c("Summer", "Winter")))

plTW<-lm(log(Pmin)~log(Cmin), data=subset(All_field2, Soil=="Plesne" & Origin=="Transplanted" & Period=="Winter"))
summary(plTW)

lbplTWa<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 2, Origin = factor("Transplanted", levels = c("Native", "Transplanted")),
                    Soil = factor("Plesne", levels = c("Plesne", "Certovo")),
                    Period = factor("Winter", levels = c("Summer", "Winter")))
lbplTWb<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 1.2, Origin = factor("Transplanted", levels = c("Native", "Transplanted")),
                    Soil = factor("Plesne", levels = c("Plesne", "Certovo")),
                    Period = factor("Winter", levels = c("Summer", "Winter")))

##Certovo
ctNS<-lm(log(Pmin)~log(Cmin), data=subset(All_field2, Soil=="Certovo" & Origin=="Native" & Period=="Summer"))
summary(ctNS)

lbctNSa<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 2, Origin = factor("Native", levels = c("Native", "Transplanted")),
                    Soil = factor("Certovo", levels = c("Plesne", "Certovo")),
                    Period = factor("Summer", levels = c("Summer", "Winter")))
lbctNSb<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 1.2, Origin = factor("Native", levels = c("Native", "Transplanted")),
                    Soil = factor("Certovo", levels = c("Plesne", "Certovo")),
                    Period = factor("Summer", levels = c("Summer", "Winter")))

ctNW<-lm(log(Pmin)~log(Cmin), data=subset(All_field2, Soil=="Certovo" & Origin=="Native" & Period=="Winter"))
summary(ctNW)

lbctNWa<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 2, Origin = factor("Native", levels = c("Native", "Transplanted")),
                    Soil = factor("Certovo", levels = c("Plesne", "Certovo")),
                    Period = factor("Winter", levels = c("Summer", "Winter")))
lbctNWb<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 1.2, Origin = factor("Native", levels = c("Native", "Transplanted")),
                    Soil = factor("Certovo", levels = c("Plesne", "Certovo")),
                    Period = factor("Winter", levels = c("Summer", "Winter")))

ctTS<-lm(log(Pmin)~log(Cmin), data=subset(All_field2, Soil=="Certovo" & Origin=="Transplanted" & Period=="Summer"))
summary(ctTS)

lbctTSa<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 2, Origin = factor("Transplanted", levels = c("Native", "Transplanted")),
                    Soil = factor("Certovo", levels = c("Plesne", "Certovo")),
                    Period = factor("Summer", levels = c("Summer", "Winter")))
lbctTSb<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 1.2, Origin = factor("Transplanted", levels = c("Native", "Transplanted")),
                    Soil = factor("Certovo", levels = c("Plesne", "Certovo")),
                    Period = factor("Summer", levels = c("Summer", "Winter")))

ctTW<-lm(log(Pmin)~log(Cmin), data=subset(All_field2, Soil=="Certovo" & Origin=="Transplanted" & Period=="Winter"))
summary(ctTW)

lbctTWa<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 2, Origin = factor("Transplanted", levels = c("Native", "Transplanted")),
                    Soil = factor("Certovo", levels = c("Plesne", "Certovo")),
                    Period = factor("Winter", levels = c("Summer", "Winter")))
lbctTWb<-data.frame(lab="Text",
                    Cmin = 100*1000, Pmin = 1.2, Origin = factor("Transplanted", levels = c("Native", "Transplanted")),
                    Soil = factor("Certovo", levels = c("Plesne", "Certovo")),
                    Period = factor("Winter", levels = c("Summer", "Winter")))

write.xlsx(Cmin_field, file = c("./Field_data_raw/Kumulativni.xlsx"))
resp_all$resp<-resp_all$resp/44*1000
resp_all2<-resp_all
resp_all2[resp_all2$Origin=="Control", "horizon"]<-"Control"
resp_all2[resp_all2$Origin=="Control", "Origin"]<-"Native"

write.xlsx(resp_all2, file = c("./Field_data_raw/Rychlosti.xlsx"))
