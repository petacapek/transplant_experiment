#LIBRARIES
library(dplyr)
library(ggplot2)
library(FME)
library(DEoptim)
library(segmented)
library(lme4)

#Reading data
#Respirations
#resp<-read.xlsx(xlsxFile = c("resp_teren_181031.xlsx"))
resp<-read.csv(file = "resp_field_181031.csv")[, -8]
resp$id<-as.character(resp$id)
resp$date_time<-as.character(resp$date_time)
resp$time.tms<-as.character(resp$time.tms)
resp$date<-as.character(resp$date)

resp[resp$id=="PL_C1", "id"]<-c("PLC_1")
resp[resp$id=="PL_C2", "id"]<-c("PLC_2")
resp[resp$id=="PL_C3", "id"]<-c("PLC_3")
resp[resp$id=="PL_C4", "id"]<-c("PLC_4")
resp[resp$id=="PL_C5", "id"]<-c("PLC_5")

#Environmental data
#env<-read.xlsx(xlsxFile = c("tms_raw_181031.xlsx"), sheet = 2)
env<-read.csv("tms_raw_181031.csv")
summary(env)

#Moisture calibration data
#wcal<-read.xlsx(xlsxFile = c("kalibrace_vlhkomery.xlsx"), sheet = 2)
wcal<-read.csv("tms_calibration.csv")[-41, ]
wcal$Soil<-as.character(wcal$Soil)
summary(wcal)

ggplot(wcal, aes(raw, W))+geom_point(cex=6, aes(colour=Soil))+xlab("Raw tms data")+
  ylab("Relative moisture (%)")

#PL
ggplot(subset(wcal[-37, ], Soil=="PL1" | Soil=="PL2"), aes(raw, W))+geom_point(cex=6, aes(colour=Soil))+xlab("Raw tms data")+
  ylab("Relative moisture (%)")+stat_function(fun = function(x){29.70*x/(24.90+x)}-28.99)

library(minpack.lm)
PLw.cal<-coef(nlsLM(W~a*raw/(b+raw)-c, subset(wcal[-37, ], Soil=="PL1" | Soil=="PL2"),
            start=list(a=1.8, b=1200, c=0.8), control = c(maxiter=1e6)))

#CT
ggplot(subset(wcal, Soil=="CT2"), aes(raw, W))+geom_point(cex=6, aes(colour=Soil))+xlab("Raw tms data")+
  ylab("Relative moisture (%)")+stat_function(fun = function(x){20.60*x/(29.95+x)}-19.78)+xlim(720,3500)
  
library(minpack.lm)
CTw.cal<-coef(nlsLM(W~a*raw/(b+raw)-c, subset(wcal, Soil=="CT2"),
              start=list(a=1.8, b=1200, c=0.8), control = c(maxiter=1e6)))

#Merging
##O horizon
###PL
####average of tms data for each time point and 
PLOenv<-data.frame(time=rep(as.character(env[,c("time")], times=5)),
                   Tsurface=c(as.numeric(env[,c("pl1_t2")]),
                          as.numeric(env[,c("pl2_t2")]),
                          as.numeric(env[,c("pl3_t2")]),
                          as.numeric(env[,c("pl4_t2")]),
                          as.numeric(env[,c("pl5_t2")])),
                   Tsoil=c(as.numeric(env[,c("pl1_t1")]),
                              as.numeric(env[,c("pl2_t1")]),
                              as.numeric(env[,c("pl3_t1")]),
                              as.numeric(env[,c("pl4_t1")]),
                              as.numeric(env[,c("pl5_t1")])),
                   Tair=c(as.numeric(env[,c("pl1_t3")]),
                           as.numeric(env[,c("pl2_t3")]),
                           as.numeric(env[,c("pl3_t3")]),
                           as.numeric(env[,c("pl4_t3")]),
                           as.numeric(env[,c("pl5_t3")])),
                   Moisture=c(as.numeric(env[,c("pl1_m")]),
                              as.numeric(env[,c("pl2_m")]),
                              as.numeric(env[,c("pl3_m")]),
                              as.numeric(env[,c("pl4_m")]),
                              as.numeric(env[,c("pl5_m")])),
                   Block=c(rep(1, times=nrow(env)),
                           rep(2, times=nrow(env)),
                           rep(3, times=nrow(env)),
                           rep(4, times=nrow(env)),
                           rep(5, times=nrow(env))))

ggplot(PLOenv[c(1:500), ], aes(as.numeric(time), Tsurface))+geom_line()+
  geom_line(data=PLOenv[c(1:500), ], aes(as.numeric(time), Tsoil), colour="red")+
  geom_line(data=PLOenv[c(1:500), ], aes(as.numeric(time), Tair), colour="blue")

#calculating the difference between surface and soil temperature
PLOenv$Tresistance<-with(PLOenv, abs(Tsurface-Tsoil))

#Moisture (relative amount in water in the whole soil sample)
PLOenv$W<-with(PLOenv, PLw.cal[1]*Moisture/(PLw.cal[2]+Moisture)-PLw.cal[3])
#replace negative values with the observed minimum moisture
PLOenv[(PLOenv$W<0 & !is.na(PLOenv$W)), "W"]<-min(PLOenv[(PLOenv$W>0 & !is.na(PLOenv$W)), "W"])

summary(PLOenv)

#extract O horizon data from transplants on PL
PLOs<-resp[(resp$id=="PLO_1" |
              resp$id=="PLO_2" |
              resp$id=="PLO_3" |
              resp$id=="PLO_4" |
              resp$id=="PLO_5" |
              resp$id=="CTO_1" |
              resp$id=="CTO_2" |
              resp$id=="CTO_3" |
              resp$id=="CTO_4" |
              resp$id=="CTO_5"), c("id", "date_time", "resp")]

#add the block assigment
PLOs$Block<-substr(PLOs$id, 5, 5)

#merge
PLO<-merge(PLOenv, PLOs, by.x=c("time", "Block"), by.y=c("date_time", "Block"), all.y = T)
#add catchment-specific label
PLO$Soil<-ifelse((PLO$id=="PLO_1" |
                    PLO$id=="PLO_2" |
                    PLO$id=="PLO_3" |
                    PLO$id=="PLO_4" |
                    PLO$id=="PLO_5"), "Plesne", "Certovo")
#add transplant related label
PLO$Origin<-ifelse(PLO$Soil=="Plesne", "Native", "Transplanted")
#add horizon information
PLO$horizon<-c("Litter")

###CT
####average of tms data for each time point and 
CTOenv<-data.frame(time=rep(as.character(env[,c("time")], times=5)),
                   Tsurface=c(as.numeric(env[,c("ct1_t2")]),
                              as.numeric(env[,c("ct2_t2")]),
                              as.numeric(env[,c("ct3_t2")]),
                              as.numeric(env[,c("ct4_t2")]),
                              as.numeric(env[,c("ct5_t2")])),
                   Tsoil=c(as.numeric(env[,c("ct1_t1")]),
                           as.numeric(env[,c("ct2_t1")]),
                           as.numeric(env[,c("ct3_t1")]),
                           as.numeric(env[,c("ct4_t1")]),
                           as.numeric(env[,c("ct5_t1")])),
                   Tair=c(as.numeric(env[,c("ct1_t3")]),
                          as.numeric(env[,c("ct2_t3")]),
                          as.numeric(env[,c("ct3_t3")]),
                          as.numeric(env[,c("ct4_t3")]),
                          as.numeric(env[,c("ct5_t3")])),
                   Moisture=c(as.numeric(env[,c("ct1_m")]),
                              as.numeric(env[,c("ct2_m")]),
                              as.numeric(env[,c("ct3_m")]),
                              as.numeric(env[,c("ct4_m")]),
                              as.numeric(env[,c("ct5_m")])),
                   Block=c(rep(6, times=nrow(env)),
                           rep(7, times=nrow(env)),
                           rep(8, times=nrow(env)),
                           rep(9, times=nrow(env)),
                           rep(10, times=nrow(env))))

ggplot(CTOenv[c(1:500), ], aes(as.numeric(time), Tsurface))+geom_line()+
  geom_line(data=CTOenv[c(1:500), ], aes(as.numeric(time), Tsoil), colour="red")+
  geom_line(data=CTOenv[c(1:500), ], aes(as.numeric(time), Tair), colour="blue")

#calculating the difference between surface and soil temperature
CTOenv$Tresistance<-with(CTOenv, abs(Tsurface-Tsoil))

#Moisture (relative amount in water in the whole soil sample)
CTOenv$W<-with(CTOenv, CTw.cal[1]*Moisture/(CTw.cal[2]+Moisture)-CTw.cal[3])
#replace negative values with observed minimum
CTOenv[(CTOenv$W<0 & !is.na(CTOenv$W)), "W"]<-min(CTOenv[(CTOenv$W>0 & !is.na(CTOenv$W)), "W"])

summary(CTOenv)

#extract O horizon data from transplants on CT
CTOs<-resp[(resp$id=="PLO_6" |
              resp$id=="PLO_7" |
              resp$id=="PLO_8" |
              resp$id=="PLO_9" |
              resp$id=="PLO_10" |
              resp$id=="CTO_6" |
              resp$id=="CTO_7" |
              resp$id=="CTO_8" |
              resp$id=="CTO_9" |
              resp$id=="CTO_10"), c("id", "date_time", "resp")]

#add the block assigment
CTOs$Block<-substr(CTOs$id, 5, 6)

#merge
CTO<-merge(CTOenv, CTOs, by.x=c("time", "Block"), by.y=c("date_time", "Block"), all.y = T)

#add catchment-specific label
CTO$Soil<-ifelse((CTO$id=="PLO_6" |
                    CTO$id=="PLO_7" |
                    CTO$id=="PLO_8" |
                    CTO$id=="PLO_9" |
                    CTO$id=="PLO_10"), "Plesne", "Certovo")
#add transplant related label
CTO$Origin<-ifelse(CTO$Soil=="Certovo", "Native", "Transplanted")
#add horizon information
CTO$horizon<-c("Litter")

##A horizon
###PL
#extract A horizon data from transplants on PL
PLAs<-resp[(resp$id=="PLA_1" |
              resp$id=="PLA_2" |
              resp$id=="PLA_3" |
              resp$id=="PLA_4" |
              resp$id=="PLA_5" |
              resp$id=="CTA_1" |
              resp$id=="CTA_2" |
              resp$id=="CTA_3" |
              resp$id=="CTA_4" |
              resp$id=="CTA_5"), c("id", "date_time", "resp")]
#add the block assigment
PLAs$Block<-substr(PLAs$id, 5, 5)

#merge
PLA<-merge(PLOenv, PLAs, by.x=c("time", "Block"), by.y=c("date_time", "Block"), all.y = T)
#add catchment-specific label
PLA$Soil<-ifelse((PLA$id=="PLA_1" |
                    PLA$id=="PLA_2" |
                    PLA$id=="PLA_3" |
                    PLA$id=="PLA_4" |
                    PLA$id=="PLA_5"), "Plesne", "Certovo")
#add transplant related label
PLA$Origin<-ifelse(PLA$Soil=="Plesne", "Native", "Transplanted")
#add horizon information
PLA$horizon<-c("Organic soil")

###CT
#extract A horizon data from transplants on CT
CTAs<-resp[(resp$id=="PLA_6" |
              resp$id=="PLA_7" |
              resp$id=="PLA_8" |
              resp$id=="PLA_9" |
              resp$id=="PLA_10" |
              resp$id=="CTA_6" |
              resp$id=="CTA_7" |
              resp$id=="CTA_8" |
              resp$id=="CTA_9" |
              resp$id=="CTA_10"), c("id", "date_time", "resp")]

#add the block assigment
CTAs$Block<-substr(CTAs$id, 5, 6)

#merge
CTA<-merge(CTOenv, CTAs, by.x=c("time", "Block"), by.y=c("date_time", "Block"), all.y = T)
#add catchment-specific label
CTA$Soil<-ifelse((CTA$id=="PLA_6" |
                    CTA$id=="PLA_7" |
                    CTA$id=="PLA_8" |
                    CTA$id=="PLA_9" |
                    CTA$id=="PLA_10"), "Plesne", "Certovo")
#add transplant related label
CTA$Origin<-ifelse(CTA$Soil=="Certovo", "Native", "Transplanted")
#add horizon information
CTA$horizon<-c("Organic soil")

###Controls - Plesne
PLCs<-resp[(resp$id=="PLC_1" |
              resp$id=="PLC_2" |
              resp$id=="PLC_3" |
              resp$id=="PLC_4" |
              resp$id=="PLC_5"), c("id", "date_time", "resp")]

#add the block assigment
PLCs$Block<-substr(PLCs$id, 5, 5)

#merge
PLC<-merge(PLOenv, PLCs, by.x=c("time", "Block"), by.y=c("date_time", "Block"))
#add catchment-specific label
PLC$Soil<-c("Plesne")
#add transplant related label
PLC$Origin<-c("Control")
#add horizon information
PLC$horizon<-c("Organic soil")

###Controls - Certovo
CTCs<-resp[(resp$id=="CTC_1" |
              resp$id=="CTC_2" |
              resp$id=="CTC_3" |
              resp$id=="CTC_4" |
              resp$id=="CTC_5"), c("id", "date_time", "resp")]

#add the block assigment
CTCs$Block<-vector("character", length = nrow(CTCs))

for(i in 1:nrow(CTCs)){
  
  if(CTCs[i, "id"]=="CTC_1"){
    CTCs[i, "Block"]<-c(6)
  }else{
    if(CTCs[i, "id"]=="CTC_2"){
      CTCs[i, "Block"]<-c(7)
    }else{
      if(CTCs[i, "id"]=="CTC_3"){
        CTCs[i, "Block"]<-c(8)
      }else{
        if(CTCs[i, "id"]=="CTC_4"){
          CTCs[i, "Block"]<-c(9)
        }else{
          CTCs[i, "Block"]<-c(10)
        }
      }
    }
  }
}

#merge
CTC<-merge(CTOenv, CTCs, by.x=c("time", "Block"), by.y=c("date_time", "Block"))
#add catchment-specific label
CTC$Soil<-c("Certovo")
#add transplant related label
CTC$Origin<-c("Control")
#add horizon information
CTC$horizon<-c("Organic soil")


#####Bind all together
resp_all<-rbind(PLO, CTO, PLA, CTA, PLC, CTC)

#rename the horizon in controls
resp_all[resp_all$Origin=="Control", "horizon"]<-c("Litter")

#Correct measured respiration rate by factual volume/area 
resp_all$resp_corr<-numeric(length = nrow(resp_all))

for(i in 1:nrow(resp_all)){
  if(resp_all$Soil[i]=="Plesne" & resp_all$horizon[i]=="Litter"){
    resp_all$resp_corr[i]<-67/78*resp_all$resp[i]
  }else{
    if(resp_all$Soil[i]=="Plesne" & resp_all$horizon[i]=="Organic soil"){
      resp_all$resp_corr[i]<-111/78*resp_all$resp[i]
    }else{
      if(resp_all$Soil[i]=="Certovo" & resp_all$horizon[i]=="Litter"){
        resp_all$resp_corr[i]<-67/78*resp_all$resp[i]
      }else{
        resp_all$resp_corr[i]<-111/78*resp_all$resp[i]
      }
    }
  }
}

resp_all[resp_all$Origin=="Control", "resp_corr"]<-resp_all[resp_all$Origin=="Control", "resp"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ANALYSIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ggplot(resp_all, aes(Tair, resp_corr))+geom_point(cex=6, aes(colour=Origin), alpha=0.5)+
  facet_wrap(horizon~Block, scales="free")

ggplot(resp_all[(resp_all$Soil=="Certovo" & resp_all$Origin=="Native"), ], aes(Tair, resp_corr))+
  geom_point(cex=6, aes(colour=horizon), alpha=0.5)+
  facet_wrap(~Block, scales="free")

#mark the outliers
resp_all$outliers<-character(length = nrow(resp_all))

resp_all[(resp_all$Block==8 & resp_all$horizon=="Litter" &
            resp_all$Origin=="Native" & resp_all$resp>1), ]

resp_all$outliers<-ifelse(substr(resp_all[,"time"], 1,10)=="2018.10.30", "YES", "NO")

ggplot(resp_all[resp_all$outliers=="NO", ], aes(Tsurface, resp_corr))+geom_point(cex=6, aes(colour=Origin), alpha=0.5)+
  facet_wrap(horizon~Block, scales="free")

#Averages
Env.av<-resp_all %>% group_by(Soil, Origin) %>% summarise(Tsoil=mean(Tsoil, na.rm = T),
                                                  Tsurface=mean(Tsurface, na.rm = T),
                                                  Tair=mean(Tair, na.rm = T),
                                                  W=mean(W, na.rm = T))


#Heterotrophic respiration rate of transplant samples is recalculated per 
#volume of soil column in the transplant and mass of the soil in the transplant.
#Add mass of the soil
resp_all$Soil_mass<-numeric(length = nrow(resp_all))#in grams

for(i in 1:nrow(resp_all)){
  if(resp_all$Soil[i]=="Plesne" & resp_all$horizon[i]=="Litter"){
    resp_all$Soil_mass[i]<-100*0.280
  }else{
    if(resp_all$Soil[i]=="Plesne" & resp_all$horizon[i]=="Organic soil"){
      resp_all$Soil_mass[i]<-100*0.280+500*0.423
    }else{
      if(resp_all$Soil[i]=="Certovo" & resp_all$horizon[i]=="Litter"){
        resp_all$Soil_mass[i]<-100*0.302
      }else{
        resp_all$Soil_mass[i]<-100*0.302+500*0.345
      }
    }
  }
}

#controls have mass of soil equal to volume of litter horizon times its bulk density
resp_all[resp_all$Origin=="Control", "Soil_mass"]<-ifelse(resp_all[resp_all$Origin=="Control", "Soil"]=="Plesne",
                                                          4.5*78*0.108, 5.5*78*0.091)

#Add the surface area of the tranplanted soil column
resp_all$Surface<-ifelse(resp_all$horizon=="Litter", 67, 111)
resp_all[resp_all$Origin=="Control", "Surface"]<-78

#Calculating respiration rate per soil mass 
resp_all$resp_mass<-with(resp_all, Surface/1e4*resp_corr/Soil_mass*60*60)#in umol/g/h

#Add volume of soil column in the transplant in cm3
resp_all$Volume<-numeric(length = nrow(resp_all))

for(i in 1:nrow(resp_all)){
  if(resp_all$Soil[i]=="Plesne" & resp_all$horizon[i]=="Litter"){
    resp_all$Volume[i]<-4.5*67
  }else{
    if(resp_all$Soil[i]=="Plesne" & resp_all$horizon[i]=="Organic soil"){
      resp_all$Volume[i]<-4.5*67+111*10.5
    }else{
      if(resp_all$Soil[i]=="Certovo" & resp_all$horizon[i]=="Litter"){
        resp_all$Volume[i]<-5.5*67
      }else{
        resp_all$Volume[i]<-5.5*67+10.5*111
      }
    }
  }
}

#Assuming only the presence of litter in controls - colar depth around 4 cm (4.5 cm Plesne, 5.5 cm Certovo)
resp_all[resp_all$Origin=="Control", "Volume"]<-ifelse(resp_all[resp_all$Origin=="Control", "Soil"]=="Plesne",
                                                       4.5*78, 5.5*78)

#Calculating respiration rate per volume of the soil
resp_all$resp_vol<-with(resp_all, Surface/1e4*resp_corr/Volume*60*60*1000)#in nmol/cm3/h

#Show averages
resp_all[resp_all$outliers=="NO", ] %>% group_by(Soil, Origin, horizon) %>% summarise(HR=mean(resp_corr, na.rm = T),
                                                           HR_mass=mean(resp_mass, na.rm = T),
                                                           HR_vol=mean(resp_vol, na.rm = T))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################Statistics###############################################
#correspondence function
correspondence<-function(obs, pred, N){
  SSres=sum(((obs-pred)^2), na.rm = T) 
  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T)
  ll=-sum(((obs-pred)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2)-length(obs[!is.na(obs)])*log(sd(obs, na.rm = T)^2*mean(obs, na.rm = T))/2
  R2<-1-SSres/SStot
  AIC<-2*N-2*ll
  
  out<-c(R2=R2, AIC=AIC, ll=ll)
  return(out)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Soil warming is associated with soil drying:

ggplot(resp_all, aes(Tair, W))+geom_point(aes(colour=horizon))+facet_grid(Soil~Origin)+
  stat_smooth(method = "lm", se=F)+
  xlab("Temperature (°C)")+
  ylab("Gravimetric water content (% of fresh weight)")+
  theme_bw()+
  theme(axis.line.x = element_blank(),
        axis.line.y =element_blank(),
        plot.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size=14,colour="black"),
        axis.text.x = element_text(size=14,colour="black"),
        axis.title = element_text(size=14,colour="black"),
        legend.text=element_text(size=14,colour="black"),
        legend.title=element_blank(),
        legend.justification=c(1,0), legend.position=c("top"),
        legend.key.size = unit(1,'lines'),
        legend.background=element_rect(fill=NA, colour=NA),
        strip.text = element_text(size=14),
        strip.background = element_rect(colour="black", fill = "white"),
        panel.spacing = unit(0.25, "lines"))

#The effect of temperature on heterotrophic soil respiration can be confounded
#by increased soil drying 
ggplot(resp_all[resp_all$outliers=="NO", ], aes(W, resp_corr))+geom_point(aes(colour=horizon))+
  facet_wrap(Soil~Origin, scales = "free")+
  stat_smooth(method = "lm", se=F)+
  ylab(expression(paste("Respiration rate (", mu, "mol ", m^{-2}~h^{-1}, ")")))+
  xlab("Gravimetric water content (% of fresh weight)")+
  theme_bw()+
  theme(axis.line.x = element_blank(),
        axis.line.y =element_blank(),
        plot.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size=14,colour="black"),
        axis.text.x = element_text(size=14,colour="black"),
        axis.title = element_text(size=14,colour="black"),
        legend.text=element_text(size=14,colour="black"),
        legend.title=element_blank(),
        legend.justification=c(1,0), legend.position=c("top"),
        legend.key.size = unit(1,'lines'),
        legend.background=element_rect(fill=NA, colour=NA),
        strip.text = element_text(size=14),
        strip.background = element_rect(colour="black", fill = "white"),
        panel.spacing = unit(0.25, "lines"))
#The figure suggests unimodal relationship as has been suggested previously

#Thus, temperature sensitivity is estimated by two different non-linear functions:
#Arrhenius equation with and without accounting for soil moisture effect.
#The equations are fitted for each experimental treatment separately and their goodness 
#of fit is compared.

#First, gravimetric water content is recalculated to volumetric water content cm3/cm3 (theta)
resp_all$theta<-with(resp_all, W/(1-W)*Soil_mass/Volume)

#Then, soil porosity of soil column (phi - cm3/cm3) is added 
resp_all$phi<-numeric(length = nrow(resp_all))

for(i in 1:nrow(resp_all)){
  if(resp_all$Soil[i]=="Plesne" & resp_all$horizon[i]=="Litter"){
    resp_all$phi[i]<-0.108
  }else{
    if(resp_all$Soil[i]=="Plesne" & resp_all$horizon[i]=="Organic soil"){
      resp_all$phi[i]<-0.246
    }else{
      if(resp_all$Soil[i]=="Certovo" & resp_all$horizon[i]=="Litter"){
        resp_all$phi[i]<-0.102
      }else{
        resp_all$phi[i]<-0.194
      }
    }
  }
}

#The volumetric water content higher than total porosity (3.3% of observations) 
#is adjusted to maximum, i.e. total porosity
for(i in 1:nrow(resp_all)){
  if(resp_all$theta[i]/resp_all$phi[i]>1 & !is.na(resp_all$theta[i])){
    resp_all$theta[i]<-resp_all$phi[i]
  }else{
    resp_all$theta[i]<-resp_all$theta[i]
  }
}

ggplot(resp_all[(resp_all$outliers=="NO"), ], 
       aes(theta/phi, resp_corr))+geom_point(aes(colour=horizon))+
  facet_wrap(Soil~Origin+Block, scales = "free")+
  stat_smooth(method = "lm", se=F, formula = y~poly(x,2))+
  ylab(expression(paste("Respiration rate (", mu, "mol ", m^{-2}~h^{-1}, ")")))+
  xlab(expression(paste("Volumetric water content (", cm^{3}~cm^{-3} ,")")))+
  theme_bw()+
  theme(axis.line.x = element_blank(),
        axis.line.y =element_blank(),
        plot.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size=14,colour="black"),
        axis.text.x = element_text(size=14,colour="black"),
        axis.title = element_text(size=14,colour="black"),
        legend.text=element_text(size=14,colour="black"),
        legend.title=element_blank(),
        legend.justification=c(1,0), legend.position=c("top"),
        legend.key.size = unit(1,'lines'),
        legend.background=element_rect(fill=NA, colour=NA),
        strip.text = element_text(size=14),
        strip.background = element_rect(colour="black", fill = "white"),
        panel.spacing = unit(0.25, "lines"))

#Data_frame is exported
write.csv(resp_all, "merged_data.csv")


#In the following function, three different linear models are defined and fitted to the data
#First model do not assume effect of water content, the second two do. 
#First model assumes linear relationship between respiration rate and relative water content
#The second model assumes polynomial linear relationship .
#fit of the models are compared for each soil core separately


#Read the function
source("./R_functions/TW.R")

TW_out<-TW(resp_all)

#Results
summary(TW_out$R2$TempWlin-TW_out$R2$Temp>0)#36/50=72%
summary(TW_out$R2$TempWpoly-TW_out$R2$Temp>0)#41/50=82%
summary(TW_out$R2$TempWpoly-TW_out$R2$TempWlin>0)#34/50=68%

summary(TW_out$AIC$TempWlin-TW_out$AIC$Temp<0)#31/50=62%
summary(TW_out$AIC$TempWpoly-TW_out$AIC$Temp<0)#33/50=60%
summary(TW_out$AIC$TempWpoly-TW_out$AIC$TempWlin<0)#27/50=54%

for(i in 1:50){
  print(TW_out$anova_comp[[i]])
}
for(i in 1:50){
  print(TW_out$AIC_comp[[i]])
}

#Now, the intercept and slope of increase/decrease of soil moisture 24 hours before the measurement 
#is added to the data frame as an explanatory variable
#slope and intercept have units of moisture per g dry weight and will be recalculated later 
#24 hours is selected based on the model excercise of Evans et al. (2016) SBB
#At the same time, mean temperature of the soil surface during preceding month
#is added to the data frame

#Three new columns are dedined - moisture inetrcept (Win) and slope (Wsl), and median surface 
#temperature (Tm)
resp_all$Win<-numeric(length = nrow(resp_all))
resp_all$Wsl<-numeric(length = nrow(resp_all))
resp_all$Tm<-numeric(length = nrow(resp_all))

#1. Plesne
##Litter (PLO1 - 5, CTO1 - 5)
##Organic soil (PLA1 - 5, CTA1 - 5)
##Controls (PLC1 - 5)

#2. Certovo
##Litter (PLO6 - 10, CTO6 - 10)
##Organic soil (PLA6 - 10, CTA6 - 10)
##Controls (PLC6 - 10)

for(i in 1:nrow(resp_all)){
  tryCatch({
    
  if(as.factor(resp_all$Block[i])=="1" |
     as.factor(resp_all$Block[i])=="2" |
     as.factor(resp_all$Block[i])=="3" |
     as.factor(resp_all$Block[i])=="4" |
     as.factor(resp_all$Block[i])=="5"){
    
    #First, find the date and block in PLOenv dataframe
    p<-which(PLOenv$Block==resp_all[i, "Block"] &
               as.character(PLOenv$time)==as.character(resp_all[i, "time"]))
    #Extract moisture data from previous day
    if(p-96>0){
      Wrunning<-as.numeric(PLOenv[c((p-96):p), "W"])
    }else{
      Wrunning<-as.numeric(PLOenv[c(1:p), "W"])
    }
    
    #Run linear regression and add the results to Win and Wsl
    resp_all[i, "Win"]<-coef(lm(Wrunning~seq(0, 24, by = 0.25)))[1]
    resp_all[i, "Wsl"]<-coef(lm(Wrunning~seq(0, 24, by = 0.25)))[2]
    #Extract Tsurface data from previous month and add the median to Tm
    if(p-2880>0){
      resp_all[i, "Tm"]<-median(as.numeric(PLOenv[c((p-2880):p), "Tsurface"]), na.rm=T)
    }else{
      resp_all[i, "Tm"]<-median(as.numeric(PLOenv[c(1:p), "Tsurface"]), na.rm=T)
    }
    
    print(i)
  }else{
    p<-which(CTOenv$Block==resp_all[i, "Block"] &
               as.character(CTOenv$time)==as.character(resp_all[i, "time"]))
    #Extract moisture data from previous day
    if(p-96>0){
      Wrunning<-as.numeric(PLOenv[c((p-96):p), "W"])
    }else{
      Wrunning<-as.numeric(PLOenv[c(1:p), "W"])
    }
    
    #Run linear regression and add the results to Win and Wsl
    resp_all[i, "Win"]<-coef(lm(Wrunning~seq(0, 24, by = 0.25)))[1]
    resp_all[i, "Wsl"]<-coef(lm(Wrunning~seq(0, 24, by = 0.25)))[2]
    #Extract Tsurface data from previous month and add the median to Tm
    if(p-2880>0){
      resp_all[i, "Tm"]<-median(as.numeric(CTOenv[c((p-2880):p), "Tsurface"]), na.rm=T)
    }else{
      resp_all[i, "Tm"]<-median(as.numeric(CTOenv[c(1:p), "Tsurface"]), na.rm=T)
    }
    print(i)
  }
    }, error=function(e){cat("Something is wrong with", print(i))})
}


#Rows 324 - 328 finished with error - will end with NA
resp_all[c(324:328), "Win"]<-NA
resp_all[c(324:328), "Wsl"]<-NA
resp_all[c(324:328), "Tm"]<-NA

summary(resp_all)

#The gravimetric water content intercept is recalculated to volumetric water content cm3/cm3 (theta)
resp_all$theta_in<-with(resp_all, Win/(1-Win)*Soil_mass/Volume)
#and normalized to pore space
resp_all$theta_rel_in<-with(resp_all, theta_in/phi)

#The slope of change of gravimetric water content is recalculated to change in relative volumetric
#water content change (cm3/cm3 of total pore space) in % per h
resp_all$theta_rel_sl<-with(resp_all, Wsl/(1-Wsl)*Soil_mass/Volume/phi*100)

#Define the position of transplants/controls within catchments
resp_all$Catchment<-ifelse(resp_all$Block<=5, "Plesne", "Certovo")
#Mark controls as they are bit different
resp_all$Control<-ifelse(resp_all$Origin=="Control", "TRUE", "FALSE")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#basic model
#Using soil surface temperature as an explanatory variable
resp_all$Temp<-1/(resp_all$Tsurface+273.15)/8.314

m0<-lmer(log(resp_corr)~Temp+poly(theta_rel, 2, raw = TRUE)+(1|Block), 
         resp_all[resp_all$outliers=="NO", ])

summary(m0)
correspondence(obs = resp_all[(!is.na(resp_all$theta_rel) & resp_all$outliers=="NO"), "resp_corr"], 
               pred = exp(predict(m0, new.data = resp_all[(!is.na(resp_all$theta_rel) & resp_all$outliers=="NO"),])),
               N = 6)

##Testing hypotheses
###1. Are controls higher?
m1a<-update(m0, .~.+Control)
summary(m1a)

correspondence(obs = resp_all[(!is.na(resp_all$theta_rel) & resp_all$outliers=="NO"), "resp_corr"], 
               pred = exp(predict(m1a, new.data = resp_all[(!is.na(resp_all$theta_rel) & resp_all$outliers=="NO"),])),
               N = 7)
anova(m0, m1a)

###2. Is there any difference between the measurements taken at different catchments?
m1b<-update(m0, .~.+Catchment)
summary(m1b)
anova(m1b)

###3. Is there any difference between the horizons?
m1c<-update(m0, .~.+horizon)
summary(m1c)
anova(m1c)


###################################Controls separately######################################
mc0<-lmer(log(resp_corr)~Temp+poly(theta_rel, 2, raw = TRUE)+(1|Block), 
          resp_all[(resp_all$outliers=="NO" & !is.na(resp_all$theta_rel)), ], 
          subset = Control=="TRUE")

summary(mc0)

#It seems that the moisture doesn't play a role - test its absence and linear term
mc0a<-lmer(log(resp_corr)~Temp+theta_rel+(1|Block), 
           resp_all[(resp_all$outliers=="NO" & !is.na(resp_all$theta_rel)), ], 
           subset = Control=="TRUE")
summary(mc0a)

mc0b<-lmer(log(resp_corr)~Temp+(1|Block), 
           resp_all[(resp_all$outliers=="NO" & !is.na(resp_all$theta_rel)), ], 
           subset = Control=="TRUE")
summary(mc0b)

anova(mc0, mc0b)
#Test indicates no effect of moisture on respiration rate of controls
#Refit the null model
mc0<-lmer(log(resp_corr)~Temp+(1|Block), 
          resp_all[(resp_all$outliers=="NO"), ], 
          subset = Control=="TRUE")

correspondence(obs = resp_all[(resp_all$Control=="TRUE" & resp_all$outliers=="NO"), "resp_corr"], 
               pred = exp(predict(mc0, new.data = resp_all[(resp_all$Control=="TRUE" & resp_all$outliers=="NO"), "resp_corr"])),
               N = 2+2)

#Is there a difference between the catchments?
mc1<-update(mc0, .~.+Catchment)
summary(mc1)
anova(mc1)
anova(mc1, mc0)

correspondence(obs = resp_all[(resp_all$Control=="TRUE" & resp_all$outliers=="NO"), "resp_corr"], 
               pred = exp(predict(mc1, new.data = resp_all[(resp_all$Control=="TRUE" & resp_all$outliers=="NO"), "resp_corr"])),
               N = 3+2)

#Is there a difference in the temperature sensitivity between the catchments?
mc2<-update(mc0, .~.:Catchment)
summary(mc2)
anova(mc2)
anova(mc0, mc2)

correspondence(obs = resp_all[(resp_all$Control=="TRUE" & resp_all$outliers=="NO"), "resp_corr"], 
               pred = exp(predict(mc2, new.data = resp_all[(resp_all$Control=="TRUE" & resp_all$outliers=="NO"), "resp_corr"])),
               N = 3+2)

#Is there an effect of changing soil moisture change?
mc3a<-update(mc0, .~.+theta_rel_sl)
summary(mc3a)
anova(mc3a)
anova(mc3a, mc0)

mc3b<-update(mc0, .~.+theta_rel_in)
summary(mc3b)
anova(mc3b)
anova(mc3b, mc0)

#Does the changing soil moisture change temperature sensitivity?
mc4a<-update(mc0, .~.:theta_rel_sl)
summary(mc4a)
anova(mc4a)
anova(mc4a, mc0)

mc4b<-update(mc0, .~.:theta_rel_in)
summary(mc4b)
anova(mc4b)
anova(mc4b, mc0)

#Does the previous temperature regime change temperature sensitivity?
mc5<-update(mc0, .~.:Tm)
summary(mc5)
anova(mc5)
anova(mc5, mc0)
