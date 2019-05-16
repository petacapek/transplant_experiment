TW<-function(d){
  #reading segmented and bbmle library
  library(segmented)
  library(bbmle)
  #defining Arrhenius plot temperature
  d$Temp_surface<-1/(d$Tsurface+273.15)/8.314
  
  #creating list that will store anova comparisons
  anova_list<-vector("list", length = length(unique(d$id)))
  names(anova_list)<-unique(d$id)
  #creating list that will store AIC comparisons
  AIC_list<-vector("list", length = length(unique(d$id)))
  names(AIC_list)<-unique(d$id)
  #creating data frame tha will store R2 values of models
  R2_data<-data.frame(id=unique(d$id), Temp=numeric(length = length(unique(d$id))),
                      TempWlin=numeric(length = length(unique(d$id))),
                      TempWpoly=numeric(length = length(unique(d$id)))
                      )
  #creating data frame tha will store AIC values of models
  AIC_data<-data.frame(id=unique(d$id), Temp=numeric(length = length(unique(d$id))),
                       TempWlin=numeric(length = length(unique(d$id))),
                       TempWpoly=numeric(length = length(unique(d$id))))
  
  #run the analysis
  for(i in unique(d$id)){
    #Arrhenius relationship
    lmTemp<-lm(log(resp_corr)~Temp_surface, d[(d$id==i & !is.na(d$theta_rel)),])
    #Linear relationship
    lmTempWlin<-lm(log(resp_corr)~Temp_surface+theta_rel, d[(d$id==i & !is.na(d$theta_rel)),])
    #Polynomial relationship
    lmTempWpoly<-lm(log(resp_corr)~Temp_surface+poly(theta_rel,2), d[(d$id==i & !is.na(d$theta_rel)),])
    
    #Store the results of anova comparison
    anova_list[[i]]<-anova(lmTemp, lmTempWlin, lmTempWpoly)
    #Store the results of AIC comparison
    AIC_list[[i]]<-AICtab(lmTemp, lmTempWlin, lmTempWpoly, base=T, sort=T, weights=T)
    
    #Store goodness of fit - R2
    R2_data[R2_data$id==i, "Temp"]<-summary(lmTemp)$adj.r.squared
    R2_data[R2_data$id==i, "TempWlin"]<-summary(lmTempWlin)$adj.r.squared
    R2_data[R2_data$id==i, "TempWpoly"]<-summary(lmTempWpoly)$adj.r.squared
    
    
    #Store goodness of fit - AIC
    AIC_data[AIC_data$id==i, "Temp"]<-AIC(lmTemp)
    AIC_data[AIC_data$id==i, "TempWlin"]<-AIC(lmTempWlin)
    AIC_data[AIC_data$id==i, "TempWpoly"]<-AIC(lmTempWpoly)
    
  }
  
  return(list(R2=R2_data, AIC=AIC_data, anova_comp=anova_list, AIC_comp=AIC_list))

}


        