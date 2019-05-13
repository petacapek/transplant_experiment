TW<-function(d){
  #reading segmented library
  library(segmented)
  #defining Arrhenius plot temperature
  d$Temp_air<-1/(d$Tair+273.15)/8.314
  
  #creating list that will store anova comparisons
  anova_list<-vector("list", length = length(unique(d$id)))
  names(anova_list)<-unique(d$id)
  #creating list that will store AIC comparisons
  AIC_list<-vector("list", length = length(unique(d$id)))
  names(AIC_list)<-unique(d$id)
  #creating data frame tha will store R2 values of models
  R2_data<-data.frame(id=unique(d$id), Temp=numeric(length = length(unique(d$id))),
                      TempWpoly=numeric(length = length(unique(d$id))),
                      TempWpiece=numeric(length = length(unique(d$id))))
  #creating data frame tha will store AIC values of models
  AIC_data<-data.frame(id=unique(d$id), Temp=numeric(length = length(unique(d$id))),
                      TempWpoly=numeric(length = length(unique(d$id))),
                      TempWpiece=numeric(length = length(unique(d$id))))
  
  #run the analysis
  for(i in unique(data$id)){
    #Arrhenius relationship
    lmA<-lm(log(resp_corr)~Temp_air, d[(d$id==i & !is.na(d$theta_rel)),])
  }

}


        