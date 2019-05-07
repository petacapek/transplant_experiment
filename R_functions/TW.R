TW<-function(data){
  #The temperature equation is defined here
  T_fun<-function(Temp, A, Ea){
    Tout<-A*exp(-Ea/8.314/(Temp+273.15))
    return(Tout)
  }
  #The temperature - moisture equation is defined here
  TW_fun<-function(Temp, theta, phi, A, Ea, K, Opt, b){
    TWout<-A*exp(-Ea/8.314/(Temp+273.15))*ifelse(theta<Opt*phi, (K+Opt*phi)/(K+theta)*(theta/Opt/phi), ((phi-theta)/(phi-Opt*phi))^b)
    return(TWout)
  }
  
  #First, both equations are fitted to the data across all replicates
  #The estimation function (estim), which returns equation parameters and their uncertainty is defined
  estim<-function(d){#d denotes any data input
    #First, the cost functions for T_fun and TW_fun are defined
    T_fun_cost<-function(x){
      #define parameters and their names
      pars<-x
      names(pars)<-c("A", "Ea")
      #measured respiration rate
      resp<-as.numeric(d[(!is.na(d$theta) & d$outliers=="NO"), "resp_corr"])
      #Temperature
      Temp<-as.numeric(d[(!is.na(d$theta) & d$outliers=="NO"), "Tair"])
      #Moisture
      #theta<-as.numeric(d[(!is.na(d$theta) & d$outliers=="NO"), "theta"])
      #porosity
      #phi<-as.numeric(d[(!is.na(d$theta) & d$outliers=="NO"), "phi"])
      yhat<-T_fun(Temp=Temp, A=pars[["A"]], Ea=pars[["Ea"]])
      
      return(sum((resp-yhat)^2))
    }
    
    TW_fun_cost<-function(x){
      #define parameters and their names
      pars<-x
      names(pars)<-c("A", "Ea", "K", "Opt", "b")
      #measured respiration rate
      resp<-as.numeric(d[(!is.na(d$theta) & d$outliers=="NO"), "resp_corr"])
      #Temperature
      Temp<-as.numeric(d[(!is.na(d$theta) & d$outliers=="NO"), "Tair"])
      #Moisture
      theta<-as.numeric(d[(!is.na(d$theta) & d$outliers=="NO"), "theta"])
      #porosity
      phi<-as.numeric(d[(!is.na(d$theta) & d$outliers=="NO"), "phi"])
      yhat<-TW_fun(Temp=Temp, theta=theta, phi=phi, 
                   A=pars[["A"]], Ea=pars[["Ea"]], K=pars[["K"]], Opt=pars[["Opt"]], b=pars[["b"]])
      
      return(sum((resp-yhat)^2))
    }
    
    #Both cost functions are minimized 
    #The minimization is done in two consecutive steps
    #First, approximate model parameters are estimated using MCMC algorithm
    #These parameters are confined by differential evolution algorithm
    
    #1. Temperature equation
    T_par_mcmc<-modMCMC(f=T_fun_cost, p=c(A=1e11, Ea=65000), 
                        lower=c(A=1, Ea=1000),
                        upper=c(A=1e30, Ea=150000), niter=50000)
    
    #lower and upper limits for parameters estimates are extracted
    T_pl<-summary(T_par_mcmc)["min",]
    T_pu<-summary(T_par_mcmc)["max",]
    
    #these limits are used to find global optimum by DEoptim algorithm
    T_opt_par<-DEoptim(fn=T_fun_cost, lower=T_pl, upper=T_pu, 
                     control = c(itermax = 10000, steptol = 50, reltol = 1e-8, 
                                 trace=FALSE, strategy=3, NP=250))
    
    #2. Temperature-Moisture equation
    TW_par_mcmc<-modMCMC(f=TW_fun_cost, p=c(A=1e11, Ea=65000, K=0.1, Opt=0.65, b=0.75), 
                        lower=c(A=1, Ea=1000, K=1e-5, Opt=0.25, b=0),
                        upper=c(A=1e30, Ea=150000, K=1e5, Opt=0.95, b=1.7), niter=50000)
    
    #lower and upper limits for parameters estimates are extracted
    TW_pl<-summary(TW_par_mcmc)["min",]
    TW_pu<-summary(TW_par_mcmc)["max",]
    
    #these limits are used to find global optimum by DEoptim algorithm
    TW_opt_par<-DEoptim(fn=TW_fun_cost, lower=TW_pl, upper=TW_pu, 
                       control = c(itermax = 10000, steptol = 50, reltol = 1e-8, 
                                   trace=FALSE, strategy=3, NP=250))
    
    #output of the function
    estim_out<-list(T_pars=T_opt_par$optim$bestmem,
                    TW_pars=TW_opt_par$optim$bestmem,
                    T_par_mcmc=T_par_mcmc,
                    TW_par_mcmc=TW_par_mcmc)
    
    return(estim_out)
  }
  
  #The function is run across all replicates
  Across_pars<-estim(data)
  
  #In the following function, the goodness of fit is evaluated
  goodness<-function(d, p){
    
    #1. Temperature equation
    #define parameters and their names
    T_pars<-p$T_pars
    names(T_pars)<-c("A", "Ea")
    #measured respiration rate
    resp<-as.numeric(d[(!is.na(d$theta) & d$outliers=="NO"), "resp_corr"])
    #Temperature
    Temp<-as.numeric(d[(!is.na(d$theta) & d$outliers=="NO"), "Tair"])
    #Moisture
    theta<-as.numeric(d[(!is.na(d$theta) & d$outliers=="NO"), "theta"])
    #porosity
    phi<-as.numeric(d[(!is.na(d$theta) & d$outliers=="NO"), "phi"])
    
    #Temperature equation
    T_yhat<-T_fun(Temp=Temp, A=T_pars[["A"]], Ea=T_pars[["Ea"]])
    
    #log likelihood is calculated
    T_ll=-sum((resp-T_yhat)^2)/2/sd(resp)^2
    #R2 is calculated
    T_R2<-1-sum((resp-T_yhat)^2)/sum((resp-mean(resp))^2)
    T_N=2
    T_AIC<-2*T_N-2*T_ll
    
    #2. Temperature-Moisture equation
    #define parameters and their names
    TW_pars<-p$TW_pars
    names(TW_pars)<-c("A", "Ea", "K", "Opt", "b")
    
    #Temperature-Moisture equation
    TW_yhat<-TW_fun(Temp=Temp, theta=theta, phi=phi, 
                    A=TW_pars[["A"]], Ea=TW_pars[["Ea"]], K=TW_pars[["K"]], Opt=TW_pars[["Opt"]], b=TW_pars[["b"]])
    
    #log likelihood is calculated
    TW_ll=-sum((resp-TW_yhat)^2)/2/sd(resp)^2
    #R2 is calculated
    TW_R2<-1-sum((resp-TW_yhat)^2)/sum((resp-mean(resp))^2)
    TW_N=5
    TW_AIC<-2*TW_N-2*TW_ll
    
    #Create output
    goodness_out<-list(Gfit=data.frame(Equation=c("Temperature", "Temperature-Moisture"),
                                       R2=c(T_R2, TW_R2),
                                       ll=c(T_ll, TW_ll),
                                       N=c(2,5),
                                       AIC=c(T_AIC, TW_AIC)),
                       Yhat=data.frame(Measured=resp, Tpred=T_yhat, TWpred=TW_yhat))
    return(goodness_out)
    
  }
  
  Across_goodness<-goodness(data, Across_pars)
  
  return(Across_goodness)
}


        