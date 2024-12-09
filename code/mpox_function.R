#Solve ODE viral dynamic########################################################
Mpoxfun<-function(pars){
  r <- as.numeric(pars[1])
  delta <- as.numeric(pars[2])
  beta <- as.numeric(pars[3])
  v <- as.numeric(pars[4])
  derivs<-function(time,y,pars){
    with(as.list(c(pars,y)),{
      dTa<--beta*Ta*V
      dV<-r*Ta*V-delta*V
      
      return(list(c(dTa,dV)))
    })
  }
  y<-c(Ta=1,V=v)
  
  times<-c(seq(Tmin,Tmax,step_size))
  out<-ode(y=y,parms=pars,times=times,func=derivs)
  out2<-cbind(time=out[,1],aV=((log10(out[,3]))))
  as.data.frame(out2)
}


Mpoxfun_tau<-function(pars){
  r <- as.numeric(pars[1])
  delta <- as.numeric(pars[2])
  beta <- as.numeric(pars[3])
  v <- as.numeric(pars[4])
  derivs<-function(time,y,pars){
    with(as.list(c(pars,y)),{
      dTa<--beta*Ta*V
      dV<-r*Ta*V-delta*V
      
      return(list(c(dTa,dV)))
    })
  }
  
  y<-c(Ta=1,V=v)
  
  times<-c(seq(Tmin,Tmax,step_size))
  times_pre<-c(seq(Tmin,Tmax_pre,-step_size))
  
  out_after<-ode(y=y,parms=pars,times=times,func=derivs)
  out_before<-ode(y=y,parms=pars,times=times_pre,func=derivs)
  #out<-lsoda(y=y,parms=pars,times=times,func=derivs,rtol=0.00004,atol=0.00000000000001)
  out<-rbind(out_before[-1,],out_after)
  out2<-cbind(Day=out[,1],aV=((log10(out[,3]))))
  as.data.frame(out2)
}


#Sample from estimated parameters###############################################
#sample_pars_pop <- function(pop,g){
  
  #pops_mean <- c("r_mean","delta_mean","beta_mean")
  #pops_sd <- c("r_sd","delta_sd","beta_sd")
  
  #pars <- matrix(0, num, length(pops_mean))
  
  #for (i in 1:length(pops_mean)) {
  #  mean_par <- pop[g,pops_mean[i]]
  #  sd_par <- pop[g,pops_sd[i]]
  #  pars[, i] <- exp(rnorm(num, mean=log(mean_par), sd=sd_par))

  #}  

#  return(pars)
#}
sample_pars_pop <- function(pop, num, rectum, saliva){
  
  inds_mean <- which(rownames(pop) %in% c("r_pop","delta_pop","beta_pop","v_pop"))
  inds_sd <- which(rownames(pop) %in% c("omega_r","omega_delta","omega_beta","omega_v"))
  
  pars <- matrix(0, num, length(inds_mean))
  
  for (i in 1:length(inds_mean)) {
    mean_par <- pop$value[inds_mean[i]]
    sd_par <- pop$value[inds_sd[i]]
    
    if (i != 1) {
      beta_rectum <- pop$value[inds_mean[i]+1]
      beta_saliva <- pop$value[inds_mean[i]+2]
      meanlog = log(mean_par) + rectum*beta_rectum + saliva*beta_saliva
      pars[, i] <- exp(rnorm(num, mean=meanlog, sd=sd_par))
      
    #} else if (i==2) {
    #  mean_delta_age <- pop$value[inds_mean[i]+1]
    #  mean_delta_vac <- pop$value[inds_mean[i]+2]
    #  meanlog = log(mean_par)+age*mean_delta_age+vac*mean_delta_vac
    #  pars[, i] <- exp(rnorm(num, mean=meanlog, sd=sd_par))
      
    } else { 
      pars[, i] <- exp(rnorm(num, mean=log(mean_par), sd=sd_par))
      
    }
    
  }  
  
  #mean_r <- pop$value[inds_mean[1]]+vac*mean_r_vac
  #mean_beta <- pop$value[inds_mean[3]]
  #sd_r <- pop$value[inds_sd[1]]
  #sd_beta <- pop$value[inds_sd[3]]
  #corr <- pop$value[12]
  #cov_r_beta <- corr*sd_r*sd_beta
  #sigma <- rbind(c(sd_r^2, cov_r_beta), c(cov_r_beta, sd_beta^2))s
  #mu <- c(log(mean_r), log(mean_beta))
  #res <- mvrnorm(n=num, mu=mu, Sigma=sigma)
  #pars[, c(1, 3)] <- exp(res)
  return(pars)
}




#Solve ODE for sampled parameters###############################################
run_ODE_pop <- function(pars){
  
  total_VL <- matrix(NA,nrow=length(seq(Tmax_pre,Tmax,step_size)),ncol=num)
  #total_VL <- matrix(NA,nrow=length(seq(Tmax_pre,Tmax,step_size)),ncol=num)
  
  for(i in 1:num){
    out <- Mpoxfun_pre(as.numeric(pars[i, ]))
    total_VL[,i] <- out$aV
  }
  return(total_VL)
}


#Individual fitting plot########################################################
ind_fit_plt<-function(Est){
  
  Fit <- list()
  for(i in 1:nrow(Est)){
    #infect <- Est$Infect[i]
    #doi <- Est$DoI[i]
    pars <- c(beta=Est$beta_SAEM[i],
              phi=Est$phi_SAEM[i],
              rho=Est$rho_SAEM[i],
              delta=Est$delta_SAEM[i],
              pi=Est$pi_SAEM[i],
              tau=Est$tau_SAEM[i],
              imtau=Est$imtau_SAEM[i],
              m=Est$m_SAEM[i])
    fitted <- Mpoxfun(pars)
    #fitted_pre <- Covfun_pre(pars)
    #d1 <- rbind(fitted,fitted_pre[-1,])
    index <- which(fitted$aV < -2)[1]
    # Check if an index was found, and then replace the value
    if (!is.na(index)) {
      fitted$aV[index] <- -2
    } 
    
    d1 <- fitted
    Code <- Est$id[i]
    #gender <- Est$Gender[i]
    #booster <- Est$X_1st_Booster_vaccine_type[i]
    #age_cat <- Est$age_cat[i]
    original_subset <- original %>%
      filter(Code == Est$id[i]) %>%
      mutate(Day = Days_from_Onset + round(pars[6])) 
    
    #S <- 100 ###100 repeats
    #P <- matrix(NA,nrow=(Tmax+1),ncol=S) #pre-sym and after-sym
    #P <- matrix(NA,nrow=(Tmax+1),ncol=S)
    
    #for(j in 1:S){
    #  pars <- c(beta=Simulated$beta[j+S*(i-1)],
    #            phi=Simulated$phi[j+S*(i-1)],
    #            rho=Simulated$rho[j+S*(i-1)],
    #            delta=Simulated$delta[j+S*(i-1)],
    #            pi=Simulated$pi[j+S*(i-1)],
    #            tau=Simulated$tau[j+S*(i-1)],
    #            imtau=Simulated$imtau[j+S*(i-1)],
    #            m=Simulated$m[j+S*(i-1)])
    #  out  <- Covfun1(pars)
      #out_pre  <- Covfun_pre(pars)
    #  P[,j] <- c(out$aV)
    #}
    
    #Min95  <- apply(P,1,function(x){quantile(x,0.025,na.rm = TRUE)})
    #Min90  <- apply(P,1,function(x){quantile(x,0.05,na.rm = TRUE)})
    #Max90  <- apply(P,1,function(x){quantile(x,0.95,na.rm = TRUE)})
    #Max95  <- apply(P,1,function(x){quantile(x,0.975,na.rm = TRUE)})
    
    fit <- cbind(d1,Min95,Min90,Max90,Max95,Code)
    fit$Min95[fit$Min95 < -2] <- -2
    fit$Min90[fit$Min90 < -2] <- -2
    fit <- merge(fit, original_subset[,c("Code","Day","VL","Age","Vaccination","Severity")], by=c("Code", "Day"),all=TRUE)
    
    Fit[[i]] <- data.frame(fit)
    
    
  }
  
  ind_fit <- map_df(Fit, ~as.data.frame(.x))
  
  #ind_fit <- merge(ind_fit_unlist, original[,c(1,6,16)], by=c("Code", "Day"),all=TRUE) 
  
  ind_fit$Code <- as.factor(ind_fit$Code)
  #ind_fit <- ind_fit %>% mutate(aV_adjust = ifelse((Day == -1) & (aV < -2 | is.na(aV)) , -1.89, aV))
  #colnames(ind_fit)[11:14] <- c("WT IgG","BA1 IgG","WT IgA","BA1 IgA")
  return(ind_fit)
}














#Calculate false negative rate##################################################
cal_false_neg <- function(value) {
  apply(pred_VL, 1, function(row) mean(row < value))
}


##Calculate confidence interval for false-negative rate#########################
cal_ci <- function(proportions, n) {
  z <- 1.96 # z-score for 95% confidence interval
  se <- sqrt(proportions * (1 - proportions) / n) # Standard error
  lower <- proportions - z * se
  upper <- proportions + z * se
  list(lower = lower, upper = upper)
}


##Calculate false-negative and confidence interval##############################
cal_fn_and_cis <- function(values, names,measure_VL) {
  n <- num
  results_list <- lapply(seq_along(values), function(i) {
    value <- values[i]
    name <- names[i]
    fn <- apply(measure_VL, 1, function(row) mean(row < value))
    #ci <- cal_ci(fn, n)
    data.frame(
      fn = fn)
    #CI_Lower = ci$lower,
    #CI_Upper = ci$upper)
  })
  
  results_df <- do.call(cbind, results_list)
  
  colnames(results_df) <- unlist(lapply(names, function(name) {
    c(name  #figure3/4
      #paste0("FN_", name) #figure2
      #paste0("LowerCI_", name),
      #paste0("UpperCI_", name)
    )
  }))
  
  results_df
}


#Simulation to generate false negative plot######
simulation_false_neg<-function(fn_name){
  pars <- sample_pars_pop(pop,fn_name)
  pred_VL <- run_ODE_pop(pars)
  
  measure_error <- matrix(rnorm(num * length(times), mean = 0, sd = 1)*2.04,  #a=2.04
                          nrow = length(times),ncol = num)
  
  measure_VL <- pred_VL + measure_error
  
  results <- cal_fn_and_cis(DL_values, DL_names,measure_VL)
  
  df_false_neg <- data.frame(times,results) %>% rename_with(~ gsub("^X", "", .), starts_with("X")) %>%
    pivot_longer(cols=-1, names_to = "DL", values_to = "FN") %>%
    mutate(FN = ifelse(times < 9.9 & fn_name == 1,1,FN))##false-neg equal to 1 before symptomatic (skin sample)
  return(df_false_neg)
}

#Calculate probability density function of illness onset at time t##############
cal_ct <- function(fn_fig3){
  ct_plt <- list()
  
  for (g in 1:(length(DL_names)+1)){
    
    if(g < (length(DL_names)+1)){
      
      P <- as.vector(t(subset(fn_fig3, DL == DL_names[g])[3]))
      c <- rep(0,Max_t/dt) # frequency of illness onset (or observed incubation period using time since immigration as a time scale)
      
      for (t in 1:(Max_t/dt)){  # incubation period
        s = k+t*dt;
        f <- rep(0, s/dt) # pdf of incubation period
        L <- rep(0, s/dt) # prob. that patients haven't developed symptoms; survival probability
        j <- rep(0, k/dt) # density of incubating population at infection-age tau at the timing of immigration (t=0); initial age distribution
        
        for (i in 1:(s/dt)){  # incubation period and survival prob.
          f[i]=plnorm(i*dt, myu, sigma)-plnorm((i-1)*dt, myu, sigma)
          L[i]=1-plnorm(i*dt, myu, sigma)
        }
        
        for (tau in 1:(k/dt)){  # the distribution of infected day
          j[tau]=exp(-r*tau*dt)
        }
        j=j/sum(j)
        
        
        for (tau in 1:(k/dt)){  # the distribution of infected day
          j[tau]=j[tau]*L[tau]
        }
        
        prob_1 <- 1-sum(j)
        
        for (tau in 1:(k/dt)){  # the distribution of infected day
          j[tau]=j[tau]*P[tau]
        }
        
        prob_2 <- 1-sum(j)
        #j=j/sum(j)
        
        
        for (tau in t:(s/dt-1)){ 
          c[t] = c[t] + f[tau]*j[tau-t+1]/L[tau-t+1]
        }
      }
      
      c=c/sum(c)
      C=cumsum(c)
      sumc=sum(c)
      f=f/sum(f)
      F=cumsum(f)
      x=seq(1*dt, Max_t, dt)
      y=seq(1*dt, k, dt)
      z=seq(1*dt, k+t*dt,dt)
      
      
      ct_plt[[g]] <- data.frame(times=x,ct=c,DL=DL_names[g],cumc=C,sumc=sumc,prob1=prob_1,prob2=prob_2-prob_1,prob3=1-prob_2)
    } else {
      for (i in 1:(Max_t/dt)){  # incubation period and survival prob.
        c[i]=plnorm(i*dt, myu, sigma)-plnorm((i-1)*dt, myu, sigma)
        #c[i]=dlnorm(i*dt, myu, sigma)
        c=c/sum(c)
        x=seq(1*dt, Max_t, dt)
        ct_plt[[g]] <- data.frame(times=x,ct=c,DL="No tests",cumc=C,sumc=sumc,prob1=0,prob2=0,prob3=1)
      }
    }
  }
  ct_plt_bind <- map_df(ct_plt, ~as.data.frame(.x))
  return(ct_plt_bind)
}


#Calculate 70th,80th,95th percentiles of post entry incubation period###########
cal_ct_tile <- function(ct_plt_bind){
  ct_sub <- ct_plt_bind %>% group_split(DL)
  tile_sub <- list()
  for (i in 1:5){
    if (i<5){
      a<-min(which(ct_sub[[i]][,4] >= 0.7))/10
      b<-min(which(ct_sub[[i]][,4] >= 0.8))/10
      c<-min(which(ct_sub[[i]][,4] >= 0.95))/10
    } else{
      a <- qlnorm(0.7, myu, sigma)
      b <- qlnorm(0.8, myu, sigma)
      c <- qlnorm(0.95, myu, sigma)
    }
    tile_sub[[i]] <- data.frame(tile=c("70%","80%","95%"),duration=c(a,b,c),DL=ct_sub[[i]]$DL[1])
  }
  tile_plt <- map_df(tile_sub, ~as.data.frame(.x))
}


#Theme for ggplot2######
mpox_theme <- function(){
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        #legend.position='none',
        axis.title.y = element_text(size=11,family="sans"),
        axis.title.x = element_text(size=11,family="sans")
  )
}
