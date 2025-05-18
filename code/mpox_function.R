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


Mpoxfun_pre<-function(pars,incu){
  r <- as.numeric(pars[1])
  delta <- as.numeric(pars[2])
  beta <- as.numeric(pars[3])
  v <- as.numeric(pars[4])
  
  if (incu == 0){
    tau <- 8
  } else {
    tau <- round(as.numeric(pars[5]),0)
  }
  
  derivs<-function(time,y,pars){
    with(as.list(c(pars,y)),{
      dTa<--beta*Ta*V
      dV<-r*Ta*V-delta*V
      
      return(list(c(dTa,dV)))
    })
  }
  
  y<-c(Ta=1,V=v)
  
  times<-c(seq(Tmin,(Tmax-tau),step_size))
  times_pre<-c(seq(Tmin,-tau,-step_size))
  
  out_after<-ode(y=y,parms=pars,times=times,func=derivs)
  out_before<-ode(y=y,parms=pars,times=times_pre,func=derivs)
  #out<-lsoda(y=y,parms=pars,times=times,func=derivs,rtol=0.00004,atol=0.00000000000001)
  out<-rbind(out_before[nrow(out_before):2,],out_after)
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
  
  pars <- matrix(0, num, (length(inds_mean)+1))
  
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
      pars[, 5] <- round(rlnormTrunc(num, meanlog=1.917, sdlog=0.592, min=1, max=26),0) #tau:incubation period

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
run_ODE_pop <- function(pars,incu){
  
  #total_VL <- matrix(NA,nrow=length(seq(Tmin,Tmax,step_size)),ncol=num)
  total_VL <- matrix(NA,nrow=length(seq(Tmin,Tmax,step_size)),ncol=num)
  
  for(i in 1:num){
    out <- Mpoxfun_pre(as.numeric(pars[i, ]),incu)
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
    pars <- c(r=Est$r_SAEM[i],
              delta=Est$delta_SAEM[i],
              beta=Est$beta_SAEM[i],
              v=Est$v_SAEM[i],
              tau=8)
    fitted <- Mpoxfun_pre(pars,incu=0)
    #fitted_pre <- Covfun_pre(pars)
    #d1 <- rbind(fitted,fitted_pre[-1,])
    #index <- which(fitted$aV < -2)[1]
    # Check if an index was found, and then replace the value
    #if (!is.na(index)) {
    #  fitted$aV[index] <- -2
    #} 
    
    #d1 <- fitted
    ID_site <- Est$id[i]
    #gender <- Est$Gender[i]
    #booster <- Est$X_1st_Booster_vaccine_type[i]
    #age_cat <- Est$age_cat[i]
    original_subset <- original_ind %>%
      filter(ID_site == Est$id[i]) %>%
      mutate(Day = `Days.post.symptoms.onset`) 
    
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
    
    fit <- cbind(fitted,ID_site,t(pars))
    #fit$Min95[fit$Min95 < -2] <- -2
    #fit$Min90[fit$Min90 < -2] <- -2
    fit <- merge(fit, original_subset[,c("ID_site","VL","censor","Day")], by=c("ID_site", "Day"),all=TRUE)
    
    Fit[[i]] <- data.frame(fit)
    
    
  }
  
  ind_fit <- map_df(Fit, ~as.data.frame(.x)) %>%
    separate(ID_site, into = c("ID", "site"), sep = "[, _/]",remove = FALSE) %>%
    mutate(ID_site = as.factor(ID_site),
           ID = as.factor(ID),
           site = as.factor(site),
           censor = as.factor(censor))
  
  #ind_fit <- merge(ind_fit_unlist, original[,c(1,6,16)], by=c("Code", "Day"),all=TRUE) 
  
  #ind_fit$Code <- as.factor(ind_fit$Code)
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

#fn_site <- fn_site[[1]]
#Simulation to generate false negative plot######
simulation_false_neg<-function(fn_site){
  rectum <- as.numeric(fn_site[2])
  saliva <- as.numeric(fn_site[3])
  
  pars <- sample_pars_pop(pop,num,rectum,saliva)
  pred_VL <- run_ODE_pop(pars,incu=1)
  
  measure_error <- matrix(rnorm(num * nrow(pred_VL), mean = 0, sd = 1)*1.4,  #a=1.4
                          nrow = nrow(pred_VL), ncol = num)
  
  measure_VL <- (pred_VL + measure_error) %>%
    replace(is.na(.), -11) #replace NA to very low value
  
  results <- cal_fn_and_cis(DL_values, DL_names, measure_VL)
  
  
  df_false_neg <- data.frame(times,results) %>% rename_with(~ gsub("^X", "", .), starts_with("X")) %>%
    pivot_longer(cols=-1, names_to = "DL", values_to = "FN") #%>%
    #mutate(FN = ifelse(times < 9.9 & fn_name == 1,1,FN))##false-neg equal to 1 before symptomatic (skin sample)
  return(df_false_neg)
  #return(measure_VL)
}

#Calculate probability density function of illness onset at time t##############
cal_ct <- function(fn_fig3){
  ct_plt <- list()
  
  for (g in 1:(length(DL_names)+1)){
    
    if(g <= length(DL_names)){
      
      P <- as.vector(t(subset(fn_fig3, DL == DL_names[g])[3]))
      c <- rep(0,Max_t/dt) # frequency of illness onset (or observed incubation period using time since immigration as a time scale)
      
      for (t in 1:(Max_t/dt)){  # incubation period
        s = k+t*dt;
        f <- rep(0, s/dt) # pdf of incubation period
        L <- rep(0, s/dt) # prob. that patients haven't developed symptoms; survival probability
        j <- rep(0, k/dt) # density of incubating population at infection-age tau at the timing of immigration (t=0); initial age distribution
        
        for (i in 1:(s/dt)){  # incubation period and survival prob.
          f[i]=plnorm(i*dt, myu, sigma)-plnorm((i-1)*dt, myu, sigma) #cumulative probability density function
          L[i]=1-plnorm(i*dt, myu, sigma)
        }
        
        for (tau in 1:(k/dt)){  
          j[tau]=exp(-r*tau*dt)
        }
        j=j/sum(j)
        
        
        for (tau in 1:(k/dt)){  
          j[tau]=j[tau]*L[tau]
        }
        
        prob_1 <- 1-sum(j)
        
        for (tau in 1:(k/dt)){  
          j[tau]=j[tau]*P[tau]
        }
        
        prob_2 <- 1-sum(j)
        #j=j/sum(j)
        
        
        for (tau in t:(s/dt-1)){ 
          c[t] = c[t] + f[tau]*j[tau-t+1]/L[tau-t+1]
        }
      }
      
      c=c/sum(c)  #should be c*step_size
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
        C=cumsum(c)
        x=seq(1*dt, Max_t, dt)
        ct_plt[[g]] <- data.frame(times=x,ct=c,DL="No tests",cumc=C,sumc=sumc,prob1=0,prob2=0,prob3=1)
     }
    }
  }
  ct_plt_bind <- map_df(ct_plt, ~as.data.frame(.x))
  return(ct_plt_bind)
}

#ct_plt_bind <- ct_plt_bind[[1]]
#Calculate 70th,80th,95th percentiles of post entry incubation period###########
cal_ct_tile <- function(ct_plt_bind){
  ct_sub <- ct_plt_bind %>% group_split(DL)
  tile_sub <- list()
  for (i in 1:5){
    if (i<5){
      a<-min(which(ct_sub[[i]][,4] >= 0.7))*step_size
      b<-min(which(ct_sub[[i]][,4] >= 0.8))*step_size
      c<-min(which(ct_sub[[i]][,4] >= 0.95))*step_size
      median_value<-min(which(ct_sub[[i]][,4] >= 0.5))*step_size
      mean_value<-sum(ct_sub[[i]]$times * ct_sub[[i]]$ct)
      variance_value<-sum((ct_sub[[i]]$times^2) * ct_sub[[i]]$ct) - mean_value^2
      
    } else{
      a <- qlnorm(0.7, myu, sigma)
      b <- qlnorm(0.8, myu, sigma)
      c <- qlnorm(0.95, myu, sigma)
      median_value <- qlnorm(0.5, myu, sigma)
    }
    tile_sub[[i]] <- data.frame(tile=c("70%","80%","95%"),duration=c(a,b,c),median=median_value, mean=mean_value,variance=variance_value,DL=ct_sub[[i]]$DL[1])
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
