library(deSolve)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggh4x)
library(ggsci)
library(viridis)

rm(list=ls())


setwd("C:/Users/yuqian.wang/NTU_Sherry/8Mpox/submit")
source("mpox_function.R")
# Setting ######################################################################

Tmin <- 0
Tmax <- 28
Tmax_pre <- -8

step_size <- 0.1
times<-c(seq(Tmin,Tmax,step_size))

DL <- log10(10^3) # Detection limit

# Figure 1 & Figure S1 Population fit ##########################################
#old
#pop <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/8Mpox/Monolix/r04_all_05/IndividualParameters/estimatedIndividualParameters.txt", row.names = 1)
#sim_pop <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/8Mpox/Monolix/r04_all_05/IndividualParameters/simulatedIndividualParameters.txt")  
#sim_ind <- split(sim_pop, f = sim_pop$id)

#new
pop <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/8Mpox/Monolix/r10_pre_sym_1point_nontau_inits_limit3_deltabetaV/populationParameters.txt", row.names = 1)
#sim_pop <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/8Mpox/Monolix/r04_all_05/IndividualParameters/simulatedIndividualParameters.txt")  
#sim_ind <- split(sim_pop, f = sim_pop$id)


#path <- "C:/Users/yuqian.wang/NTU_Sherry/8Mpox/mpox-VL.xlsx"
#mysheets <- excel_sheets(path)
#original <- lapply(mysheets, function(x) read_excel(path, sheet = x)) 

original <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/8Mpox/submit/data/combine_VL_pre_and_sym_1point_limit3.csv")
original$site <- factor(original$site, levels=c('Rectum', 'Saliva','Oropharynx'))
original_site <- split(original, f = original$site)


df_cov <- data.frame(group = c("Rectum", "Saliva", "Oropharynx"),
                     rectum = c(1, 0, 0),
                     saliva = c(0, 1, 0))

num = 10000
Fit <- list()
ind_fit <- list()

for (g in 1:3){
  group <- df_cov$group[g]
  rectum <- df_cov$rectum[g]
  saliva <- df_cov$saliva[g]
  par <- c(r=pop["r_pop","value"],
           #r=pop["r_pop","value"]*exp(rectum*pop["beta_r_site_Rectum", "value"])*exp(saliva*pop["beta_r_site_Saliva", "value"]),
           #delta=pop["delta_pop","value"],
           delta=pop["delta_pop","value"]*exp(rectum*pop["beta_delta_site_Rectum", "value"])*exp(saliva*pop["beta_delta_site_Saliva", "value"]),
           #beta=pop["beta_pop","value"],
           beta=pop["beta_pop","value"]*exp(rectum*pop["beta_beta_site_Rectum", "value"])*exp(saliva*pop["beta_beta_site_Saliva", "value"]),
           #v=pop["v_pop","value"])
           v=pop["v_pop","value"]*exp(rectum*pop["beta_v_site_Rectum", "value"])*exp(saliva*pop["beta_v_site_Saliva", "value"]))
  #par <- c(r=pop["gamma_pop","value"],
  #         delta=pop["delta_pop","value"],
  #         beta=pop["beta1_pop","value"]*(10^-5),
  #         v=0.01)
  
  #n=227
  #par<-c(
  #  r=original_ind[n,"r_SAEM"],
  #  delta=original_ind[n,"delta_SAEM"],
  #  beta=original_ind[n,"beta_SAEM"],
  #  v=original_ind[n,"v_SAEM"])
  
  #best_fit <- Mpoxfun(par)
  best_fit <- Mpoxfun_pre(par)
  #plot(best_fit)
  
  #pars <- sim_ind[[g]][,3:5]
  pars <- sample_pars_pop(pop, num, rectum, saliva)
  total_VL <- run_ODE_pop(pars)
  
  MeanVL <- apply(total_VL,1,function(x){quantile(x,0.5,na.rm=T)})
  Min95  <- apply(total_VL,1,function(x){quantile(x,0.025,na.rm=T)})
  Max95  <- apply(total_VL,1,function(x){quantile(x,0.975,na.rm=T)})
  Min50  <- apply(total_VL,1,function(x){quantile(x,0.25,na.rm=T)})
  Max50  <- apply(total_VL,1,function(x){quantile(x,0.75,na.rm=T)})
  
  Fit[[g]] <- cbind(best_fit,MeanVL,Min95,Max95,Min50,Max50,group)
  Fit[[g]] <- data.frame(Fit[[g]])
  colnames(Fit[[g]]) <- c("time","best_fit","MeanVL","Min95","Max95","Min50","Max50","Sample")
  Fit[[g]]$Max95[Fit[[g]]$Max95 > 10] <- 10
  Fit[[g]]$Max50[Fit[[g]]$Max50 > 10] <- 10
  Fit[[g]]$Min95[Fit[[g]]$Min95 < -2] <- -2
  Fit[[g]]$Min50[Fit[[g]]$Min50 < -2] <- -2
  
  ind_df <- original_site[[g]] %>% 
    mutate(time = Days.post.symptoms.onset)
           #censor = if_else(VL %in% 3, 1, 0))
  
  ind_fit[[g]] <- merge(Fit[[g]], ind_df[,c("time","VL","censor")], by=c("time"),all=TRUE)
}

combine_pop <- map_df(ind_fit, ~as.data.frame(.x))

combine_pop$Sample <- factor(combine_pop$Sample, levels=c('Rectum', 'Saliva','Oropharynx'))
 
#Calculate viral shedding duration
shedding <- combine_pop %>%
  filter(time > 3) %>%
  group_by(Sample) %>%
  mutate(
    shedding_mean = ifelse(MeanVL <= 3, time, NA),
    shedding_best_fit = ifelse(best_fit <= 3, time, NA)
  ) %>%
  summarize(
    shedding_mean = min(shedding_mean, na.rm = TRUE),
    shedding_best_fit = min(shedding_best_fit, na.rm = TRUE)
  ) %>%
  ungroup()

##Figure S1 Estimated viral load curve##########################################

ggplot(data=combine_pop) +
  #geom_point(aes(x=time,y=VL,colour = cut(VL, c(-Inf,3,Inf))),size=0.5, shape=16, stroke = 3) +
  #geom_point(aes(x=time,y=VL,colour = Sample),shape=16,alpha=0.5) +
  geom_jitter(aes(x = time, y = VL, colour = Sample), height = 0.3, width = 0.3, size = 0.8, shape = 16)+
  #geom_jitter(aes(x=time,y=VL,colour = Sample)) +
  #geom_line(aes(x=time,y=best_fit),lwd=1, color ="#7FA2C5") +
  geom_line(aes(x=time,y=MeanVL,colour=Sample),lwd=1) +
  #geom_line(aes(x=time,y=MeanVL),lwd=1, color ="green") +
  geom_hline(yintercept=DL, linetype="dashed", color = "red",alpha=0.7) +
  #geom_ribbon(aes(x=time,ymin=Min95,ymax=Max95),fill="#7FA2C5",alpha=0.2) +
  geom_ribbon(aes(x=time,ymin=Min95,ymax=Max95,fill=Sample),alpha=0.1) +
  geom_ribbon(aes(x=time,ymin=Min50,ymax=Max50,fill=Sample),alpha=0.3) +
  #geom_ribbon(aes(x=time,ymin=Min50,ymax=Max50),fill="#7FA2C5",alpha=0.4) +
  #geom_ribbon(aes(x=time,ymin=Min50,ymax=Max50,fill=Sample),alpha=0.1) +
  xlab("Time after symptom onset (Days)") +
  ylab("Viral RNA load/n(copies/ml)")  +
  #scale_x_continuous(breaks=seq(-8,28,by=4),labels = expression(-8,-4,0,4,8,12,16,20,24,28),limits=c(-9,28)) +
  scale_x_continuous(breaks=seq(-8,28,by=4),labels = expression(-8,-4,0,4,8,12,16,20,24,28),limits=c(-8,28)) +
  scale_y_continuous(breaks=seq(-2,10,by=2),labels = expression(10^-2,10^0,10^2,10^4,10^6,10^8,10^10),limits=c(-2,10.5)) +
  #scale_color_manual(
  #  values = c("(-Inf,3]" = "#cc718b",
  #             "(3, Inf]" = "#7FA2C5"))#+
  facet_wrap(vars(Sample))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  mpox_theme()

ggsave("figure/Figure1.png", width = 8, height = 2.5,bg = "white")

#Figure Individual fit#########################################################
original_ind <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/8Mpox/Monolix/r10_pre_sym_1point_nontau_inits_RDeltaBetaV/IndividualParameters/estimatedIndividualParameters.txt") #%>%
  #mutate(Severity = case_when(intubated_ever == 1|icu_ever == 1|death == 1  ~ "Critical",
  #                            supp_o2_ever == 0 ~ "Mild",
  #                           supp_o2_ever == 1 ~ "Severe"),
  #       Vaccination = case_when(vaccinated == 0 ~ "No",
  #                               vaccinated == 1 ~ "Yes"))
colnames(original)[c(1,16)] <- c("Code","VL")

#target-cell limited model
Est <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/8Mpox/Monolix/r10_pre_sym_1point_nontau_inits_limit3_deltabetaV/IndividualParameters/estimatedIndividualParameters.txt", sep = ",", comment.char = "", header = T)
Simulated <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/8Mpox/Monolix/r10_pre_sym_1point_nontau_inits_limit3_deltabetaV/IndividualParameters/simulatedIndividualParameters.txt", sep = ",", comment.char = "", header = T)



#Generate figures
ind_fit <- list()
ind_fit <- ind_fit_plt1(Est,Simulated)

pdf(paste0("C:/Users/yuqian.wang/NTU_Sherry/2trial_design/5code/github/Figure/r12", ".pdf"), 11, 10)
for (i in seq(1, length(unique(ind_fit$Code)), 49)) {
  print(
    ggplot(ind_fit[ind_fit$Code %in% levels(ind_fit$Code)[i:(i+48)],]) +
      geom_point(aes(x=Day,y=VL,colour = cut(VL, c(-Inf,-1.89,Inf))),size=0.5, shape=16, stroke = 3) +
      #geom_point(aes(x=Day,y=VL,colour = cut(VL, c(-Inf, 17, 19, Inf))),color="#cc718b",size=0.5, shape=16, stroke = 3) +
      geom_line(aes(x=Day,y=aV),lwd=1, color ="#7FA2C5") +
      geom_ribbon(aes(x=Day,ymin=Min90,ymax=Max90), fill="#7FA2C5", alpha=0.2) +
      geom_text(aes(x=23,y=11,label = paste(
        ifelse(!is.na(Age), paste(Age,","), ""),
        ifelse(!is.na(Vaccination), paste(Vaccination,","), ""),
        ifelse(!is.na(Severity), Severity, ""))), size=3) + 
      facet_wrap(vars(Code), ncol=7, nrow=7)+
      xlab("Day after symptom onset") +
      ylab("Viral RNA load/n(copies/ml)")  +
      scale_x_continuous(breaks=seq(-10,40,by=10),labels = expression(-10,0,10,20,30,40),limits=c(-5,41)) +
      scale_y_continuous(breaks=seq(-2,10,by=3),labels = expression(10^-2,10^1,10^4,10^7,10^10),limits=c(-3,12)) +
      scale_color_manual(#name = "qsec",
        values = c("(-Inf,-1.89]" = "#cc718b",
                   "(-1.89, Inf]" = "#7FA2C5"),
        labels = c("<= 17", "17 < qsec <= 19", "> 19"))+
      theme(axis.text = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position='none',
            axis.title.y = element_text(size=11,family="sans"),
            axis.title.x = element_text(size=11,family="sans")))
  
}
dev.off()





##Figure 1 Estimated viral load of different sites##############################

ggplot(data=combine_pop) +
  geom_line(aes(x=time,y=best_fit,color =Sample),lwd=1) +
  geom_hline(yintercept=DL, linetype="dashed", color = "red") +
  geom_ribbon(aes(x=time,ymin=Min95,ymax=Max95,fill=Sample),alpha=0.4) +
  xlab("Time after infection (Days)") +
  ylab("Viral RNA load/n(copies/ml)")  +
  scale_x_continuous(breaks=seq(0,40,by=4),labels = expression(0,4,8,12,16,20,24,28,32,36,40),limits=c(-1,40)) +
  scale_y_continuous(breaks=seq(-4,10,by=2),labels = expression(10^-4,10^-2,10^0,10^2,10^4,10^6,10^8,10^10),limits=c(-3,10)) +
  mpox_theme()+
  theme(legend.position = "right")

ggsave("plot/Figure1.png", width = 7, height = 4.5,bg = "white")

write.csv(combine_pop,"output/VLpop.csv")

# Figure 2 false-negative #######################################################

num = 10000
Tmin <- 0
Tmax <- 37
step_size <- 1
times <- c(seq(Tmin,Tmax,step_size))

# Define values and their corresponding names
DL_values <- c(log10(10),log10(250),log10(1000))
DL_names <- c("10", "250", "1000")

fn_list <- list()

fn_name <- list(Skin = 1,
                Rectum = 2,
                Saliva=3,
                Oropharynx=4)

fn_list <- map(fn_name, simulation_false_neg)

combine_fn <- map_df(fn_list, ~as.data.frame(.x), .id = "Sample") %>%
  mutate(Sample = ifelse(Sample == "Skin", "Skin lesion", Sample))
  
combine_fn$Sample <- factor(combine_fn$Sample, levels=c('Skin lesion', 'Rectum', 'Saliva','Oropharynx'))
combine_fn$DL <- factor(combine_fn$DL, levels=c('10', '250', '1000'))

#calculate min false-negative rate
min <- combine_fn %>%
  #filter(time>10 & best_fit<=3) %>%
  group_by(Sample,DL) %>%
  slice(which.min(FN))


ggplot(data=combine_fn) +
  geom_step(aes(x=times,y=FN,linetype=DL,colour=Sample),lwd=0.7) +
  facet_wrap(vars(Sample))+
  xlab("Time after infection (Days)") +
  ylab("False-negative rate")  +
  scale_x_continuous(breaks=seq(0,40,by=4),labels = expression(0,4,8,12,16,20,24,28,32,36,40),limits=c(-1,40)) +
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  scale_linetype_manual(values=c("longdash","dotted","solid")) +  
  labs(linetype = "Detection limit\n(copies/mL)")+
  guides(color = "none")+  
  mpox_theme()+
  theme(legend.position = "right")

ggsave("plot/Figure2.png", width = 7, height = 4.2,bg = "white")

# Figure 3 & Figure S2 importation incubation ##################################

dt=0.1; # time clock
k=17; # maximum day of incubation period of mpox
Max_t=20; # maximum day of the simulation
R0=1.3 # reproduction number for mpox
T=8.7 # mean of serial interval for mpox
r=(R0-1)/T # growth rate incase the serial interval follows exponential distribution for mpox
#r=log(R0)/T # growth rate incase the serial interval follows rectangular distribution for mpox
myu=log(9.9) # a parameter for incubation period distribution (lognormal)
sigma=0.3 # a parameter for incubation period distribution (lognormal)

num = 10000
Tmin <- 0
Tmax <- 37
step_size <- 0.1
times <- c(seq(Tmin,Tmax,step_size))

DL_values <- c(log10(10),log10(250),log10(1000),1000000)
DL_names <- c("10", "250", "1000","HS")


fn_name <- list(Skin = 1,
                Rectum = 2,
                Saliva=3,
                Oropharynx=4)

fn_fig3 <- list()
fn_fig3 <- map(fn_name, simulation_false_neg)


##Figure S2 post-entry incubation period distribution###########################

ct_plt_bind <- list()

ct_plt_bind <- map(fn_fig3,cal_ct)
combine_ct <- map_df(ct_plt_bind, ~as.data.frame(.x),.id = "Sample") %>%
  mutate(Sample = ifelse(Sample == "Skin", "Skin lesion", Sample),
         DL= case_when(DL == "10" ~ "HS+PCR1",
                       DL == "250" ~ "HS+PCR2",
                       DL == "1000" ~ "HS+PCR3",
                       DL == "HS" ~ "HS",
                       DL == "No tests" ~ "No tests"))

combine_ct$Sample <- factor(combine_ct$Sample, levels=c('Skin lesion', 'Rectum', 'Saliva','Oropharynx'))
combine_ct$DL <- factor(combine_ct$DL, levels=c('HS+PCR1', 'HS+PCR2', 'HS+PCR3','HS','No tests'))

combine_ct_sub <- subset(combine_ct, !((Sample == 'Skin lesion' & DL %in% c('HS+PCR1', 'HS+PCR2', 'HS+PCR3'))| DL== 'No tests'))

ggplot(data=combine_ct_sub) +
  geom_line(aes(x=times,y=ct,linetype=DL,colour=Sample,lwd=DL)) +
  xlab("Time after immigration (days)") +
  ylab("Probability of illness onset")+
  scale_x_continuous(breaks=seq(0,20,by=4),labels = expression(0,4,8,12,16,20),limits=c(0,20)) +
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  scale_linewidth_manual(values=c(0.4,0.7,0.5,0.6,0.3)) +
  scale_linetype_manual(values=c("solid","longdash","dotted","dotdash","dashed")) +
  facet_wrap(vars(Sample))+
  labs(linetype = "Tests",linewidth = "Tests")+
  guides(color = "none")+ 
  mpox_theme()+
  theme(legend.position = "right")

ggsave("plot/FigureS2.png", width = 7.5, height = 4.7,bg = "white")


## Figure 3 bar chart of effectiveness of health screening and PCR##############

combine_ct_sub <- subset(combine_ct, Sample != 'Skin lesion')

pie_plt <- combine_ct_sub[combine_ct_sub$DL != "No tests",c(1,4,7,8,9)] %>% 
  group_by(Sample,DL) %>% 
  slice(1) %>%
  pivot_longer(-c(1:2), names_to = "type",values_to = "prob") %>%
  mutate(type = case_when(type == "prob1"  ~ "Health Screening",
                          type == "prob2"  ~ "PCR",
                          type == "prob3"  ~ "Undetected"),
         labels = scales::percent(prob %>% round(2)))

pie_plt$type <- factor(pie_plt$type, levels=c('Undetected','PCR', 'Health Screening'))

ggplot(pie_plt, aes(x = DL, y = prob, fill=type)) +
  geom_bar(stat = "identity", position="stack") +
  facet_wrap(vars(Sample))+
  scale_fill_brewer(palette="Blues")+
  labs(fill="Infected traveler") +
  ylab("Infected traveler")+
  geom_text(aes(label = labels), size=3.5,
            position = position_stack(vjust = 0.5))+
  mpox_theme()+
  theme(legend.position = "right")

ggsave("plot/Figure3.png", width = 10, height = 3.1,bg = "white")

#Figure 4 70th, 80th, 95th percentile of post-entry incubation period ##########

tile_plt <- list()

tile_plt <- map(ct_plt_bind,cal_ct_tile)
combine_tile <- map_df(tile_plt, ~as.data.frame(.x),.id = "Sample")


combine_tile <- combine_tile %>%
  mutate(DL = case_when(DL == "10"  ~ "HS+PCR1", 
                        DL == "250" ~ "HS+PCR2",
                        DL == "1000" ~ "HS+PCR3",
                        DL == "HS" ~ "HS",
                        DL == "No tests" ~ "No tests"))

combine_tile_sub <- subset(combine_tile, Sample != 'Skin')
combine_tile_sub$Sample <- factor(combine_tile_sub$Sample, levels=c('Rectum', 'Saliva','Oropharynx'))
combine_tile_sub$DL <- factor(combine_tile_sub$DL, levels=c('HS+PCR1', 'HS+PCR2', 'HS+PCR3','HS','No tests'))

ggplot(data=combine_tile_sub, aes(x=DL, y=duration)) +
  geom_bar(aes(fill=DL),stat="identity",alpha=0.6,width = 0.65)+
  scale_y_continuous(breaks=seq(0,15,by=5),labels = expression(0,5,10,15),limits=c(0,16.3)) +
  facet_grid2(vars(Sample),vars(tile), axes = "y",remove_labels = "x")+
  scale_fill_brewer(palette = "Dark2")+
  ylab("Time after immigration (Days)")  +
  theme_bw()+
  mpox_theme()

ggsave("plot/Figure4.png", width = 10, height = 4,bg = "white")

write.csv(combine_tile_sub, "output/percentile_incubation.csv")
