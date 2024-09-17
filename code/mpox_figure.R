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


setwd("C:/Users/yuqian.wang/NTU_Sherry/8Mpox/code")
source("mpox_function.R")
# Setting ######################################################################

Tmin <- 0
Tmax <- 47
#Tmax_pre <- 5

step_size <- 1
times<-c(seq(Tmin,Tmax,step_size))

DL <- log10(10^3) # Detection limit

# Figure 1 & Figure S1 Population fit ##########################################

pop <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/8Mpox/Monolix/r04_all_05/IndividualParameters/estimatedIndividualParameters.txt", row.names = 1)
sim_pop <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/8Mpox/Monolix/r04_all_05/IndividualParameters/simulatedIndividualParameters.txt")  
sim_ind <- split(sim_pop, f = sim_pop$id)

path <- "C:/Users/yuqian.wang/NTU_Sherry/8Mpox/mpox-VL.xlsx"
mysheets <- excel_sheets(path)

original <- lapply(mysheets, function(x) read_excel(path, sheet = x)) 

df_cov <- data.frame(group = c("Skin lesion", "Rectum", "Saliva", "Oropharynx"))

num = 1000
Fit <- list()
ind_fit <- list()

for (g in 1:4){
  group <- df_cov$group[g]
  par <- c(r=pop[g,"r_mean"],
           delta=pop[g,"delta_mean"],
           beta=pop[g,"beta_mean"])
  best_fit <- Covfun(par)
  pars <- sim_ind[[g]][,3:5]
  total_VL <- run_ODE_pop(pars)
  
  MeanVL <- apply(total_VL,1,function(x){quantile(x,0.5,na.rm=T)})
  Min95  <- apply(total_VL,1,function(x){quantile(x,0.025,na.rm=T)})
  Max95  <- apply(total_VL,1,function(x){quantile(x,0.975,na.rm=T)})
  Min50  <- apply(total_VL,1,function(x){quantile(x,0.25,na.rm=T)})
  Max50  <- apply(total_VL,1,function(x){quantile(x,0.75,na.rm=T)})
  
  Fit[[g]] <- cbind(best_fit,MeanVL,Min95,Max95,Min50,Max50,group)
  Fit[[g]] <- data.frame(Fit[[g]])
  colnames(Fit[[g]]) <- c("time","best_fit","MeanVL","Min95","Max95","Min50","Max50","Sample")
  Fit[[g]]$Min95[Fit[[g]]$Min95 < -3] <- -3
  Fit[[g]]$Min50[Fit[[g]]$Min50 < -3] <- -3
  
  ind_df <- original[[g]] %>% 
    mutate(time = Days + 8.3,
           censor = if_else(VL %in% 3, 1, 0))
  
  ind_fit[[g]] <- merge(Fit[[g]], ind_df[,c(3,4,5)], by=c("time"),all=TRUE)
}

combine_pop <- map_df(ind_fit, ~as.data.frame(.x))

combine_pop$Sample <- factor(combine_pop$Sample, levels=c('Skin lesion', 'Rectum', 'Saliva','Oropharynx'))
 
#Calculate viral shedding duration
max <- combine_pop %>%
  filter(time>10 & best_fit<=3) %>%
  group_by(Sample) %>%
  slice(which.max(best_fit))

##Figure S1 Estimated viral load curve##########################################

ggplot(data=combine_pop) +
  geom_point(aes(x=time,y=VL,colour = cut(VL, c(-Inf,3,Inf))),size=0.5, shape=16, stroke = 3) +
  geom_line(aes(x=time,y=best_fit),lwd=1, color ="#7FA2C5") +
  geom_hline(yintercept=DL, linetype="dashed", color = "red") +
  geom_ribbon(aes(x=time,ymin=Min95,ymax=Max95),fill="#7FA2C5",alpha=0.2) +
  geom_ribbon(aes(x=time,ymin=Min50,ymax=Max50),fill="#7FA2C5",alpha=0.4) +
  xlab("Time after infection (Days)") +
  ylab("Viral RNA load/n(copies/ml)")  +
  scale_x_continuous(breaks=seq(0,40,by=4),labels = expression(0,4,8,12,16,20,24,28,32,36,40),limits=c(-1,40)) +
  scale_y_continuous(breaks=seq(-4,10,by=2),labels = expression(10^-4,10^-2,10^0,10^2,10^4,10^6,10^8,10^10),limits=c(-3,10)) +
  scale_color_manual(
    values = c("(-Inf,3]" = "#cc718b",
               "(3, Inf]" = "#7FA2C5"))+
  facet_wrap(vars(Sample))+
  mpox_theme()

ggsave("plot/FigureS1.png", width = 7, height = 4.5,bg = "white")

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
