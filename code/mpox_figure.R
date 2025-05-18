library(deSolve)
library(EnvStats)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggh4x)
library(ggsci)
library(viridis)
library(GGally)
library(cowplot)

rm(list=ls())

# Set your path to the project folder before running the script
# Example: 
MPOXpath <- "~/path_to_project/"


source(paste0(MPOXpath,"Mpox/code/mpox_function.R"))
# Setting ######################################################################

Tmin <- 0
Tmax <- 36
Tmax_pre <- -8

step_size <- 0.1
times<-c(seq(Tmax_pre,Tmax,step_size))

DL <- log10(10^3) # Detection limit

# Figure 1 & Figure S1 Population fit ##########################################
pop <- read.csv(paste0(MPOXpath,"Mpox/Monolix/r10_pre_sym_1point_nontau_inits_limit3_deltabetaV/populationParameters.txt"), row.names = 1)


original <- read.csv(paste0(MPOXpath,"Mpox/data/combine_VL_pre_and_sym_1point_limit3.csv"))
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
           delta=pop["delta_pop","value"]*exp(rectum*pop["beta_delta_site_Rectum", "value"])*exp(saliva*pop["beta_delta_site_Saliva", "value"]),
           beta=pop["beta_pop","value"]*exp(rectum*pop["beta_beta_site_Rectum", "value"])*exp(saliva*pop["beta_beta_site_Saliva", "value"]),
           v=pop["v_pop","value"]*exp(rectum*pop["beta_v_site_Rectum", "value"])*exp(saliva*pop["beta_v_site_Saliva", "value"]),
           tau=8) #mean of incubation period

  best_fit <- Mpoxfun_pre(par,incu=0)
  
  pars <- sample_pars_pop(pop, num, rectum, saliva)
  total_VL <- run_ODE_pop(pars,incu=0)
  
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
  filter(time > 5) %>%
  group_by(Sample) %>%
  mutate(
    shedding_mean = ifelse(MeanVL <= 3, time, NA),
    shedding_best_fit = ifelse(best_fit <= 3, time, NA)
  ) %>%
  summarize(
    #shedding_mean = min(shedding_mean, na.rm = TRUE),
    shedding_best_fit = min(shedding_best_fit, na.rm = TRUE)
  ) %>%
  ungroup()

##Figure S1 Estimated viral load curve##########################################

ggplot(data=combine_pop) +
  geom_point(aes(x=time,y=VL,colour = Sample),shape=16,alpha=0.5) +
  geom_line(aes(x=time,y=MeanVL,colour=Sample),lwd=1) +
  geom_hline(yintercept=DL, linetype="dashed", color = "red",alpha=0.7) +
  geom_ribbon(aes(x=time,ymin=Min95,ymax=Max95,fill=Sample),alpha=0.1) +
  geom_ribbon(aes(x=time,ymin=Min50,ymax=Max50,fill=Sample),alpha=0.3) +
  xlab("Time after symptom onset (Days)") +
  ylab("Viral load (copies/ml)")  +
  scale_x_continuous(breaks=seq(-8,28,by=4),labels = expression(-8,-4,0,4,8,12,16,20,24,28),limits=c(-8,28)) +
  scale_y_continuous(breaks=seq(-2,10,by=2),labels = expression(10^-2,10^0,10^2,10^4,10^6,10^8,10^10),limits=c(-2,10.5)) +
  facet_wrap(vars(Sample))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  mpox_theme()

ggsave("figure/Figure1.png", width = 8, height = 2.5,bg = "white")


#Figure Individual fit#########################################################
original_ind <- read.csv(paste0(MPOXpath,"Mpox/data/combine_VL_pre_and_sym_1point_limit3.csv")) 
#colnames(original)[c(1,16)] <- c("Code","VL")

#individual fit from Monolix
Est <- read.csv(paste0(MPOXpath,"Mpox/Monolix/r10_pre_sym_1point_nontau_inits_limit3_deltabetaV/IndividualParameters/estimatedIndividualParameters.txt"), sep = ",", comment.char = "", header = T)

#Generate figures
ind_fit <- list()
ind_fit <- ind_fit_plt(Est) %>%
  mutate(ID = as.character(ID))

ind_fit$site <- factor(ind_fit$site, levels=c('Rectum', 'Saliva','Oropharynx'))

ind_fit_remove_censorsite <- ind_fit %>%
  group_by(ID_site) %>%
  #mutate(ID = as.character(ID)) %>%
  filter(any(censor == 0)) %>%
  ungroup() 

ind_fit_sub <- ind_fit_remove_censorsite %>%
  mutate(ID = ifelse(
    grepl("^P[1-9]$", ID),
    sprintf("P0%d", as.numeric(sub("P", "", ID))),
    ID))

  ggplot(ind_fit_sub) +
      geom_jitter(aes(x=Day,y=VL,colour=site,shape=censor),size=2,height=0.1,width = 0.2) +
      geom_line(aes(x=Day,y=aV,colour=site),lwd=0.7) +
      facet_wrap(vars(ID), ncol=8, nrow=10)+
      xlab("Day after symptom onset") +
      ylab("Viral RNA load/n(copies/ml)")  +
      scale_x_continuous(breaks=seq(0,30,by=10),labels = expression(0,10,20,30),limits=c(-8,30)) +
      scale_y_continuous(breaks=seq(0,9,by=3),labels = expression(10^0,10^3,10^6,10^9),limits=c(-1,11)) +
      scale_shape_manual(values = c("1" = 1, "0" = 16))+
      scale_color_brewer(palette = "Dark2")+
      guides(shape = "none")+
      mpox_theme()+
      theme(legend.position='bottom')
  
ggsave(paste0(MPOXpath,"figure/FigureS1_remove_censorsite.png"), width = 10, height = 12,bg = "white")

#Correlation oro and saliva####################################################

# Function to generate a ggpairs plot for a given variable prefix
colnames(Est)[2]<-"gamma_SAEM"
Est$`log10(V)_SAEM` <- log10(Est$v_SAEM) #for visualization we use log10 transform here

cor_plot <- function(variable) {
  test <- Est[,c("id","gamma_SAEM","delta_SAEM","beta_SAEM","v_SAEM")] %>%
    separate(id, into = c("ID", "site"), sep = "[, _/]") %>%
    select(ID, site, starts_with(variable)) %>%
    mutate(site = case_when(site == "Rectum" ~ "rectum",
                            site == "Saliva" ~ "saliva",
                            site == "Oropharynx" ~ "oropharynx")) %>%
    pivot_wider(
      names_from = site,
      values_from = starts_with(variable),
      names_glue = paste0(variable,"_{site}")
      )
  
  # Create the plot
  plt <-ggpairs(
    test[,-1]#, # Exclude ID
    #upper = list(continuous = wrap("cor", method = "spearman", size = 4))
    #lower = list(continuous = "points"),
    #diag = list(continuous = "barDiag")
  ) + mpox_theme()+
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)  # Adjust angle and size
    )
    
  
  ggmatrix_gtable <- ggmatrix_gtable(plt)
}

# Generate plots for each variable
plot_r <- cor_plot("gamma")
plot_beta <- cor_plot("beta")
plot_delta <- cor_plot("delta")
plot_v <- cor_plot("v")


# Combine the plots in a 2x2 grid
combined_plot <- plot_grid(
  plot_r, plot_beta,
  plot_delta, plot_v,
  labels = c("A", "B", "C", "D"),
  ncol = 2
)
# Display the combined plot
print(combined_plot)

ggsave(paste0(MPOXpath,"figure/Figure_v_cor.png"), width = 10, height = 8,bg = "white")


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

ggsave(paste0(MPOXpath,"plot/Figure1.png"), width = 7, height = 4.5,bg = "white")

write.csv(combine_pop,paste0(MPOXpath,"output/VLpop.csv"))

# Figure 2 false-negative #######################################################
num = 10000

Tmin <- 0
Tmax <- 30
#Tmax_pre <- -8

step_size <- 1
#times<-c(seq(Tmax_pre,Tmax,step_size))
times<-c(seq(Tmin,Tmax,step_size)) #start from day 0 of infection

# Define values and their corresponding names
DL_values <- c(log10(10),log10(250),log10(1000))
DL_names <- c("10", "250", "1000")

fn_list <- list()

fn_site <- list(c("Rectum", 1,0),
                c("Saliva", 0,1),
                c("Oropharynx",0,0))

combine_fn <- map(fn_site, simulation_false_neg) %>%
  set_names(map_chr(fn_site, ~ .x[1])) %>%
  map_df(~as.data.frame(.x), .id = "Sample")


combine_fn$Sample <- factor(combine_fn$Sample, levels=c('Rectum', 'Saliva','Oropharynx'))
combine_fn$DL <- factor(combine_fn$DL, levels=c('10', '250', '1000'))

#calculate min false-negative rate
min <- combine_fn %>%
  #filter(time>10 & best_fit<=3) %>%
  group_by(Sample,DL) %>%
  slice(which.min(FN))


ggplot(data=combine_fn) +
  geom_step(aes(x=times,y=FN,linetype=DL,colour=Sample),lwd=0.7) +
  facet_wrap(vars(Sample))+
  xlab("Time after infection (Days)")+
  ylab("False-negative rate")+
  ylim(0,1)+
  scale_x_continuous(breaks=seq(0,32,by=4),labels = expression(0,4,8,12,16,20,24,28,32),limits=c(-1,32))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  scale_linetype_manual(values=c("longdash","dotted","solid"))+
  labs(linetype = "Detection limit\n(copies/mL)")+
  guides(color = "none")+  
  mpox_theme()+
  theme(legend.position = "right")


ggsave(paste0(MPOXpath,"figure/Figure_neg_10k.png"), width = 8, height = 2.5,bg = "white")

##visualize simulation data (with tau)################################

total_VL <- pred_VL + measure_error

fit <- cbind(times,MeanVL,Min95,Max95,Min50,Max50)


ggplot(fit) +
  geom_line(aes(x=times,y=MeanVL),lwd=1) +
  geom_ribbon(aes(x=times,ymin=Min95,ymax=Max95),alpha=0.1) +
  geom_ribbon(aes(x=times,ymin=Min50,ymax=Max50),alpha=0.3) +
  xlab("Time after infection (Days)") +
  ylab("Viral RNA load/n(copies/ml)")  +
  mpox_theme()


# Figure 3 & Figure S2 importation incubation ##################################

dt=0.1; # time clock
k=14; # maximum duration of travel (set 3 weeks here; can change to 2weeks/4weeks)
Max_t=30; # maximum day of the simulation
R0=1.5 # reproduction number for mpox
R0=1 # endemic state
T=8.7 # mean of serial interval for mpox
r=(R0-1)/T # growth rate incase the serial interval follows exponential distribution for mpox
#r=log(R0)/T # growth rate incase the serial interval follows rectangular distribution for mpox
myu=1.917 # a parameter for incubation period distribution (lognormal)
sigma=0.592 # a parameter for incubation period distribution (lognormal)

num = 1000
Tmin <- 0
Tmax <- k+Max_t
step_size <- 0.1
times<-c(seq(Tmin,Tmax,step_size))

DL_values <- c(log10(10),log10(250),log10(1000),1000000)
DL_names <- c("10", "250", "1000","HS")

fn_fig3 <- list()

fn_site <- list(c("Rectum", 1,0),
                c("Saliva", 0,1),
                c("Oropharynx",0,0))

fn_fig3 <- map(fn_site, simulation_false_neg) %>%
  set_names(map_chr(fn_site, ~ .x[1]))  

combine_fig3 <- map_df(fn_fig3, ~as.data.frame(.x), .id = "Sample")


##Figure S2 post-entry incubation period distribution###########################


ct_plt_bind <- list()

ct_plt_bind <- map(fn_fig3,cal_ct) 

combine_ct_plt1 <- ct_plt_bind %>%
  map_df(~as.data.frame(.x),.id = "Sample") %>%
  mutate(DL= case_when(DL == "10" ~ "HS+PCR1",
                       DL == "250" ~ "HS+PCR2",
                       DL == "1000" ~ "HS+PCR3",
                       DL == "HS" ~ "HS",
                       DL == "No tests" ~ "No tests"),
         Sample = factor(Sample, levels = c('Rectum', 'Saliva', 'Oropharynx')),
         DL = factor(DL, levels = c('HS+PCR1', 'HS+PCR2', 'HS+PCR3', 'HS', 'No tests'))) %>%
  subset(!(DL== 'No tests'))

combine_ct_plt$group <- "Epidemic"
combine_ct_plt1$group <- "Endemic"
#ombine_ct_sub2$group <- "exp2"

combine_sub2<-as.data.frame(rbind(combine_ct_plt,combine_ct_plt1)) %>%
  mutate(group = factor(group, levels = c('Epidemic','Endemic')))


ggplot(data=combine_ct_plt1) +
  geom_line(aes(x=times,y=ct,linetype=DL,colour=Sample,lwd=DL)) +
  xlab("Time after immigration (days)") +
  ylab("Probability of illness onset")+
  scale_x_continuous(breaks=seq(0,28,by=4),labels = expression(0,4,8,12,16,20,24,28),limits=c(0,28)) +
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  scale_linewidth_manual(values=c(0.4,0.7,0.5,0.6,0.3)) +
  scale_linetype_manual(values=c("solid","longdash","dotted","dotdash","dashed")) +
  facet_wrap(vars(Sample))+
  labs(linetype = "Tests",linewidth = "Tests")+
  guides(color = "none")+ 
  mpox_theme()+
  theme(legend.position = "right")

ggsave(paste0(MPOXpath,"figure/FigureS21.png"), width = 7, height = 3.8,bg = "white")

## Figure 3 bar chart of effectiveness of health screening and PCR##############
pie_plt <- combine_ct_plt[combine_ct_plt$DL != "No tests",c(1,4,7,8,9)] %>% 
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
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 20, vjust = 0.75),
        axis.title.x = element_blank())

ggsave(paste0(MPOXpath,"figure/Figure3_uni.png"), width = 8, height = 2.5,bg = "white")

#Figure 4 70th, 80th, 95th percentile of post-entry incubation period ##########

tile_plt <- list()

tile_plt <- map(ct_plt_bind,cal_ct_tile)
combine_tile <- map_df(tile_plt, ~as.data.frame(.x),.id = "Sample") %>%
  mutate(DL = case_when(DL == "10"  ~ "HS+PCR1", 
                        DL == "250" ~ "HS+PCR2",
                        DL == "1000" ~ "HS+PCR3",
                        DL == "HS" ~ "HS",
                        DL == "No tests" ~ "No tests"))

combine_tile$Sample <- factor(combine_tile$Sample, levels=c('Rectum', 'Saliva','Oropharynx'))
combine_tile$DL <- factor(combine_tile$DL, levels=c('HS+PCR1', 'HS+PCR2', 'HS+PCR3','HS','No tests'))

ggplot(data=combine_tile, aes(x=DL, y=duration)) +
  geom_bar(aes(fill=DL),stat="identity",alpha=0.6,width = 0.65)+
  facet_grid2(vars(tile), vars(Sample), axes = "y",remove_labels = "x")+
  scale_fill_brewer(palette = "Dark2")+
  ylab("Time after immigration (Days)")  +
  theme_bw()+
  mpox_theme()+
  theme(legend.position="none",
        axis.text.x = element_text(angle = 20, vjust = 0.75),
        axis.title.x = element_blank())

ggsave(paste0(MPOXpath,"figure/Figure4_uni.png"), width = 7, height = 3, bg = "white")

write.csv(combine_tile, paste0(MPOXpath,"output/percentile_incubation_exp.csv"))
