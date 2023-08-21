rm(list = ls())
setwd("/Users/lingfengluo/Dropbox (University of Michigan)/Lingfeng Research/R Tutorial Package/surtvep/")

library(devtools)

########################################################################
######fit based on kidney/lung/breast cancer data:
########################################################################

################################################################################################################################################
#get time-to-event data from kidney (z_kidney, event_kidney, and time_kidney)
model.kidney  <- coxtp(z = z_kidney, event = event_kidney, time=time_kidney, lambda = seq(0.1,1,0.1) * ((dim(z)[1])^(1/3)),
                       threads = 3, parallel = TRUE)
IC.kidney <- IC(model.kidney)


beta_cancer_kidney_cancerall        <- get.tvcoef(IC.kidney$model.TIC)
colnames(beta_cancer_kidney_cancerall) <- c("age_50","age_50_59","age_70", 
                                            "race_black", "race_hispanic", "race_asian", "race_other", 
                                            "yod_2009_2012","yod_2013_2017",
                                            "distant", "regional_both", "regional_ext", "regional_lymph", "stage_unknown")
beta_cancer_kidney_cancerall$time   <- time_kidney
beta_cancer_kidney_cancerall$constant   <- 1

confint_kidney <- confint(IC.kidney$model.TIC)
beta_cancer_kidney_cancerall$distant_low <- confint_kidney$tvef$distant[,2]
beta_cancer_kidney_cancerall$distant_up <- confint_kidney$tvef$distant[,3]
beta_cancer_kidney_cancerall$regional_both_low <- confint_kidney$tvef$regional_both[,2]
beta_cancer_kidney_cancerall$regional_both_up <- confint_kidney$tvef$regional_both[,3]
beta_cancer_kidney_cancerall$regional_ext_low <- confint_kidney$tvef$regional_ext[,2]
beta_cancer_kidney_cancerall$regional_ext_up <- confint_kidney$tvef$regional_ext[,3]
beta_cancer_kidney_cancerall$regional_ext_low <- confint_kidney$tvef$regional_ext[,2]
beta_cancer_kidney_cancerall$regional_ext_up <- confint_kidney$tvef$regional_ext[,3]
beta_cancer_kidney_cancerall$regional_lymph_low <- confint_kidney$tvef$regional_lymph[,2]
beta_cancer_kidney_cancerall$regional_lymph_up <- confint_kidney$tvef$regional_lymph[,3]
################################################################################################################################################

################################################################################################################################################
#get time-to-event data from lung (z_lung, event_lung, and time_lung)
model.lung  <- coxtp(z = z_lung, event = event_lung, time=time_lung, lambda = seq(0.1,1,0.1) * ((dim(z)[1])^(1/3)),
                     threads = 3, parallel = TRUE)
IC.lung <- IC(model.lung)


beta_cancer_lung_cancerall        <- get.tvcoef(IC.lung$model.TIC)
colnames(beta_cancer_lung_cancerall) <- c("age_50","age_50_59","age_70", 
                                          "race_black", "race_hispanic", "race_asian", "race_other", 
                                          "yod_2009_2012","yod_2013_2017",
                                          "distant", "regional_both", "regional_ext", "regional_lymph", "stage_unknown")
beta_cancer_lung_cancerall$time   <- time_lung
beta_cancer_lung_cancerall$constant   <- 1

confint_lung <- confint(IC.lung$model.TIC)
beta_cancer_lung_cancerall$distant_low <- confint_lung$tvef$distant[,2]
beta_cancer_lung_cancerall$distant_up <- confint_lung$tvef$distant[,3]
beta_cancer_lung_cancerall$regional_both_low <- confint_lung$tvef$regional_both[,2]
beta_cancer_lung_cancerall$regional_both_up <- confint_lung$tvef$regional_both[,3]
beta_cancer_lung_cancerall$regional_ext_low <- confint_lung$tvef$regional_ext[,2]
beta_cancer_lung_cancerall$regional_ext_up <- confint_lung$tvef$regional_ext[,3]
beta_cancer_lung_cancerall$regional_ext_low <- confint_lung$tvef$regional_ext[,2]
beta_cancer_lung_cancerall$regional_ext_up <- confint_lung$tvef$regional_ext[,3]
beta_cancer_lung_cancerall$regional_lymph_low <- confint_lung$tvef$regional_lymph[,2]
beta_cancer_lung_cancerall$regional_lymph_up <- confint_lung$tvef$regional_lymph[,3]
################################################################################################################################################


################################################################################################################################################
#get time-to-event data from breast (z_breast, event_breast, and time_breast)
model.breast  <- coxtp(z = z_breast, event = event_breast, time=time_breast, lambda = seq(0.1,1,0.1) * ((dim(z)[1])^(1/3)),
                       threads = 3, parallel = TRUE)
IC.breast <- IC(model.breast)


beta_cancer_breast_cancerall        <- get.tvcoef(IC.breast$model.TIC)
colnames(beta_cancer_breast_cancerall) <- c("age_50","age_50_59","age_70", 
                                            "race_black", "race_hispanic", "race_asian", "race_other", 
                                            "yod_2009_2012","yod_2013_2017",
                                            "distant", "regional_both", "regional_ext", "regional_lymph", "stage_unknown")
beta_cancer_breast_cancerall$time   <- time_breast
beta_cancer_breast_cancerall$constant   <- 1

confint_breast <- confint(IC.breast$model.TIC)
beta_cancer_breast_cancerall$distant_low <- confint_breast$tvef$distant[,2]
beta_cancer_breast_cancerall$distant_up <- confint_breast$tvef$distant[,3]
beta_cancer_breast_cancerall$regional_both_low <- confint_breast$tvef$regional_both[,2]
beta_cancer_breast_cancerall$regional_both_up <- confint_breast$tvef$regional_both[,3]
beta_cancer_breast_cancerall$regional_ext_low <- confint_breast$tvef$regional_ext[,2]
beta_cancer_breast_cancerall$regional_ext_up <- confint_breast$tvef$regional_ext[,3]
beta_cancer_breast_cancerall$regional_ext_low <- confint_breast$tvef$regional_ext[,2]
beta_cancer_breast_cancerall$regional_ext_up <- confint_breast$tvef$regional_ext[,3]
beta_cancer_breast_cancerall$regional_lymph_low <- confint_breast$tvef$regional_lymph[,2]
beta_cancer_breast_cancerall$regional_lymph_up <- confint_breast$tvef$regional_lymph[,3]


########################################################################
######generate plots:
########################################################################

pcancer_beta_cancer_breast_cancerall = ggplot(data=beta_cancer_breast_cancerall, aes(x=time)) +
  geom_line(aes(y=distant, color='Distant', linetype = 'Distant'), size = 0.9) +
  geom_ribbon(aes(ymin=distant_low, ymax = distant_up),fill = "red", alpha = 0.2) +
  
  geom_line(aes(y=regional_both, color='Regional both', linetype = 'Regional both'), size = 0.9) +
  geom_ribbon(aes(ymin=regional_both_low, ymax = regional_both_up), fill = "dark green",alpha = 0.2) +
  
  geom_line(aes(y=regional_ext, color='Regional extend', linetype = 'Regional extend'), size = 0.9) +
  geom_ribbon(aes(ymin=regional_ext_low, ymax = regional_ext_up),fill = "orange", alpha = 0.2) +
  
  geom_line(aes(y=regional_lymph, color='Regional lymph', linetype = 'Regional lymph'), size = 0.9) +
  geom_ribbon(aes(ymin=regional_lymph_low, ymax = regional_lymph_up),fill = "blue", alpha = 0.2) +
  
  geom_line(aes(y=constant, color='Localized', linetype = 'Localized'), size = 0.9) +
  
  ggtitle("Cancer Death, Breast") +
  scale_x_continuous(name='Years since diagnosis', limits=c(0,12), breaks = c(seq(0,12,2))) +
  scale_y_continuous(name='Hazard Ratio (log-scale)', limits=c(-1,4.2), breaks=log(c(0.33,0.5,1,2,5,10,20,40,80)), labels = c(0.33,0.5,1,2,5,10,20,40,80)) +
  
  theme_bw() +  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name='Stage', values = c('Distant'='red', 'Localized'='black',
                                              'Regional both' = 'dark green', 'Regional extend'='orange', 'Regional lymph' = 'blue')) +  
  scale_linetype_manual(name='Stage', values = c('Distant'= 2, 'Localized'= 1, 'Regional both'= 6, "Regional extend" = 11, 'Regional lymph' = 8)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(text= element_text(size=14)) + theme(axis.text= element_text(size=14)) +
  theme(axis.title.y = element_text(margin= margin(t=0, r=10, b=0, l=0))) +
  theme(legend.title = element_text(size=18), legend.text = element_text(size=16), legend.key.width = unit(2, 'cm')) +
  theme(plot.title = element_text(size=14))
pcancer_beta_cancer_breast_cancerall


pcancer_beta_cancer_kidney_cancerall = ggplot(data=beta_cancer_kidney_cancerall, aes(x=time)) +
  geom_line(aes(y=distant, color='Distant', linetype = 'Distant'), size = 0.9) +
  geom_ribbon(aes(ymin=distant_low, ymax = distant_up),fill = "red", alpha = 0.2) +
  
  geom_line(aes(y=regional_both, color='Regional both', linetype = 'Regional both'), size = 0.9) +
  geom_ribbon(aes(ymin=regional_both_low, ymax = regional_both_up), fill = "dark green",alpha = 0.2) +
  
  geom_line(aes(y=regional_ext, color='Regional extend', linetype = 'Regional extend'), size = 0.9) +
  geom_ribbon(aes(ymin=regional_ext_low, ymax = regional_ext_up),fill = "orange", alpha = 0.2) +
  
  geom_line(aes(y=regional_lymph, color='Regional lymph', linetype = 'Regional lymph'), size = 0.9) +
  geom_ribbon(aes(ymin=regional_lymph_low, ymax = regional_lymph_up),fill = "blue", alpha = 0.2) +
  
  geom_line(aes(y=constant, color='Localized', linetype = 'Localized'), size = 0.9) +
  
  ggtitle("Cancer Death, Kidney") +
  scale_x_continuous(name='Years since diagnosis', limits=c(0,12), breaks = c(seq(0,12,2))) +
  scale_y_continuous(name='Hazard Ratio (log-scale)', limits=c(-1,4.2), breaks=log(c(0.33,0.5,1,2,5,10,20,40,80)), labels = c(0.33,0.5,1,2,5,10,20,40,80)) +
  
  theme_bw() +  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name='Stage', values = c('Distant'='red', 'Localized'='black',
                                              'Regional both' = 'dark green', 'Regional extend'='orange', 'Regional lymph' = 'blue')) +  
  scale_linetype_manual(name='Stage', values = c('Distant'= 2, 'Localized'= 1, 'Regional both'= 6, "Regional extend" = 11, 'Regional lymph' = 8)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(text= element_text(size=14)) + theme(axis.text= element_text(size=14)) +
  theme(axis.title.y = element_text(margin= margin(t=0, r=10, b=0, l=0))) +
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14)) +
  theme(plot.title = element_text(size=14))
pcancer_beta_cancer_kidney_cancerall

pcancer_beta_cancer_lung_cancerall = ggplot(data=beta_cancer_lung_cancerall, aes(x=time)) +
  geom_line(aes(y=distant, color='Distant', linetype = 'Distant'), size = 0.9) +
  geom_ribbon(aes(ymin=distant_low, ymax = distant_up),fill = "red", alpha = 0.2) +
  
  geom_line(aes(y=regional_both, color='Regional both', linetype = 'Regional both'), size = 0.9) +
  geom_ribbon(aes(ymin=regional_both_low, ymax = regional_both_up), fill = "dark green",alpha = 0.2) +
  
  geom_line(aes(y=regional_ext, color='Regional extend', linetype = 'Regional extend'), size = 0.9) +
  geom_ribbon(aes(ymin=regional_ext_low, ymax = regional_ext_up),fill = "orange", alpha = 0.2) +
  
  geom_line(aes(y=regional_lymph, color='Regional lymph', linetype = 'Regional lymph'), size = 0.9) +
  geom_ribbon(aes(ymin=regional_lymph_low, ymax = regional_lymph_up),fill = "blue", alpha = 0.2) +
  
  geom_line(aes(y=constant, color='Localized', linetype = 'Localized'), size = 0.9) +
  
  ggtitle("Cancer Death, Lung") +
  scale_x_continuous(name='Years since diagnosis', limits=c(0,12), breaks = c(seq(0,12,2))) +
  scale_y_continuous(name='Hazard Ratio (log-scale)', limits=c(-1,4.2), breaks=log(c(0.33,0.5,1,2,5,10,20,40,80)), labels = c(0.33,0.5,1,2,5,10,20,40,80)) +
  
  theme_bw() +  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name='Stage', values = c('Distant'='red', 'Localized'='black',
                                              'Regional both' = 'dark green', 'Regional extend'='orange', 'Regional lymph' = 'blue')) +  
  scale_linetype_manual(name='Stage', values = c('Distant'= 2, 'Localized'= 1, 'Regional both'= 6, "Regional extend" = 11, 'Regional lymph' = 8)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(text= element_text(size=14)) + theme(axis.text= element_text(size=14)) +
  theme(axis.title.y = element_text(margin= margin(t=0, r=10, b=0, l=0))) +
  theme(legend.title = element_text(size=16), legend.text = element_text(size=14)) +
  theme(plot.title = element_text(size=14))
pcancer_beta_cancer_lung_cancerall


plot <- ggarrange(pcancer_beta_cancer_kidney_cancerall+ rremove("ylab") + rremove("xlab"), 
                  pcancer_beta_cancer_lung_cancerall+ rremove("ylab") + rremove("xlab"),
                  pcancer_beta_cancer_breast_cancerall + rremove("ylab") + rremove("xlab"), 
                  ncol = 3, nrow = 1,
                  common.legend = TRUE)
plot <- annotate_figure(plot,
                        left = textGrob("Hazard Ratio (log-scale)",  rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                        bottom = textGrob("Years since diagnosis", gp = gpar(cex = 1.3)))

plot

png(filename=paste0("hr_appendix_stage_joss.png"), height=6, width=13, res=600, units="in")
plot
dev.off()

