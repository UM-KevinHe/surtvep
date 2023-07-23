rm(list = ls())
setwd("/Users/lingfengluo/Dropbox (University of Michigan)/Lingfeng Research/R Tutorial Package/surtvep/")

library(devtools)

load_all()


###
data("ExampleData")
z     <- ExampleData$x
time  <- ExampleData$time
event <- ExampleData$event

###
fit.tv        <- coxtv(z = z, event = event, time=time)
fit.penalize  <- coxtp(z = z, event = event, time=time, lambda = c(1,10,100))

IC <- IC(fit.penalize)

### 
plot.tv <- plot(fit.tv, ylim = c(-2,2), title = "coxtv")
plot.tv
plot.tp <- plot(IC$model.mAIC, ylim = c(-2,2), title = "coxtp")
plot.tp

library(grid)
library(ggpubr)
plot_all <- ggarrange(plot.tv + rremove("xlab") + rremove("ylab"),
                      plot.tp + rremove("xlab") + rremove("ylab"),
                      nrow = 1, ncol = 2,
                      common.legend = TRUE, legend = "right")
# plot_all2 <- annotate_figure(plot_all,
#                             top = text_grob("White, Male", face = "bold", size = 14),
#                             bottom = textGrob("time", rot = 0, vjust = 0, hjust = 0.5, gp = gpar(cex = 1.1)),
#                             left = textGrob("coefficient", rot = 90, vjust = 1,gp = gpar(cex = 1.1)))
plot_all
###
# png(file = "/Users/lingfengluo/Dropbox (University of Michigan)/Lingfeng Research/R Tutorial Package/tutorialPaper_jooss/surtvep-jooss-main/coxtv.png",   # The directory you want to save the file in
#     units="in", width=8, height=5, res=400)
# plot.tv
# dev.off()
# 
# 
# png(file = "/Users/lingfengluo/Dropbox (University of Michigan)/Lingfeng Research/R Tutorial Package/tutorialPaper_jooss/surtvep-jooss-main/coxtp.png",   # The directory you want to save the file in
#     units="in", width=8, height=5, res=400)
# plot.tp
# dev.off()



png(file = "/Users/lingfengluo/Dropbox (University of Michigan)/Lingfeng Research/R Tutorial Package/tutorialPaper_jooss/surtvep-jooss-main/coxtv_tp.png",   # The directory you want to save the file in
    units="in", width=16, height=5, res=400)
plot_all
dev.off()
