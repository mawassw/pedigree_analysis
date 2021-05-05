
##############plots####################
require(ggplot2)
require(tidyr)

pedigrees$pre_size <- as.numeric(pedigrees$pre_size)
pedigrees$post_size <- as.numeric(pedigrees$post_size)

pedigrees_test <- pedigrees[c(1,3,6,7)]

tiff("plot4.tiff", units ="in", width = 12, height = 6, res = 300)
pedigrees_test %>%
  gather("Type", "Value",-region) %>%
  ggplot(aes(region, Value, fill = Type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() + ylab("Individuals") + coord_flip() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,500000))
dev.off()
#pedigree depth
tiff("plot2.tiff", units ="in", width = 8, height = 5, res = 300)
ggplot(data = pedigrees, aes(x=region, y=maxpeddepth)) + 
  geom_bar(stat="identity") + ylab("Pedigree depth(equ. N of generations)")+theme_bw() + coord_flip() +
  scale_y_continuous(expand = c(0, 0), limits = c(0,11), breaks = seq(0,11,1))
dev.off()
#percent founders
tiff("plot3.tiff", units ="in", width = 8, height = 5, res = 300)
ggplot(data = pedigrees[-c(6,13),], aes(x=region, y=missing_parents)) + 
  geom_bar(stat="identity") +ylab("N of individuals with unknown parentage") + theme_bw() + coord_flip() +
  scale_y_continuous(expand = c(0, 0), limits = c(0,480), breaks = seq(0,450,50))
dev.off()
#processing time
tiff("plot5.tiff", units ="in", width = 8, height = 5, res = 300)
ggplot(pedigrees, aes(x=post_size, y=processtime)) + 
  geom_area()+ ylab("processing time (m)") + xlab("Pedigree size (ind)")+ #geom_smooth(method=lm) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,510000)) + 
  geom_label_repel(aes(label = region), 
                   box.padding = 0.35, 
                   point.padding = 0.5, 
                   segment.color = 'grey50')+ theme_bw()
dev.off()

#plot
pedigrees_test_1 <- pedigrees[c(1,10,11)]

tiff("plot6.tiff", units ="in", width = 8, height = 10, res = 300)
ggplot(dat, aes(x=region, y=fertility)) +
  geom_boxplot() + ylab("Fertility")+
  scale_y_continuous(expand = c(0,0), limits = c(0,20), breaks = seq(0,20,5)) +theme_bw() 
dev.off()



#plot correlation between pedigree variables
library(corrplot)
tiff("corrplot.tiff", units = "in", width = 10, height = 10, res = 300)
corrplot.mixed(cor(df_precision[,-c(1,2,5:13)]),
               lower = "number", 
               upper = "circle",
               tl.col = "black",
               tl.pos = "lt")
dev.off()
#############plot simple power analyses
power_gaspesie <- cbind.data.frame(read.table(file = "C:/Users/Crampton-Mawass/Desktop/Work/Pedigree/summaries/power_ped_gaspesie.txt"), "region" = factor(rep("Gaspésie", 51)))
power_bois_franc <- cbind.data.frame(read.table(file = "C:/Users/Crampton-Mawass/Desktop/Work/Pedigree/summaries/power_ped_basstlaurent.txt"), "region" = factor(rep("Bas-St-Laurent", 51)))
power_laurentides <- cbind.data.frame(read.table(file = "C:/Users/Crampton-Mawass/Desktop/Work/Pedigree/summaries/power_ped_laurentides.txt"), "region" = factor(rep("Laurentides", 51)))
power_bois_franc <- cbind.data.frame(read.table(file = "C:/Users/Crampton-Mawass/Desktop/Work/Pedigree/summaries/power_ped_bois_franc.txt"), "region" = factor(rep("Bois-Francs", 51)))
power_charlevoix <- cbind.data.frame(read.table(file = "C:/Users/Crampton-Mawass/Desktop/Work/Pedigree/summaries/power_ped_charlevoix.txt"), "region" = factor(rep("Charlevoix", 51)))
power_cote_beaupre <- cbind.data.frame(read.table(file = "C:/Users/Crampton-Mawass/Desktop/Work/Pedigree/summaries/power_ped_cote_beaupre.txt"), "region" = factor(rep("Côte-de-Beaupré", 51)))

power_all <- rbind(power_gaspesie, power_bois_franc,power_laurentides,power_bois_franc, power_charlevoix, power_cote_beaupre)
power_all$region <- factor(power_all$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"))

tiff("result_power.tiff", units ="in", width = 6, height = 5, res = 300)
ggplot(power_all, aes(x=VA, y=power, group = region, color=region)) + 
  geom_point(alpha=0.2)+geom_line() + ylab("Power") +xlab("Additive genetic variance")+geom_vline(xintercept = 0.05, linetype="dashed")+
  theme_bw()
dev.off()
########
#plot MCMC

tiff("result_mcmc_1st.tiff", units ="in", width = 5, height = 7, res = 300)
p1 <- ggplot(df, aes(x=region, y=additive)) + 
  geom_boxplot() + ylab("Estimated Additive Genetic Variance") +xlab("Region")+
  geom_hline(aes(yintercept = additive), data.frame(Va = c(0.1, 0.3), additive = c(0.1, 0.3)),linetype="dashed")+
  ggtitle("1st Scenario using MCMC")+facet_grid(Va ~ .)+stat_summary(fun="mean", color="black", shape=15)+
  theme_bw()+theme(legend.position = "none",axis.text.x=element_text(angle=80, vjust=0.5, hjust = 0.5),plot.title = element_text(hjust = 0.5))

p1 <- p1+geom_text(data=annotation, aes( x=x, y=y, label=label), 
                   color="red", 
                   size=4, angle=0, fontface="bold")
p1
dev.off()

###plotting model results#####
###1st scenario, simple architecture with Va= 0.1 & 0.3
#plot MCMC

#plot type 1 - not used
#tiff("model_1stscenario/result_mcmc.tiff", units ="in", width = 5, height = 7, res = 300)
#p1 <- ggplot(df, aes(x=region, y=additive)) + 
  #geom_boxplot() + ylab("Estimated Additive Genetic Variance") +xlab("Region")+
  #geom_hline(aes(yintercept = additive), data.frame(Va = c(0.1, 0.3), additive = c(0.1, 0.3)),linetype="dashed")+
  #ggtitle("1st Scenario using MCMC")+facet_grid(Va ~ .)+stat_summary(fun="mean", color="black", shape=15)+
  #theme_bw()+theme(legend.position = "none",axis.text.x=element_text(angle=80, vjust=0.5, hjust = 0.5),plot.title = element_text(hjust = 0.5))

#p1 <- p1+geom_text(data=annotation, aes( x=x, y=y, label=label), 
                   #color="red", 
                   #size=4, angle=0, fontface="bold")
#p1
#dev.off()

#plot type 2 - used
tiff("model_1stscenario/result_mcmc.tiff.tiff", units="in", width=6, height=7, res=300)
ggplot(df_mcmc, aes(x = region, y = Estimate)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper)) + 
  geom_hline(aes(yintercept = Estimate), data.frame(Va = c(0.1, 0.3), Estimate = c(0.1, 0.3)),linetype="dashed") +
  scale_x_discrete(limits = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix","Côte-de-Beaupré","Bois-Francs")) +
  labs(x = "Region",
       y = "Estimated Additive Genetic Variance (Estimate +/- 95% CIs)", color = "Region") +geom_hline(yintercept = 0)+
  theme_bw() + facet_grid(Va ~ .)+coord_flip()+theme(legend.position = "none",axis.text.x=element_text(angle=80, vjust=0.5, hjust = 0.5),plot.title = element_text(hjust = 0.5))
dev.off()

#####2nd scenario, complex architecture (Vm = 0.1) with Va= 0.1 & 0.3
#plot MCMC

#plot type 1 - not used
#tiff("model_complex/result_mcmc_complex.tiff", units ="in", width = 5, height = 7, res = 300)
#p1 <- ggplot(df, aes(x=region, y=additive)) + 
  #geom_boxplot() + ylab("Estimated Additive Genetic Variance") +xlab("Region")+
  #geom_hline(aes(yintercept = additive), data.frame(Va = c(0.1, 0.3), additive = c(0.1, 0.3)),linetype="dashed")+
  #ggtitle("2nd Scenario using MCMC (proper param.)")+facet_grid(Va ~ .)+stat_summary(fun="mean", color="black", shape=15)+
  #theme_bw()+theme(legend.position = "none",axis.text.x=element_text(angle=80, vjust=0.5, hjust = 0.5),plot.title = element_text(hjust = 0.5))

#p1 <- p1+geom_text(data=annotation, aes( x=x, y=y, label=label), 
                   #color="red", 
                   #size=4, angle=0, fontface="bold")
#p1
#dev.off()
#Variable Va
tiff("model_complex/result_mcmc_complex.tiff", units="in", width=6, height=7, res=300)
ggplot(df_mcmc_complex, aes(x = region, y = Estimate, group=Var)) +
  geom_pointrange(position=position_dodge(width=0.5),aes(ymin = Lower, ymax = Upper, color=Var)) + 
  geom_hline(aes(yintercept = Estimate), data.frame(Va = c(0.1,0.3, 0.3), Estimate = c(0.1, 0.1, 0.3)),linetype="dashed", col=c("black","cyan3","lightcoral")) +
  scale_x_discrete(limits = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix","Côte-de-Beaupré","Bois-Francs")) +
  labs(x = "Region",
       y = "Estimated Variance (Estimate +/- 95% CIs)", color = "Variance type") +geom_hline(yintercept = 0)+
  theme_bw() + facet_grid(Va ~ .)+coord_flip()+theme(legend.position = "bottom",axis.text.x=element_text(angle=80, vjust=0.5, hjust = 0.5),plot.title = element_text(hjust = 0.5))
dev.off()
#Variable Vm
tiff("model_complex/result_mcmc_complex_2.tiff", units="in", width=6, height=7, res=300)
ggplot(df_mcmc_complex_2, aes(x = region, y = Estimate, group=Var)) +
  geom_pointrange(position=position_dodge(width=0.5),aes(ymin = Lower, ymax = Upper, color=Var)) + 
  geom_hline(aes(yintercept = Estimate), data.frame(Vm = c(0.3,0.3, 0.5,0.5), Estimate = c(0.1,0.3,0.1, 0.5)),linetype="dashed",col=c("lightcoral","cyan3","lightcoral","cyan3")) +
  scale_x_discrete(limits = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix","Côte-de-Beaupré","Bois-Francs")) +
  labs(x = "Region",
       y = "Estimated Variance (Estimate +/- 95% CIs)", color = "Variance type") +geom_hline(yintercept = 0)+scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5))+
  theme_bw() + facet_grid(Vm ~ .)+coord_flip()+theme(legend.position = "bottom",axis.text.x=element_text(angle=80, vjust=0.5, hjust = 0.5),plot.title = element_text(hjust = 0.5))
dev.off()

#####3rd scenario, simple architecture with Va= 0.1 & 0.3 and r=0.1 & 0.3
#plot MCMC areas
library(ggpubr)
#for r=0.1
mcmc_area_0.1 <- mcmc_areas(df_0.1, prob = 0.95)+xlab("Estimated additive variance")+geom_vline(xintercept = 0.1, linetype="dashed")
mcmc_area_0.3 <- mcmc_areas(df_0.3,prob = 0.95)+xlab("Estimated additive variance")+geom_vline(xintercept = 0.3, linetype="dashed")

tiff("model_bivariate/result_mcmc_bivariate_areas.tiff", units ="in", width = 8, height = 5, res = 300)
ggarrange(mcmc_area_0.1, mcmc_area_0.3)
dev.off()
#for r=0.5
mcmc_area_0.1_1 <- mcmc_areas(df_0.1_1, prob = 0.95)+xlab("Estimated additive variance")+geom_vline(xintercept = 0.1, linetype="dashed")
mcmc_area_0.3_1 <- mcmc_areas(df_0.3_1,prob = 0.95)+xlab("Estimated additive variance")+geom_vline(xintercept = 0.3, linetype="dashed")

tiff("model_bivariate/result_mcmc_bivariate_areas_1.tiff", units ="in", width = 8, height = 5, res = 300)
ggarrange(mcmc_area_0.1_1, mcmc_area_0.3_1)
dev.off()

#plot genetic covariance
#for r=0.1
mcmc_area_r_0.1 <- mcmc_areas(df_r_0.1, prob = 0.95)+xlab("Estimated additive covariance")+geom_vline(xintercept = 0.01, linetype="dashed")
mcmc_area_r_0.1_1 <- mcmc_areas(df_r_0.1_1,prob = 0.95)+xlab("Estimated additive covariance")+geom_vline(xintercept = 0.0173, linetype="dashed")

tiff("model_bivariate/result_mcmc_bivariate_areas_r.tiff", units ="in", width = 8, height = 5, res = 300)
ggarrange(mcmc_area_r_0.1, mcmc_area_r_0.1_1)
dev.off()
#for r=0.5
mcmc_area_r_0.5 <- mcmc_areas(df_r_0.5, prob = 0.95)+xlab("Estimated additive covariance")+geom_vline(xintercept = 0.05, linetype="dashed")
mcmc_area_r_0.5_1 <- mcmc_areas(df_r_0.5_1,prob = 0.95)+xlab("Estimated additive covariance")+geom_vline(xintercept = 0.0866, linetype="dashed")

tiff("model_bivariate/result_mcmc_bivariate_areas_r_1.tiff", units ="in", width = 8, height = 5, res = 300)
ggarrange(mcmc_area_r_0.5, mcmc_area_r_0.5_1)
dev.off()
