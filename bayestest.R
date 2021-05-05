library(bayestestR)
library(Rmisc)
bayesttest <- function(posterior, V=0){
  describe_posterior(posterior,
                     centrality = "MAP",
                     ci = 0.95,
                     ci_method = "HDI",
                     test = c("pd","p_map","rope","equitest"),
                     rope_range = c(V*-0.1, V*0.1),
                     rope_ci = 0.95)
}

#1st scenario
#VA=0.1
bayest_gaspesie <- bayesttest(mcmc_gaspesie, V=0.1)
bayest_laurentides <- bayesttest(mcmc_laurentides, V=0.1)
bayest_basstlaurent <- bayesttest(mcmc_basstlaurent, V=0.1)
bayest_bois_franc <- bayesttest(mcmc_bois_franc, V=0.1)
bayest_cote_beaupre <- bayesttest(mcmc_cote_beaupre, V=0.1)
bayest_charlevoix <- bayesttest(mcmc_charlevoix, V=0.1)
#VA=0.3
bayest_gaspesie_2 <- bayesttest(mcmc_gaspesie_2, V=0.3)
bayest_laurentides_2 <- bayesttest(mcmc_laurentides_2, V=0.3)
bayest_basstlaurent_2 <- bayesttest(mcmc_basstlaurent_2, V=0.3)
bayest_bois_franc_2 <- bayesttest(mcmc_bois_franc_2, V=0.3)
bayest_cote_beaupre_2 <- bayesttest(mcmc_cote_beaupre_2, V=0.3)
bayest_charlevoix_2 <- bayesttest(mcmc_charlevoix_2, V=0.3)

bayest_results <- data.frame("Va" = as.factor(c(rep(0.1,30),rep(0.3,30))),
                              "region" = c(rep("Gaspésie",5),rep("Bas-St-Laurent",5),rep("Laurentides",5),
                                          rep("Charlevoix",5),rep("Côte-de-Beaupré",5),rep("Bois-Francs",5),
                                          rep("Gaspésie",5),rep("Bas-St-Laurent",5),rep("Laurentides",5),
                                          rep("Charlevoix",5),rep("Côte-de-Beaupré",5),rep("Bois-Francs",5)),
                             "pd"= c(bayest_gaspesie$pd,bayest_laurentides$pd,bayest_basstlaurent$pd,
                                     bayest_bois_franc$pd,bayest_cote_beaupre$pd,bayest_charlevoix$pd,
                                     bayest_gaspesie_2$pd,bayest_laurentides_2$pd,bayest_basstlaurent_2$pd,
                                     bayest_bois_franc_2$pd,bayest_cote_beaupre_2$pd,bayest_charlevoix_2$pd),
                             "ROPE"= c(bayest_gaspesie$ROPE_Percentage,bayest_laurentides$ROPE_Percentage,bayest_basstlaurent$ROPE_Percentage,
                                       bayest_bois_franc$ROPE_Percentage,bayest_cote_beaupre$ROPE_Percentage,bayest_charlevoix$ROPE_Percentage,
                                       bayest_gaspesie_2$ROPE_Percentage,bayest_laurentides_2$ROPE_Percentage,bayest_basstlaurent_2$ROPE_Percentage,
                                       bayest_bois_franc_2$ROPE_Percentage,bayest_cote_beaupre_2$ROPE_Percentage,bayest_charlevoix_2$ROPE_Percentage))

bayest_results$region<-as.factor(bayest_results$region)
bayest_results$region <- factor(bayest_results$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"))

summary_bayestest <- cbind(summarySE(bayest_results, measurevar = "ROPE", groupvars = c("region","Va"))[,c(1,2,4)],
                           "pd"=summarySE(bayest_results, measurevar = "pd", groupvars = c("region","Va"))[,c(4)])

#plot
ggplot(summary_bayestest, aes(x=region)) + 
 geom_point(aes(y=pd*100, shape = "pd",col = "pd"), size=4) + 
  geom_point(aes(y=ROPE*100, shape = "ROPE",col = "ROPE"),size=4) +
  xlab("Region")+
  scale_y_continuous(
    
    # Features of the first axis
    name = "Probability of direction",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*1, name="ROPE")
  )+
  facet_grid(Va ~ .)+labs(shape = "Bayesian index", color = "Bayesian index") + 
  theme_bw()+theme(axis.text.x=element_text(angle=80, vjust=0.5, hjust = 0.5))

