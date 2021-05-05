

################## All Pedigrees ###################
##pedigree statistics
pedigrees <- read.table(file = "C:/Users/walid/OneDrive/Bureau/Work/Pedigree/data/pedigrees.txt")

pedigrees <- pedigrees[c(1,2,3,4,6,10),]

pedigrees$region<- factor(c("Bas-St-Laurent","Bois-Francs", "Charlevoix", "Côte-de-Beaupré","Gaspésie","Laurentides"))

pedigrees$region <- factor(pedigrees$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"))

library(readxl)

pedigrees_all <- as.data.frame(read_excel("C:/Users/walid/OneDrive/Bureau/Work/Pedigree/pedigrees.xlsx"))

names(pedigrees_all)[c(1,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20)] <- c("Pedigree","Fullsibs",
                                                                        "Matsibs","Patsibs","Matgrdmothers",
                                                                        "Matgrdfathers","Patgrdmothers","Patgrdfathers",
                                                                        "Maxdepth","Meanmatsibsize","Meanpatsibsize",
                                                                        "Meaninbreeding","Meancompleteness","MeanExpGenealdepth",
                                                                        "VarExpGenealdepth","Missingparentage")
library(corrplot)
M <- cor(pedigrees_all[,-c(1,6,7,9,11)])
tiff("corr.tiff", units ="in", width = 12, height = 12, res = 300)
corrplot(M, method = "number", type = "upper")
dev.off()
#####################

#index_f <- which(saguenay$pere %in% saguenay$ind)
#fathers <- cbind.data.frame("ind"=saguenay[-index_f,]$pere, "father"=rep(0L,length(saguenay[-index_f,]$pere)),"mother"=rep(0L,length(saguenay[-index_f,]$pere)),"sex"=rep(1L,length(saguenay[-index_f,]$pere)),"yob"=rep(0L,length(saguenay[-index_f,]$pere)))
#fathers <- fathers[!fathers$ind ==0,]
#fathers <- fathers[!duplicated(fathers$ind),]

#index_m <- which(saguenay$mere %in% saguenay$ind)
#mothers <- cbind.data.frame("ind"=saguenay[-index_m,]$mere, "father"=rep(0L,length(saguenay[-index_m,]$mere)),"mother"=rep(0L,length(saguenay[-index_m,]$mere)),"sex"=rep(2L,length(saguenay[-index_m,]$mere)),"yob"=rep(0L,length(saguenay[-index_m,]$mere)))
#mothers <- mothers[!mothers$ind ==0,]
#mothers <- mothers[!duplicated(mothers$ind),]

#ped_sag <- rbind(fathers,mothers,ped_sag)
#ped_sag <- ped_sag[!duplicated(ped_sag$ind),]
#ped_sag <- ped_sag[!is.na(ped_sag$ind),]
#save(ped_sag, file = "C:/Users/walid/OneDrive/Bureau/Work/Pedigree/ped_sag.rda")

########## first and last cohort assignment ########

#best function

cohort_parent <- function(df){
  df$firstcohort <- rep(0,nrow(df))
  df$lastcohort <- rep(0,nrow(df))
  
  for (i in 1:nrow(df)) {
    if (df$id[i] %in% df$sire) {
      yobs <- na.omit(df[df$sire == df$id[i],]$yob)
      df$firstcohort[i] <- min(yobs)
      df$lastcohort[i] <- max(yobs)
    }
    else if (df$id[i] %in% df$dam) {
      yobs <- na.omit(df[df$dam == df$id[i],]$yob)
      df$firstcohort[i] <- min(yobs)
      df$lastcohort[i] <- max(yobs)
    }
  }
  df[df == 0] <- NA
  return(df)
}

test<-cohort_parent(test)
#########founder########
ped$founder <- rep(0L, nrow(ped))
ped[is.na(ped$dam) == T & is.na(ped$datn) == T,]$founder <- 1

###################@ package PEDINF @#######################
library("pedinf")
library("tidyverse")

pedinf_dat <- ped_gaspesie[,c(1,2,3,4)]
#pedinf_dat[pedinf_dat == 0] <- NA

prep_pedinf <- function(pedigree){
  pedigree$id <- as.numeric(as.character(pedigree$id))
  pedigree$dam <- as.numeric(as.character(pedigree$dam))
  pedigree$sire <- as.numeric(as.character(pedigree$sire))
  pedigree$male <- rep(FALSE, nrow(pedigree))
  pedigree[pedigree$sex == 0,]$male <- TRUE
  return(pedigree)
}
#paternal lineages size and id

assign_paternal_lineage <- function(pedinf){
  test <- with(pedinf, load_individuals(pid = id,
                                       pid_mom = rep(0L, nrow(pedinf)),
                                       pid_dad = sire,
                                       is_male = male,
                                       error_on_pid_not_found = TRUE,
                                       progress = TRUE))
  pedigree_pedinf <- build_pedigrees_recursive(test, progress = TRUE) 

  test1 <- pedinf %>% mutate(ped_size = get_pedigree_size_from_pid(test, id)) 

  return(test1$ped_size)
}
#maternal lineages size and id
assign_maternal_lineage <- function(pedinf){
  test <- with(pedinf, load_individuals(pid = id,
                                        pid_mom = dam,
                                        pid_dad = rep(0L, nrow(pedinf_dat)),
                                        is_male = male,
                                        error_on_pid_not_found = TRUE,
                                        progress = TRUE))
  pedigree_pedinf <- build_pedigrees_recursive(test, progress = TRUE) 
  
  test1 <- pedinf %>% mutate(ped_size = get_pedigree_size_from_pid(test, id)) 
  
  return(test1$ped_size)
}
##########################assigning error rates#############################
#Pedigrees to be eliminated due to low pedigree depth (estrie, outaouais, quebec)

#input pruned pedigrees
load(file = "C:/Users/walid/OneDrive/Bureau/Work/Pedigree/data/ped_basstlaurent.rda")
load(file = "C:/Users/walid/OneDrive/Bureau/Work/Pedigree/data/ped_bois_franc.rda")
load(file = "C:/Users/walid/OneDrive/Bureau/Work/Pedigree/data/ped_charlevoix.rda")
load(file = "C:/Users/walid/OneDrive/Bureau/Work/Pedigree/data/ped_cote_beaupre.rda")
load(file = "C:/Users/walid/OneDrive/Bureau/Work/Pedigree/data/ped_gaspesie.rda")
load(file = "C:/Users/walid/OneDrive/Bureau/Work/Pedigree/data/ped_laurentides.rda")


ped_list<- list(ped_basstlaurent, ped_bois_franc,ped_charlevoix,
                ped_cote_beaupre,ped_cotesud,ped_gaspesie,
                ped_ile_montreal,ped_ilemad,ped_lanaudiere,
                ped_laurentides,ped_mauricie,
                ped_richelieu,ped_saguenay)

write.table(orderPed(ped_basstlaurent[,c(1,3,2)]), file = "ped_basstlaurent.txt", row.names = FALSE)
write.table(orderPed(ped_bois_franc[,c(1,3,2)]), file = "ped_bois_franc.txt", row.names = FALSE)
write.table(orderPed(ped_charlevoix[,c(1,3,2)]), file = "ped_charlevoix.txt", row.names = FALSE)
write.table(orderPed(ped_cote_beaupre[,c(1,3,2)]), file = "ped_cote_beaupre.txt", row.names = FALSE)
write.table(orderPed(ped_gaspesie[,c(1,3,2)]), file = "ped_gaspesie.txt", row.names = FALSE)
write.table(orderPed(ped_laurentides[,c(1,3,2)]), file = "ped_laurentides.txt", row.names = FALSE)

#(mtDNA/mother: 0.00383 [s.e.=0.00009]; Ychr/father(including EPP 0.41%): 0.0080 [s.e.=0.0007])
####
#individuals with patlineage size and matlineage size of 1 are changed to NA since it corresponds to them being the head of the line
ped_gaspesie[ped_gaspesie$patlineage_size == 1,]$patlineage_size <- NA
ped_gaspesie[ped_gaspesie$matlineage_size == 1,]$matlineage_size <- NA

#function assign error rates
assign_error_rates <- function(pedigree, fathererror = 0, mothererror = 0){
  pedigree$fatherErrorProb <- numeric(length = nrow(pedigree))
  pedigree$motherErrorProb <- numeric(length = nrow(pedigree))
  pedigree[which(is.na(pedigree$sire) == T),]$fatherErrorProb <- 1
  pedigree[which(is.na(pedigree$sire) == F),]$fatherErrorProb <- as.numeric(fathererror)
  pedigree[which(is.na(pedigree$dam) == T),]$motherErrorProb <- 1
  pedigree[which(is.na(pedigree$dam) == F),]$motherErrorProb <- as.numeric(mothererror)
  return(pedigree)
}
#experiment with different rates of genealogical errors
ped_gaspesie <- assign_error_rates(ped_gaspesie, fathererror = 0.0080, mothererror = 0.00383)

#function assign founders
assign_founder <- function(pedigree){
  pedigree$founder <- rep(0L, nrow(pedigree))
  pedigree[which(is.na(pedigree$dam) == T & is.na(pedigree$sire) == T & is.na(pedigree$yob) == T),]$founder <- 1
  return(pedigree)
}

ped_gaspesie <- assign_founder(ped_gaspesie)
#function to assign sampling status of records
sampled <- function(ped){
  ped$sampled <- as.numeric(rep(1,nrow(ped)))
  return(ped)
}

ped_gaspesie <- sampled(ped_gaspesie)
##########censoring#############
#censoring using paternal and maternal lineage size
#lineage size if theoretically set at a minimum of 10 for both sides of the family

#censor <- function(ped, lineagesize){
#  if(is.numeric(lineagesize) == F){
#    print("lineage size must be numeric")
#  }
#  for (l in lineagesize){
#    m <- as.numeric(as.character(ped[which(ped$patlineage_size == l & ped$matlineage_size == l),]$dam))
#    p <- as.numeric(as.character(ped[which(ped$patlineage_size == l & ped$matlineage_size == l),]$sire))
#    parents <- c(m,p)
#    ped <- ped[which(ped$patlineage_size != l & ped$matlineage_size != l),]
#    ped <- subset(ped, id != parents)
#  }
#  colnames <- colnames(ped)
#  ped <- fixPedigree(ped, dat = ped)
#  ped <- ped[,-c(2,3)]
#  colnames(ped)<-colnames
#  return(ped)
#}

#ped_gaspesie <- censor(ped_gaspesie,seq(2,9,1))

####
filter <- function(ped, lineagesize){
  if(is.numeric(lineagesize) == F){
    print("lineage size must be numeric")
  }
  df <- ped[which(ped$patlineage_size >= lineagesize | ped$matlineage_size >= lineagesize),c(1,2,3)]
  df <- orderPed(fixPedigree(df))
  df <- merge(df, ped[,-c(2,3)], by="id")
  return(df)
}

ped_gaspesie <- filter(ped_gaspesie, 10)
###
########################@ Pedantics @#####################
require("pedantics")
require("MCMCglmm")

gaspesiePed <- fixPedigree(ped_gaspesie[,c(1,2,3)]) #extract pedigree with id, dam and sire only

###########@ Visualizing @########
tiff("gaspesie_ped.tiff", units ="in", width = 10, height = 5, res = 300)
drawPedigree(gaspesiePed) #maternal links in red and paternal links in blue
dev.off()

#better visualisation was achieved using the program PedigreeViewer

write.table(orderPed(gaspesiePed[,c(1,3,2)]), file = "ped_gaspesie_2.txt", row.names = FALSE)#to use with pedigree viewer
##############@ Pedigree Summary @###########
pedsummary <- function(ped){
  p <- fixPedigree(ped[,c(1,2,3)])
  capture.output(pedStatSummary(pedigreeStats(p, graphicalReport = 'n', lowMem = TRUE)), file = paste(paste(deparse(substitute(ped))),"_summary.txt", sep = ""))
}

pedsummary(ped_gaspesie)

pedStatSummary(pedigreeStats(gaspesiePed, graphicalReport = 'n')) #summary of pedigree statistics

##############@ Example: Gaspesie @####################
#using value of 0 instead of NA for years works with pedantics functions

ped_gaspesie[is.na(ped_gaspesie$yob),]$yob <- 0
ped_gaspesie[is.na(ped_gaspesie$firstcohort),]$firstcohort <- 0
ped_gaspesie[is.na(ped_gaspesie$lastcohort),]$lastcohort <- 0

#pedantics function outputing "true" plausible pedigree based on supplied parameters
fullPedigree<-rpederr(ped_gaspesie[,c(1,2,3)],founders=ped_gaspesie$founder,
                      sex=ped_gaspesie$sex,cohort=ped_gaspesie$yob,
                      sireE=ped_gaspesie$fatherErrorProb,
                      damE=ped_gaspesie$motherErrorProb,first=ped_gaspesie$firstcohort,
                      last=ped_gaspesie$lastcohort,samp=ped_gaspesie$sampled)$truePedigree

fullPedigree<-orderPed(fixPedigree(fullPedigree))

########@ generating phenotypes @###########
#univariate/1st scenario

simPhen <- phensim(fullPedigree, randomA=diag(c(0.3,0.1),2), randomE=0.6,parentalA = c("d","m"), returnAllEffects = TRUE)$allEffects #generate phenotype with preset heritability
simPhen <- merge(ped, simPhen, by = "id")


#####MCMCglmm
colnames(simPhen)[1] <- "animal"
prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002)))

model <- MCMCglmm(trait_1 ~ 1,random = ~animal, rcov = ~units, family = "gaussian",data = simPhen, prior = prior, pedigree = gaspesiePed, nitt = 100000,burnin = 10000, thin = 10)

model2 <- MCMCglmm(trait_1 ~ 1,random = ~animal, rcov = ~units, family = "gaussian",data = simPhen, prior = prior, pedigree = gaspesiePed, nitt = 100000,burnin = 10000, thin = 10)

model3 <- MCMCglmm(trait_1 ~ 1,random = ~animal, rcov = ~units, family = "gaussian",data = simPhen, prior = prior, pedigree = gaspesiePed, nitt = 100000,burnin = 10000, thin = 10)

summary(model)
posterior.mode(model$VCV[,1]) #mode of additive genetic variance Va
HPDinterval(model$VCV[,1]) #confidence intervals at 95% of Va

#test with small sample size
#individuals chosen randomly from dataset
test <- merge(ped_gaspesie, simPhen, by = "id")

test_1 <- test[is.na(test$dam) == F, ]
test_1 <- test_1[sample(nrow(test_1), 144), ]

colnames(test_1)[1] <- "animal"

ped <- prunePed(gaspesiePed,test_1$animal)

#naive prior
prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002)))

model <- MCMCglmm(trait_1 ~ 1,random = ~animal,rcov = ~units, family = "gaussian",data = test_1, prior = prior, pedigree = ped, nitt = 4500000,burnin = 500000, thin = 500)

summary(model)
posterior.mode(model$VCV)
plot(model$VCV)

#informative priors
p.var <- var(test_1$trait_1, na.rm = TRUE)

prior_1 <- list(R = list(V=p.var * 0.1, nu=1), G = list(G1 = list(V=p.var * 0.9, nu=1)))

model_infprior <- MCMCglmm(trait_1 ~ 1,random = ~animal,rcov = ~units, family = "gaussian",data = test_1, prior = prior_1, pedigree = ped, nitt = 4500000,burnin = 500000, thin = 500)

summary(model_infprior)
posterior.mode(model_infprior$VCV)
plot(model_infprior$VCV)

herit <- model_infprior$VCV[,1]/(model_infprior$VCV[,1]+model_infprior$VCV[,2])
posterior.mode(herit)

#######
#random individuals from later cohort
test_1 <- test[is.na(test$dam) == F, ]

test_1 <- test_1[which(test_1$yob %in% 1830:1849), ]

test_1 <- test_1[sample(nrow(test_1), 144), ]

colnames(test_1)[1] <- "animal"

ped <- prunePed(gaspesiePed,test_1$animal)
#naive prior
prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002)))

model <- MCMCglmm(trait_1 ~ 1,random = ~animal,rcov = ~units, family = "gaussian",data = test_1, prior = prior, pedigree = ped, nitt = 4500000,burnin = 500000, thin = 500)

summary(model)
posterior.mode(model$VCV)
plot(model$VCV)

#informative priors
p.var <- var(test_1$trait_1, na.rm = TRUE)

prior_1 <- list(R = list(V=p.var * 0.9, nu=1), G = list(G1 = list(V=p.var * 0.1, nu=1)))

model_infprior <- MCMCglmm(trait_1 ~ 1,random = ~animal,rcov = ~units, family = "gaussian",data = test_1, prior = prior_1, pedigree = ped, nitt = 4500000,burnin = 500000, thin = 500)

summary(model_infprior)
posterior.mode(model_infprior$VCV)
plot(model_infprior$VCV)
#######
#univariate complex architecture/2nd scenario
randomA<-diag(c(0.3,0.1),2)
randomE<-0.6
parentalA<-c("d","m")
parentalE<-c("d")


simPhen <- phensim(fullPedigree, randomA=randomA, randomE=randomE,
                   parentalA=parentalA)$phenotypes

colnames(simPhen)[1] <- "animal"

#MCMCglmm
prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002),G2 = list(V=1, nu=0.002)))

model <- MCMCglmm(trait_1 ~ 1,random = ~animal + dam,family = "gaussian",data = simPhen, prior = prior, pedigree = gaspesiePed, nitt = 450000,burnin = 50000, thin = 50)

#####two correlated traits
#genetic correlation
r <- x/sqrt(a*b) #formula for genetic correlation
Va1<-0.3 #genetic variance of trait 1
Va2<-0.1 #genetic variance of trait 2; fixed at 0.1
Ve1<-1-Va1 #environmental variance of trait 1
Ve2<-1-Va2 #environmental variance of trait 2
rg<-0.5 #expected genetic correlation between trait 1 and trait 2
re<-0.3 #expected environmental correlation between trait 1 and trait 2; fixed at 0.3
xg <- rg*sqrt(Va1*Va2) #formula for genetic covariance
xe <- re*sqrt(Ve1*Ve2) #formula for environmental covariance

randomA<-matrix(c(Va1,xg,xg,Va2),2)
randomE<-matrix(c(Ve1,xe,xe,Ve2),2)
parentalA<-c("d","d")
parentalE<-c("d","m")


simPhen <- phensim(fullPedigree, traits =2, randomA=randomA, 
                   randomE=randomE)$phenotypes

simPhen <- merge(ped_gaspesie, simPhen, by = "id")
dat <- subset(simPhen, yob != 0 & firstcohort != 0)
ped <- orderPed(prunePed(ped_gaspesie[,c(1,2,3)], dat$id))


##############################
## we’ll do plenty of replicates of each scenario:
n<-500

## here is a decent range of heritabilities to simulate:
heritabilities<-seq(0,0.4,by=0.02)
## we will need to keep track of the p-values resulting
## from each replicate of each scenario, and for each
## scenario, we’ll want to keep track of the proportion
## of significant tests:
p_values<-matrix(NA,n,length(heritabilities))
power<-array(dim=length(heritabilities))
## here we’ll loop through all the scenarios and replicates,
## keeping track of the p-value at the end of each iteration
## through the work flow, and stopping to calculate the
## proportion of significant tests each time we get through
## all the replicates for each heritability
for(x in 1:length(heritabilities)) {
  for(y in 1:n) {
    simPhen<-phensim(fullPedigree,randomA=heritabilities[x],
                     randomE=1-heritabilities[x])$phenotypes
    simPhen$mum<-ped_gaspesie$dam[match(simPhen$id,ped_gaspesie$id)]
    simPhen$mumPhenotype<-simPhen$trait_1[match(simPhen$mum,
                                                as.numeric(as.character(simPhen$id)))]
    p_values[y,x]<-anova(lm(simPhen$trait_1~simPhen$mumPhenotype))[1,5]
  }
  power[x]<-table(p_values[,x]<0.05)["TRUE"]/n
}
power_results<-cbind.data.frame(heritabilities,power)

###########interpreting models#########
#function to fix colnames
colnames_fun <- function(dat){
  colnames(dat)<- paste("sim",seq(1:5), sep = "")
  return(dat)
}
##1st scenario, simple architecture with Va= 0.1 & 0.3
#Va=0.1
#main mcmc dataframes
mcmc_gaspesie <- colnames_fun(readRDS(file = "model_1stscenario/mcmc_sim_ped_gaspesie_0.1.rds"))
mcmc_laurentides <- colnames_fun(readRDS(file = "model_1stscenario/mcmc_sim_ped_laurentides_0.1.rds"))
mcmc_basstlaurent <- colnames_fun(readRDS(file = "model_1stscenario/mcmc_sim_ped_basstlaurent_0.1.rds"))
mcmc_bois_franc <- colnames_fun(readRDS(file = "model_1stscenario/mcmc_sim_ped_bois_franc_0.1.rds"))
mcmc_cote_beaupre <- colnames_fun(readRDS(file = "model_1stscenario/mcmc_sim_ped_cote_beaupre_0.1.rds"))
mcmc_charlevoix <- colnames_fun(readRDS(file = "model_1stscenario/mcmc_sim_ped_charlevoix_0.1.rds"))

#do not use if using all iterations
#mode
mcmc_gaspesie_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie)))
mcmc_laurentides_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides)))
mcmc_basstlaurent_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent)))
mcmc_bois_franc_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc)))
mcmc_cote_beaupre_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre)))
mcmc_charlevoix_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix)))
#mean
mcmc_gaspesie_mean <- as.vector(colMeans(as.mcmc(mcmc_gaspesie)))
mcmc_laurentides_mean <- as.vector(colMeans(as.mcmc(mcmc_laurentides)))
mcmc_basstlaurent_mean <- as.vector(colMeans(as.mcmc(mcmc_basstlaurent)))
mcmc_bois_franc_mean <- as.vector(colMeans(as.mcmc(mcmc_bois_franc)))
mcmc_cote_beaupre_mean <- as.vector(colMeans(as.mcmc(mcmc_cote_beaupre)))
mcmc_charlevoix_mean <- as.vector(colMeans(as.mcmc(mcmc_charlevoix)))
#Va=0.3
#main mcmc dataframes
mcmc_gaspesie_2 <- colnames_fun(readRDS(file = "model_1stscenario/mcmc_sim_ped_gaspesie_0.3.rds"))
mcmc_laurentides_2 <- colnames_fun(readRDS(file = "model_1stscenario/mcmc_sim_ped_laurentides_0.3.rds"))
mcmc_basstlaurent_2 <- colnames_fun(readRDS(file = "model_1stscenario/mcmc_sim_ped_basstlaurent_0.3.rds"))
mcmc_bois_franc_2 <- colnames_fun(readRDS(file = "model_1stscenario/mcmc_sim_ped_bois_franc_0.3.rds"))
mcmc_cote_beaupre_2 <- colnames_fun(readRDS(file = "model_1stscenario/mcmc_sim_ped_cote_beaupre_0.3.rds"))
mcmc_charlevoix_2 <- colnames_fun(readRDS(file = "model_1stscenario/mcmc_sim_ped_charlevoix_0.3.rds"))

#do not use if using all iterations
#mode
mcmc_gaspesie_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_2)))
mcmc_laurentides_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_2)))
mcmc_basstlaurent_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_2)))
mcmc_bois_franc_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_2)))
mcmc_cote_beaupre_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_2)))
mcmc_charlevoix_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_2)))
#mean
mcmc_gaspesie_2_mean <- as.vector(colMeans(as.mcmc(mcmc_gaspesie_2)))
mcmc_laurentides_2_mean <- as.vector(colMeans(as.mcmc(mcmc_laurentides_2)))
mcmc_basstlaurent_2_mean <- as.vector(colMeans(as.mcmc(mcmc_basstlaurent_2)))
mcmc_bois_franc_2_mean <- as.vector(colMeans(as.mcmc(mcmc_bois_franc_2)))
mcmc_cote_beaupre_2_mean <- as.vector(colMeans(as.mcmc(mcmc_cote_beaupre_2)))
mcmc_charlevoix_2_mean <- as.vector(colMeans(as.mcmc(mcmc_charlevoix_2)))

#if using all iterations from simulations
df<-data.frame("Va" = as.factor(c(rep(0.1,6000),rep(0.3,6000))),
               "region" = c(rep("Gaspésie",1000),rep("Bas-St-Laurent",1000),rep("Laurentides",1000),
                            rep("Charlevoix",1000),rep("Côte-de-Beaupré",1000),rep("Bois-Francs",1000),
                            rep("Gaspésie",1000),rep("Laurentides",1000),rep("Bas-St-Laurent",1000),
                            rep("Charlevoix",1000),rep("Côte-de-Beaupré",1000),rep("Bois-Francs",1000)),
               "additive" = c(mcmc_gaspesie,mcmc_basstlaurent,mcmc_laurentides,mcmc_charlevoix,mcmc_cote_beaupre,mcmc_bois_franc,
                              mcmc_gaspesie_2,mcmc_basstlaurent_2,mcmc_laurentides_2,mcmc_charlevoix_2,mcmc_cote_beaupre_2,mcmc_bois_franc_2))

df$region<-as.factor(df$region)
df$region <- factor(df$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix","Côte-de-Beaupré","Bois-Franc"))

#if using just posterior mode from simulations
df<-data.frame("Va" = as.factor(c(rep(0.1,30),rep(0.3,30))),
               "region" = c(rep("Gaspésie",5),rep("Bas-St-Laurent",5),rep("Laurentides",5),
                            rep("Charlevoix",5),rep("Côte-de-Beaupré",5),rep("Bois-Francs",5),
                            rep("Gaspésie",5),rep("Laurentides",5),rep("Bas-St-Laurent",5),
                            rep("Charlevoix",5),rep("Côte-de-Beaupré",5),rep("Bois-Francs",5)),
               "additive" = c(mcmc_gaspesie_mode,mcmc_basstlaurent_mode,mcmc_laurentides_mode,mcmc_charlevoix_mode,mcmc_cote_beaupre_mode,
                              mcmc_bois_franc_mode,mcmc_gaspesie_2_mode,mcmc_basstlaurent_2_mode,mcmc_laurentides_2_mode,mcmc_charlevoix_2_mode,
                              mcmc_cote_beaupre_2_mode,mcmc_bois_franc_2_mode))

df$region<-as.factor(df$region)
df$region <- factor(df$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"))

#calculating precision based on posterior mode

prec_0.1_mcmc <- sqrt((df[which(df$Va == 0.1),]$additive - 0.1)^2)
prec_0.3_mcmc <- sqrt((df[which(df$Va == 0.3),]$additive - 0.3)^2)

n <- 1
m <- n+4
prec_0.1 <- 0
for (x in 1:6) {
  prec_0.1[x] <- mean(prec_0.1_mcmc[n:m])
  n <- n+5
  m <- m+5
}
prec_0.1 <- round(prec_0.1,digits = 3)

n <- 1
m <- n+4
prec_0.3 <- 0
for (x in 1:6) {
  prec_0.3[x] <- mean(prec_0.3_mcmc[n:m])
  n <- n+5
  m <- m+5
}
prec_0.3 <- round(prec_0.3,digits = 3)

########summary
#Va=0.1
mcmc_gaspesie_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_gaspesie))))
mcmc_laurentides_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_laurentides))))
mcmc_basstlaurent_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent))))
mcmc_bois_franc_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_bois_franc))))
mcmc_cote_beaupre_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre))))
mcmc_charlevoix_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_charlevoix))))

mcmc_gaspesie_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie)))[1:5])
mcmc_laurentides_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides)))[1:5])
mcmc_basstlaurent_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent)))[1:5])
mcmc_bois_franc_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc)))[1:5])
mcmc_cote_beaupre_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre)))[1:5])
mcmc_charlevoix_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix)))[1:5])

mcmc_gaspesie_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie)))[6:10])
mcmc_laurentides_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides)))[6:10])
mcmc_basstlaurent_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent)))[6:10])
mcmc_bois_franc_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc)))[6:10])
mcmc_cote_beaupre_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre)))[6:10])
mcmc_charlevoix_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix)))[6:10])
#Va=0.3
mcmc_gaspesie_2_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_2))))
mcmc_laurentides_2_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_laurentides_2))))
mcmc_basstlaurent_2_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_2))))
mcmc_bois_franc_2_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_2))))
mcmc_cote_beaupre_2_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_2))))
mcmc_charlevoix_2_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_2))))

mcmc_gaspesie_2_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_2)))[1:5])
mcmc_laurentides_2_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_2)))[1:5])
mcmc_basstlaurent_2_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_2)))[1:5])
mcmc_bois_franc_2_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_2)))[1:5])
mcmc_cote_beaupre_2_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_2)))[1:5])
mcmc_charlevoix_2_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_2)))[1:5])

mcmc_gaspesie_2_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_2)))[6:10])
mcmc_laurentides_2_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_2)))[6:10])
mcmc_basstlaurent_2_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_2)))[6:10])
mcmc_bois_franc_2_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_2)))[6:10])
mcmc_cote_beaupre_2_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_2)))[6:10])
mcmc_charlevoix_2_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_2)))[6:10])
#dataframe
df_mcmc <- data.frame(region = c("Gaspésie", "Bas-St-Laurent","Laurentides",
                                 "Charlevoix", "Côte-de-Beaupré", "Bois-Francs",
                                 "Gaspésie", "Bas-St-Laurent","Laurentides",
                                 "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"),
                               Estimate = c(mcmc_gaspesie_mean,
                                            mcmc_basstlaurent_mean,
                                            mcmc_laurentides_mean,
                                            mcmc_charlevoix_mean,
                                            mcmc_cote_beaupre_mean,
                                            mcmc_bois_franc_mean,
                                            mcmc_gaspesie_2_mean,
                                            mcmc_basstlaurent_2_mean,
                                            mcmc_laurentides_2_mean,
                                            mcmc_charlevoix_2_mean,
                                            mcmc_cote_beaupre_2_mean,
                                            mcmc_bois_franc_2_mean),
                               Lower = c(mcmc_gaspesie_lower,
                                         mcmc_basstlaurent_lower,
                                         mcmc_laurentides_lower,
                                         mcmc_charlevoix_lower,
                                         mcmc_cote_beaupre_lower,
                                         mcmc_bois_franc_lower,
                                         mcmc_gaspesie_2_lower,
                                         mcmc_basstlaurent_2_lower,
                                         mcmc_laurentides_2_lower,
                                         mcmc_charlevoix_2_lower,
                                         mcmc_cote_beaupre_2_lower,
                                         mcmc_bois_franc_2_lower),
                               Upper = c(mcmc_gaspesie_upper,
                                         mcmc_basstlaurent_upper,
                                         mcmc_laurentides_upper,
                                         mcmc_charlevoix_upper,
                                         mcmc_cote_beaupre_upper,
                                         mcmc_bois_franc_upper,
                                         mcmc_gaspesie_2_upper,
                                         mcmc_basstlaurent_2_upper,
                                         mcmc_laurentides_2_upper,
                                         mcmc_charlevoix_2_upper,
                                         mcmc_cote_beaupre_2_upper,
                                         mcmc_bois_franc_2_upper),
                                Va = as.factor(c(rep(0.1,6),rep(0.3,6))))

df_mcmc$region<-as.factor(df_mcmc$region)
df_mcmc$region <- factor(df_mcmc$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix","Côte-de-Beaupré","Bois-Francs"))

#######
#####2nd scenario, complex architecture including maternal effects (Vm=0.1) with Va= 0.1 & 0.3
#Va=0.1 
mcmc_gaspesie_complex <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_gaspesie_diag(c(0.1, 0.1), 2).rds"))
mcmc_laurentides_complex <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_laurentides_diag(c(0.1, 0.1), 2).rds"))
mcmc_basstlaurent_complex <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_basstlaurent_diag(c(0.1, 0.1), 2).rds"))
mcmc_bois_franc_complex <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_bois_franc_diag(c(0.1, 0.1), 2).rds"))
mcmc_cote_beaupre_complex <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_cote_beaupre_diag(c(0.1, 0.1), 2).rds"))
mcmc_charlevoix_complex <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_charlevoix_diag(c(0.1, 0.1), 2).rds"))
#& Vm=0.1
mcmc_gaspesie_complex_dam <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_gaspesie_diag(c(0.1, 0.1), 2).rds"))
mcmc_laurentides_complex_dam <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_laurentides_diag(c(0.1, 0.1), 2).rds"))
mcmc_basstlaurent_complex_dam <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_basstlaurent_diag(c(0.1, 0.1), 2).rds"))
mcmc_bois_franc_complex_dam <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_bois_franc_diag(c(0.1, 0.1), 2).rds"))
mcmc_cote_beaupre_complex_dam <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_cote_beaupre_diag(c(0.1, 0.1), 2).rds"))
mcmc_charlevoix_complex_dam <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_charlevoix_diag(c(0.1, 0.1), 2).rds"))

#do not use if using all iterations
#VA
mcmc_gaspesie_complex_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex)))
mcmc_laurentides_complex_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex)))
mcmc_basstlaurent_complex_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex)))
mcmc_bois_franc_complex_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex)))
mcmc_cote_beaupre_complex_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex)))
mcmc_charlevoix_complex_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex)))
#VM
mcmc_gaspesie_complex_dam_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex_dam)))
mcmc_laurentides_complex_dam_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex_dam)))
mcmc_basstlaurent_complex_dam_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex_dam)))
mcmc_bois_franc_complex_dam_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex_dam)))
mcmc_cote_beaupre_complex_dam_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex_dam)))
mcmc_charlevoix_complex_dam_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex_dam)))
#VA
mcmc_gaspesie_complex_mean <- as.vector(colMeans(as.mcmc(mcmc_gaspesie_complex)))
mcmc_laurentides_complex_mean <- as.vector(colMeans(as.mcmc(mcmc_laurentides_complex)))
mcmc_basstlaurent_complex_mean <- as.vector(colMeans(as.mcmc(mcmc_basstlaurent_complex)))
mcmc_bois_franc_complex_mean <- as.vector(colMeans(as.mcmc(mcmc_bois_franc_complex)))
mcmc_cote_beaupre_complex_mean <- as.vector(colMeans(as.mcmc(mcmc_cote_beaupre_complex)))
mcmc_charlevoix_complex_mean <- as.vector(colMeans(as.mcmc(mcmc_charlevoix_complex)))
#VM
mcmc_gaspesie_complex_dam_mean <- as.vector(colMeans(as.mcmc(mcmc_gaspesie_complex_dam)))
mcmc_laurentides_complex_dam_mean <- as.vector(colMeans(as.mcmc(mcmc_laurentides_complex_dam)))
mcmc_basstlaurent_complex_dam_mean <- as.vector(colMeans(as.mcmc(mcmc_basstlaurent_complex_dam)))
mcmc_bois_franc_complex_dam_mean <- as.vector(colMeans(as.mcmc(mcmc_bois_franc_complex_dam)))
mcmc_cote_beaupre_complex_dam_mean <- as.vector(colMeans(as.mcmc(mcmc_cote_beaupre_complex_dam)))
mcmc_charlevoix_complex_dam_mean <- as.vector(colMeans(as.mcmc(mcmc_charlevoix_complex_dam)))

#Va=0.3 
mcmc_gaspesie_complex_2 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_gaspesie_diag(c(0.3, 0.1), 2).rds"))
mcmc_laurentides_complex_2 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_laurentides_diag(c(0.3, 0.1), 2).rds"))
mcmc_basstlaurent_complex_2 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_basstlaurent_diag(c(0.3, 0.1), 2).rds"))
mcmc_bois_franc_complex_2 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_bois_franc_diag(c(0.3, 0.1), 2).rds"))
mcmc_cote_beaupre_complex_2 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_cote_beaupre_diag(c(0.3, 0.1), 2).rds"))
mcmc_charlevoix_complex_2 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_charlevoix_diag(c(0.3, 0.1), 2).rds"))
#& Vm=0.1
mcmc_gaspesie_complex_2_dam <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_gaspesie_diag(c(0.3, 0.1), 2).rds"))
mcmc_laurentides_complex_2_dam <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_laurentides_diag(c(0.3, 0.1), 2).rds"))
mcmc_basstlaurent_complex_2_dam <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_basstlaurent_diag(c(0.3, 0.1), 2).rds"))
mcmc_bois_franc_complex_2_dam <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_bois_franc_diag(c(0.3, 0.1), 2).rds"))
mcmc_cote_beaupre_complex_2_dam <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_cote_beaupre_diag(c(0.3, 0.1), 2).rds"))
mcmc_charlevoix_complex_2_dam <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_charlevoix_diag(c(0.3, 0.1), 2).rds"))

#do not use if using all iterations
#VA
mcmc_gaspesie_complex_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex_2)))
mcmc_laurentides_complex_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex_2)))
mcmc_basstlaurent_complex_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex_2)))
mcmc_bois_franc_complex_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex_2)))
mcmc_cote_beaupre_complex_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex_2)))
mcmc_charlevoix_complex_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex_2)))
#VM
mcmc_gaspesie_complex_2_dam_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex_2_dam)))
mcmc_laurentides_complex_2_dam_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex_2_dam)))
mcmc_basstlaurent_complex_2_dam_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex_2_dam)))
mcmc_bois_franc_complex_2_dam_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex_2_dam)))
mcmc_cote_beaupre_complex_2_dam_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex_2_dam)))
mcmc_charlevoix_complex_2_dam_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex_2_dam)))
#VA
mcmc_gaspesie_complex_2_mean <- as.vector(colMeans(as.mcmc(mcmc_gaspesie_complex_2)))
mcmc_laurentides_complex_2_mean <- as.vector(colMeans(as.mcmc(mcmc_laurentides_complex_2)))
mcmc_basstlaurent_complex_2_mean <- as.vector(colMeans(as.mcmc(mcmc_basstlaurent_complex_2)))
mcmc_bois_franc_complex_2_mean <- as.vector(colMeans(as.mcmc(mcmc_bois_franc_complex_2)))
mcmc_cote_beaupre_complex_2_mean <- as.vector(colMeans(as.mcmc(mcmc_cote_beaupre_complex_2)))
mcmc_charlevoix_complex_2_mean <- as.vector(colMeans(as.mcmc(mcmc_charlevoix_complex_2)))
#VM
mcmc_gaspesie_complex_2_dam_mean <- as.vector(colMeans(as.mcmc(mcmc_gaspesie_complex_2_dam)))
mcmc_laurentides_complex_2_dam_mean <- as.vector(colMeans(as.mcmc(mcmc_laurentides_complex_2_dam)))
mcmc_basstlaurent_complex_2_dam_mean <- as.vector(colMeans(as.mcmc(mcmc_basstlaurent_complex_2_dam)))
mcmc_bois_franc_complex_2_dam_mean <- as.vector(colMeans(as.mcmc(mcmc_bois_franc_complex_2_dam)))
mcmc_cote_beaupre_complex_2_dam_mean <- as.vector(colMeans(as.mcmc(mcmc_cote_beaupre_complex_2_dam)))
mcmc_charlevoix_complex_2_dam_mean <- as.vector(colMeans(as.mcmc(mcmc_charlevoix_complex_2_dam)))

#if using all iterations from simulations
df_complex<-data.frame("Va" = as.factor(c(rep(0.1,6000),rep(0.3,6000))),
               "region" = c(rep("Gaspésie",1000),rep("Bas-St-Laurent",1000),rep("Laurentides",1000),
                            rep("Charlevoix",1000),rep("Côte-de-Beaupré",1000),rep("Bois-Francs",1000),
                            rep("Gaspésie",1000),rep("Laurentides",1000),rep("Bas-St-Laurent",1000),
                            rep("Charlevoix",1000),rep("Côte-de-Beaupré",1000),rep("Bois-Francs",1000)),
               "additive" = c(mcmc_gaspesie_complex,mcmc_basstlaurent_complex,mcmc_laurentides_complex,mcmc_charlevoix_complex,mcmc_cote_beaupre_complex,mcmc_bois_franc_complex,
                              mcmc_gaspesie_complex_2,mcmc_basstlaurent_complex_2,mcmc_laurentides_complex_2,mcmc_charlevoix_complex_2,mcmc_cote_beaupre_complex_2,mcmc_bois_franc_complex_2))

df$region<-as.factor(df$region)
df$region <- factor(df$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix","Côte-de-Beaupré","Bois-Franc"))

#if using just posterior mode from simulations
df_complex<-data.frame("Va" = as.factor(c(rep(0.1,30),rep(0.3,30))),
               "Vm" = as.factor(c(rep(0.1,60))),
               "region" = c(rep("Gaspésie",5),rep("Bas-St-Laurent",5),rep("Laurentides",5),
                            rep("Charlevoix",5),rep("Côte-de-Beaupré",5),rep("Bois-Francs",5),
                            rep("Gaspésie",5),rep("Laurentides",5),rep("Bas-St-Laurent",5),
                            rep("Charlevoix",5),rep("Côte-de-Beaupré",5),rep("Bois-Francs",5)),
               "additive" = c(mcmc_gaspesie_complex_mode,mcmc_basstlaurent_complex_mode,mcmc_laurentides_complex_mode,mcmc_charlevoix_complex_mode,
                              mcmc_cote_beaupre_complex_mode,mcmc_bois_franc_complex_mode,mcmc_gaspesie_complex_2_mode,mcmc_basstlaurent_complex_2_mode,
                              mcmc_laurentides_complex_2_mode,mcmc_charlevoix_complex_2_mode,mcmc_cote_beaupre_complex_2_mode,mcmc_bois_franc_complex_2_mode),
               "maternal" = c(mcmc_gaspesie_complex_dam_mode,mcmc_basstlaurent_complex_dam_mode,mcmc_laurentides_complex_dam_mode,mcmc_charlevoix_complex_dam_mode,
                              mcmc_cote_beaupre_complex_dam_mode,mcmc_bois_franc_complex_dam_mode,mcmc_gaspesie_complex_2_dam_mode,mcmc_basstlaurent_complex_2_dam_mode,
                              mcmc_laurentides_complex_2_dam_mode,mcmc_charlevoix_complex_2_dam_mode,mcmc_cote_beaupre_complex_2_dam_mode,mcmc_bois_franc_complex_2_dam_mode))

df_complex$region<-as.factor(df_complex$region)
df_complex$region <- factor(df_complex$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"))

prec_0.1_mcmc_complex <- sqrt((df_complex[which(df_complex$Va == 0.1),]$additive - 0.1)^2)
prec_0.3_mcmc_complex <- sqrt((df_complex[which(df_complex$Va == 0.3),]$additive - 0.3)^2)

n <- 1
m <- n+4
prec_0.1_complex <- 0
for (x in 1:6) {
  prec_0.1_complex[x] <- mean(prec_0.1_mcmc_complex[n:m])
  n <- n+5
  m <- m+5
}
prec_0.1_complex <- round(prec_0.1_complex,digits = 3)

n <- 1
m <- n+4
prec_0.3_complex <- 0
for (x in 1:6) {
  prec_0.3_complex[x] <- mean(prec_0.3_mcmc_complex[n:m])
  n <- n+5
  m <- m+5
}
prec_0.3_complex <- round(prec_0.3_complex,digits = 3)

########summary
#Va=0.1
mcmc_gaspesie_complex_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex))))
mcmc_laurentides_complex_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex))))
mcmc_basstlaurent_complex_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex))))
mcmc_bois_franc_complex_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex))))
mcmc_cote_beaupre_complex_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex))))
mcmc_charlevoix_complex_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex))))

mcmc_gaspesie_complex_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex)))[1:5])
mcmc_laurentides_complex_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex)))[1:5])
mcmc_basstlaurent_complex_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex)))[1:5])
mcmc_bois_franc_complex_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex)))[1:5])
mcmc_cote_beaupre_complex_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex)))[1:5])
mcmc_charlevoix_complex_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex)))[1:5])

mcmc_gaspesie_complex_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex)))[6:10])
mcmc_laurentides_complex_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex)))[6:10])
mcmc_basstlaurent_complex_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex)))[6:10])
mcmc_bois_franc_complex_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex)))[6:10])
mcmc_cote_beaupre_complex_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex)))[6:10])
mcmc_charlevoix_complex_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex)))[6:10])
#Vm=0.1
mcmc_gaspesie_complex_mean_dam <- mean(as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex_dam))))
mcmc_laurentides_complex_mean_dam <- mean(as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex_dam))))
mcmc_basstlaurent_complex_mean_dam <- mean(as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex_dam))))
mcmc_bois_franc_complex_mean_dam <- mean(as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex_dam))))
mcmc_cote_beaupre_complex_mean_dam <- mean(as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex_dam))))
mcmc_charlevoix_complex_mean_dam <- mean(as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex_dam))))

mcmc_gaspesie_complex_lower_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex_dam)))[1:5])
mcmc_laurentides_complex_lower_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex_dam)))[1:5])
mcmc_basstlaurent_complex_lower_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex_dam)))[1:5])
mcmc_bois_franc_complex_lower_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex_dam)))[1:5])
mcmc_cote_beaupre_complex_lower_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex_dam)))[1:5])
mcmc_charlevoix_complex_lower_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex_dam)))[1:5])

mcmc_gaspesie_complex_upper_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex_dam)))[6:10])
mcmc_laurentides_complex_upper_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex_dam)))[6:10])
mcmc_basstlaurent_complex_upper_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex_dam)))[6:10])
mcmc_bois_franc_complex_upper_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex_dam)))[6:10])
mcmc_cote_beaupre_complex_upper_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex_dam)))[6:10])
mcmc_charlevoix_complex_upper_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex_dam)))[6:10])
#Va=0.3
mcmc_gaspesie_complex_2_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex_2))))
mcmc_laurentides_complex_2_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex_2))))
mcmc_basstlaurent_complex_2_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex_2))))
mcmc_bois_franc_complex_2_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex_2))))
mcmc_cote_beaupre_complex_2_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex_2))))
mcmc_charlevoix_complex_2_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex_2))))

mcmc_gaspesie_complex_2_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex_2)))[1:5])
mcmc_laurentides_complex_2_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex_2)))[1:5])
mcmc_basstlaurent_complex_2_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex_2)))[1:5])
mcmc_bois_franc_complex_2_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex_2)))[1:5])
mcmc_cote_beaupre_complex_2_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex_2)))[1:5])
mcmc_charlevoix_complex_2_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex_2)))[1:5])

mcmc_gaspesie_complex_2_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex_2)))[6:10])
mcmc_laurentides_complex_2_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex_2)))[6:10])
mcmc_basstlaurent_complex_2_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex_2)))[6:10])
mcmc_bois_franc_complex_2_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex_2)))[6:10])
mcmc_cote_beaupre_complex_2_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex_2)))[6:10])
mcmc_charlevoix_complex_2_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex_2)))[6:10])
#Vm=0.1
mcmc_gaspesie_complex_2_mean_dam <- mean(as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex_2_dam))))
mcmc_laurentides_complex_2_mean_dam <- mean(as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex_2_dam))))
mcmc_basstlaurent_complex_2_mean_dam <- mean(as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex_2_dam))))
mcmc_bois_franc_complex_2_mean_dam <- mean(as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex_2_dam))))
mcmc_cote_beaupre_complex_2_mean_dam <- mean(as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex_2_dam))))
mcmc_charlevoix_complex_2_mean_dam <- mean(as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex_2_dam))))

mcmc_gaspesie_complex_2_lower_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex_2_dam)))[1:5])
mcmc_laurentides_complex_2_lower_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex_2_dam)))[1:5])
mcmc_basstlaurent_complex_2_lower_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex_2_dam)))[1:5])
mcmc_bois_franc_complex_2_lower_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex_2_dam)))[1:5])
mcmc_cote_beaupre_complex_2_lower_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex_2_dam)))[1:5])
mcmc_charlevoix_complex_2_lower_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex_2_dam)))[1:5])

mcmc_gaspesie_complex_2_upper_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex_2_dam)))[6:10])
mcmc_laurentides_complex_2_upper_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex_2_dam)))[6:10])
mcmc_basstlaurent_complex_2_upper_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex_2_dam)))[6:10])
mcmc_bois_franc_complex_2_upper_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex_2_dam)))[6:10])
mcmc_cote_beaupre_complex_2_upper_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex_2_dam)))[6:10])
mcmc_charlevoix_complex_2_upper_dam <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex_2_dam)))[6:10])
#dataframe
df_mcmc_complex <- data.frame(region = rep(c("Gaspésie", "Bas-St-Laurent","Laurentides",
                                 "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"),4),
                      Estimate = c(mcmc_gaspesie_complex_mean,
                                   mcmc_basstlaurent_complex_mean,
                                   mcmc_laurentides_complex_mean,
                                   mcmc_charlevoix_complex_mean,
                                   mcmc_cote_beaupre_complex_mean,
                                   mcmc_bois_franc_complex_mean,
                                   mcmc_gaspesie_complex_2_mean,
                                   mcmc_basstlaurent_complex_2_mean,
                                   mcmc_laurentides_complex_2_mean,
                                   mcmc_charlevoix_complex_2_mean,
                                   mcmc_cote_beaupre_complex_2_mean,
                                   mcmc_bois_franc_complex_2_mean,
                                   mcmc_gaspesie_complex_mean_dam,
                                   mcmc_basstlaurent_complex_mean_dam,
                                   mcmc_laurentides_complex_mean_dam,
                                   mcmc_charlevoix_complex_mean_dam,
                                   mcmc_cote_beaupre_complex_mean_dam,
                                   mcmc_bois_franc_complex_mean_dam,
                                   mcmc_gaspesie_complex_2_mean_dam,
                                   mcmc_basstlaurent_complex_2_mean_dam,
                                   mcmc_laurentides_complex_2_mean_dam,
                                   mcmc_charlevoix_complex_2_mean_dam,
                                   mcmc_cote_beaupre_complex_2_mean_dam,
                                   mcmc_bois_franc_complex_2_mean_dam),
                      Lower = c(mcmc_gaspesie_complex_lower,
                                mcmc_basstlaurent_complex_lower,
                                mcmc_laurentides_complex_lower,
                                mcmc_charlevoix_complex_lower,
                                mcmc_cote_beaupre_complex_lower,
                                mcmc_bois_franc_complex_lower,
                                mcmc_gaspesie_complex_2_lower,
                                mcmc_basstlaurent_complex_2_lower,
                                mcmc_laurentides_complex_2_lower,
                                mcmc_charlevoix_complex_2_lower,
                                mcmc_cote_beaupre_complex_2_lower,
                                mcmc_bois_franc_complex_2_lower,
                                mcmc_gaspesie_complex_lower_dam,
                                mcmc_basstlaurent_complex_lower_dam,
                                mcmc_laurentides_complex_lower_dam,
                                mcmc_charlevoix_complex_lower_dam,
                                mcmc_cote_beaupre_complex_lower_dam,
                                mcmc_bois_franc_complex_lower_dam,
                                mcmc_gaspesie_complex_2_lower_dam,
                                mcmc_basstlaurent_complex_2_lower_dam,
                                mcmc_laurentides_complex_2_lower_dam,
                                mcmc_charlevoix_complex_2_lower_dam,
                                mcmc_cote_beaupre_complex_2_lower_dam,
                                mcmc_bois_franc_complex_2_lower_dam),
                      Upper = c(mcmc_gaspesie_complex_upper,
                                mcmc_basstlaurent_complex_upper,
                                mcmc_laurentides_complex_upper,
                                mcmc_charlevoix_complex_upper,
                                mcmc_cote_beaupre_complex_upper,
                                mcmc_bois_franc_complex_upper,
                                mcmc_gaspesie_complex_2_upper,
                                mcmc_basstlaurent_complex_2_upper,
                                mcmc_laurentides_complex_2_upper,
                                mcmc_charlevoix_complex_2_upper,
                                mcmc_cote_beaupre_complex_2_upper,
                                mcmc_bois_franc_complex_2_upper,
                                mcmc_gaspesie_complex_upper_dam,
                                mcmc_basstlaurent_complex_upper_dam,
                                mcmc_laurentides_complex_upper_dam,
                                mcmc_charlevoix_complex_upper_dam,
                                mcmc_cote_beaupre_complex_upper_dam,
                                mcmc_bois_franc_complex_upper_dam,
                                mcmc_gaspesie_complex_2_upper_dam,
                                mcmc_basstlaurent_complex_2_upper_dam,
                                mcmc_laurentides_complex_2_upper_dam,
                                mcmc_charlevoix_complex_2_upper_dam,
                                mcmc_cote_beaupre_complex_2_upper_dam,
                                mcmc_bois_franc_complex_2_upper_dam),
                      Va = as.factor(rep(c(rep(0.1,6),rep(0.3,6)),2)),
                      Var = as.factor(rep(c(rep("Additive",12),rep("Maternal",12)))))

df_mcmc_complex$region<-as.factor(df_mcmc_complex$region)
df_mcmc_complex$region <- factor(df_mcmc_complex$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix","Côte-de-Beaupré","Bois-Francs"))

#Va=0.1 with increased Vm
#Va=0.1 
mcmc_gaspesie_complex_03 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_gaspesie_diag(c(0.1, 0.3), 2).rds"))
mcmc_laurentides_complex_03 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_laurentides_diag(c(0.1, 0.3), 2).rds"))
mcmc_basstlaurent_complex_03 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_basstlaurent_diag(c(0.1, 0.3), 2).rds"))
mcmc_bois_franc_complex_03 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_bois_franc_diag(c(0.1, 0.3), 2).rds"))
mcmc_cote_beaupre_complex_03 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_cote_beaupre_diag(c(0.1, 0.3), 2).rds"))
mcmc_charlevoix_complex_03 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_charlevoix_diag(c(0.1, 0.3), 2).rds"))
#& Vm=0.3
mcmc_gaspesie_complex_dam_03 <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_gaspesie_diag(c(0.1, 0.3), 2).rds"))
mcmc_laurentides_complex_dam_03 <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_laurentides_diag(c(0.1, 0.3), 2).rds"))
mcmc_basstlaurent_complex_dam_03 <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_basstlaurent_diag(c(0.1, 0.3), 2).rds"))
mcmc_bois_franc_complex_dam_03 <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_bois_franc_diag(c(0.1, 0.3), 2).rds"))
mcmc_cote_beaupre_complex_dam_03 <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_cote_beaupre_diag(c(0.1, 0.3), 2).rds"))
mcmc_charlevoix_complex_dam_03 <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_charlevoix_diag(c(0.1, 0.3), 2).rds"))

#VA
mcmc_gaspesie_complex_mode_03 <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex_03)))
mcmc_laurentides_complex_mode_03 <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex_03)))
mcmc_basstlaurent_complex_mode_03 <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex_03)))
mcmc_bois_franc_complex_mode_03 <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex_03)))
mcmc_cote_beaupre_complex_mode_03 <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex_03)))
mcmc_charlevoix_complex_mode_03 <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex_03)))
#VM
mcmc_gaspesie_complex_dam_mode_03 <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex_dam_03)))
mcmc_laurentides_complex_dam_mode_03 <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex_dam_03)))
mcmc_basstlaurent_complex_dam_mode_03 <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex_dam_03)))
mcmc_bois_franc_complex_dam_mode_03 <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex_dam_03)))
mcmc_cote_beaupre_complex_dam_mode_03 <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex_dam_03)))
mcmc_charlevoix_complex_dam_mode_03 <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex_dam_03)))
#VA
mcmc_gaspesie_complex_mean_03 <- as.vector(colMeans(as.mcmc(mcmc_gaspesie_complex_03)))
mcmc_laurentides_complex_mean_03 <- as.vector(colMeans(as.mcmc(mcmc_laurentides_complex_03)))
mcmc_basstlaurent_complex_mean_03 <- as.vector(colMeans(as.mcmc(mcmc_basstlaurent_complex_03)))
mcmc_bois_franc_complex_mean_03 <- as.vector(colMeans(as.mcmc(mcmc_bois_franc_complex_03)))
mcmc_cote_beaupre_complex_mean_03 <- as.vector(colMeans(as.mcmc(mcmc_cote_beaupre_complex_03)))
mcmc_charlevoix_complex_mean_03 <- as.vector(colMeans(as.mcmc(mcmc_charlevoix_complex_03)))
#VM
mcmc_gaspesie_complex_dam_mean_03 <- as.vector(colMeans(as.mcmc(mcmc_gaspesie_complex_dam_03)))
mcmc_laurentides_complex_dam_mean_03 <- as.vector(colMeans(as.mcmc(mcmc_laurentides_complex_dam_03)))
mcmc_basstlaurent_complex_dam_mean_03 <- as.vector(colMeans(as.mcmc(mcmc_basstlaurent_complex_dam_03)))
mcmc_bois_franc_complex_dam_mean_03 <- as.vector(colMeans(as.mcmc(mcmc_bois_franc_complex_dam_03)))
mcmc_cote_beaupre_complex_dam_mean_03 <- as.vector(colMeans(as.mcmc(mcmc_cote_beaupre_complex_dam_03)))
mcmc_charlevoix_complex_dam_mean_03 <- as.vector(colMeans(as.mcmc(mcmc_charlevoix_complex_dam_03)))

#Va=0.1
mcmc_gaspesie_complex_05 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_gaspesie_diag(c(0.1, 0.5), 2).rds"))
mcmc_laurentides_complex_05 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_laurentides_diag(c(0.1, 0.5), 2).rds"))
mcmc_basstlaurent_complex_05 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_basstlaurent_diag(c(0.1, 0.5), 2).rds"))
mcmc_bois_franc_complex_05 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_bois_franc_diag(c(0.1, 0.5), 2).rds"))
mcmc_cote_beaupre_complex_05 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_cote_beaupre_diag(c(0.1, 0.5), 2).rds"))
mcmc_charlevoix_complex_05 <- colnames_fun(readRDS(file = "model_complex/mcmc_simped_charlevoix_diag(c(0.1, 0.5), 2).rds"))
#& Vm=0.5
mcmc_gaspesie_complex_dam_05 <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_gaspesie_diag(c(0.1, 0.5), 2).rds"))
mcmc_laurentides_complex_dam_05 <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_laurentides_diag(c(0.1, 0.5), 2).rds"))
mcmc_basstlaurent_complex_dam_05 <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_basstlaurent_diag(c(0.1, 0.5), 2).rds"))
mcmc_bois_franc_complex_dam_05 <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_bois_franc_diag(c(0.1, 0.5), 2).rds"))
mcmc_cote_beaupre_complex_dam_05 <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_cote_beaupre_diag(c(0.1, 0.5), 2).rds"))
mcmc_charlevoix_complex_dam_05 <- colnames_fun(readRDS(file = "model_complex/mcmc_sim_damped_charlevoix_diag(c(0.1, 0.5), 2).rds"))
#VA
mcmc_gaspesie_complex_mode_05 <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex_05)))
mcmc_laurentides_complex_mode_05 <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex_05)))
mcmc_basstlaurent_complex_mode_05 <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex_05)))
mcmc_bois_franc_complex_mode_05 <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex_05)))
mcmc_cote_beaupre_complex_mode_05 <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex_05)))
mcmc_charlevoix_complex_mode_05 <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex_05)))
#VM
mcmc_gaspesie_complex_dam_mode_05 <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex_dam_05)))
mcmc_laurentides_complex_dam_mode_05 <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex_dam_05)))
mcmc_basstlaurent_complex_dam_mode_05 <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex_dam_05)))
mcmc_bois_franc_complex_dam_mode_05 <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex_dam_05)))
mcmc_cote_beaupre_complex_dam_mode_05 <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex_dam_05)))
mcmc_charlevoix_complex_dam_mode_05 <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex_dam_05)))
#VA
mcmc_gaspesie_complex_mean_05 <- as.vector(colMeans(as.mcmc(mcmc_gaspesie_complex_05)))
mcmc_laurentides_complex_mean_05 <- as.vector(colMeans(as.mcmc(mcmc_laurentides_complex_05)))
mcmc_basstlaurent_complex_mean_05 <- as.vector(colMeans(as.mcmc(mcmc_basstlaurent_complex_05)))
mcmc_bois_franc_complex_mean_05 <- as.vector(colMeans(as.mcmc(mcmc_bois_franc_complex_05)))
mcmc_cote_beaupre_complex_mean_05 <- as.vector(colMeans(as.mcmc(mcmc_cote_beaupre_complex_05)))
mcmc_charlevoix_complex_mean_05 <- as.vector(colMeans(as.mcmc(mcmc_charlevoix_complex_05)))
#VM
mcmc_gaspesie_complex_dam_mean_05 <- as.vector(colMeans(as.mcmc(mcmc_gaspesie_complex_dam_05)))
mcmc_laurentides_complex_dam_mean_05 <- as.vector(colMeans(as.mcmc(mcmc_laurentides_complex_dam_05)))
mcmc_basstlaurent_complex_dam_mean_05 <- as.vector(colMeans(as.mcmc(mcmc_basstlaurent_complex_dam_05)))
mcmc_bois_franc_complex_dam_mean_05 <- as.vector(colMeans(as.mcmc(mcmc_bois_franc_complex_dam_05)))
mcmc_cote_beaupre_complex_dam_mean_05 <- as.vector(colMeans(as.mcmc(mcmc_cote_beaupre_complex_dam_05)))
mcmc_charlevoix_complex_dam_mean_05 <- as.vector(colMeans(as.mcmc(mcmc_charlevoix_complex_dam_05)))
#if using just posterior mode from simulations
df_complex_m<-data.frame("Va" = as.factor(rep(0.1,60)),
                       "Vm" = as.factor(c(rep(0.3,30),rep(0.5,30))),
                       "region" = c(rep("Gaspésie",5),rep("Bas-St-Laurent",5),rep("Laurentides",5),
                                    rep("Charlevoix",5),rep("Côte-de-Beaupré",5),rep("Bois-Francs",5),
                                    rep("Gaspésie",5),rep("Laurentides",5),rep("Bas-St-Laurent",5),
                                    rep("Charlevoix",5),rep("Côte-de-Beaupré",5),rep("Bois-Francs",5)),
                       "additive" = c(mcmc_gaspesie_complex_mode_03,mcmc_basstlaurent_complex_mode_03,mcmc_laurentides_complex_mode_03,mcmc_charlevoix_complex_mode_03,
                                      mcmc_cote_beaupre_complex_mode_03,mcmc_bois_franc_complex_mode_03,mcmc_gaspesie_complex_mode_05,mcmc_basstlaurent_complex_mode_05,
                                      mcmc_laurentides_complex_mode_05,mcmc_charlevoix_complex_mode_05,mcmc_cote_beaupre_complex_mode_05,mcmc_bois_franc_complex_mode_05),
                       "maternal" = c(mcmc_gaspesie_complex_dam_mode_03,mcmc_basstlaurent_complex_dam_mode_03,mcmc_laurentides_complex_dam_mode_03,mcmc_charlevoix_complex_dam_mode_03,
                                      mcmc_cote_beaupre_complex_dam_mode_03,mcmc_bois_franc_complex_dam_mode_03,mcmc_gaspesie_complex_dam_mode_05,mcmc_basstlaurent_complex_dam_mode_05,
                                      mcmc_laurentides_complex_dam_mode_05,mcmc_charlevoix_complex_dam_mode_05,mcmc_cote_beaupre_complex_dam_mode_05,mcmc_bois_franc_complex_dam_mode_05))

df_complex_m$region<-as.factor(df_complex_m$region)
df_complex_m$region <- factor(df_complex_m$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"))

prec_mcmc_complex_m <- sqrt((df_complex_m[which(df_complex_m$Va == 0.1),]$additive - 0.1)^2)

n <- 1
m <- n+4
prec_complex_m <- 0
for (x in 1:12) {
  prec_complex_m[x] <- mean(prec_mcmc_complex_m[n:m])
  n <- n+5
  m <- m+5
}
prec_complex_m <- round(prec_complex_m,digits = 3)

#Va=0.1
mcmc_gaspesie_complex_03_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex_03))))
mcmc_laurentides_complex_03_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex_03))))
mcmc_basstlaurent_complex_03_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex_03))))
mcmc_bois_franc_complex_03_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex_03))))
mcmc_cote_beaupre_complex_03_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex_03))))
mcmc_charlevoix_complex_03_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex_03))))

mcmc_gaspesie_complex_03_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex_03)))[1:5])
mcmc_laurentides_complex_03_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex_03)))[1:5])
mcmc_basstlaurent_complex_03_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex_03)))[1:5])
mcmc_bois_franc_complex_03_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex_03)))[1:5])
mcmc_cote_beaupre_complex_03_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex_03)))[1:5])
mcmc_charlevoix_complex_03_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex_03)))[1:5])

mcmc_gaspesie_complex_03_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex_03)))[6:10])
mcmc_laurentides_complex_03_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex_03)))[6:10])
mcmc_basstlaurent_complex_03_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex_03)))[6:10])
mcmc_bois_franc_complex_03_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex_03)))[6:10])
mcmc_cote_beaupre_complex_03_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex_03)))[6:10])
mcmc_charlevoix_complex_03_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex_03)))[6:10])
#Vm=0.3
mcmc_gaspesie_complex_mean_dam_03 <- mean(as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex_dam_03))))
mcmc_laurentides_complex_mean_dam_03 <- mean(as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex_dam_03))))
mcmc_basstlaurent_complex_mean_dam_03 <- mean(as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex_dam_03))))
mcmc_bois_franc_complex_mean_dam_03 <- mean(as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex_dam_03))))
mcmc_cote_beaupre_complex_mean_dam_03 <- mean(as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex_dam_03))))
mcmc_charlevoix_complex_mean_dam_03 <- mean(as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex_dam_03))))

mcmc_gaspesie_complex_lower_dam_03 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex_dam_03)))[1:5])
mcmc_laurentides_complex_lower_dam_03 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex_dam_03)))[1:5])
mcmc_basstlaurent_complex_lower_dam_03 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex_dam_03)))[1:5])
mcmc_bois_franc_complex_lower_dam_03 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex_dam_03)))[1:5])
mcmc_cote_beaupre_complex_lower_dam_03 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex_dam_03)))[1:5])
mcmc_charlevoix_complex_lower_dam_03 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex_dam_03)))[1:5])

mcmc_gaspesie_complex_upper_dam_03 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex_dam_03)))[6:10])
mcmc_laurentides_complex_upper_dam_03 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex_dam_03)))[6:10])
mcmc_basstlaurent_complex_upper_dam_03 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex_dam_03)))[6:10])
mcmc_bois_franc_complex_upper_dam_03 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex_dam_03)))[6:10])
mcmc_cote_beaupre_complex_upper_dam_03 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex_dam_03)))[6:10])
mcmc_charlevoix_complex_upper_dam_03 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex_dam_03)))[6:10])
#Va=0.1
mcmc_gaspesie_complex_05_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex_05))))
mcmc_laurentides_complex_05_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex_05))))
mcmc_basstlaurent_complex_05_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex_05))))
mcmc_bois_franc_complex_05_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex_05))))
mcmc_cote_beaupre_complex_05_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex_05))))
mcmc_charlevoix_complex_05_mean <- mean(as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex_05))))

mcmc_gaspesie_complex_05_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex_05)))[1:5])
mcmc_laurentides_complex_05_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex_05)))[1:5])
mcmc_basstlaurent_complex_05_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex_05)))[1:5])
mcmc_bois_franc_complex_05_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex_05)))[1:5])
mcmc_cote_beaupre_complex_05_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex_05)))[1:5])
mcmc_charlevoix_complex_05_lower <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex_05)))[1:5])

mcmc_gaspesie_complex_05_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex_05)))[6:10])
mcmc_laurentides_complex_05_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex_05)))[6:10])
mcmc_basstlaurent_complex_05_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex_05)))[6:10])
mcmc_bois_franc_complex_05_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex_05)))[6:10])
mcmc_cote_beaupre_complex_05_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex_05)))[6:10])
mcmc_charlevoix_complex_05_upper <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex_05)))[6:10])
#Vm=0.5
mcmc_gaspesie_complex_mean_dam_05 <- mean(as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_complex_dam_05))))
mcmc_laurentides_complex_mean_dam_05 <- mean(as.vector(posterior.mode(as.mcmc(mcmc_laurentides_complex_dam_05))))
mcmc_basstlaurent_complex_mean_dam_05 <- mean(as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_complex_dam_05))))
mcmc_bois_franc_complex_mean_dam_05 <- mean(as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_complex_dam_05))))
mcmc_cote_beaupre_complex_mean_dam_05 <- mean(as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_complex_dam_05))))
mcmc_charlevoix_complex_mean_dam_05 <- mean(as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_complex_dam_05))))

mcmc_gaspesie_complex_lower_dam_05 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex_dam_05)))[1:5])
mcmc_laurentides_complex_lower_dam_05 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex_dam_05)))[1:5])
mcmc_basstlaurent_complex_lower_dam_05 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex_dam_05)))[1:5])
mcmc_bois_franc_complex_lower_dam_05 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex_dam_05)))[1:5])
mcmc_cote_beaupre_complex_lower_dam_05 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex_dam_05)))[1:5])
mcmc_charlevoix_complex_lower_dam_05 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex_dam_05)))[1:5])

mcmc_gaspesie_complex_upper_dam_05 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_gaspesie_complex_dam_05)))[6:10])
mcmc_laurentides_complex_upper_dam_05 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_laurentides_complex_dam_05)))[6:10])
mcmc_basstlaurent_complex_upper_dam_05 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_basstlaurent_complex_dam_05)))[6:10])
mcmc_bois_franc_complex_upper_dam_05 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_bois_franc_complex_dam_05)))[6:10])
mcmc_cote_beaupre_complex_upper_dam_05 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_cote_beaupre_complex_dam_05)))[6:10])
mcmc_charlevoix_complex_upper_dam_05 <- mean(as.vector(HPDinterval(as.mcmc(mcmc_charlevoix_complex_dam_05)))[6:10])
#dataframe
df_mcmc_complex_2 <- data.frame(region = rep(c("Gaspésie", "Bas-St-Laurent","Laurentides",
                                             "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"),4),
                              Estimate = c(mcmc_gaspesie_complex_03_mean,
                                           mcmc_basstlaurent_complex_03_mean,
                                           mcmc_laurentides_complex_03_mean,
                                           mcmc_charlevoix_complex_03_mean,
                                           mcmc_cote_beaupre_complex_03_mean,
                                           mcmc_bois_franc_complex_03_mean,
                                           mcmc_gaspesie_complex_05_mean,
                                           mcmc_basstlaurent_complex_05_mean,
                                           mcmc_laurentides_complex_05_mean,
                                           mcmc_charlevoix_complex_05_mean,
                                           mcmc_cote_beaupre_complex_05_mean,
                                           mcmc_bois_franc_complex_05_mean,
                                           mcmc_gaspesie_complex_mean_dam_03,
                                           mcmc_basstlaurent_complex_mean_dam_03,
                                           mcmc_laurentides_complex_mean_dam_03,
                                           mcmc_charlevoix_complex_mean_dam_03,
                                           mcmc_cote_beaupre_complex_mean_dam_03,
                                           mcmc_bois_franc_complex_mean_dam_03,
                                           mcmc_gaspesie_complex_mean_dam_05,
                                           mcmc_basstlaurent_complex_mean_dam_05,
                                           mcmc_laurentides_complex_mean_dam_05,
                                           mcmc_charlevoix_complex_mean_dam_05,
                                           mcmc_cote_beaupre_complex_mean_dam_05,
                                           mcmc_bois_franc_complex_mean_dam_05),
                              Lower = c(mcmc_gaspesie_complex_03_lower,
                                        mcmc_basstlaurent_complex_03_lower,
                                        mcmc_laurentides_complex_03_lower,
                                        mcmc_charlevoix_complex_03_lower,
                                        mcmc_cote_beaupre_complex_03_lower,
                                        mcmc_bois_franc_complex_03_lower,
                                        mcmc_gaspesie_complex_05_lower,
                                        mcmc_basstlaurent_complex_05_lower,
                                        mcmc_laurentides_complex_05_lower,
                                        mcmc_charlevoix_complex_05_lower,
                                        mcmc_cote_beaupre_complex_05_lower,
                                        mcmc_bois_franc_complex_05_lower,
                                        mcmc_gaspesie_complex_lower_dam_03,
                                        mcmc_basstlaurent_complex_lower_dam_03,
                                        mcmc_laurentides_complex_lower_dam_03,
                                        mcmc_charlevoix_complex_lower_dam_03,
                                        mcmc_cote_beaupre_complex_lower_dam_03,
                                        mcmc_bois_franc_complex_lower_dam_03,
                                        mcmc_gaspesie_complex_lower_dam_05,
                                        mcmc_basstlaurent_complex_lower_dam_05,
                                        mcmc_laurentides_complex_lower_dam_05,
                                        mcmc_charlevoix_complex_lower_dam_05,
                                        mcmc_cote_beaupre_complex_lower_dam_05,
                                        mcmc_bois_franc_complex_lower_dam_05),
                              Upper = c(mcmc_gaspesie_complex_03_upper,
                                        mcmc_basstlaurent_complex_03_upper,
                                        mcmc_laurentides_complex_03_upper,
                                        mcmc_charlevoix_complex_03_upper,
                                        mcmc_cote_beaupre_complex_03_upper,
                                        mcmc_bois_franc_complex_03_upper,
                                        mcmc_gaspesie_complex_05_upper,
                                        mcmc_basstlaurent_complex_05_upper,
                                        mcmc_laurentides_complex_05_upper,
                                        mcmc_charlevoix_complex_05_upper,
                                        mcmc_cote_beaupre_complex_05_upper,
                                        mcmc_bois_franc_complex_05_upper,
                                        mcmc_gaspesie_complex_upper_dam_03,
                                        mcmc_basstlaurent_complex_upper_dam_03,
                                        mcmc_laurentides_complex_upper_dam_03,
                                        mcmc_charlevoix_complex_upper_dam_03,
                                        mcmc_cote_beaupre_complex_upper_dam_03,
                                        mcmc_bois_franc_complex_upper_dam_03,
                                        mcmc_gaspesie_complex_upper_dam_05,
                                        mcmc_basstlaurent_complex_upper_dam_05,
                                        mcmc_laurentides_complex_upper_dam_05,
                                        mcmc_charlevoix_complex_upper_dam_05,
                                        mcmc_cote_beaupre_complex_upper_dam_05,
                                        mcmc_bois_franc_complex_upper_dam_05),
                              Vm = as.factor(rep(c(rep(0.3,6),rep(0.5,6)),2)),
                              Var = as.factor(rep(c(rep("Additive",12),rep("Maternal",12)))))

df_mcmc_complex_2$region<-as.factor(df_mcmc_complex_2$region)
df_mcmc_complex_2$region <- factor(df_mcmc_complex_2$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix","Côte-de-Beaupré","Bois-Francs"))

#######
#####3rd scenario, simple architecture with genetic correlation with Va= 0.1 & 0.3 and r= 0.1 & 0.3#####
#Va=0.1 & r=0.1
mcmc_gaspesie_bivariate <- readRDS(file = "model_bivariate/mcmc_sim_ped_gaspesie_matrix(c(0.1, 0.01, 0.01, 0.1), 2).rds")
mcmc_laurentides_bivariate <- readRDS(file = "model_bivariate/mcmc_sim_ped_laurentides_matrix(c(0.1, 0.01, 0.01, 0.1), 2).rds")
mcmc_basstlaurent_bivariate <- readRDS(file = "model_bivariate/mcmc_sim_ped_basstlaurent_matrix(c(0.1, 0.01, 0.01, 0.1), 2).rds")
mcmc_bois_franc_bivariate <- readRDS(file = "model_bivariate/mcmc_sim_ped_bois_franc_matrix(c(0.1, 0.01, 0.01, 0.1), 2).rds")
mcmc_cote_beaupre_bivariate <- readRDS(file = "model_bivariate/mcmc_sim_ped_cote_beaupre_matrix(c(0.1, 0.01, 0.01, 0.1), 2).rds")
mcmc_charlevoix_bivariate <- readRDS(file = "model_bivariate/mcmc_sim_ped_charlevoix_matrix(c(0.1, 0.01, 0.01, 0.1), 2).rds")
#Va=0.1 & r=0.5
mcmc_gaspesie_bivariate_1 <- readRDS(file = "model_bivariate/mcmc_sim_ped_gaspesie_matrix(c(0.1, 0.05, 0.05, 0.1), 2).rds")
mcmc_laurentides_bivariate_1 <- readRDS(file = "model_bivariate/mcmc_sim_ped_laurentides_matrix(c(0.1, 0.05, 0.05, 0.1), 2).rds")
mcmc_basstlaurent_bivariate_1 <- readRDS(file = "model_bivariate/mcmc_sim_ped_basstlaurent_matrix(c(0.1, 0.05, 0.05, 0.1), 2).rds")
mcmc_bois_franc_bivariate_1 <- readRDS(file = "model_bivariate/mcmc_sim_ped_bois_franc_matrix(c(0.1, 0.05, 0.05, 0.1), 2).rds")
mcmc_cote_beaupre_bivariate_1 <- readRDS(file = "model_bivariate/mcmc_sim_ped_cote_beaupre_matrix(c(0.1, 0.05, 0.05, 0.1), 2).rds")
mcmc_charlevoix_bivariate_1 <- readRDS(file = "model_bivariate/mcmc_sim_ped_charlevoix_matrix(c(0.1, 0.05, 0.05, 0.1), 2).rds")

#do not use if using all iterations
mcmc_gaspesie_bivariate_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_bivariate[,2])))
mcmc_laurentides_bivariate_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_bivariate[,2])))
mcmc_basstlaurent_bivariate_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_bivariate[,2])))
mcmc_bois_franc_bivariate_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_bivariate[,2])))
mcmc_cote_beaupre_bivariate_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_bivariate[,2])))
mcmc_charlevoix_bivariate_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_bivariate[,2])))

mcmc_gaspesie_bivariate_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_bivariate_1[,2])))
mcmc_laurentides_bivariate_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_bivariate_1[,2])))
mcmc_basstlaurent_bivariate_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_bivariate_1[,2])))
mcmc_bois_franc_bivariate_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_bivariate_1[,2])))
mcmc_cote_beaupre_bivariate_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_bivariate_1[,2])))
mcmc_charlevoix_bivariate_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_bivariate_1[,2])))

#Va=0.3 & r=0.1
mcmc_gaspesie_bivariate_2 <- readRDS(file = "model_bivariate/mcmc_sim_ped_gaspesie_matrix(c(0.3, 0.0173, 0.0173, 0.1), 2).rds")
mcmc_laurentides_bivariate_2 <- readRDS(file = "model_bivariate/mcmc_sim_ped_laurentides_matrix(c(0.3, 0.0173, 0.0173, 0.1), 2).rds")
mcmc_basstlaurent_bivariate_2 <- readRDS(file = "model_bivariate/mcmc_sim_ped_basstlaurent_matrix(c(0.3, 0.0173, 0.0173, 0.1), 2).rds")
mcmc_bois_franc_bivariate_2 <- readRDS(file = "model_bivariate/mcmc_sim_ped_bois_franc_matrix(c(0.3, 0.0173, 0.0173, 0.1), 2).rds")
mcmc_cote_beaupre_bivariate_2 <- readRDS(file = "model_bivariate/mcmc_sim_ped_cote_beaupre_matrix(c(0.3, 0.0173, 0.0173, 0.1), 2).rds")
mcmc_charlevoix_bivariate_2 <- readRDS(file = "model_bivariate/mcmc_sim_ped_charlevoix_matrix(c(0.3, 0.0173, 0.0173, 0.1), 2).rds")

#Va=0.3 & r=0.5
mcmc_gaspesie_bivariate_2_1 <- readRDS(file = "model_bivariate/mcmc_sim_ped_gaspesie_matrix(c(0.3, 0.0866, 0.0866, 0.1), 2).rds")
mcmc_laurentides_bivariate_2_1 <- readRDS(file = "model_bivariate/mcmc_sim_ped_laurentides_matrix(c(0.3, 0.0866, 0.0866, 0.1), 2).rds")
mcmc_basstlaurent_bivariate_2_1 <- readRDS(file = "model_bivariate/mcmc_sim_ped_basstlaurent_matrix(c(0.3, 0.0866, 0.0866, 0.1), 2).rds")
mcmc_bois_franc_bivariate_2_1 <- readRDS(file = "model_bivariate/mcmc_sim_ped_bois_franc_matrix(c(0.3, 0.0866, 0.0866, 0.1), 2).rds")
mcmc_cote_beaupre_bivariate_2_1 <- readRDS(file = "model_bivariate/mcmc_sim_ped_cote_beaupre_matrix(c(0.3, 0.0866, 0.0866, 0.1), 2).rds")
mcmc_charlevoix_bivariate_2_1 <- readRDS(file = "model_bivariate/mcmc_sim_ped_charlevoix_matrix(c(0.3, 0.0866, 0.0866, 0.1), 2).rds")

#do not use if using all iterations
mcmc_gaspesie_bivariate_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_bivariate_2[,2])))
mcmc_laurentides_bivariate_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_bivariate_2[,2])))
mcmc_basstlaurent_bivariate_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_bivariate_2[,2])))
mcmc_bois_franc_bivariate_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_bivariate_2[,2])))
mcmc_cote_beaupre_bivariate_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_bivariate_2[,2])))
mcmc_charlevoix_bivariate_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_bivariate_2[,2])))

mcmc_gaspesie_bivariate_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_bivariate_2_1[,2])))
mcmc_laurentides_bivariate_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_bivariate_2_1[,2])))
mcmc_basstlaurent_bivariate_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_bivariate_2_1[,2])))
mcmc_bois_franc_bivariate_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_bivariate_2_1[,2])))
mcmc_cote_beaupre_bivariate_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_bivariate_2_1[,2])))
mcmc_charlevoix_bivariate_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_bivariate_2_1[,2])))
#for r=0.1
#if using all iterations from simulations
df_0.1<-cbind.data.frame("Gaspésie"= mcmc_gaspesie_bivariate[,2],"Bas-St-Laurent"=mcmc_basstlaurent_bivariate[,2],
                         "Laurentides"=mcmc_laurentides_bivariate[,2],"Charlevoix"=mcmc_charlevoix_bivariate[,2],
                         "Côte-de-Beaupré"=mcmc_cote_beaupre_bivariate[,2],"Bois-Francs"=mcmc_bois_franc_bivariate[,2])
df_0.3 <- cbind.data.frame("Gaspésie"= mcmc_gaspesie_bivariate_2[,2],"Bas-St-Laurent"=mcmc_basstlaurent_bivariate_2[,2],
                           "Laurentides"=mcmc_laurentides_bivariate_2[,2],"Charlevoix"=mcmc_charlevoix_bivariate_2[,2],
                           "Côte-de-Beaupré"=mcmc_cote_beaupre_bivariate_2[,2],"Bois-Francs"=mcmc_bois_franc_bivariate_2[,2])

#if using just posterior mode from simulations

df_bivariate<-data.frame("Va" = as.factor(c(rep(0.1,6),rep(0.3,6))),
               "region" = c(rep("Gaspésie",1),rep("Bas-St-Laurent",1),rep("Laurentides",1),
                            rep("Charlevoix",1),rep("Côte-de-Beaupré",1),rep("Bois-Francs",1),
                            rep("Gaspésie",1),rep("Laurentides",1),rep("Bas-St-Laurent",1),
                            rep("Charlevoix",1),rep("Côte-de-Beaupré",1),rep("Bois-Francs",1)),
               "additive" = c(mcmc_gaspesie_bivariate_mode,mcmc_basstlaurent_bivariate_mode,mcmc_laurentides_bivariate_mode,mcmc_charlevoix_bivariate_mode,mcmc_cote_beaupre_bivariate_mode,mcmc_bois_franc_bivariate_mode,
                              mcmc_gaspesie_bivariate_2_mode,mcmc_basstlaurent_bivariate_2_mode,mcmc_laurentides_bivariate_2_mode,mcmc_charlevoix_bivariate_2_mode,mcmc_cote_beaupre_bivariate_2_mode,mcmc_bois_franc_bivariate_2_mode))

df_bivariate$region<-as.factor(df_bivariate$region)
df_bivariate$region <- factor(df_bivariate$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"))

prec_0.1_mcmc_bivariate <- sqrt((df_bivariate[which(df_bivariate$Va == 0.1),]$additive - 0.1)^2)
prec_0.3_mcmc_bivariate <- sqrt((df_bivariate[which(df_bivariate$Va == 0.3),]$additive - 0.3)^2)

prec_0.1_bivariate <- round(prec_0.1_mcmc_bivariate,digits = 3)
prec_0.3_bivariate <- round(prec_0.3_mcmc_bivariate,digits = 3)

#for r=0.5
#if using all iterations from simulations
df_0.1_1<-cbind.data.frame("Gaspésie"= mcmc_gaspesie_bivariate_1[,2],"Bas-St-Laurent"=mcmc_basstlaurent_bivariate_1[,2],
                         "Laurentides"=mcmc_laurentides_bivariate_1[,2],"Charlevoix"=mcmc_charlevoix_bivariate_1[,2],
                         "Côte-de-Beaupré"=mcmc_cote_beaupre_bivariate_1[,2],"Bois-Francs"=mcmc_bois_franc_bivariate_1[,2])
df_0.3_1 <- cbind.data.frame("Gaspésie"= mcmc_gaspesie_bivariate_2_1[,2],"Bas-St-Laurent"=mcmc_basstlaurent_bivariate_2_1[,2],
                           "Laurentides"=mcmc_laurentides_bivariate_2_1[,2],"Charlevoix"=mcmc_charlevoix_bivariate_2_1[,2],
                           "Côte-de-Beaupré"=mcmc_cote_beaupre_bivariate_2_1[,2],"Bois-Francs"=mcmc_bois_franc_bivariate_2_1[,2])
#if using just posterior mode from simulations
df_bivariate_1<-data.frame("Va" = as.factor(c(rep(0.1,6),rep(0.3,6))),
               "region" = c(rep("Gaspésie",1),rep("Bas-St-Laurent",1),rep("Laurentides",1),
                            rep("Charlevoix",1),rep("Côte-de-Beaupré",1),rep("Bois-Francs",1),
                            rep("Gaspésie",1),rep("Laurentides",1),rep("Bas-St-Laurent",1),
                            rep("Charlevoix",1),rep("Côte-de-Beaupré",1),rep("Bois-Francs",1)),
               "additive" = c(mcmc_gaspesie_bivariate_1_mode,mcmc_basstlaurent_bivariate_1_mode,mcmc_laurentides_bivariate_1_mode,mcmc_charlevoix_bivariate_1_mode,mcmc_cote_beaupre_bivariate_1_mode,mcmc_bois_franc_bivariate_1_mode,
                              mcmc_gaspesie_bivariate_2_1_mode,mcmc_basstlaurent_bivariate_2_1_mode,mcmc_laurentides_bivariate_2_1_mode,mcmc_charlevoix_bivariate_2_1_mode,mcmc_cote_beaupre_bivariate_2_1_mode,mcmc_bois_franc_bivariate_2_1_mode))

df_bivariate_1$region<-as.factor(df_bivariate_1$region)
df_bivariate_1$region <- factor(df_bivariate_1$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"))

prec_0.1_mcmc_bivariate_1 <- ((df_bivariate_1[which(df_bivariate_1$Va == 0.1),]$additive - 0.1)^2)
prec_0.3_mcmc_bivariate_1 <- ((df_bivariate_1[which(df_bivariate_1$Va == 0.3),]$additive - 0.3)^2)

prec_0.1_bivariate_1 <- round(prec_0.1_mcmc_bivariate_1,digits = 3)
prec_0.3_bivariate_1 <- round(prec_0.3_mcmc_bivariate_1,digits = 3)


#####genetic correlation
#r=0.1
mcmc_gaspesie_bivariate_cov_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_bivariate[,1])))
mcmc_laurentides_bivariate_cov_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_bivariate[,1])))
mcmc_basstlaurent_bivariate_cov_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_bivariate[,1])))
mcmc_bois_franc_bivariate_cov_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_bivariate[,1])))
mcmc_cote_beaupre_bivariate_cov_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_bivariate[,1])))
mcmc_charlevoix_bivariate_cov_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_bivariate[,1])))

mcmc_gaspesie_bivariate_cov_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_bivariate_2[,1])))
mcmc_laurentides_bivariate_cov_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_bivariate_2[,1])))
mcmc_basstlaurent_bivariate_cov_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_bivariate_2[,1])))
mcmc_bois_franc_bivariate_cov_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_bivariate_2[,1])))
mcmc_cote_beaupre_bivariate_cov_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_bivariate_2[,1])))
mcmc_charlevoix_bivariate_cov_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_bivariate_2[,1])))

mcmc_gaspesie_bivariate_corr_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_bivariate[,1])/sqrt(as.mcmc(mcmc_gaspesie_bivariate[,2])*rep(0.1,1000))))
mcmc_laurentides_bivariate_corr_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_bivariate[,1])/sqrt(as.mcmc(mcmc_laurentides_bivariate[,2])*rep(0.1,1000))))
mcmc_basstlaurent_bivariate_corr_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_bivariate[,1])/sqrt(as.mcmc(mcmc_basstlaurent_bivariate[,2])*rep(0.1,1000))))
mcmc_bois_franc_bivariate_corr_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_bivariate[,1])/sqrt(as.mcmc(mcmc_bois_franc_bivariate[,2])*rep(0.1,1000))))
mcmc_cote_beaupre_bivariate_corr_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_bivariate[,1])/sqrt(as.mcmc(mcmc_cote_beaupre_bivariate[,2])*rep(0.1,1000))))
mcmc_charlevoix_bivariate_corr_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_bivariate[,1])/sqrt(as.mcmc(mcmc_charlevoix_bivariate[,2])*rep(0.1,1000))))

mcmc_gaspesie_bivariate_corr_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_bivariate_2[,1])/sqrt(as.mcmc(mcmc_gaspesie_bivariate_2[,2])*rep(0.1,1000))))
mcmc_laurentides_bivariate_corr_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_bivariate_2[,1])/sqrt(as.mcmc(mcmc_laurentides_bivariate_2[,2])*rep(0.1,1000))))
mcmc_basstlaurent_bivariate_corr_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_bivariate_2[,1])/sqrt(as.mcmc(mcmc_basstlaurent_bivariate_2[,2])*rep(0.1,1000))))
mcmc_bois_franc_bivariate_corr_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_bivariate_2[,1])/sqrt(as.mcmc(mcmc_bois_franc_bivariate_2[,2])*rep(0.1,1000))))
mcmc_cote_beaupre_bivariate_corr_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_bivariate_2[,1])/sqrt(as.mcmc(mcmc_cote_beaupre_bivariate_2[,2])*rep(0.1,1000))))
mcmc_charlevoix_bivariate_corr_2_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_bivariate_2[,1])/sqrt(as.mcmc(mcmc_charlevoix_bivariate_2[,2])*rep(0.1,1000))))
#r=0.5
mcmc_gaspesie_bivariate_cov_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_bivariate_1[,1])))
mcmc_laurentides_bivariate_cov_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_bivariate_1[,1])))
mcmc_basstlaurent_bivariate_cov_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_bivariate_1[,1])))
mcmc_bois_franc_bivariate_cov_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_bivariate_1[,1])))
mcmc_cote_beaupre_bivariate_cov_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_bivariate_1[,1])))
mcmc_charlevoix_bivariate_cov_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_bivariate_1[,1])))

mcmc_gaspesie_bivariate_cov_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_bivariate_2_1[,1])))
mcmc_laurentides_bivariate_cov_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_bivariate_2_1[,1])))
mcmc_basstlaurent_bivariate_cov_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_bivariate_2_1[,1])))
mcmc_bois_franc_bivariate_cov_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_bivariate_2_1[,1])))
mcmc_cote_beaupre_bivariate_cov_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_bivariate_2_1[,1])))
mcmc_charlevoix_bivariate_cov_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_bivariate_2_1[,1])))

mcmc_gaspesie_bivariate_corr_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_bivariate_1[,1])/sqrt(as.mcmc(mcmc_gaspesie_bivariate_1[,2])*rep(0.1,1000))))
mcmc_laurentides_bivariate_corr_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_bivariate_1[,1])/sqrt(as.mcmc(mcmc_laurentides_bivariate_1[,2])*rep(0.1,1000))))
mcmc_basstlaurent_bivariate_corr_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_bivariate_1[,1])/sqrt(as.mcmc(mcmc_basstlaurent_bivariate_1[,2])*rep(0.1,1000))))
mcmc_bois_franc_bivariate_corr_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_bivariate_1[,1])/sqrt(as.mcmc(mcmc_bois_franc_bivariate_1[,2])*rep(0.1,1000))))
mcmc_cote_beaupre_bivariate_corr_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_bivariate_1[,1])/sqrt(as.mcmc(mcmc_cote_beaupre_bivariate_1[,2])*rep(0.1,1000))))
mcmc_charlevoix_bivariate_corr_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_bivariate_1[,1])/sqrt(as.mcmc(mcmc_charlevoix_bivariate_1[,2])*rep(0.1,1000))))

mcmc_gaspesie_bivariate_corr_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_gaspesie_bivariate_2[,1])/sqrt(as.mcmc(mcmc_gaspesie_bivariate_2_1[,2])*rep(0.1,1000))))
mcmc_laurentides_bivariate_corr_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_laurentides_bivariate_2[,1])/sqrt(as.mcmc(mcmc_laurentides_bivariate_2_1[,2])*rep(0.1,1000))))
mcmc_basstlaurent_bivariate_corr_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_basstlaurent_bivariate_2[,1])/sqrt(as.mcmc(mcmc_basstlaurent_bivariate_2_1[,2])*rep(0.1,1000))))
mcmc_bois_franc_bivariate_corr_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_bois_franc_bivariate_2[,1])/sqrt(as.mcmc(mcmc_bois_franc_bivariate_2_1[,2])*rep(0.1,1000))))
mcmc_cote_beaupre_bivariate_corr_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_cote_beaupre_bivariate_2[,1])/sqrt(as.mcmc(mcmc_cote_beaupre_bivariate_2_1[,2])*rep(0.1,1000))))
mcmc_charlevoix_bivariate_corr_2_1_mode <- as.vector(posterior.mode(as.mcmc(mcmc_charlevoix_bivariate_2[,1])/sqrt(as.mcmc(mcmc_charlevoix_bivariate_2_1[,2])*rep(0.1,1000))))

#if using just posterior mode from simulations
#for r=0.1
df_corr<-data.frame("cov" = as.factor(c(rep(0.01,6),rep(0.0173,6))),
               "region" = c(rep("Gaspésie",1),rep("Bas-St-Laurent",1),rep("Laurentides",1),
                            rep("Charlevoix",1),rep("Côte-de-Beaupré",1),rep("Bois-Francs",1),
                            rep("Gaspésie",1),rep("Bas-St-Laurent",1),rep("Laurentides",1),
                            rep("Charlevoix",1),rep("Côte-de-Beaupré",1),rep("Bois-Francs",1)),
               "additive" = c(mcmc_gaspesie_bivariate_cov_mode,mcmc_basstlaurent_bivariate_cov_mode,
                              mcmc_laurentides_bivariate_cov_mode,mcmc_charlevoix_bivariate_cov_mode,
                              mcmc_cote_beaupre_bivariate_cov_mode,mcmc_bois_franc_bivariate_cov_mode,
                              mcmc_gaspesie_bivariate_cov_2_mode,mcmc_basstlaurent_bivariate_cov_2_mode,
                              mcmc_laurentides_bivariate_cov_2_mode,mcmc_charlevoix_bivariate_cov_2_mode,
                              mcmc_cote_beaupre_bivariate_cov_2_mode,mcmc_bois_franc_bivariate_cov_2_mode))

df_corr$region<-as.factor(df_corr$region)
df_corr$region <- factor(df_corr$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"))

prec_0.1_mcmc_bivariate_corr <- c(sqrt((df_corr[which(df_corr$cov == 0.01),]$additive - 0.01)^2), sqrt((df_corr[which(df_corr$cov == 0.0173),]$additive - 0.0173)^2))

prec_0.1_bivariate_corr <- round(prec_0.1_mcmc_bivariate_corr,digits = 3)

#for r=0.5
df_corr_1<-data.frame("cov" = as.factor(c(rep(0.05,6),rep(0.0866,6))),
                    "region" = c(rep("Gaspésie",1),rep("Bas-St-Laurent",1),rep("Laurentides",1),
                                 rep("Charlevoix",1),rep("Côte-de-Beaupré",1),rep("Bois-Francs",1),
                                 rep("Gaspésie",1),rep("Bas-St-Laurent",1),rep("Laurentides",1),
                                 rep("Charlevoix",1),rep("Côte-de-Beaupré",1),rep("Bois-Francs",1)),
                    "additive" = c(mcmc_gaspesie_bivariate_cov_1_mode,mcmc_basstlaurent_bivariate_cov_1_mode,
                                   mcmc_laurentides_bivariate_cov_1_mode,mcmc_charlevoix_bivariate_cov_1_mode,
                                   mcmc_cote_beaupre_bivariate_cov_1_mode,mcmc_bois_franc_bivariate_cov_1_mode,
                                   mcmc_gaspesie_bivariate_cov_2_1_mode,mcmc_basstlaurent_bivariate_cov_2_1_mode,
                                   mcmc_laurentides_bivariate_cov_2_1_mode,mcmc_charlevoix_bivariate_cov_2_1_mode,
                                   mcmc_cote_beaupre_bivariate_cov_2_1_mode,mcmc_bois_franc_bivariate_cov_2_1_mode))

df_corr_1$region<-as.factor(df_corr_1$region)
df_corr_1$region <- factor(df_corr_1$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"))

prec_0.5_mcmc_bivariate_corr <- c(sqrt((df_corr_1[which(df_corr_1$cov == 0.05),]$additive - 0.05)^2), sqrt((df_corr_1[which(df_corr_1$cov == 0.0866),]$additive - 0.0866)^2))

prec_0.5_bivariate_corr <- round(prec_0.5_mcmc_bivariate_corr,digits = 3)


#for r=0.1
#if using all iterations from simulations
df_r_0.1<-cbind.data.frame("Gaspésie"= mcmc_gaspesie_bivariate[,1],"Bas-St-Laurent"=mcmc_basstlaurent_bivariate[,1],
                         "Laurentides"=mcmc_laurentides_bivariate[,1],"Charlevoix"=mcmc_charlevoix_bivariate[,1],
                         "Côte-de-Beaupré"=mcmc_cote_beaupre_bivariate[,1],"Bois-Francs"=mcmc_bois_franc_bivariate[,1])
df_r_0.5 <- cbind.data.frame("Gaspésie"= mcmc_gaspesie_bivariate_1[,1],"Bas-St-Laurent"=mcmc_basstlaurent_bivariate_1[,1],
                           "Laurentides"=mcmc_laurentides_bivariate_1[,1],"Charlevoix"=mcmc_charlevoix_bivariate_1[,1],
                           "Côte-de-Beaupré"=mcmc_cote_beaupre_bivariate_1[,1],"Bois-Francs"=mcmc_bois_franc_bivariate_1[,1])
#for r=0.5
#if using all iterations from simulations
df_r_0.1_1<-cbind.data.frame("Gaspésie"= mcmc_gaspesie_bivariate_1[,1],"Bas-St-Laurent"=mcmc_basstlaurent_bivariate_1[,1],
                           "Laurentides"=mcmc_laurentides_bivariate_1[,1],"Charlevoix"=mcmc_charlevoix_bivariate_1[,1],
                           "Côte-de-Beaupré"=mcmc_cote_beaupre_bivariate_1[,1],"Bois-Francs"=mcmc_bois_franc_bivariate_1[,1])
df_r_0.5_1 <- cbind.data.frame("Gaspésie"= mcmc_gaspesie_bivariate_2_1[,1],"Bas-St-Laurent"=mcmc_basstlaurent_bivariate_2_1[,1],
                             "Laurentides"=mcmc_laurentides_bivariate_2_1[,1],"Charlevoix"=mcmc_charlevoix_bivariate_2_1[,1],
                             "Côte-de-Beaupré"=mcmc_cote_beaupre_bivariate_2_1[,1],"Bois-Francs"=mcmc_bois_franc_bivariate_2_1[,1])

######breeding values########
#Va=0.1
#model
mcmc_gaspesie_bv <- readRDS(file = "model_bv/mcmc_bv_sim_ped_gaspesie_0.1.rds")
mcmc_laurentides_bv <- readRDS(file = "model_bv/mcmc_bv_sim_ped_laurentides_0.1.rds")
mcmc_basstlaurent_bv <- readRDS(file = "model_bv/mcmc_bv_sim_ped_basstlaurent_0.1.rds")
mcmc_bois_franc_bv <- readRDS(file = "model_bv/mcmc_bv_sim_ped_bois_franc_0.1.rds")
mcmc_cote_beaupre_bv <- readRDS(file = "model_bv/mcmc_bv_sim_ped_cote_beaupre_0.1.rds")
mcmc_charlevoix_bv <- readRDS(file = "model_bv/mcmc_bv_sim_ped_charlevoix_0.1.rds")
#sim BV
gaspesie_bv <- read.table(file = "model_bv/eff_sim_ped_gaspesie_0.1.txt")
laurentides_bv <- read.table(file = "model_bv/eff_sim_ped_laurentides_0.1.txt")
basstlaurent_bv <- read.table(file = "model_bv/eff_sim_ped_basstlaurent_0.1.txt")
bois_franc_bv <- read.table(file = "model_bv/eff_sim_ped_bois_franc_0.1.txt")
cote_beaupre_bv <- read.table(file = "model_bv/eff_sim_ped_cote_beaupre_0.1.txt")
charlevoix_bv <- read.table(file = "model_bv/eff_sim_ped_charlevoix_0.1.txt")

#Va=0.3
mcmc_gaspesie_bv_2 <- readRDS(file = "model_bv/mcmc_bv_sim_ped_gaspesie_0.3.rds")
mcmc_laurentides_bv_2 <- readRDS(file = "model_bv/mcmc_bv_sim_ped_laurentides_0.3.rds")
mcmc_basstlaurent_bv_2 <- readRDS(file = "model_bv/mcmc_bv_sim_ped_basstlaurent_0.3.rds")
mcmc_bois_franc_bv_2 <- readRDS(file = "model_bv/mcmc_bv_sim_ped_bois_franc_0.3.rds")
mcmc_cote_beaupre_bv_2 <- readRDS(file = "model_bv/mcmc_bv_sim_ped_cote_beaupre_0.3.rds")
mcmc_charlevoix_bv_2 <- readRDS(file = "model_bv/mcmc_bv_sim_ped_charlevoix_0.3.rds")
#sim BV
gaspesie_bv_2 <- read.table(file = "model_bv/eff_sim_ped_gaspesie_0.3.txt")
laurentides_bv_2 <- read.table(file = "model_bv/eff_sim_ped_laurentides_0.3.txt")
basstlaurent_bv_2 <- read.table(file = "model_bv/eff_sim_ped_basstlaurent_0.3.txt")
bois_franc_bv_2 <- read.table(file = "model_bv/eff_sim_ped_bois_franc_0.3.txt")
cote_beaupre_bv_2 <- read.table(file = "model_bv/eff_sim_ped_cote_beaupre_0.3.txt")
charlevoix_bv_2 <- read.table(file = "model_bv/eff_sim_ped_charlevoix_0.3.txt")

#BV extract function
BVcorr <- function(Dataset, initBV, ModelName){
  nitt <- nrow(ModelName$Sol) #number of iterations within the model
  colnames(Dataset)[1] <- "animal" #fix first column name to animal
  colnames(initBV)[1] <- "animal" #fix first column name to aninal
  Dataset <- merge(Dataset, initBV, by.x="animal") #merge initial BV with main dataset
  EBVs <- cbind(Dataset["animal"],Dataset["a_tr1"]) #create dataframe with animal and initial BVs
  for(i in 1:nitt){
    BV <- ModelName$Sol[i,] #group all residuals by iteration (model row)
    BV <- BV[substr(names(BV),1,7)=="animal."] #substract the wanted breeding values per individual 
    #change AFR in line 7 depending on the model
    names(BV)<-substr(names(BV),8,15) #substract the names of the indvidual
    #change the margins 17,22 in line 9 to the margins that match individual's ID
    
    EBVs<-merge(EBVs,BV,by.x="animal",by.y=0) #merge and match breeding values to their individual in the dataset
    colnames(EBVs)[dim(EBVs)[2]]<-paste("BreedingValue",i,sep = "_") #name the column
  }
  cor.df <- data.frame(cor(EBVs[,2],EBVs[,-c(1,2)]))
  cor.df <- transform(cor.df, mean=apply(cor.df, 1, mean))
  cor.df <- transform(cor.df, sd=apply(cor.df, 1, sd))
  cor.df <- transform(cor.df, se=apply(cor.df, 1,function(x) x[length(x)]/sqrt(length(x))))
  return(as.numeric(cor.df[c(length(cor.df)-2,length(cor.df)-1,length(cor.df))])) #return calculate mean, sd, and se of corr coef between intial BV and each iteration of EBVs
}

#Va=0.1
my.list1 <- list(ped_gaspesie, ped_basstlaurent,ped_laurentides,ped_charlevoix,ped_cote_beaupre,ped_bois_franc)
my.list2 <- list(gaspesie_bv,basstlaurent_bv,laurentides_bv,charlevoix_bv,cote_beaupre_bv,bois_franc_bv)
my.list3 <- list(mcmc_gaspesie_bv, mcmc_basstlaurent_bv, mcmc_laurentides_bv, mcmc_charlevoix_bv, mcmc_cote_beaupre_bv, mcmc_bois_franc_bv)

corr <- matrix(numeric(0),3,6)
colnames(corr) <- c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix", "Côte-de-Beaupré", "Bois-Francs")
for (i in 1:6) {
  corr[,i] <- BVcorr(data.frame(my.list1[i]), data.frame(my.list2[i]), my.list3[[i]])
}

df_corr<- data.frame(Region= c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"),
                     Corr= as.vector(corr[1,]),
                     sd= as.vector(corr[2,]),
                     se= as.vector(corr[3,]))

df_corr$Region<-as.factor(df_corr$Region)
df_corr$Region <- factor(df_corr$Region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix","Côte-de-Beaupré","Bois-Francs"))
#Va=0.3
my.list1 <- list(ped_gaspesie, ped_basstlaurent,ped_laurentides,ped_charlevoix,ped_cote_beaupre,ped_bois_franc)
my.list2 <- list(gaspesie_bv_2,basstlaurent_bv_2,laurentides_bv_2,charlevoix_bv_2,cote_beaupre_bv_2,bois_franc_bv_2)
my.list3 <- list(mcmc_gaspesie_bv_2, mcmc_basstlaurent_bv_2, mcmc_laurentides_bv_2, mcmc_charlevoix_bv_2, mcmc_cote_beaupre_bv_2, mcmc_bois_franc_bv_2)

corr_2 <- matrix(numeric(0),3,6)
colnames(corr_2) <- c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix", "Côte-de-Beaupré", "Bois-Francs")
for (i in 1:6) {
  corr_2[,i] <- BVcorr(data.frame(my.list1[i]), data.frame(my.list2[i]), my.list3[[i]])
}

df_corr_2<- data.frame(Region= c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"),
                     Corr= as.vector(corr_2[1,]),
                     sd= as.vector(corr_2[2,]),
                     se= as.vector(corr_2[3,]))

df_corr_2$Region<-as.factor(df_corr_2$Region)
df_corr_2$Region <- factor(df_corr_2$Region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix","Côte-de-Beaupré","Bois-Francs"))

df_corr_all <- rbind.data.frame(df_corr, df_corr_2)
df_corr_all$Va<- as.factor(c(rep(0.1,6), rep(0.3,6)))

########precision#########
precision <-data.frame("Va" = as.factor(c(rep(0.1,6),rep(0.3,6),
                                          rep(0.1,6),rep(0.3,6),
                                          rep(0.1,6),rep(0.3,6),
                                          rep(0.1,6),rep(0.3,6))),
                       "region" = c(rep(c("Gaspésie","Bas-St-Laurent","Laurentides","Charlevoix","Côte-de-Beaupré","Bois-Francs"),8)),
                       "precision" = c(prec_0.1, prec_0.3, prec_0.1_complex, prec_0.3_complex, prec_0.1_bivariate,
                                       prec_0.3_bivariate, prec_0.1_bivariate_1, prec_0.3_bivariate_1))

df_corr_1$region<-as.factor(df_corr_1$region)
df_corr_1$region <- factor(df_corr_1$region, levels = c("Gaspésie", "Bas-St-Laurent","Laurentides", "Charlevoix", "Côte-de-Beaupré", "Bois-Francs"))

colnames(pedigrees_all)[1] <- "region"
df_precision <- merge(precision, pedigrees_all, by="region")

summary_precision <- summarySE(df_precision, measurevar = "precision", groupvars = c("region","Va"))

pedigrees_all$precision <-  summary_precision$precision

plot(co)