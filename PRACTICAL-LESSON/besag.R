### --- 0. Loading the packages --- ####
library(maptools)
library(rgdal)
library(spdep)
library(lattice)
library(latticeExtra)
library(viridis)
library(gridExtra)
library(RColorBrewer)
library(INLA)

#BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE)
library(Rgraphviz)
library(graph)
library(R2OpenBUGS)
library(coda)

### --- 1. Loading the data --- ####
#Dataset
data <- readRDS("data/processed/london_suic/london_suicides.RDS")
Nareas <- length(data[,1])

#Spatial polygon data frame
london.gen <- readOGR("data/processed/london_suic", "LDNSuicides")


### --- 2. Checking if the data of the sp and the data.frame match --- ####
# The order of the areas needs to be the same between 
# the data and the spatial polygon object obtained importing 
# the shapefile, so we re-order the data.
data$NAME == london.gen$NAME #Not the same
boroughs <- london.gen
data.boroughs <- attr(boroughs, "data")
order <- match(data.boroughs$NAME,data$NAME)
data <- data[order,]
data$NAME == london.gen$NAME

### ----- 2.1. Plotting the data --- ####
london.gen$SMR_raw <- data$y/data$E 
SMR_raw.cutoff<- c(0.6, 0.9, 1.0, 1.1,  1.8)
SMR_raw_disc = cut(london.gen$SMR_raw,
                   breaks         = SMR_raw.cutoff,
                   include.lowest = TRUE)

london.gen$SMR_raw_disc <- SMR_raw_disc

png("images/SMR_raw.png", width = 1000, height = 800, res = 150)
spplot(london.gen,
       c("SMR_raw_disc"),
       col.regions = brewer.pal(9,'Blues')[c(2,4,6,8)],
       main        = "SMR raw",
       par.settings =
         list(axis.line = list(col =  'transparent')))

dev.off()

# Covariates
london.gen$x1 <- data$x1
london.gen$x2 <- data$x2

png("images/covs.png", width = 1200, height = 800, res = 150)
grid.arrange(
  spplot(london.gen, c("x1"),
         main = c("index of social deprivation"),
         #col.regions = rev(viridis_pal(option = "B")(101)),
         col.regions = colorRampPalette(brewer.pal(9,'Reds'))(101),
         cuts        = 100,
         colorkey=list(space="bottom", space = "bottom"),
         par.settings =
           list(axis.line = list(col =  'transparent',
                                 legend.ticks = 'black'))),
  spplot(london.gen, c("x2"),
         main = c("index of social fragmentation"),
         col.regions = colorRampPalette(brewer.pal(9,'Reds'))(101),
         cuts        = 100,
         colorkey=list(space="bottom", space = "bottom"),
         par.settings =
           list(axis.line = list(col =  'transparent',
                                 legend.ticks = 'black'))),
  ncol = 2)
dev.off()


### --- 3. Defining neighbor relation --- ####
temp <- poly2nb(london.gen)



### ----- 3.1. Plotting the neighbors --- ####
plot_map_neig <- function(neig)
{
  plot(london.gen)
  plot(london.gen[neig, ], border="white", 
       col="red", add=TRUE)
  
  plot(london.gen[temp[[neig]], ], 
       border="white", 
       col="pink", add=TRUE)
  london.gen[temp[[neig]],]$NAME
}

plot_map_neig(30)




### --- 4. Preparing data to fit the model using OPENBUGS --- ####
### ----- 4.1. Creating openbugs structure of the graph --- ####
london.gen_bugs <- nb2WB(temp)


### ----- 4.2. Preparing the data for bugs --- ####
data_bugs <- list(adj = london.gen_bugs$adj, 
                   w   = london.gen_bugs$weights, 
                   num = london.gen_bugs$num,
                   x1  = data$x1,
                   x2  = data$x2,
                   O   = data$y, 
                   E   = data$E, 
                   n   = length(data$y))

### ----- 4.3. Initial values for OPENbugs --- ####
ini <- function(){list(beta0 = runif(1, -10, 10), 
                       beta1 = rnorm(1, 0, 1),
                       sdsp  = runif(1, 0 , 10))}


### ----- 4.4. Parameters --- ####
parameters<- c("beta0", "beta1", "sdsp", "sp", "smr_prob", "smr_R")

#parameters<- c("beta0", "beta1", "R", "prob", "sdsp", "sp")


### ----- 4.5. Defining the model in BUGS code --- ####
model.file <- file.path(tempdir(), "model.txt") 
cat("model{ 
  for(i in 1:n){ 
    O[i]~ dpois(mu[i]) 
    log(mu[i]) <- log(E[i]) + beta0 + beta1*x1[i] + sp[i] 
    smr_R[i] <- exp(beta0 + beta1*x1[i] + sp[i]) 
    smr_prob[i] <- step(smr_R[i]-1)
  }
  
  sp[1:n] ~ car.normal(adj[], w[], num[], precsp) 
  beta0 ~ dflat() 
  beta1 ~ dnorm(0, 0.0001)
  precsp <- pow(sdsp,-2) 
  sdsp ~ dunif(0,10)
  }", file = model.file)


### ----- 4.6. Calling bug from R --- ####
mod1 <- bugs(data       = data_bugs, 
             inits      = ini, 
             parameters = parameters, 
             model.file = model.file, 
             n.chains   = 3, 
             n.iter     = 10000, 
             n.burnin   = 2000, 
             n.thin     = 3, #
             debug      = TRUE, 
             DIC        = TRUE,
             codaPkg    = TRUE)


### ----- 4.7. Checking chains --- ####

### ----- 4.8. Posterior distributions of the fixed effects --- ####
mod1_mcmc <- read.bugs(mod1)
saveRDS(mod1_mcmc, "rds/model1_smr.RDS")
mod1_mcmc <- readRDS("rds/model1_smr.RDS")
summary(mod1_mcmc)[[2]] %>% round(., 3)%>% .[1:2,]

mod1_mcmc
### --- chains --- ####
png("images/chains_smr.png", width = 1000, height = 600, res = 120)
lattice::xyplot(mod1_mcmc[,c(1,2)])
dev.off()

### --- Posterior distributions --- ####
png("images/chains_dens_smr.png", width = 1400, height = 700, res = 120)
plot(mod1_mcmc[,1:2])
dev.off()

#Thre chains work, so, we can compute the P(beta1>0)
(mod1_mcmc[[1]][,"beta0"] > 0) %>% table(.) %>%  c(.)/length(mod1_mcmc[[1]][,"beta0"])
(mod1_mcmc[[1]][,"beta1"] > 0) %>% table(.) %>%  c(.)/length(mod1_mcmc[[1]][,"beta1"])


png("images/autocor.png", width = 1400, height = 1000, res = 120)
acfplot(mod1_mcmc[,1:2], lag.max=20)
dev.off()


### ----- 4.8. Posterior distributions of the hyperpars --- ####
### --- chains --- ####
png("images/chains_dens_smr_hyper.png", width = 1000, height = 400, res = 120)
plot(mod1_mcmc[,4])
dev.off()

summary(mod1_mcmc)[[2]] %>% round(., 3)%>% .[4,]


png("images/autocor.png", width = 1400, height = 1000, res = 120)
acfplot(mod1_mcmc[,1:2], lag.max=20)
dev.off()

### ----- 4.9. Posterior distribution of the random effects --- ####
summary(mod1_mcmc)[[1]][rownames(summary(mod1_mcmc)[[1]])]

summary(mod1_mcmc)[[1]] %>% 
  as.data.frame(.) %>%
  cbind(rowname = rownames(.), .) %>% 
  dplyr::filter(stringr::str_detect(rowname, "^sp")) %>%
  dplyr::select(Mean, SD) -> sp_post



london.gen$SPmean <- sp_post$Mean
london.gen$SPsd <- sp_post$SD

png("images/random_effect_spatial_besag.png", width = 1500, height = 700, 
    res = 100)
grid.arrange(
  spplot(london.gen, c("SPmean"),
         main = c("Mean posterior of S"),
         #col.regions = rev(viridis_pal(option = "B")(101)),
         col.regions = colorRampPalette(brewer.pal(9,'Blues'))(101),
         cuts        = 100,
         colorkey=list(space = "bottom", space = "bottom"),
         par.settings =
           list(axis.line = list(col =  'transparent',
                                 legend.ticks = 'black'))),
  spplot(london.gen, c("SPsd"),
         main = c("Sd posterior of S"),
         col.regions = colorRampPalette(brewer.pal(9,'Blues'))(101),
         cuts        = 100,
         colorkey=list(space="bottom", space = "bottom"),
         par.settings =
           list(axis.line = list(col =  'transparent',
                                 legend.ticks = 'black'))),
  ncol = 2)
dev.off()



### ----- 4.9. Posterior distribution of suicides mortality --- ####
summary(mod1_mcmc)[[1]] %>% 
  as.data.frame(.) %>%
  cbind(rowname = rownames(.), .) %>% 
  dplyr::filter(stringr::str_detect(rowname, "^smr_R")) %>%
  dplyr::select(Mean, SD) -> smr_R_post

summary(mod1_mcmc)[[1]] %>% 
  as.data.frame(.) %>%
  cbind(rowname = rownames(.), .) %>% 
  dplyr::filter(stringr::str_detect(rowname, "^smr_prob")) %>%
  dplyr::select(Mean, SD) -> smr_prob


london.gen$SMR_mean <- smr_R_post$Mean  # mean
london.gen$SMR_sd <- smr_R_post$SD #s
london.gen$SMR_p1 <- smr_prob$Mean # probability to be greater than 1

png("images/SMR_besag.png", width = 1500, height = 600, 
    res = 100)
grid.arrange(spplot(london.gen,
                    c("SMR_mean"),
                    col.regions = colorRampPalette(brewer.pal(9,'Blues'))(101),
                    cuts         = 100,
                    main        = "SMR mean ",
                    colorkey=list(space="bottom"),
                    par.settings =
                      list(axis.line = list(col =  'transparent'))),
             spplot(london.gen,
                    c("SMR_sd"),
                    col.regions = colorRampPalette(brewer.pal(9,'Blues'))(101),
                    cuts         = 100,
                    main        = "SMR sd ",
                    colorkey=list(space="bottom"),
                    par.settings =
                      list(axis.line = list(col =  'transparent'))), ncol = 2)

dev.off()


### ----- 4.9.1. Posterior distribution of suicides SMR with cutoff--- ####
## Also, the probability for SMR to be greater than 1.
SMR.cutoff<- c(0.6, 0.9, 1.0, 1.1,  1.8)
SMR_p1.cutoff <- c(0,0.2,0.8,1)

SMR_disc = cut(london.gen$SMR_mean,
               breaks         = SMR.cutoff,
               include.lowest = TRUE)

SMR_p1_disc = cut(london.gen$SMR_p1,
                  breaks         = SMR_p1.cutoff,
                  include.lowest = TRUE)


london.gen$SMR_disc <- SMR_disc
london.gen$SMR_p1_disc <- SMR_p1_disc

png("images/SMR_besag_disc.png", width = 1500, height = 600, res = 150)
grid.arrange(spplot(london.gen,
                    c("SMR_disc"),
                    col.regions = brewer.pal(9,'Blues')[c(2,4,6,8)],
                    main        = "SMR ",
                    par.settings =
                      list(axis.line = list(col =  'transparent'))),
             spplot(london.gen,
                    c("SMR_p1_disc"),
                    col.regions = brewer.pal(9,'Blues')[c(3,6,9)],
                    main        = "p(SMR > 1) ",
                    par.settings =
                      list(axis.line = list(col =  'transparent'))), ncol = 2)
dev.off()











### --- 5. Fitting the model using INLA bym effect --- ####
#This create a file called ``LDN.graph'' with the graph for INLA
nb2INLA("data/processed/london_suic/LDN.graph", temp)

### ----- 5.1. Plotting the generated graph --- ####
H <- inla.read.graph(filename="data/processed/london_suic/LDN.graph")
image(inla.graph2matrix(H),xlab="",ylab="")

### ----- 5.2. Adding ids for the random effects --- ####
S <- U <- seq(1,32)
data <- cbind(data, S, U)


### ----- 5.3. Formula --- ####
formula <- y ~ 1 + f(S, 
                     model       = "besag", 
                     graph       = H,
                     scale.model = TRUE,
                     hyper       = 
                       list(prec = list(prior="loggamma",param = c(1,0.001)))) +
  f(U, 
    model       = "iid",
    hyper       = 
      list(prec = list(prior="loggamma",param = c(1,0.001))))

### ----- 5.4. Model --- ####
mod.suicides <- inla(formula,
                     family          = "poisson",
                     data            = data,
                     E               = E,
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.predictor = list(compute=TRUE, cdf=c(log(1))))

summary(mod.suicides)


### ----- 5.5. Posterior distribution of the random effects --- ####
london.gen$SPmean <- round(mod.suicides$summary.random$S[["mean"]], 4)
london.gen$SPsd <- round(mod.suicides$summary.random$S[["sd"]],5)

#png("images/random_effect_spatial_bym.png", width = 1500, height = 700, 
#    res = 100)
grid.arrange(
  spplot(london.gen, c("SPmean"),
         main = c("Mean posterior of S"),
         #col.regions = rev(viridis_pal(option = "B")(101)),
         col.regions = colorRampPalette(brewer.pal(9,'Blues'))(101),
         cuts        = 100,
         colorkey=list(space="bottom", space = "bottom"),
         par.settings =
           list(axis.line = list(col =  'transparent',
                                 legend.ticks = 'black'))),
  spplot(london.gen, c("SPsd"),
         main = c("Sd posterior of S"),
         col.regions = colorRampPalette(brewer.pal(9,'Blues'))(101),
         cuts        = 100,
         colorkey=list(space="bottom", space = "bottom"),
         par.settings =
           list(axis.line = list(col =  'transparent',
                                 legend.ticks = 'black'))),
  ncol = 2)
#dev.off()

### ----- 5.6. Posterior distribution of suicides mortality --- ####
london.gen$SMR_mean <- mod.suicides$summary.fitted.values$mean # mean
london.gen$SMR_sd <- mod.suicides$summary.fitted.values$sd #s
london.gen$SMR_median <- mod.suicides$summary.fitted.values$`0.5quant` # median
london.gen$SMR_q025 <- mod.suicides$summary.fitted.values$`0.025quant` # quantile
london.gen$SMR_q975 <- mod.suicides$summary.fitted.values$`0.975quant` # quantile
london.gen$SMR_p1 <- 1 - mod.suicides$summary.fitted.values$`1 cdf` # probability to be greater than 1

#png("images/SMR_bym.png", width = 1500, height = 600, 
#    res = 100)
grid.arrange(spplot(london.gen,
                    c("SMR_mean"),
                    col.regions = colorRampPalette(brewer.pal(9,'Blues'))(101),
                    cuts         = 100,
                    main        = "SMR mean ",
                    colorkey=list(space="bottom"),
                    par.settings =
                      list(axis.line = list(col =  'transparent'))),
             spplot(london.gen,
                    c("SMR_sd"),
                    col.regions = colorRampPalette(brewer.pal(9,'Blues'))(101),
                    cuts         = 100,
                    main        = "SMR sd ",
                    colorkey=list(space="bottom"),
                    par.settings =
                      list(axis.line = list(col =  'transparent'))), ncol = 2)

#dev.off()


### ----- 5.7. Posterior distribution of suicides SMR with cutoff--- ####
## Also, the probability for SMR to be greater than 1.
SMR.cutoff<- c(0.6, 0.9, 1.0, 1.1,  1.8)
SMR_p1.cutoff <- c(0,0.2,0.8,1)

SMR_disc = cut(london.gen$SMR_mean,
               breaks         = SMR.cutoff,
               include.lowest = TRUE)

SMR_p1_disc = cut(london.gen$SMR_p1,
                  breaks         = SMR_p1.cutoff,
                  include.lowest = TRUE)


london.gen$SMR_disc <- SMR_disc
london.gen$SMR_p1_disc <- SMR_p1_disc

#png("images/SMR_bym_disc.png", width = 1500, height = 600, res = 150)
grid.arrange(spplot(london.gen,
                    c("SMR_disc"),
                    col.regions = brewer.pal(9,'Blues')[c(2,4,6,8)],
                    main        = "SMR ",
                    par.settings =
                      list(axis.line = list(col =  'transparent'))),
             spplot(london.gen,
                    c("SMR_p1_disc"),
                    col.regions = brewer.pal(9,'Blues')[c(3,6,9)],
                    main        = "p(SMR > 1) ",
                    par.settings =
                      list(axis.line = list(col =  'transparent'))), ncol = 2)

#dev.off()

### --- 5. Excercise --- ###
#Now, add the covariates (deprivation - x1 and social fragmentation - x2) and 
#repeat the steps