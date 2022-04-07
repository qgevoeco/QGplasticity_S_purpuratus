# Be sure to set the working directory to wherever data and models are stored
# setwd()


# Load a few packages that will be used below:
library(MCMCglmm)
library(lmerTest) #library(lme4)
library(Rmisc)
library(ggplot2)
library(gridExtra)

# use `remotes` package to install `wolakR` for posterior summary functions:
## remotes::install_github("matthewwolak/wolakR")
library(wolakR)
library(wesanderson)  #<-- contains color palette to use

# Load data
parN <- read.table("data_parN.txt", header = TRUE)
parU <- read.table("data_parU.txt", header = TRUE)

  # Set factors
  parN <- within(parN, {
    treat_dev <- as.factor(treat_dev)
    animal <- as.factor(animal)
    cultureFac <- as.factor(cultureFac)
    Dam <- as.factor(Dam)
    Sire <- as.factor(Sire)
    cross <- as.factor(cross)
  })

  parU <- within(parU, {
    treat_dev <- as.factor(treat_dev)
    animal <- as.factor(animal)
    cultureFac <- as.factor(cultureFac)
    Dam <- as.factor(Dam)
    Sire <- as.factor(Sire)
    cross <- as.factor(cross)
  })

# Percent abnormalities
percAb_parN <- read.table("percAb_parN.txt", header = TRUE)
  percAb_parN <- within(percAb_parN, {
    culture <- as.factor(culture)
    treat_dev <- as.factor(treat_dev)
    blockFac <- as.factor(blockFac)
  })
  
percAb_parU <- read.table("percAb_parU.txt", header = TRUE)
  percAb_parU <- within(percAb_parU, {
    culture <- as.factor(culture)
    treat_dev <- as.factor(treat_dev)
    blockFac <- as.factor(blockFac)
  })

# Egg diameter data
eggs <- read.table("data_eggs.txt", header = TRUE)

  # Set factors
  eggs <- within(eggs, {
    treat_adult <- as.factor(treat_adult)
    tube <- as.factor(tube)
    blockFac <- as.factor(block)
  })
  



# Load the models
## See file "analyses_QGplasticity_S_purpuratus.R" for code to generate models
## Assumes models are in the directory that contains this GitHub repository

### Spicule Length ######################
load(file = "../spiMod_parN_nit2200k.rdata")
  # replace 1 copy of covariance with the correlation
  spiMod_parN$VCV[, 2] <- posterior.cor(spiMod_parN$VCV[, 1:4])[, 2]
load(file = "../spiMod_parU_nit2200k.rdata")
  # replace 1 copy of covariance with the correlation
  spiMod_parU$VCV[, 2] <- posterior.cor(spiMod_parU$VCV[, 1:4])[, 2]


### Body Length #########################
load(file = paste0("../bodMod_parN_nit2130k.rdata"))
  # replace 1 copy of covariance with the correlation
  bodMod_parN$VCV[, 2] <- posterior.cor(bodMod_parN$VCV[, 1:4])[, 2]
load(file = paste0("../bodMod_parU_nit2130k.rdata"))
  # replace 1 copy of covariance with the correlation
  bodMod_parU$VCV[, 2] <- posterior.cor(bodMod_parU$VCV[, 1:4])[, 2]



# Create a few colors for plots
clP16 <- wes_palette("GrandBudapest2", n = 16, type = "continuous") 
clPurp <- clP16[c(1:6, 15:16)]
  class(clPurp) <- class(clP16)
  attr(clPurp, "name") <- "S. purpuratus"

  NNcl <- clPurp[2]  #<-- "Pink"
  NUcl <- clPurp[6]  #<-- "Light Blue"
  UNcl <- clPurp[4]  #<-- "Lavender"
  UUcl <- clPurp[8]  #<-- "Blue"


################################################################################

# 				RESULTS	

################################################################################

##########################################
#####   Eggs Diameter               ######
##########################################
# Egg diameter
glm_egg <- lmer(Average ~ treat_adult + (1 | blockFac) + (1 | blockFac:tube),
  data = eggs)
anova(glm_egg)


# Egg size predictor of body size
eggsSumm <- summarySE(eggs, measurevar = "Average", groupvars = c("tube"),
  na.rm = TRUE) 
  eggsSumm$Dam <- eggsSumm$tube
  eggsSumm$Dam <- sub("^", "dam_", eggsSumm$Dam ) #specifies dam for each number


PspiSumm <- summarySE(rbind(cbind(parN, treat_adult = "N"),
                            cbind(parU, treat_adult = "U")),
  measurevar = "Length.spi",
  groupvars = c("treat_adult", "treat_dev", "Dam"),
  na.rm = TRUE)

PbodSumm=summarySE(rbind(cbind(parN, treat_adult = "N"),
                            cbind(parU, treat_adult = "U")),
  measurevar = "Length.bod",
  groupvars = c("treat_adult", "treat_dev", "Dam"),
  na.rm = TRUE)

summPspi <- merge(PspiSumm, eggsSumm, by = "Dam")
summPspi$treat <- paste0(summPspi$treat_adult, summPspi$treat_dev)

summPbod <- merge(PbodSumm, eggsSumm, by = "Dam")
summPbod$treat <- paste0(summPbod$treat_adult, summPbod$treat_dev)


eggspicule_lm <- lm(Length.spi ~ Average, data = summPspi)
summary(eggspicule_lm)

eggbod_lm <- lm(Length.bod ~ Average, data = summPbod)
summary(eggbod_lm)










##########################################
##### Quantitative Genetic Models   ######
##########################################

# Tables summarizing posterior distributions of parental treatment fixed effects
##XXX Note, models setup WITHOUT an intercept
## `postTable()` lists posterior mode, mean, and HPD interval
postTable(spiMod_parN$Sol[, 1:2])
postTable(spiMod_parU$Sol[, 1:2])

## Compare larval spicule lengths (parent N vs. U) for upwelling development evironment
### what is probability that UU-NU difference is > 0.6 (6 microns)
mean((spiMod_parU$Sol[, "treat_devU"] - spiMod_parN$Sol[, "treat_devU"]) > 0.6)

postTable(bodMod_parN$Sol[, 1:2])
postTable(bodMod_parU$Sol[, 1:2])






# Spicule Length
## Calculate IAs and h2s
spih2 <- as.mcmc({
      cbind(spiParNdevN_h2 = with(spiMod_parN,
      		VCV[, 1] / rowSums(VCV[, c(1, 5:9)])),
            spiParNdevU_h2 = with(spiMod_parN,
            	VCV[, 4] / rowSums(VCV[, c(4, 5:8, 10)])),
	    spiParUdevN_h2 = with(spiMod_parU,
	    	VCV[, 1] / rowSums(VCV[, c(1, 5:9)])),
	    spiParUdevU_h2 = with(spiMod_parU,
	    	VCV[, 4] / rowSums(VCV[, c(4, 5:8, 10)])))
	})
## IAs use fixed effect treatment estimate as mean
### doing so marginalizes over other fixed effects (measurer)
#### include measurer in predicted value at average measurer covariate (0 by
#### centering and hence 0*estimated effect size)
spiIa <- as.mcmc({
      cbind(spiParNdevN_Ia = with(spiMod_parN, VCV[, 1] / (Sol[, "treat_devN"]^2)),
            spiParNdevU_Ia = with(spiMod_parN, VCV[, 4] / (Sol[, "treat_devU"]^2)),
	    spiParUdevN_Ia = with(spiMod_parU, VCV[, 1] / (Sol[, "treat_devN"]^2)),
	    spiParUdevU_Ia = with(spiMod_parU, VCV[, 4] / (Sol[, "treat_devU"]^2)))
	})







# Body Length
## Calculate CVAs and h2s
bodh2 <- as.mcmc({
      cbind(bodParNdevN_h2 = with(bodMod_parN,
      		VCV[, 1] / rowSums(VCV[, c(1, 5:9)])),
            bodParNdevU_h2 = with(bodMod_parN,
            	VCV[, 4] / rowSums(VCV[, c(4, 5:8, 10)])),
	    bodParUdevN_h2 = with(bodMod_parU,
	    	VCV[, 1] / rowSums(VCV[, c(1, 5:9)])),
	    bodParUdevU_h2 = with(bodMod_parU,
	    	VCV[, 4] / rowSums(VCV[, c(4, 5:8, 10)])))
	})
## IAs use fixed effect treatment estimate as mean
### doing so marginalizes over other fixed effects (measurer)
#### include measurer in predicted value at average measurer covariate (0 by
#### centering and hence 0*estimated effect size)
bodIa <- as.mcmc({
      cbind(bodParNdevN_Ia = with(bodMod_parN, VCV[, 1] / (Sol[, "treat_devN"]^2)),
            bodParNdevU_Ia = with(bodMod_parN, VCV[, 4] / (Sol[, "treat_devU"]^2)),
	    bodParUdevN_Ia = with(bodMod_parU, VCV[, 1] / (Sol[, "treat_devN"]^2)),
	    bodParUdevU_Ia = with(bodMod_parU, VCV[, 4] / (Sol[, "treat_devU"]^2)))
	})





#################################
# Parameter Estimate Differences
#################################

# Spicule length
spiNUminNN_VAdiff <- spiMod_parN$VCV[, "treat_devU:treat_devU.animal"] -
  spiMod_parN$VCV[, "treat_devN:treat_devN.animal"]
spiUUminUN_VAdiff <- spiMod_parU$VCV[, "treat_devU:treat_devU.animal"] -
  spiMod_parU$VCV[, "treat_devN:treat_devN.animal"]

  postTable(spiNUminNN_VAdiff)
    mean(spiNUminNN_VAdiff > 0.1)
  postTable(spiUUminUN_VAdiff)
    mean(spiUUminUN_VAdiff > 0.1)


spiNUminNN_Iadiff <- spiIa[, "spiParNdevU_Ia"] - spiIa[, "spiParNdevN_Ia"]
spiUUminUN_Iadiff <- spiIa[, "spiParUdevU_Ia"] - spiIa[, "spiParUdevN_Ia"]
  postTable(spiNUminNN_Iadiff)
    mean(spiNUminNN_Iadiff > 0.0025) 
  postTable(spiUUminUN_Iadiff)
    mean(spiUUminUN_Iadiff > 0.0025) 


spiNUminUU_Iadiff <- spiIa[, "spiParNdevU_Ia"] - spiIa[, "spiParUdevU_Ia"]
  postTable(spiNUminUU_Iadiff)
    mean(spiNUminUU_Iadiff > 0.0025) 


##############################################
# Body size differences 
## in cross-environment genetic correlation
##############################################

# How much of posterior density is negative in Non-Upwelling
mean(bodMod_parN$VCV[, 2] < 0)

# How much of posterior density is positive in Upwelling 
mean(bodMod_parU$VCV[, 2] > 0)





#############################
# Variance Component Tables
#############################
## `postTable()` lists posterior mode, mean, and HPD interval

# Spicule Length
postTable(spiMod_parN$VCV)
postTable(spiMod_parU$VCV)
  postTable(spih2)
  postTable(spiIa)
  

# Body Length
postTable(bodMod_parN$VCV)
postTable(bodMod_parU$VCV)
  postTable(bodh2)
  postTable(bodIa)












#########################
# Simulate Priors
#########################
  
set.seed(10001)
Nsim <- 10000    # Number of prior distribution samples
# Simulate the prior for additive genetic (co)variances
## nu=k+1
simAddGenPrior <- rpeIW(Nsim, V = diag(2)*0.02, nu = 3,
		alpha.mu = rep(0, 2), alpha.V = diag(2)*1000)

# Simulate a 1-dimension prior for all univariate varcomps that aren't VA
simUniPrior <- replicate(ncol(spiMod_parN$VCV) - 6,
  rpeIW(Nsim, V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1)*1000))

# Simulate 2-dimension residual covariance prior BUT DROP THE COVARIANCES
simResidPrior <- rIW(V = diag(2), nu = 2, n = Nsim)

# Create a prior distribution for each of the results
simPriors <- as.mcmc({cbind(
  VA1 = simAddGenPrior[, 1], VA2 = simAddGenPrior[, 4],
  COVA12 = simAddGenPrior[, 2],
  V1 = simUniPrior[, 1], V2 = simUniPrior[, 2],
    V3 = simUniPrior[, 3], V4 = simUniPrior[, 4],
  VR1 = simResidPrior[, 1], VR2 = simResidPrior[, 4])
  })


      
# Now create some combined statistic priors
h2Prior <- as.mcmc({
  simPriors[, "VA1"] / rowSums(simPriors[, c("VA1",
                                                "V1", "V2", "V3", "V4", "VR1")])
})

IaPrior <- as.mcmc(simPriors[, "VA1"] / (rnorm(nrow(simPriors), 0, sqrt(1e10))^2))

 
 
 
 
 
 
 
 
 
 
 
    
################################################################################

# 				FIGURES

################################################################################



####################

# Figure 2   

####################
tiff("../Fig2.tiff", width = 8, height = 5, units = "in",
  res = 500, compression = "jpeg")          

## Model Treatment Mean Reaction Norms
par(mfrow = c(1, 2), mar = c(5, 5, 2, 0.5), cex.lab = 1.5, cex.axis = 1.25)
# SPICULE LENGTH
## Note multiplying all by 10 to convert to micrometers
### (millimeter scale of response already multiplied by 100)
tmpSpiModLst <- as.mcmc(10 * cbind(NN = spiMod_parN$Sol[, "treat_devN"],
      UN = spiMod_parU$Sol[, "treat_devN"],
      NU = spiMod_parN$Sol[, "treat_devU"],
      UU = spiMod_parU$Sol[, "treat_devU"]))
xaxs <- c(0.9,1.1, 2.9,3.1)
plot(x = rep(xaxs, 2), y = apply(tmpSpiModLst, MARGIN = 2, FUN = range),
  axes = FALSE, type = "n",
  xlim = c(0.5, 3.5), ylim = c(80, 132.5),
  xlab = "Larval Environment", ylab = ~Spicule~Length~(mu*m))
 tmpSpiModMean <- apply(tmpSpiModLst, MARGIN = 2, FUN = mean, na.rm = TRUE)
 tmpSpiModCI <- HPDinterval(tmpSpiModLst)
 # lines/reaction norms
 lines(x = c(0.9, 2.9), y = tmpSpiModMean[c(1,3)], lwd = 3)
 lines(x = c(1.1, 3.1), y = tmpSpiModMean[c(2,4)], lwd = 3)
 # CIs
 arrows(x0 = xaxs,
   y0 = tmpSpiModCI[, "lower"], y1 = tmpSpiModCI[, "upper"],
   length = 0.075, angle = 90, code = 3, col = c(NNcl, UNcl, NUcl, UUcl), lwd = 3)   
 # means
 points(x = xaxs, y = tmpSpiModMean,
   pch = c(21, 24, 21, 24), bg = c(NNcl, UNcl, NUcl, UUcl), cex = 2.2, lwd = 2)
 # modes
# points(x = xaxs, y = posterior.mode(tmpSpiModLst),
#   pch = 8, col = c(NNcl, UNcl, NUcl, UUcl), cex = 2)
   
axis(1, at = c(1,3), labels = c("Non-Upwelling", "Upwelling"))
axis(2, at = seq(80, 130, 10))

mtext(text = expression((bolditalic(a))),
  side = 3, line = -0.4, adj = -0.3, cex = 1.3)

# Add Parent and Larvae Treatment letters/labels next to points
text(x = xaxs + rep(c(-1, 1) * 0.275, 2), y = tmpSpiModMean,
  labels = names(tmpSpiModMean))



# BODY LENGTH
## Note multiplying all by 10 to convert to micrometers
### (millimeter scale of response already multiplied by 100)
tmpBodModLst <- as.mcmc(10 * cbind(NN = bodMod_parN$Sol[, "treat_devN"],
      UN = bodMod_parU$Sol[, "treat_devN"],
      NU = bodMod_parN$Sol[, "treat_devU"],
      UU = bodMod_parU$Sol[, "treat_devU"]))
xaxs <- c(0.9,1.1, 2.9,3.1)
plot(x = rep(xaxs, 2), y = apply(tmpBodModLst, MARGIN = 2, FUN = range),
  axes = FALSE, type = "n",
  xlim = c(0.5, 3.5), ylim = c(125, 157.5),
  xlab = "Larval Environment", ylab = ~Body~Length~(mu*m))
 tmpBodModMean <- apply(tmpBodModLst, MARGIN = 2, FUN = mean, na.rm = TRUE)
 tmpBodModCI <- HPDinterval(tmpBodModLst)
 # lines/reaction norms
 lines(x = c(0.9, 2.9), y = tmpBodModMean[c(1,3)], lwd = 3)
 lines(x = c(1.1, 3.1), y = tmpBodModMean[c(2,4)], lwd = 3)
 # Std. Dev.
 arrows(x0 = xaxs,
   y0 = tmpBodModCI[, "lower"], y1 = tmpBodModCI[, "upper"],
   length = 0.075, angle = 90, code = 3, col = c(NNcl, UNcl, NUcl, UUcl), lwd = 3)   
 # means
 points(x = xaxs, y = tmpBodModMean,
   pch = c(21, 24, 21, 24), bg = c(NNcl, UNcl, NUcl, UUcl), cex = 2.2, lwd = 2)
 # modes
# points(x = xaxs, y = posterior.mode(tmpBodModLst),
#    pch = 8, col = c(NNcl, UNcl, NUcl, UUcl), cex = 2)
   
axis(1, at = c(1,3), labels = c("Non-Upwelling", "Upwelling"))
axis(2, at = seq(125, 155, 5))
  
mtext(text = expression((bolditalic(b))),
  side = 3, line = -0.4, adj = -0.3, cex = 1.3)
  
# Add Parent and Larvae Treatment letters/labels next to points
text(x = xaxs + rep(c(-1, 1) * 0.275, 2), y = tmpBodModMean,
  labels = names(tmpBodModMean))


dev.off()











######################################

# spicule Length
## Fgure 3

######################################
tmpModN <- spiMod_parN
tmpModU <- spiMod_parU
tmph2 <- spih2
tmpIa <- spiIa


tiff("../Fig3.tiff", width = 9, height = 10, units = "in",
  res = 500, compression = "jpeg")          
par(mfcol = c(4, 3), oma = c(1, 7, 4, 0),
  mar = c(4.5,4.6,3,1.5), mgp = c(2.5, 1, 0), cex.lab = 1.55, cex.axis = 1.25)    

lettLine <- 1.5  #<-- vertical line placement of panel letter
lettAdj <- -0.3 #<-- horizonal adjustment of panel letter
#################
# VAs
#################
genXlim1in <- c(0, 2.0)
genXlim2in <- c(0, 3.5)
genYlim1in <- c(0, 5)
genYlim2in <- c(0, 1.2)
 ################### 
 # Parent N  #####
  NVA1 <- postPlot(tmpModN$VCV[, 1], plotHist = TRUE,
	xlim = genXlim1in, ylim = genYlim1in,
	main = "",
	xlab = "", ylab = "Density",
	histcol = NNcl)
	# estimate prior densities at posterior histogram midpoints
	poh <- NVA1$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(a))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise ylabel increases 1 line each plot??
  NVA2 <- postPlot(tmpModN$VCV[, 4], plotHist = TRUE,
	xlim = genXlim2in, ylim = genYlim2in,
	main = "",
	xlab = "", ylab = "Density",
	histcol = NUcl)
	# estimate prior densities at posterior histogram midpoints
	poh <- NVA2$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(b))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
 ################### 
 # Parent U  #####
par(mgp = c(2.5, 1, 0))
  UVA1 <- postPlot(tmpModU$VCV[, 1], plotHist = TRUE,
	xlim = genXlim1in, ylim = genYlim1in,
	main = "",
	xlab = "", ylab = "Density",
	histcol = UNcl)
	# estimate prior densities at posterior histogram midpoints
	poh <- UVA1$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(c))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)

par(mgp = c(2.5, 1, 0))
  UVA2 <- postPlot(tmpModU$VCV[, 4], plotHist = TRUE,
	xlim = genXlim2in, ylim = genYlim2in,
	main = "",
	xlab = expression(V[ A ]), ylab = "Density",
	histcol = UUcl)
	# estimate prior densities at posterior histogram midpoints
	poh <- UVA2$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(d))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
mtext(text = "Additive genetic variance", #    expression(V[ A ])
     outer = TRUE, side = 3, line = 1, adj = 0.05, cex = 1.2)

#################
# h2s
#################
ylimin1 <- ylimin2 <- c(0, 11)
 ################### 
 # Parent N  #####
  postPlot(tmph2[, 1], #denscol = clPurp[1],
	main = "",
	xlab = "", ylab = "",
	histbreaks = 50, histcol = NNcl,
	prior = h2Prior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = c(0, 1), ylim = ylimin1)
   mtext(text = expression((bolditalic(e))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
  # Density approximation makes prior go past h2=1: COVER THIS UP
  lines(x = c(1.0, 1.2), y = c(0, 0), lwd = 8, col = "white", lend = 2)

  postPlot(tmph2[, 2], #denscol = clPurp[1],
	main = "",
	xlab = "", ylab = "",
	histbreaks = 50, histcol = NUcl,
	prior = h2Prior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = c(0, 1), ylim = ylimin1)
   mtext(text = expression((bolditalic(f))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
  # Density approximation makes prior go past h2=1: COVER THIS UP
  lines(x = c(1.0, 1.2), y = c(0, 0), lwd = 8, col = "white", lend = 2)
 ################### 
 # Parent U  #####
  postPlot(tmph2[, 3], #denscol = clPurp[1],
	main = "",
	xlab = "", ylab = "",
	histbreaks = 50, histcol = UNcl,
	prior = h2Prior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = c(0, 1), ylim = ylimin2)
   mtext(text = expression((bolditalic(g))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
  # Density approximation makes prior go past h2=1: COVER THIS UP
  lines(x = c(1.0, 1.2), y = c(0, 0), lwd = 8, col = "white", lend = 2)
   
par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise xlab floats each time
  postPlot(tmph2[, 4], #denscol = clPurp[1],
	main = "",
	xlab = expression("h"^2), ylab = "",
	histbreaks = 50, histcol = UUcl,
	prior = h2Prior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = c(0, 1), ylim = ylimin2)
   mtext(text = expression((bolditalic(h))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
  # Density approximation makes prior go past h2=1: COVER THIS UP
  lines(x = c(1.0, 1.2), y = c(0, 0), lwd = 8, col = "white", lend = 2)
mtext(text = "heritability", # expression(paste("heritability"[   ]))
     outer = TRUE, side = 3, line = 1, adj = 0.5, cex = 1.2)

#################
# evolvabilities (Ia)
#################
xlimIn1 <- c(0, 0.015)
xlimIn2 <- c(0, 0.05)
ylimIn1 <- c(0, 650)
ylimIn2 <- c(0, 150)
 ################### 
 # Parent N  #####
  postPlot(tmpIa[, 1], #denscol = clPurp[1],
	main = "",
	xlab = "", ylab = "",
	histbreaks = 50, histcol = NNcl,
	prior = IaPrior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = xlimIn1, ylim = ylimIn1)
   mtext(text = expression((bolditalic(i))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)

  postPlot(tmpIa[, 2], #denscol = clPurp[1],
	main = "",
	xlab = "", ylab = "",
	histbreaks = 50, histcol = NUcl,
	prior = IaPrior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = xlimIn2, ylim = ylimIn2)
   mtext(text = expression((bolditalic(j))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
 ################### 
 # Parent U  #####
  postPlot(tmpIa[, 3], #denscol = clPurp[1],
	main = "",
	xlab = "", ylab = "",
	histbreaks = 50, histcol = UNcl,
	prior = IaPrior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = xlimIn1, ylim = ylimIn1)
   mtext(text = expression((bolditalic(k))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
   

par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise xlab floats each time
  postPlot(tmpIa[, 4], #denscol = clPurp[1],
	main = "",
	xlab = expression("I"[ A ]), ylab = "",
	histbreaks = 50, histcol = UUcl,
	prior = IaPrior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = xlimIn2, ylim = ylimIn2)
   mtext(text = expression((bolditalic(l))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
mtext(text = "Evolvability", # expression(paste("Evolvability (I"[ A ],")"))
     outer = TRUE, side = 3, line = 1, adj = 0.9, cex = 1.2)


mtext(text = "Parent Upwelling",
     outer = TRUE, side = 2, line = 5, adj = 0.2, cex = 1.5)
mtext(text = "Parent Non-Upwelling",
     outer = TRUE, side = 2, line = 5, adj = 0.85, cex = 1.5)
  mtext(text = "Larval Upwelling        Larval Non-Upwelling",
     outer = TRUE, side = 2, line = 2.5, adj = 0.93, cex = 1.2)
  mtext(text = "Larval Upwelling        Larval Non-Upwelling",
     outer = TRUE, side = 2, line = 2.5, adj = 0.1, cex = 1.2)

dev.off()



























######################################
# body Length
## Fgure 4
######################################
tmpModN <- bodMod_parN
tmpModU <- bodMod_parU
tmph2 <- bodh2
tmpIa <- bodIa

tiff("../Fig4.tiff", width = 9, height = 10, units = "in",
  res = 500, compression = "jpeg")          
par(mfcol = c(4, 3), oma = c(1, 7, 4, 0),
  mar = c(4.5,4.6,3,1.5), mgp = c(2.5, 1, 0), cex.lab = 1.55, cex.axis = 1.25)    

lettLine <- 1.5  #<-- vertical line placement of panel letter
lettAdj <- -0.3 #<-- horizonal adjustment of panel letter
#################
# VAs
#################
genXlim1in <- c(0, 1.4)
genXlim2in <- c(0, 1.0)
genYlim1in <- c(0, 7.5)
genYlim2in <- c(0, 20)
 ################### 
 # Parent N  #####
  NVA1 <- postPlot(tmpModN$VCV[, 1], plotHist = TRUE,
	xlim = genXlim1in, ylim = genYlim1in,
	main = "",
	xlab = "", ylab = "Density",
	histcol = NNcl)
	# estimate prior densities at posterior histogram midpoints
	poh <- NVA1$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(a))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise ylabel increases 1 line each plot??
  NVA2 <- postPlot(tmpModN$VCV[, 4], plotHist = TRUE,
	xlim = genXlim1in, ylim = genYlim1in,
	main = "",
	xlab = "", ylab = "Density",
	histcol = NUcl)
	# estimate prior densities at posterior histogram midpoints
	poh <- NVA2$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(b))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
 ################### 
 # Parent U  #####
par(mgp = c(2.5, 1, 0))
  UVA1 <- postPlot(tmpModU$VCV[, 1], plotHist = TRUE,
	xlim = genXlim2in, ylim = genYlim2in,
	main = "",
	xlab = "", ylab = "Density",
	histcol = UNcl)
	# estimate prior densities at posterior histogram midpoints
	poh <- UVA1$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(c))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)

par(mgp = c(2.5, 1, 0))
  UVA2 <- postPlot(tmpModU$VCV[, 4], plotHist = TRUE,
	xlim = genXlim2in, ylim = genYlim2in,
	main = "",
	xlab = expression(V[ A ]), ylab = "Density",
	histcol = UUcl)
	# estimate prior densities at posterior histogram midpoints
	poh <- UVA2$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(d))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
mtext(text = "Additive genetic variance", #    expression(V[ A ])
     outer = TRUE, side = 3, line = 1, adj = 0.05, cex = 1.2)

#################
# h2s
#################
ylimin1 <- c(0, 12)
ylimin2 <- c(0, 25)
 ################### 
 # Parent N  #####
  postPlot(tmph2[, 1], #denscol = clPurp[1],
	main = "",
	xlab = "", ylab = "",
	histbreaks = 50, histcol = NNcl,
	prior = h2Prior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = c(0, 1), ylim = ylimin1)
   mtext(text = expression((bolditalic(e))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
  # Density approximation makes prior go past h2=1: COVER THIS UP
  lines(x = c(1.0, 1.2), y = c(0, 0), lwd = 8, col = "white", lend = 2)

  postPlot(tmph2[, 2], #denscol = clPurp[1],
	main = "",
	xlab = "", ylab = "",
	histbreaks = 50, histcol = NUcl,
	prior = h2Prior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = c(0, 1), ylim = ylimin1)
   mtext(text = expression((bolditalic(f))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
  # Density approximation makes prior go past h2=1: COVER THIS UP
  lines(x = c(1.0, 1.2), y = c(0, 0), lwd = 8, col = "white", lend = 2)
 ################### 
 # Parent U  #####
  postPlot(tmph2[, 3], #denscol = clPurp[1],
	main = "",
	xlab = "", ylab = "",
	histbreaks = 50, histcol = UNcl,
	prior = h2Prior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = c(0, 1), ylim = ylimin2)
   mtext(text = expression((bolditalic(g))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
  # Density approximation makes prior go past h2=1: COVER THIS UP
  lines(x = c(1.0, 1.2), y = c(0, 0), lwd = 8, col = "white", lend = 2)
   
par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise xlab floats each time
  postPlot(tmph2[, 4], #denscol = clPurp[1],
	main = "",
	xlab = expression("h"^2), ylab = "",
	histbreaks = 50, histcol = UUcl,
	prior = h2Prior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = c(0, 1), ylim = ylimin2)
   mtext(text = expression((bolditalic(h))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
  # Density approximation makes prior go past h2=1: COVER THIS UP
  lines(x = c(1.0, 1.2), y = c(0, 0), lwd = 8, col = "white", lend = 2)
mtext(text = "heritability", # expression(paste("heritability"[   ]))
     outer = TRUE, side = 3, line = 1, adj = 0.5, cex = 1.2)

#################
# evolvabilities (Ia)
#################
xlimIn1 <- c(0, 0.0065)
xlimIn2 <- c(0, 0.0065)
ylimIn1 <- c(0, 1200)
ylimIn2 <- c(0, 3800)
 ################### 
 # Parent N  #####
  postPlot(tmpIa[, 1], #denscol = clPurp[1],
	main = "",
	xlab = "", ylab = "",
	histbreaks = 50, histcol = NNcl,
	prior = IaPrior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = xlimIn1, ylim = ylimIn1)
   mtext(text = expression((bolditalic(i))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)

  postPlot(tmpIa[, 2], #denscol = clPurp[1],
	main = "",
	xlab = "", ylab = "",
	histbreaks = 50, histcol = NUcl,
	prior = IaPrior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = xlimIn2, ylim = ylimIn1)
   mtext(text = expression((bolditalic(j))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
 ################### 
 # Parent U  #####
  postPlot(tmpIa[, 3], #denscol = clPurp[1],
	main = "",
	xlab = "", ylab = "",
	histbreaks = 50, histcol = UNcl,
	prior = IaPrior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = xlimIn1, ylim = ylimIn2)
   mtext(text = expression((bolditalic(k))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
   

par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise xlab floats each time
  postPlot(tmpIa[, 4], #denscol = clPurp[1],
	main = "",
	xlab = expression("I"[ A ]), ylab = "",
	histbreaks = 50, histcol = UUcl,
	prior = IaPrior, prange = "prior", priorcol = "grey50", priorlwd = 4,
	xlim = xlimIn2, ylim = ylimIn2)
   mtext(text = expression((bolditalic(l))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
mtext(text = "Evolvability", # expression(paste("Evolvability (I"[ A ],")"))
     outer = TRUE, side = 3, line = 1, adj = 0.9, cex = 1.2)


mtext(text = "Parent Upwelling",
     outer = TRUE, side = 2, line = 5, adj = 0.2, cex = 1.5)
mtext(text = "Parent Non-Upwelling",
     outer = TRUE, side = 2, line = 5, adj = 0.85, cex = 1.5)
  mtext(text = "Larval Upwelling        Larval Non-Upwelling",
     outer = TRUE, side = 2, line = 2.5, adj = 0.93, cex = 1.2)
  mtext(text = "Larval Upwelling        Larval Non-Upwelling",
     outer = TRUE, side = 2, line = 2.5, adj = 0.1, cex = 1.2)

dev.off()






















######################################
# cross-environment correlations
## Fgure 5
######################################

# Extract random effects
## First get names of all effects from column names of the Sol object
## Then search for additive genetic values (BV) which are prefixed by "animal"
## Match with order of animals from data (parN/parU)

### Spicule length
nms_spiModN <- dimnames(spiMod_parN$Sol)[[2]]
  slnAnimNms_spiModN <- sapply(strsplit(nms_spiModN, "animal."), FUN = "[", i = 2)
  animBVind_spiModN <- match(parN$animal, slnAnimNms_spiModN)
nms_spiModU <- dimnames(spiMod_parU$Sol)[[2]]
  slnAnimNms_spiModU <- sapply(strsplit(nms_spiModU, "animal."), FUN = "[", i = 2)
  animBVind_spiModU <- match(parU$animal, slnAnimNms_spiModU)

### Body length
nms_bodModN <- dimnames(bodMod_parN$Sol)[[2]]
  slnAnimNms_bodModN <- sapply(strsplit(nms_bodModN, "animal."), FUN = "[", i = 2)
  animBVind_bodModN <- match(parN$animal, slnAnimNms_bodModN)
nms_bodModU <- dimnames(bodMod_parU$Sol)[[2]]
  slnAnimNms_bodModU <- sapply(strsplit(nms_bodModU, "animal."), FUN = "[", i = 2)
  animBVind_bodModU <- match(parU$animal, slnAnimNms_bodModU)



# Go through each MCMC sample to calculate posterior of family mean values
## calc. each family's mean add. genetic value (BV) of larvae in each treatment

### Spicule length
#### parent N
spiModN_FamTrtMeanPost <- matrix(NA, nrow = nrow(spiMod_parN$Sol), ncol = 2*20)
for(i in 1:nrow(spiMod_parN$Sol)){
  # mean (across larvae) of each family-by-treatment combination
  i_FamTrtMean <- aggregate(spiMod_parN$Sol[i, animBVind_spiModN] ~ parN$cross +
    parN$treat_dev, FUN = mean)
  if(i == 1){
    dimnames(spiModN_FamTrtMeanPost) <- list(NULL, 
    	paste(as.character(i_FamTrtMean[, 1]), as.character(i_FamTrtMean[, 2]), 
    	  sep = "_"))
  }
  spiModN_FamTrtMeanPost[i, ] <- i_FamTrtMean[, 3]
}  

#### parentU
spiModU_FamTrtMeanPost <- matrix(NA, nrow = nrow(spiMod_parU$Sol), ncol = 2*20)
for(i in 1:nrow(spiMod_parU$Sol)){
  # mean (across larvae) of each family-by-treatment combination
  i_FamTrtMean <- aggregate(spiMod_parU$Sol[i, animBVind_spiModU] ~ parU$cross +
    parU$treat_dev, FUN = mean)
  if(i == 1){
    dimnames(spiModU_FamTrtMeanPost) <- list(NULL, 
    	paste(as.character(i_FamTrtMean[, 1]), as.character(i_FamTrtMean[, 2]), 
    	  sep = "_"))
  }
  spiModU_FamTrtMeanPost[i, ] <- i_FamTrtMean[, 3]
}  



### Body length
#### parent N
bodModN_FamTrtMeanPost <- matrix(NA, nrow = nrow(bodMod_parN$Sol), ncol = 2*20)
for(i in 1:nrow(bodMod_parN$Sol)){
  # mean (across larvae) of each family-by-treatment combination
  i_FamTrtMean <- aggregate(bodMod_parN$Sol[i, animBVind_bodModN] ~ parN$cross +
    parN$treat_dev, FUN = mean)
  if(i == 1){
    dimnames(bodModN_FamTrtMeanPost) <- list(NULL, 
    	paste(as.character(i_FamTrtMean[, 1]), as.character(i_FamTrtMean[, 2]), 
    	  sep = "_"))
  }
  bodModN_FamTrtMeanPost[i, ] <- i_FamTrtMean[, 3]
}  

#### parent U
bodModU_FamTrtMeanPost <- matrix(NA, nrow = nrow(bodMod_parU$Sol), ncol = 2*20)
for(i in 1:nrow(bodMod_parU$Sol)){
  # mean (across larvae) of each family-by-treatment combination
  i_FamTrtMean <- aggregate(bodMod_parU$Sol[i, animBVind_bodModU] ~ parU$cross +
    parU$treat_dev, FUN = mean)
  if(i == 1){
    dimnames(bodModU_FamTrtMeanPost) <- list(NULL, 
    	paste(as.character(i_FamTrtMean[, 1]), as.character(i_FamTrtMean[, 2]), 
    	  sep = "_"))
  }
  bodModU_FamTrtMeanPost[i, ] <- i_FamTrtMean[, 3]
}  


# Assign mcmc class and attributes so these will be treated as posteriors
class(spiModN_FamTrtMeanPost) <- class(spiModU_FamTrtMeanPost) <-
  class(bodModN_FamTrtMeanPost) <- class(bodModU_FamTrtMeanPost) <- "mcmc"

  attr(spiModN_FamTrtMeanPost, "mcpar") <- attr(spiMod_parN$Sol, "mcpar")
  attr(spiModU_FamTrtMeanPost, "mcpar") <- attr(spiMod_parU$Sol, "mcpar")
  attr(bodModN_FamTrtMeanPost, "mcpar") <- attr(bodMod_parN$Sol, "mcpar")
  attr(bodModU_FamTrtMeanPost, "mcpar") <- attr(bodMod_parU$Sol, "mcpar")



###################     DIVERSION/ASIDE   ######################################
# Posterior (1) mode versus (2) mean of family mean additive genetic value
## What single metric better represents the posterior probability?

### First look at the probability distributions
#### plot each family and development environment probability distribution

#### Spicule Length: Parent N - Larvae N
tiff("../FigESM5_S7.tiff", width = 12, height = 10, units = "in",
  res = 500, compression = "jpeg")          
par(mfrow = c(4, 5), mar = c(6, 6, 0, 0))
  ignoreOutput <- apply(spiModN_FamTrtMeanPost[, 1:20], MARGIN = 2,
    FUN = function(x) postPlot(as.mcmc(x), xlab = "", ylab = ""))  
  mtext(text = "Density",
     outer = TRUE, side = 2, line = -2.5, adj = 0.55, cex = 2.2)
  mtext(text = "Family mean additive genetic value",
     outer = TRUE, side = 1, line = -1.3, adj = 0.5, cex = 2.2)
dev.off()
#### Spicule Length: Parent N - Larvae U
tiff("../FigESM5_S8.tiff", width = 12, height = 10, units = "in",
  res = 500, compression = "jpeg")          
par(mfrow = c(4, 5), mar = c(6, 6, 0, 0))
  ignoreOutput <- apply(spiModN_FamTrtMeanPost[, 21:40], MARGIN = 2,
    FUN = function(x) postPlot(as.mcmc(x), xlab = "", ylab = ""))
  mtext(text = "Density",
     outer = TRUE, side = 2, line = -2.5, adj = 0.55, cex = 2.2)
  mtext(text = "Family mean additive genetic value",
     outer = TRUE, side = 1, line = -1.3, adj = 0.5, cex = 2.2)
dev.off()

#### Spicule Length: Parent U - Larvae N
tiff("../FigESM5_S9.tiff", width = 12, height = 10, units = "in",
  res = 500, compression = "jpeg")          
par(mfrow = c(4, 5), mar = c(6, 6, 0, 0))
  ignoreOutput <- apply(spiModU_FamTrtMeanPost[, 1:20], MARGIN = 2,
    FUN = function(x) postPlot(as.mcmc(x), xlab = "", ylab = ""))   
  mtext(text = "Density",
     outer = TRUE, side = 2, line = -2.5, adj = 0.55, cex = 2.2)
  mtext(text = "Family mean additive genetic value",
     outer = TRUE, side = 1, line = -1.3, adj = 0.5, cex = 2.2)
dev.off()
#### Spicule Length: Parent U - Larvae U
tiff("../FigESM5_S10.tiff", width = 12, height = 10, units = "in",
  res = 500, compression = "jpeg")          
par(mfrow = c(4, 5), mar = c(6, 6, 0, 0))
  ignoreOutput <- apply(spiModU_FamTrtMeanPost[, 21:40], MARGIN = 2,
    FUN = function(x) postPlot(as.mcmc(x), xlab = "", ylab = ""))   
  mtext(text = "Density",
     outer = TRUE, side = 2, line = -2.5, adj = 0.55, cex = 2.2)
  mtext(text = "Family mean additive genetic value",
     outer = TRUE, side = 1, line = -1.3, adj = 0.5, cex = 2.2)
dev.off()



#### Body Length: Parent N - Larvae N
tiff("../FigESM5_S11.tiff", width = 12, height = 10, units = "in",
  res = 500, compression = "jpeg")          
par(mfrow = c(4, 5), mar = c(6, 6, 0, 0))
  ignoreOutput <- apply(bodModN_FamTrtMeanPost[, 1:20], MARGIN = 2,
    FUN = function(x) postPlot(as.mcmc(x), xlab = "", ylab = ""))   
  mtext(text = "Density",
     outer = TRUE, side = 2, line = -2.5, adj = 0.55, cex = 2.2)
  mtext(text = "Family mean additive genetic value",
     outer = TRUE, side = 1, line = -1.3, adj = 0.5, cex = 2.2)
dev.off()
#### Body Length: Parent N - Larvae U
tiff("../FigESM5_S12.tiff", width = 12, height = 10, units = "in",
  res = 500, compression = "jpeg")          
par(mfrow = c(4, 5), mar = c(6, 6, 0, 0))
  ignoreOutput <- apply(bodModN_FamTrtMeanPost[, 21:40], MARGIN = 2,
    FUN = function(x) postPlot(as.mcmc(x), xlab = "", ylab = ""))   
  mtext(text = "Density",
     outer = TRUE, side = 2, line = -2.5, adj = 0.55, cex = 2.2)
  mtext(text = "Family mean additive genetic value",
     outer = TRUE, side = 1, line = -1.3, adj = 0.5, cex = 2.2)
dev.off()

#### Body Length: Parent U - Larvae N
tiff("../FigESM5_S13.tiff", width = 12, height = 10, units = "in",
  res = 500, compression = "jpeg")          
par(mfrow = c(4, 5), mar = c(6, 6, 0, 0))
  ignoreOutput <- apply(bodModU_FamTrtMeanPost[, 1:20], MARGIN = 2,
    FUN = function(x) postPlot(as.mcmc(x), xlab = "", ylab = ""))   
  mtext(text = "Density",
     outer = TRUE, side = 2, line = -2.5, adj = 0.55, cex = 2.2)
  mtext(text = "Family mean additive genetic value",
     outer = TRUE, side = 1, line = -1.3, adj = 0.5, cex = 2.2)
dev.off()
#### Body Length: Parent U - Larvae U
tiff("../FigESM5_S14.tiff", width = 12, height = 10, units = "in",
  res = 500, compression = "jpeg")          
par(mfrow = c(4, 5), mar = c(6, 6, 0, 0))
  ignoreOutput <- apply(bodModU_FamTrtMeanPost[, 21:40], MARGIN = 2,
    FUN = function(x) postPlot(as.mcmc(x), xlab = "", ylab = ""))   
  mtext(text = "Density",
     outer = TRUE, side = 2, line = -2.5, adj = 0.55, cex = 2.2)
  mtext(text = "Family mean additive genetic value",
     outer = TRUE, side = 1, line = -1.3, adj = 0.5, cex = 2.2)
dev.off()

#XXX Often very leptokurtic (though Body Length parent N show less kurtosis) and
## sometimes quite skewed. When posterior distributions show skew, posterior mean
## and modes may indicate different attributes of the distribution. In this case,
## a posterior mode gives the parameter value with highest posterior probability,
## so we will use that descriptive statistic to summarize the family mean additive
## genetic value. The posterior mode also describes the posterior probability
## with the highest posterior density (i.e., a very, very narrow credible interval)
## which is often a quite high probablity as the distribution is more leptokurtic
#################    END Diversion/Aside  ######################################

# Calculate posterior mode of mean family additive genetic value (BV)
## then rank them so can assign y-axis position

## Spicule length
### parent N
spiModN_FamMeanBV_devNmode <- apply(spiModN_FamTrtMeanPost[, 1:20], MARGIN = 2,
  FUN = function(x) posterior.mode(as.mcmc(x)))
spiModN_FamMeanBV_devUmode <- apply(spiModN_FamTrtMeanPost[, 21:40], MARGIN = 2,
  FUN = function(x) posterior.mode(as.mcmc(x)))
#### (-1 so rank 1 at "top" of figure)
spiModN_FamMeanBV_devNmode_rank <- -1*rank(spiModN_FamMeanBV_devNmode) 
spiModN_FamMeanBV_devUmode_rank <- -1*rank(spiModN_FamMeanBV_devUmode)

### parent U
spiModU_FamMeanBV_devNmode <- apply(spiModU_FamTrtMeanPost[, 1:20], MARGIN = 2,
  FUN = function(x) posterior.mode(as.mcmc(x)))
spiModU_FamMeanBV_devUmode <- apply(spiModU_FamTrtMeanPost[, 21:40], MARGIN = 2,
  FUN = function(x) posterior.mode(as.mcmc(x)))
#### (-1 so rank 1 at "top" of figure)
spiModU_FamMeanBV_devNmode_rank <- -1*rank(spiModU_FamMeanBV_devNmode)
spiModU_FamMeanBV_devUmode_rank <- -1*rank(spiModU_FamMeanBV_devUmode)

## Body length
### parent N
bodModN_FamMeanBV_devNmode <- apply(bodModN_FamTrtMeanPost[, 1:20], MARGIN = 2,
  FUN = function(x) posterior.mode(as.mcmc(x)))
bodModN_FamMeanBV_devUmode <- apply(bodModN_FamTrtMeanPost[, 21:40], MARGIN = 2,
  FUN = function(x) posterior.mode(as.mcmc(x)))
#### (-1 so rank 1 at "top" of figure)
bodModN_FamMeanBV_devNmode_rank <- -1*rank(bodModN_FamMeanBV_devNmode)
bodModN_FamMeanBV_devUmode_rank <- -1*rank(bodModN_FamMeanBV_devUmode)

### parent U 
bodModU_FamMeanBV_devNmode <- apply(bodModU_FamTrtMeanPost[, 1:20], MARGIN = 2,
  FUN = function(x) posterior.mode(as.mcmc(x)))
bodModU_FamMeanBV_devUmode <- apply(bodModU_FamTrtMeanPost[, 21:40], MARGIN = 2,
  FUN = function(x) posterior.mode(as.mcmc(x)))
#### (-1 so rank 1 at "top" of figure)
bodModU_FamMeanBV_devNmode_rank <- -1*rank(bodModU_FamMeanBV_devNmode)
bodModU_FamMeanBV_devUmode_rank <- -1*rank(bodModU_FamMeanBV_devUmode)
  
 










tiff("../Fig5.tiff", width = 12, height = 7, units = "in",
  res = 500, compression = "jpeg")          
mar1 <- c(5.5, 6, 5.5, 2.3)
mar2 <- c(5.5, 3, 5.5, 5.3)

par(mfrow = c(2, 4), mar = mar1,
  cex.lab = 1.55, cex.axis = 1.25)    

lettLine <- 1.4  #<-- vertical line placement of panel letter
lettAdj <- -0.3 #<-- horizonal adjustment of panel letter
lettAdj2 <- -0.2
# Create generic vector of line colors to connect N-U larval envrionment ranks
lineCols <- rep(c("black", "grey70"), each = 10)
#################
# reaction norms
#################
##############
# SPICULE LENGTH
##############
## Parent N
plot(x = rep(c(1,3), each = 20), y = rep(-1*seq(20), 2),
  axes = FALSE, type = "n",
  xlab = "", ylab = "Spicule length\nfamily mean genetic rank")

 for(i in 1:20){
  # lines/reaction norms
     ## map line colors onto ranking order according to the N larval environment
     lineColsO <- lineCols[-1 * spiModN_FamMeanBV_devNmode_rank]
   lines(x = c(1.1, 2.9),
     y = c(spiModN_FamMeanBV_devNmode_rank[i], spiModN_FamMeanBV_devUmode_rank[i]),
     lwd = 3, col = lineColsO[i])
  # modes
   points(x = c(1.1, 2.9), y = c(spiModN_FamMeanBV_devNmode_rank[i],
       spiModN_FamMeanBV_devUmode_rank[i]),
     pch = 21, bg = c(NNcl, NUcl), cex = 2)
 }

axis(1, at = c(1.1, 2.9), labels = c("Non-Upwelling", "Upwelling"))
mtext(text = expression((bolditalic(a))),
  side = 3, line = lettLine, adj = lettAdj, cex = 1.3)

par(mar = mar2)
## Parent U
plot(x = rep(c(1,3), each = 20), y = rep(-1*seq(20), 2),
  axes = FALSE, type = "n",
  xlab = "", ylab = "")

 for(i in 1:20){
  # lines/reaction norms
     ## map line colors onto ranking order according to the N larval environment
     lineColsO <- lineCols[-1 * spiModU_FamMeanBV_devNmode_rank]
   lines(x = c(1.1, 2.9),
     y = c(spiModU_FamMeanBV_devNmode_rank[i], spiModU_FamMeanBV_devUmode_rank[i]),
     lwd = 3, col = lineColsO[i])
  # modes
   points(x = c(1.1, 2.9), y = c(spiModU_FamMeanBV_devNmode_rank[i],
       spiModU_FamMeanBV_devUmode_rank[i]),
     pch = 21, bg = c(UNcl, UUcl), cex = 2)
 }

axis(1, at = c(1.1, 2.9), labels = c("Non-Upwelling", "Upwelling"))
mtext(text = expression((bolditalic(b))),
  side = 3, line = lettLine, adj = lettAdj2, cex = 1.3)


##############
# BODY LENGTH
##############
par(mar = mar1)
## Parent N
plot(x = rep(c(1,3), each = 20), y = rep(-1*seq(20), 2),
  axes = FALSE, type = "n",
  xlab = "", ylab = "Body length\nfamily mean genetic rank")

 for(i in 1:20){
  # lines/reaction norms
     ## map line colors onto ranking order according to the N larval environment
     lineColsO <- lineCols[-1 * bodModN_FamMeanBV_devNmode_rank]
   lines(x = c(1.1, 2.9),
     y = c(bodModN_FamMeanBV_devNmode_rank[i], bodModN_FamMeanBV_devUmode_rank[i]),
     lwd = 3, col = lineColsO[i])
  # modes
   points(x = c(1.1, 2.9), y = c(bodModN_FamMeanBV_devNmode_rank[i],
       bodModN_FamMeanBV_devUmode_rank[i]),
     pch = 21, bg = c(NNcl, NUcl), cex = 2)
 }

axis(1, at = c(1.1, 2.9), labels = c("Non-Upwelling", "Upwelling"))
mtext(text = expression((bolditalic(c))),
  side = 3, line = lettLine, adj = lettAdj, cex = 1.3)


## Parent U
par(mar = mar2)
plot(x = rep(c(1,3), each = 20), y = rep(-1*seq(20), 2),
  axes = FALSE, type = "n",
  xlab = "", ylab = "")

 for(i in 1:20){
  # lines/reaction norms
     ## map line colors onto ranking order according to the N larval environment
     lineColsO <- lineCols[-1 * bodModU_FamMeanBV_devNmode_rank]
   lines(x = c(1.1, 2.9),
     y = c(bodModU_FamMeanBV_devNmode_rank[i], bodModU_FamMeanBV_devUmode_rank[i]),
     lwd = 3, col = lineColsO[i])
  # modes
   points(x = c(1.1, 2.9), y = c(bodModU_FamMeanBV_devNmode_rank[i],
       bodModU_FamMeanBV_devUmode_rank[i]),
     pch = 21, bg = c(UNcl, UUcl), cex = 2)
 }

axis(1, at = c(1.1, 2.9), labels = c("Non-Upwelling", "Upwelling"))
mtext(text = expression((bolditalic(d))),
  side = 3, line = lettLine, adj = lettAdj2, cex = 1.3)

     
##########################
# CORa posteriors
##########################
ylimIn <- c(0, 1.2)
xlabExp <- expression("Genetic correlation"[ N-U])
################### 
 ################### 
 # Spicule length
 ###################
# Parent N  #####
par(mar = mar1)
  postPlot(spiMod_parN$VCV[, 2], plotHist = TRUE,
	xlim = c(-1.05, 1.05), ylim = ylimIn,
	xlab = xlabExp, ylab = "Density",
	denscol = clPurp[3])
     # Add prior for CORRELATION
     abline(h = 0.5, col = "grey50", lwd = 4)
   mtext(text = expression((bolditalic(e))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)

# Parent U  #####
par(mar = mar2, mgp = c(3, 1, 0))
  postPlot(spiMod_parU$VCV[, 2], plotHist = TRUE,
	xlim = c(-1.05, 1.05), ylim = ylimIn,
	xlab = xlabExp, ylab = "",
	denscol = clPurp[3])
     # Add prior for CORRELATION
     abline(h = 0.5, col = "grey50", lwd = 4)
   mtext(text = expression((bolditalic(f))),
     side = 3, line = lettLine, adj = lettAdj2, cex = 1.3)
 ################### 
 # Body length
 ###################
 # Parent N  #####
par(mar = mar1, mgp = c(3, 1, 0))
  postPlot(bodMod_parN$VCV[, 2], plotHist = TRUE,
	xlim = c(-1.05, 1.05), ylim = ylimIn,
	xlab = xlabExp, ylab = "Density",
	denscol = clPurp[3])
     # Add prior for CORRELATION
     abline(h = 0.5, col = "grey50", lwd = 4)
   mtext(text = expression((bolditalic(g))),
     side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
 # Parent N  #####
par(mar = mar2, mgp = c(3, 1, 0))
  postPlot(bodMod_parU$VCV[, 2], plotHist = TRUE,
	xlim = c(-1.05, 1.05), ylim = ylimIn,
	xlab = xlabExp, ylab = "",
	denscol = clPurp[3])
     # Add prior for CORRELATION
     abline(h = 0.5, col = "grey50", lwd = 4)
   mtext(text = expression((bolditalic(h))),
     side = 3, line = lettLine, adj = lettAdj2, cex = 1.3)


mtext(text = "Larval Environment", 
     outer = TRUE, side = 3, line = -25.5, adj = 0.5, cex = 1.2)

mtext(text = "Parent Non-Upwelling",
     outer = TRUE, side = 3, line = -1.9, adj = 0.073, cex = 1.25)
mtext(text = "Parent Upwelling",
     outer = TRUE, side = 3, line = -1.9, adj = 0.34, cex = 1.25)
mtext(text = "Parent Non-Upwelling",
     outer = TRUE, side = 3, line = -1.9, adj = 0.67, cex = 1.25)
mtext(text = "Parent Upwelling",
     outer = TRUE, side = 3, line = -1.9, adj = 0.915, cex = 1.25)

dev.off()










################################################################################
################################################################################
################################################################################
################################################################################
################################################################################











################################################################################
################################################################################

#			SUPPLEMENTARY MATERIALS 

################################################################################
################################################################################


############################################
#  Figure S1: Egg diameter
############################################


tiff("../FigESM3_S1.tiff", width = 9, height = 4, units = "in",
  res = 500, compression = "jpeg")          

EggPlot <- ggplot(eggs, aes(treat_adult, 1000*Average, color = treat_adult)) +
	labs(y = Egg~Diameter~(mu*m), x = "Parental Treatment") +
	scale_color_manual(values = c("N" = NNcl, "U" = UUcl)) +
	scale_fill_manual(values=c("N" = NNcl, "U" = UUcl)) +
	geom_boxplot() + 
	geom_point() +
	theme_bw() +
	theme(legend.position = "none")

SpiPlot <- ggplot(data = summPspi,
	aes(x = 1000*Average, y = 1000*Length.spi,
	  color = treat, fill = treat)) +
	#geom_smooth(method = "lm", se=FALSE, color="black", formula=y~x)+
	geom_errorbar(aes(ymin = 1000*(Length.spi - se.x),
	  ymax = 1000*(Length.spi + se.x))) +
	geom_errorbarh(aes(xmin = 1000*(Average - se.y),
	  xmax = 1000*(Average + se.y))) +
	labs(y = Prism~Spicule~Length~(mu*m), x = Egg~Diameter~(mu*m)) +
	scale_color_manual(values = c("NN" = NNcl, "NU" = NUcl,
	                              "UN" = UNcl,"UU" = UUcl)) +
	scale_fill_manual(values = c("NN" = NNcl, "NU" = NUcl,
	                             "UN" = UNcl,"UU" = UUcl)) +
	geom_point() + theme_bw() + theme(legend.position = "none")

BodPlot <- ggplot(data = summPbod,
	aes(x = 1000*Average, y = 1000*Length.bod,
	  color = treat, fill = treat)) +
	#geom_smooth(method = "lm", se=FALSE, color="black", formula=y~x)+
	geom_errorbar(aes(ymin = 1000*(Length.bod - se.x),
	  ymax = 1000*(Length.bod + se.x))) +
	geom_errorbarh(aes(xmin = 1000*(Average - se.y),
	  xmax = 1000*(Average + se.y))) +
	labs(y = Prism~Body~Length~(mu*m), x = Egg~Diameter~(mu*m)) +
	scale_color_manual(values = c("NN" = NNcl, "NU" = NUcl,
	                              "UN" = UNcl,"UU" = UUcl)) +
	scale_fill_manual(values = c("NN" = NNcl, "NU" = NUcl,
	                             "UN" = UNcl,"UU" = UUcl)) +
	geom_point() + theme_bw() + theme(legend.position = "none")

grid.arrange(EggPlot, SpiPlot, BodPlot, nrow = 1)





dev.off()








############################################
#  Figure S2: Proportion Abnormality
############################################
# Model treatment mean reaction norms

tiff("../FigESM3_S2.tiff", width = 4, height = 5, units = "in",
  res = 500, compression = "jpeg")          

par(mar = c(5, 5, 2, 0.5), cex.lab = 1.5, cex.axis = 1.25)

tmpPercAbLst <- cbind(NN = percAb_parN$PercAb[which(percAb_parN$treat_dev == "N")],
      UN = percAb_parU$PercAb[which(percAb_parN$treat_dev == "N")],
      NU = percAb_parN$PercAb[which(percAb_parN$treat_dev == "U")],
      UU = percAb_parU$PercAb[which(percAb_parN$treat_dev == "U")])
xaxs <- c(0.9,1.1, 2.9,3.1)
plot(x = rep(xaxs, 2), y = apply(tmpPercAbLst, MARGIN = 2, FUN = range),
  axes = FALSE, type = "n",
  xlim = c(0.5, 3.5), ylim = c(0.12, 0.28),
  xlab = "Larval Environment", ylab = "Proportion Abnormal")
 tmpPercAbMean <- apply(tmpPercAbLst, MARGIN = 2, FUN = mean, na.rm = TRUE)
 tmpPercAbSD <- apply(tmpPercAbLst, MARGIN = 2, FUN = sd, na.rm = TRUE)
 tmpPercAbSE <- tmpPercAbSD / sqrt(apply(tmpPercAbLst, MARGIN = 2, FUN = length))
 # lines/reaction norms
 lines(x = c(0.9, 2.9), y = tmpPercAbMean[c(1,3)], lwd = 3)
 lines(x = c(1.1, 3.1), y = tmpPercAbMean[c(2,4)], lwd = 3)
 # Std. Dev.
 arrows(x0 = xaxs,
   y0 = tmpPercAbMean - 1*tmpPercAbSE, y1 = tmpPercAbMean + 1*tmpPercAbSE,
   length = 0.1, angle = 90, code = 3, col = c(NNcl, UNcl, NUcl, UUcl), lwd = 3)   
 # means
 points(x = xaxs, y = tmpPercAbMean,
   pch = c(21, 24, 21, 24), bg = c(NNcl, UNcl, NUcl, UUcl), cex = 2, lwd = 2)
 # modes
# points(x = xaxs, y = posterior.mode(tmpBodModLst),
#    pch = 8, col = c(NNcl, UNcl, NUcl, UUcl), cex = 2)
   
axis(1, at = c(1,3), labels = c("Non-Upwelling", "Upwelling"))
axis(2, at = seq(0.12, 0.28, 0.04), labels = seq(0.12, 0.28, 0.04))
  
#mtext(text = expression((bolditalic(c))),
#  side = 3, line = -0.4, adj = -0.2, cex = 1.3)

dev.off()





############################################

# Posterior/Prior distribution Figures
# other variance components

############################################


######################
# Spicule length
######################
tmpModN <- spiMod_parN
tmpModU <- spiMod_parU




tiff("../FigESM4_S3.tiff", width = 12, height = 5, units = "in",
  res = 500, compression = "jpeg")          
par(mfcol = c(2, 4), oma = c(1, 7, 4, 0),
  mar = c(4.5,4.6,3,1.5), mgp = c(2.5, 1, 0), cex.lab = 1.55, cex.axis = 1.25)    

lettLine <- 1.5  #<-- vertical line placement of panel letter
lettAdj <- -0.3 #<-- horizonal adjustment of panel letter
#################
# Dam variance
#################
 ################### 
 # Parent N  #####
  NVdam <- postPlot(tmpModN$VCV[, 6], plotHist = TRUE,
	main = "",
	xlab = "", ylab = "Density")
	# estimate prior densities at posterior histogram midpoints
	poh <- NVdam$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")      
     mtext(text = expression((bolditalic(a))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise ylabel increases 1 line each plot??
 ################### 
 # Parent U  #####
  UVdam <- postPlot(tmpModU$VCV[, 6], plotHist = TRUE,
	main = "",
	xlab = expression(V[ Dam ]), ylab = "Density")
	# estimate prior densities at posterior histogram midpoints
	poh <- UVdam$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(b))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)

mtext(text = "Dam variance", 
     outer = TRUE, side = 3, line = 1, adj = 0.1, cex = 1.2)

#################
# Sire variance
#################
 ################### 
 # Parent N  #####
  NVsire <- postPlot(tmpModN$VCV[, 7], plotHist = TRUE,
	xlim = c(0, 3.5),
	main = "",
	xlab = "", ylab = "")
	# estimate prior densities at posterior histogram midpoints
	poh <- NVsire$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")      
     mtext(text = expression((bolditalic(c))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise ylabel increases 1 line each plot??
 ################### 
 # Parent U  #####
  UVsire <- postPlot(tmpModU$VCV[, 7], plotHist = TRUE,
	xlim = c(0, 3.5),
	main = "",
	xlab = expression(V[ Sire ]), ylab = "")
	# estimate prior densities at posterior histogram midpoints
	poh <- UVsire$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(d))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)

mtext(text = "Sire variance", 
     outer = TRUE, side = 3, line = 1, adj = 0.39, cex = 1.2)
#################
# Culture variance
#################
 ################### 
 # Parent N  #####
  NVc <- postPlot(tmpModN$VCV[, 5], plotHist = TRUE,
	xlim = c(0, 0.24), ylim = c(0, 21),
	main = "",
	xlab = "", ylab = "")
	# estimate prior densities at posterior histogram midpoints
	poh <- NVc$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")      
     mtext(text = expression((bolditalic(e))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise ylabel increases 1 line each plot??
 ################### 
 # Parent U  #####
  UVc <- postPlot(tmpModU$VCV[, 5], plotHist = TRUE,
	xlim = c(0, 0.24), ylim = c(0, 21),
	main = "",
	xlab = expression(V[ Culture ]), ylab = "")
	# estimate prior densities at posterior histogram midpoints
	poh <- UVc$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(f))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
mtext(text = "Culture variance", 
     outer = TRUE, side = 3, line = 1, adj = 0.66, cex = 1.2)

#################
# Block variance
#################
 ################### 
 # Parent N  #####
  NVb <- postPlot(tmpModN$VCV[, "blockFac"], plotHist = TRUE, histbreaks = 200,
	xlim = c(0, 4), ylim = c(0, 8),
	main = "",
	xlab = "", ylab = "")
	# estimate prior densities at posterior histogram midpoints
	poh <- NVb$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")      
     mtext(text = expression((bolditalic(g))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise ylabel increases 1 line each plot??
 ################### 
 # Parent U  #####
  UVb <- postPlot(tmpModU$VCV[, "blockFac"], plotHist = TRUE, histbreaks = 200,
	xlim = c(0, 4), ylim = c(0, 4),
	main = "",
	xlab = expression(V[ Block ]), ylab = "")
	# estimate prior densities at posterior histogram midpoints
	poh <- UVb$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(h))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
mtext(text = "Block variance", 
     outer = TRUE, side = 3, line = 1, adj = 0.95, cex = 1.2)


mtext(text = " Parent\nUpwelling",
     outer = TRUE, side = 2, line = 3, adj = 0.2, cex = 1.3)
mtext(text = "Parent  \n Non-Upwelling",
     outer = TRUE, side = 2, line = 3, adj = 0.85, cex = 1.3)


dev.off()










######################
# Residual Variances
######################
tiff("../FigESM4_S4.tiff", width = 5, height = 10, units = "in",
  res = 500, compression = "jpeg")          
par(mfcol = c(4, 1), oma = c(1, 7, 4, 0),
  mar = c(4.5,4.6,3,1), mgp = c(2.5, 1, 0), cex.lab = 1.55, cex.axis = 1.25)    

lettLine <- 1.5  #<-- vertical line placement of panel letter
lettAdj <- -0.3 #<-- horizonal adjustment of panel letter
 ################### 
 # Parent N  #####
  NVr1 <- postPlot(tmpModN$VCV[, 9], plotHist = TRUE,
	xlim = c(0, 1.3), ylim = c(0, 5),
	main = "",
	xlab = "", ylab = "Density",
	histcol = NNcl)
     mtext(text = expression((bolditalic(a))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise ylabel increases 1 line each plot??
  NVr2 <- postPlot(tmpModN$VCV[, 10], plotHist = TRUE,
	xlim = c(0, 2), ylim = c(0, 3),
	main = "",
	xlab = "", ylab = "Density",
	histcol = NUcl)
     mtext(text = expression((bolditalic(b))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
 ################### 
 # Parent U  #####
par(mgp = c(2.5, 1, 0))
  UVr1 <- postPlot(tmpModU$VCV[, 9], plotHist = TRUE,
	xlim = c(0, 1.3), ylim = c(0, 5),
	main = "",
	xlab = "", ylab = "Density",
	histcol = UNcl)
     mtext(text = expression((bolditalic(c))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)

par(mgp = c(2.5, 1, 0))
  UVr2 <- postPlot(tmpModU$VCV[, 10], plotHist = TRUE,
	xlim = c(0, 2), ylim = c(0, 3),
	main = "",
	xlab = expression(V[ Residual ]), ylab = "Density",
	histcol = UUcl)
     mtext(text = expression((bolditalic(d))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)


mtext(text = "Parent Upwelling",
     outer = TRUE, side = 2, line = 5, adj = 0.2, cex = 1.5)
mtext(text = "Parent Non-Upwelling",
     outer = TRUE, side = 2, line = 5, adj = 0.85, cex = 1.5)
  mtext(text = "Larval Upwelling         Larval Non-Upwelling",
     outer = TRUE, side = 2, line = 2.5, adj = 0.93, cex = 1.2)
  mtext(text = "Larval Upwelling         Larval Non-Upwelling",
     outer = TRUE, side = 2, line = 2.5, adj = 0.1, cex = 1.2)

dev.off()

















######################
# Body length
######################
tmpModN <- bodMod_parN
tmpModU <- bodMod_parU




tiff("../FigESM4_S5.tiff", width = 12, height = 5, units = "in",
  res = 500, compression = "jpeg")          
par(mfcol = c(2, 4), oma = c(1, 7, 4, 0),
  mar = c(4.5,4.6,3,1.5), mgp = c(2.5, 1, 0), cex.lab = 1.55, cex.axis = 1.25)    

lettLine <- 1.5  #<-- vertical line placement of panel letter
lettAdj <- -0.3 #<-- horizonal adjustment of panel letter
#################
# Dam variance
#################
 ################### 
 # Parent N  #####
  NVdam <- postPlot(tmpModN$VCV[, 6], plotHist = TRUE,
        xlim = c(0, 2.0), ylim = c(0, 6),
	main = "",
	xlab = "", ylab = "Density")
	# estimate prior densities at posterior histogram midpoints
	poh <- NVdam$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")      
     mtext(text = expression((bolditalic(a))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise ylabel increases 1 line each plot??
 ################### 
 # Parent U  #####
  UVdam <- postPlot(tmpModU$VCV[, 6], plotHist = TRUE,
        xlim = c(0, 2.0), ylim = c(0, 9),
	main = "",
	xlab = expression(V[ Dam ]), ylab = "Density")
	# estimate prior densities at posterior histogram midpoints
	poh <- UVdam$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(b))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)

mtext(text = "Dam variance", 
     outer = TRUE, side = 3, line = 1, adj = 0.1, cex = 1.2)

#################
# Sire variance
#################
 ################### 
 # Parent N  #####
  NVsire <- postPlot(tmpModN$VCV[, 7], plotHist = TRUE,
        xlim = c(0, 0.58), ylim = c(0, 95),
	main = "",
	xlab = "", ylab = "")
	# estimate prior densities at posterior histogram midpoints
	poh <- NVsire$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")      
     mtext(text = expression((bolditalic(c))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise ylabel increases 1 line each plot??
 ################### 
 # Parent U  #####
  UVsire <- postPlot(tmpModU$VCV[, 7], plotHist = TRUE,
        xlim = c(0, 1.0), ylim = c(0, 17),
	main = "",
	xlab = expression(V[ Sire ]), ylab = "")
	# estimate prior densities at posterior histogram midpoints
	poh <- UVsire$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(d))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)

mtext(text = "Sire variance", 
     outer = TRUE, side = 3, line = 1, adj = 0.39, cex = 1.2)
#################
# Culture variance
#################
 ################### 
 # Parent N  #####
  NVc <- postPlot(tmpModN$VCV[, 5], plotHist = TRUE,
        xlim = c(0, 0.16), ylim = c(0, 30),
	main = "",
	xlab = "", ylab = "")
	# estimate prior densities at posterior histogram midpoints
	poh <- NVc$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")      
     mtext(text = expression((bolditalic(e))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise ylabel increases 1 line each plot??
 ################### 
 # Parent U  #####
  UVc <- postPlot(tmpModU$VCV[, 5], plotHist = TRUE,
        xlim = c(0, 0.3), ylim = c(0, 18),
	main = "",
	xlab = expression(V[ Culture ]), ylab = "")
	# estimate prior densities at posterior histogram midpoints
	poh <- UVc$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(f))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
mtext(text = "Culture variance", 
     outer = TRUE, side = 3, line = 1, adj = 0.66, cex = 1.2)
#################
# Block variance
#################
 ################### 
 # Parent N  #####
  NVb <- postPlot(tmpModN$VCV[, "blockFac"], plotHist = TRUE, histbreaks = 200,
	xlim = c(0, 1.5), ylim = c(0, 10),
	main = "",
	xlab = "", ylab = "")
	# estimate prior densities at posterior histogram midpoints
	poh <- NVb$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")      
     mtext(text = expression((bolditalic(g))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise ylabel increases 1 line each plot??
 ################### 
 # Parent U  #####
  UVb <- postPlot(tmpModU$VCV[, "blockFac"], plotHist = TRUE, histbreaks = 200,
	xlim = c(0, 1.5), ylim = c(0, 10),
	main = "",
	xlab = expression(V[ Block ]), ylab = "")
	# estimate prior densities at posterior histogram midpoints
	poh <- UVb$postDensity$histogram
      	prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10),
      	  2*poh$mids[1]/10)
       	prd <- df(prdx / 1000, df1 = 1, df2 = 1, ncp = (0^2)/1000)
        # scale so total of: density * width of histogram bar, summed over all bars = 1
         prda <- sum(2*poh$mids[1]/10 * prd)# prd area
         sprd <- prd / prda
        lines(sprd ~ prdx, lwd = 4, col = "grey50")
     mtext(text = expression((bolditalic(h))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
mtext(text = "Block variance", 
     outer = TRUE, side = 3, line = 1, adj = 0.95, cex = 1.2)


mtext(text = " Parent\nUpwelling",
     outer = TRUE, side = 2, line = 3, adj = 0.2, cex = 1.3)
mtext(text = "Parent  \n Non-Upwelling",
     outer = TRUE, side = 2, line = 3, adj = 0.85, cex = 1.3)



dev.off()










######################
# Residual Variances
######################
tiff("../FigESM4_S6.tiff", width = 5, height = 10, units = "in",
  res = 500, compression = "jpeg")          
par(mfcol = c(4, 1), oma = c(1, 7, 4, 0),
  mar = c(4.5,4.6,3,1), mgp = c(2.5, 1, 0), cex.lab = 1.55, cex.axis = 1.25)    

lettLine <- 1.5  #<-- vertical line placement of panel letter
lettAdj <- -0.3 #<-- horizonal adjustment of panel letter
 ################### 
 # Parent N  #####
  NVr1 <- postPlot(tmpModN$VCV[, 9], plotHist = TRUE,
	xlim = c(0, 1.0), ylim = c(0, 10),
	main = "",
	xlab = "", ylab = "Density",
	histcol = NNcl)
     mtext(text = expression((bolditalic(a))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
par(mgp = c(2.5, 1, 0))  #<--XXX Otherwise ylabel increases 1 line each plot??
  NVr2 <- postPlot(tmpModN$VCV[, 10], plotHist = TRUE,
	xlim = c(0, 0.8), ylim = c(0, 10),
	main = "",
	xlab = "", ylab = "Density",
	histcol = NUcl)
     mtext(text = expression((bolditalic(b))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)
 ################### 
 # Parent U  #####
par(mgp = c(2.5, 1, 0))
  UVr1 <- postPlot(tmpModU$VCV[, 9], plotHist = TRUE,
	xlim = c(0, 1), ylim = c(0, 10),
	main = "",
	xlab = "", ylab = "Density",
	histcol = UNcl)
     mtext(text = expression((bolditalic(c))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)

par(mgp = c(2.5, 1, 0))
  UVr2 <- postPlot(tmpModU$VCV[, 10], plotHist = TRUE,
	xlim = c(0, 0.8), ylim = c(0, 12),
	main = "",
	xlab = expression(V[ Residual ]), ylab = "Density",
	histcol = UUcl)
     mtext(text = expression((bolditalic(d))),
       side = 3, line = lettLine, adj = lettAdj, cex = 1.3)


mtext(text = "Parent Upwelling",
     outer = TRUE, side = 2, line = 5, adj = 0.2, cex = 1.5)
mtext(text = "Parent Non-Upwelling",
     outer = TRUE, side = 2, line = 5, adj = 0.85, cex = 1.5)
  mtext(text = "Larval Upwelling         Larval Non-Upwelling",
     outer = TRUE, side = 2, line = 2.5, adj = 0.93, cex = 1.2)
  mtext(text = "Larval Upwelling         Larval Non-Upwelling",
     outer = TRUE, side = 2, line = 2.5, adj = 0.1, cex = 1.2)

dev.off()


