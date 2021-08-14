rm(list = ls())

library(MCMCglmm)

# Be sure to set the working directory to wherever data and pedigree are stored
## setwd()
################################################################################
# Load data
parN <- read.table("data_parN.txt", header = TRUE)
parU <- read.table("data_parU.txt", header = TRUE)

# Set factors
parN$treat_dev <- as.factor(parN$treat_dev)
parN$animal <- as.factor(parN$animal)
parN$tubeFac <- as.factor(parN$tubeFac)
parN$Dam <- as.factor(parN$Dam)
parN$Sire <- as.factor(parN$Sire)
parN$cross <- as.factor(parN$cross)

parU$treat_dev <- as.factor(parU$treat_dev)
parU$animal <- as.factor(parU$animal)
parU$tubeFac <- as.factor(parU$tubeFac)
parU$Dam <- as.factor(parU$Dam)
parU$Sire <- as.factor(parU$Sire)
parU$cross <- as.factor(parU$cross)

# Load pedigrees
ped_parN <- read.table("pedigree_parN.txt", header = TRUE)
ped_parU <- read.table("pedigree_parU.txt", header = TRUE)

# Create A-inverse relatedness matrices
Ainv_parN <- inverseA(ped_parN)$Ainv
Ainv_parU <- inverseA(ped_parU)$Ainv


################################################################################
################################################################################
# 		MODELS
################################################################################
# Setup prior
quadPE4 <- list(R = list(V = diag(2), nu = 2), 
  G = list(G1 = list(V = diag(2)*0.02, nu = 3, alpha.mu = rep(0, 2), alpha.V = diag(2)*1000),
    G2 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*1000),
    G3 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*1000),
    G4 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*1000)))


##############################
# SPICULE LENGTH
##############################
nsample <- 2000
THIN <- 1000; BURN <- 200000;
 (NITT <- BURN + nsample*THIN)

spiMod_parN <- MCMCglmm(scLength.spi ~ 0 + treat_dev + measurer.spi,
	random = ~ us(treat_dev):animal + tubeFac + Dam + Sire,
	rcov = ~ idh(treat_dev):units, 
        ginverse = list(animal = Ainv_parN),
	prior = quadPE4,
	data = parN,
	nitt = NITT, burnin = BURN, thin = THIN,
	pr = TRUE)
save("spiMod_parN",
  file = paste0("spiMod_parN_nit", NITT/1000, "k.rdata"))

spiMod_parU <- MCMCglmm(scLength.spi ~ 0 + treat_dev + measurer.spi,
	random = ~ us(treat_dev):animal + tubeFac + Dam + Sire,
	rcov = ~ idh(treat_dev):units, 
       ginverse = list(animal = Ainv_parU),
	prior = quadPE4,
	data = parU,
	nitt = NITT, burnin = BURN, thin = THIN,
	pr = TRUE)
save("spiMod_parU",
  file = paste0("spiMod_parU_nit", NITT/1000, "k.rdata"))


##############################
# BODY LENGTH
##############################
nsample <- 2000
THIN <- 1000; BURN <- 130000;
 (NITT <- BURN + nsample*THIN)

bodMod_parN <- MCMCglmm(scLength.bod ~ 0 + treat_dev + measurer.bod,
	random = ~ us(treat_dev):animal + tubeFac + Dam + Sire,
	rcov = ~ idh(treat_dev):units, 
        ginverse = list(animal = Ainv_parN),
	prior = quadPE4,
	data = parN,
	nitt = NITT, burnin = BURN, thin = THIN,
	pr = TRUE)
save("bodMod_parN",
  file = paste0("bodMod_parN_nit", NITT/1000, "k.rdata"))


bodMod_parU <- MCMCglmm(scLength.bod ~ 0 + treat_dev + measurer.bod,
	random = ~ us(treat_dev):animal + tubeFac + Dam + Sire,
	rcov = ~ idh(treat_dev):units, 
        ginverse = list(animal = Ainv_parU),
	prior = quadPE4,
	data = parU,
	nitt = NITT, burnin = BURN, thin = THIN,
	pr = TRUE)
save("bodMod_parU",
  file = paste0("bodMod_parU_nit", NITT/1000, "k.rdata"))




