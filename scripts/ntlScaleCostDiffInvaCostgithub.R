##########################################################################################################
## national-level prediction & extrapolation of costs based on socio-economic traits
## InvaCost Workshop 2022
## November 2022 v1.1
## CJA Bradshaw + Emma Hudgins + Phil Hulme + Brian Leung + Franck Courchamp + Melina Kourantidou +
## Pierre Courtois + Anna Turbelin + Shana McDermott + Katherine Lee + Alok Bang + Thomas Bodey + Phillip Haubrock
##########################################################################################################

library(invacost)
library(lme4)
library(dismo)
library(gbm)
library(boot)
library(VIM)
library(mice)
library(performance)
library(sjPlot)
library(rworldmap)
library(rgeos)
library(SpatialEpi)
library(nlme)
library(rcompanion)

# source files
source("new_lmer_AIC_tables3.R") # change path as required
source("r.squared.R") # change path as required

AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

## functions
delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.IC(AIC.vec); wAIC.vec <- weight.IC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}


# INVACOST
data(invacost)
str(invacost)

#########################
  # cycle through countries
  dat.rel <- subset(invacost, Method_reliability == "High")
  cntryvec <- names(table(dat.rel$Official_country))
  cntry.vec <- cntryvec[-(grep("/", cntryvec))]
  cntry.vec <- cntry.vec[-1]
  lcntry <- length(cntry.vec)
  
  st.yr.vec <- seq(2000, 2018, 3)
  en.yr.vec <- c(st.yr.vec[-1], 2020) 
  iter <- 1000
  itdiv <- iter/10
  
  cntry.dam.md.boot <- cntry.mgm.md.boot <- cntry.dam.mgm.ratio.md.boot <- cntry.dam.mgm.raw.md.boot <-
    cntry.dam.mgm.ratio.r.md <- cntry.dam.gen.md.boot <- cntry.mgm.gen.md.boot <- 
    cntry.pdam.md.boot <- cntry.pmgm.md.boot <- rep(NA,lcntry)

  for (c in 1:lcntry) { # 1:lcntry
  
    cntry.dat.rel <- subset(dat.rel, Official_country == cntry.vec[c])
    dim(cntry.dat.rel)
    cntry.dat.rel$Applicable_year
  
    # observed costs only
    cntry.obs <- subset(cntry.dat.rel, Implementation == "Observed")
  
    # management costs only
    cntry.obs.mgm <- subset(cntry.obs, Type_of_cost_merged == "Management")
  
    # damage costs only
    cntry.obs.dam <- subset(cntry.obs, Type_of_cost_merged == "Damage")
  
    # create resampling loop by 3-year increment starting from 2000
    max.samp <- max(table(cntry.dat.rel$Applicable_year))
    
    cntry.mgm.per.md.mat <- cntry.mgm.per.ngen.mat <- cntry.dam.per.md.mat <- cntry.dam.per.ngen.mat <- 
      cntry.mgmPgen.mat <- cntry.damPgen.mat <- matrix(data=NA, ncol=length(st.yr.vec), nrow=iter)
    
      for (t in 1:length(st.yr.vec)) {
        for (i in 1:iter) {
          
        # management only
        cntry.mgm.per <- cntry.obs.mgm[which(cntry.obs.mgm$Applicable_year >= st.yr.vec[t] & cntry.obs.mgm$Applicable_year >= en.yr.vec[t]), ]
        cntry.mgml <- dim(cntry.mgm.per)[1]
        cntry.mgm.sub <- sample(1:cntry.mgml, max.samp, replace=T)
        cntry.mgm.resamp <- cntry.mgm.per[cntry.mgm.sub,]
        cntry.mgm.per.md.mat[i,t] <- median(cntry.mgm.resamp$Cost_estimate_per_year_2017_USD_exchange_rate, na.rm=T)
        cntry.mgm.per.ngen.mat[i,t] <- length(table(cntry.mgm.resamp$Genus))
        cntry.mgmPgen.mat[i,t] <- cntry.mgm.per.md.mat[i,t] / cntry.mgm.per.ngen.mat[i,t] 
        
        # damage only
        cntry.dam.per <- cntry.obs.dam[which(cntry.obs.dam$Applicable_year >= st.yr.vec[t] & cntry.obs.dam$Applicable_year >= en.yr.vec[t]), ]
        cntry.daml <- dim(cntry.dam.per)[1]
        cntry.dam.sub <- sample(1:cntry.daml, max.samp, replace=T)
        cntry.dam.resamp <- cntry.dam.per[cntry.dam.sub,]
        cntry.dam.per.md.mat[i,t] <- median(cntry.dam.resamp$Cost_estimate_per_year_2017_USD_exchange_rate, na.rm=T)
        cntry.dam.per.ngen.mat[i,t] <- length(table(cntry.dam.resamp$Genus))
        cntry.damPgen.mat[i,t] <- cntry.dam.per.md.mat[i,t] / cntry.dam.per.ngen.mat[i,t] 
  
        if (i %% itdiv==0) print(i) 
        
        } # end i
        
        print("###################")
        print(paste(st.yr.vec[t], "to", en.yr.vec[t], sep=" "))
        print("###################")
        
      } # end t
      
    cntry.mgm.per.ngen.md <- apply(cntry.mgm.per.ngen.mat, MARGIN=2, median, na.rm=T)
    cntry.mgm.per.ngen.lo <- apply(cntry.mgm.per.ngen.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
    cntry.mgm.per.ngen.up <- apply(cntry.mgm.per.ngen.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
    
    cntry.mgm.per.md <- apply(cntry.mgm.per.md.mat, MARGIN=2, median, na.rm=T)
    cntry.mgm.per.lo <- apply(cntry.mgm.per.md.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
    cntry.mgm.per.up <- apply(cntry.mgm.per.md.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
     
    cntry.mgmPgen.per.md <- apply(cntry.mgmPgen.mat, MARGIN=2, median, na.rm=T)
    cntry.mgmPgen.per.lo <- apply(cntry.mgmPgen.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
    cntry.mgmPgen.per.up <- apply(cntry.mgmPgen.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
     
    cntry.dam.per.ngen.md <- apply(cntry.dam.per.ngen.mat, MARGIN=2, median, na.rm=T)
    cntry.dam.per.ngen.lo <- apply(cntry.dam.per.ngen.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
    cntry.dam.per.ngen.up <- apply(cntry.dam.per.ngen.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
    
    cntry.dam.per.md <- apply(cntry.dam.per.md.mat, MARGIN=2, median, na.rm=T)
    cntry.dam.per.lo <- apply(cntry.dam.per.md.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
    cntry.dam.per.up <- apply(cntry.dam.per.md.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
     
    cntry.damPgen.per.md <- apply(cntry.damPgen.mat, MARGIN=2, median, na.rm=T)
    cntry.damPgen.per.lo <- apply(cntry.damPgen.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
    cntry.damPgen.per.up <- apply(cntry.damPgen.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
    
    cntry.dam.mgm.ratio.md <- cntry.damPgen.per.md/cntry.mgmPgen.per.md
    cntry.dam.mgm.ratio.up <- cntry.damPgen.per.lo/cntry.mgmPgen.per.lo
    cntry.dam.mgm.ratio.lo <- cntry.damPgen.per.up/cntry.mgmPgen.per.up
    
    cntry.dam.mgm.raw.md <- cntry.dam.per.md / cntry.mgm.per.md
    cntry.dam.mgm.raw.lo <- cntry.dam.per.lo / cntry.mgm.per.lo
    cntry.dam.mgm.raw.up <- cntry.dam.per.up / cntry.mgm.per.up

    cntry.pdam.md <- cntry.dam.per.md / (cntry.mgm.per.md + cntry.dam.per.md)
    cntry.pdam.lo <- cntry.dam.per.lo / (cntry.mgm.per.lo + cntry.dam.per.lo)
    cntry.pdam.up <- cntry.dam.per.up / (cntry.mgm.per.up + cntry.dam.per.up)

    cntry.pmgm.md <- cntry.mgm.per.md / (cntry.mgm.per.md + cntry.dam.per.md)
    cntry.pmgm.lo <- cntry.mgm.per.lo / (cntry.mgm.per.lo + cntry.dam.per.lo)
    cntry.pmgm.up <- cntry.mgm.per.up / (cntry.mgm.per.up + cntry.dam.per.up)
    
    cntry.dam.mgm.ratio.r <- log10(cntry.dam.mgm.ratio.md[2:length(cntry.dam.mgm.ratio.md)] / cntry.dam.mgm.ratio.md[1:(length(cntry.dam.mgm.ratio.md)-1)])
    if (length(cntry.dam.mgm.ratio.r) > 0) {
      cntry.dam.mgm.ratio.r.md[c] <- median(cntry.dam.mgm.ratio.r, na.rm=T)
    }
    
    cntry.dam.gen <- as.vector(na.omit(c(cntry.dam.per.ngen.md,cntry.dam.per.ngen.lo,cntry.dam.per.ngen.up)))
    if (length(cntry.dam.gen[cntry.dam.gen > 0]) > 0) {
      cntry.dam.gen.md.boot[c] <- mean(sample(cntry.dam.gen[cntry.dam.gen > 0], 10000, replace=T))
    }
    cntry.dam <- as.vector(na.omit(c(cntry.dam.per.md,cntry.dam.per.lo,cntry.dam.per.up)))
    if (length(cntry.dam) > 0) {
      cntry.dam.md.boot[c] <- mean(sample(cntry.dam, 10000, replace=T))
    }
    
    cntry.mgm.gen <- as.vector(na.omit(c(cntry.mgm.per.ngen.md,cntry.mgm.per.ngen.lo,cntry.mgm.per.ngen.up)))
    if (length(cntry.mgm.gen[cntry.mgm.gen > 0]) > 0) {
      cntry.mgm.gen.md.boot[c] <- mean(sample(cntry.mgm.gen[cntry.mgm.gen > 0], 10000, replace=T))
    }
    cntry.mgm <- as.vector(na.omit(c(cntry.mgm.per.md,cntry.mgm.per.lo,cntry.mgm.per.up)))
    if (length(cntry.mgm) > 0) {
      cntry.mgm.md.boot[c] <- mean(sample(cntry.mgm, 10000, replace=T))
    }
    
    cntry.dam.mgm.ratio <- as.vector(na.omit(c(cntry.dam.mgm.ratio.md, cntry.dam.mgm.ratio.up, cntry.dam.mgm.ratio.lo)))
    if (length(cntry.dam.mgm.ratio) > 0) {
      cntry.dam.mgm.ratio.md.boot[c] <- mean(sample(cntry.dam.mgm.ratio, 10000, replace=T))
    }
    
    cntry.dam.mgm.raw <- as.vector(na.omit(c(cntry.dam.mgm.raw.md, cntry.dam.mgm.raw.up, cntry.dam.mgm.raw.lo)))
    if (length(cntry.dam.mgm.raw) > 0) {
      cntry.dam.mgm.raw.md.boot[c] <- mean(sample(cntry.dam.mgm.raw, 10000, replace=T))
    }
    
    cntry.pdam <- as.vector(na.omit(c(cntry.pdam.md, cntry.pdam.up, cntry.pdam.lo)))
    if (length(cntry.pdam) > 0) {
      cntry.pdam.md.boot[c] <- mean(sample(cntry.pdam, 10000, replace=T))
    }
    cntry.pmgm <- as.vector(na.omit(c(cntry.pmgm.md, cntry.pmgm.up, cntry.pmgm.lo)))
    if (length(cntry.pmgm) > 0) {
      cntry.pmgm.md.boot[c] <- mean(sample(cntry.pmgm, 10000, replace=T))
    }
    
    print("###################")
    print(paste(cntry.vec[c], " (", round(c/lcntry*100, 0), "% complete)", sep=""))
    print("###################")
  } # end c
  
ratio.out <- data.frame(cntry.vec,cntry.dam.mgm.ratio.md.boot,cntry.dam.mgm.raw.md.boot,cntry.dam.mgm.ratio.r.md,
                        cntry.pdam.md.boot, cntry.pmgm.md.boot)
colnames(ratio.out) <- c("country","damPgen.mgmPgen","dam.mgm","ratio.r","pdam","pmgm")
write.csv(ratio.out, "ratioOut.csv")
hist(log10(ratio.out$dam.mgm),main="",xlab="log10 damage:management")
hist((ratio.out$ratio.r),main="",xlab="damage:management r")
hist(logit(ratio.out$pdam),main="",xlab="damage:management r")
hist(logit(ratio.out$pmgm),main="",xlab="damage:management r")

# dam by mgm relationship (across countries)
dam.mgm.out <- data.frame(cntry.vec, cntry.dam.md.boot, cntry.mgm.md.boot)
colnames(dam.mgm.out) <- c("country","dam","mgm")
plot(log10(dam.mgm.out$mgm), log10(dam.mgm.out$dam), pch=19)
fit.dam.mgm <- lm(log10(dam.mgm.out$dam) ~ log10(dam.mgm.out$mgm))
summary(fit.dam.mgm)
abline(fit.dam.mgm, lty=2, col="red")
linreg.ER(log10(dam.mgm.out$mgm), log10(dam.mgm.out$dam))

check_model(fit.dam.mgm)
plot_model(fit.dam.mgm, show.values=T, vline.color = "purple")

# dam by number of genera
dam.ngen.out <- data.frame(cntry.vec, cntry.dam.gen.md.boot, cntry.dam.md.boot)
colnames(dam.ngen.out) <- c("country", "ngen", "dam")
head(dam.ngen.out)
plot(log10(dam.ngen.out$ngen), log10(dam.ngen.out$dam), pch=19, xlab="log10 number of genera", ylab="log10 annual damage cost")
fit <- lm(log10(dam.ngen.out$dam) ~ log10(dam.ngen.out$ngen))
summary(fit)
abline(fit, lty=2, col="red")

# mgm by number of genera
mgm.ngen.out <- data.frame(cntry.vec, cntry.mgm.gen.md.boot, cntry.mgm.md.boot)
colnames(mgm.ngen.out) <- c("country", "ngen", "mgm")
head(mgm.ngen.out)
plot(log10(mgm.ngen.out$ngen), log10(mgm.ngen.out$mgm), pch=19, xlab="log10 number of genera", ylab="log10 annual management cost")
fit <- lm(log10(mgm.ngen.out$mgm) ~ log10(mgm.ngen.out$ngen))
summary(fit)
abline(fit, lty=2, col="red")

# import cntry.code
cont.cntry <- read.csv("continent.countryINVACOST.csv", header=T)

# merge with cntry.code
ratio.cntry <- merge(ratio.out, cont.cntry, by="country")
head(ratio.cntry)

dam.mgm.cntry.out <- merge(dam.mgm.out, cont.cntry, by="country")
table(dam.mgm.cntry.out$cont)
dam.mgm.cntry.out$reg2 <- ifelse(dam.mgm.cntry.out$cont == "NAM" | dam.mgm.cntry.out$cont == "CAR", "NAMCAR", dam.mgm.cntry.out$cont)
dam.mgm.cntry.out$reg2 <- ifelse(dam.mgm.cntry.out$cont == "EUR" | dam.mgm.cntry.out$cont == "ME", "EURME", dam.mgm.cntry.out$reg2)
dam.mgm.cntry.out$reg2 <- ifelse(dam.mgm.cntry.out$cont == "ASIA" | dam.mgm.cntry.out$cont == "OC", "ASIAOC", dam.mgm.cntry.out$reg2)
table(dam.mgm.cntry.out$reg2)

write.csv(dam.mgm.cntry.out, "damMgmOut.csv")

dam.ngen.cntry.out <- merge(dam.ngen.out, cont.cntry, by="country")
table(dam.ngen.cntry.out$cont)
dam.ngen.cntry.out$reg2 <- ifelse(dam.ngen.cntry.out$cont == "NAM" | dam.ngen.cntry.out$cont == "CAR", "NAMCAR", dam.ngen.cntry.out$cont)
dam.ngen.cntry.out$reg2 <- ifelse(dam.ngen.cntry.out$cont == "EUR" | dam.ngen.cntry.out$cont == "ME", "EURME", dam.ngen.cntry.out$reg2)
dam.ngen.cntry.out$reg2 <- ifelse(dam.ngen.cntry.out$cont == "ASIA" | dam.ngen.cntry.out$cont == "OC", "ASIAOC", dam.ngen.cntry.out$reg2)
table(dam.ngen.cntry.out$reg2)

write.csv(dam.ngen.cntry.out, "damGenOut.csv")

mgm.ngen.cntry.out <- merge(mgm.ngen.out, cont.cntry, by="country")
table(mgm.ngen.cntry.out$cont)
mgm.ngen.cntry.out$reg2 <- ifelse(mgm.ngen.cntry.out$cont == "NAM" | mgm.ngen.cntry.out$cont == "CAR", "NAMCAR", mgm.ngen.cntry.out$cont)
mgm.ngen.cntry.out$reg2 <- ifelse(mgm.ngen.cntry.out$cont == "EUR" | mgm.ngen.cntry.out$cont == "ME", "EURME", mgm.ngen.cntry.out$reg2)
mgm.ngen.cntry.out$reg2 <- ifelse(mgm.ngen.cntry.out$cont == "ASIA" | mgm.ngen.cntry.out$cont == "OC", "ASIAOC", mgm.ngen.cntry.out$reg2)
table(mgm.ngen.cntry.out$reg2)

write.csv(mgm.ngen.cntry.out, "mgmGenOut.csv")


# GDP per capita
gdp <- read.csv("GDPpc.csv")
head(gdp)
ratio.gdp <- merge(ratio.cntry, gdp, by="cntry.code", all = T)
ratio.gdp$GDPpcRecent <- as.numeric(ratio.gdp$GDPpcRecent)
head(ratio.gdp)

plot(log10(ratio.gdp$GDPpcRecent), log10(ratio.gdp$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 GDP per capita")
fit <- lm(log10(ratio.gdp$dam.mgm) ~ log10(ratio.gdp$GDPpcRecent))
summary(fit)
abline(fit,lty=2, col="red")

plot(log10(ratio.gdp$GDPpcRecent), (ratio.gdp$ratio.r), pch=19, ylab="damage:management r", xlab="log10 GDP per capita")
fit <- lm((ratio.gdp$ratio.r) ~ log10(ratio.gdp$GDPpcRecent))
summary(fit)
abline(fit,lty=2, col="red")

# by continent
table(ratio.gdp$cont)
EUR.dat <- subset(ratio.gdp, cont=="EUR")
plot(log10(EUR.dat$GDPpcRecent), log10(EUR.dat$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 GDP per capita")
fit <- lm(log10(EUR.dat$dam.mgm) ~ log10(EUR.dat$GDPpcRecent))
abline(fit,lty=2, col="red")

ASIA.dat <- subset(ratio.gdp, cont=="ASIA")
plot(log10(ASIA.dat$GDPpcRecent), log10(ASIA.dat$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 GDP per capita")
fit <- lm(log10(ASIA.dat$dam.mgm) ~ log10(ASIA.dat$GDPpcRecent))
abline(fit,lty=2, col="red")

AFR.dat <- subset(ratio.gdp, cont=="AFR")
plot(log10(AFR.dat$GDPpcRecent), log10(AFR.dat$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 GDP per capita")
fit <- lm(log10(AFR.dat$dam.mgm) ~ log10(AFR.dat$GDPpcRecent))
abline(fit,lty=2, col="red")

NAM.dat <- subset(ratio.gdp, cont=="NAM")
plot(log10(NAM.dat$GDPpcRecent), log10(NAM.dat$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 GDP per capita")
fit <- lm(log10(NAM.dat$dam.mgm) ~ log10(NAM.dat$GDPpcRecent))
abline(fit,lty=2, col="red")

SA.dat <- subset(ratio.gdp, cont=="SA")
plot(log10(SA.dat$GDPpcRecent), log10(SA.dat$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 GDP per capita")
fit <- lm(log10(SA.dat$dam.mgm) ~ log10(SA.dat$GDPpcRecent))
abline(fit,lty=2, col="red")

OC.dat <- subset(ratio.gdp, cont=="OC")
plot(log10(OC.dat$GDPpcRecent), log10(OC.dat$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 GDP per capita")
fit <- lm(log10(OC.dat$dam.mgm) ~ log10(OC.dat$GDPpcRecent))
abline(fit,lty=2, col="red")



## FAO proportion of GDP devoted to agriculture
## Agriculture, forestry, and fishing, value added (% of GDP)
# Agriculture, forestry, and fishing corresponds to ISIC divisions 1-3 and includes forestry, hunting, and fishing, as well as
# cultivation of crops and livestock production. Value added is the net output of a sector after adding up all outputs and subtracting
# intermediate inputs. It is calculated without making deductions for depreciation of fabricated assets or depletion and degradation
# of natural resources. The origin of value added is determined by the International Standard Industrial Classification (ISIC),
# revision 4. Note: For VAB countries, gross value added at factor cost is used as the denominator.
FAOag <- read.csv("faoag.csv")
VAag.all <- subset(FAOag, Item == 'Value Added (Agriculture, Forestry and Fishing)' & ElementCode == 6139)
head(VAag.all)

# cycle through countries for most recent value
VAag.curr <- VAag.all[1,]
cntry.vec <- names(table(VAag.all$Area))
for (c in 1:length(cntry.vec)) {
  cntry.it <- subset(VAag.all, Area == cntry.vec[c])
  VAag.curr <- rbind(VAag.curr, cntry.it[which(cntry.it$Year == max(cntry.it$Year, na.rm=T)),])
}
VAag.curr <- VAag.curr[-1,]
head(VAag.curr)
VAag1 <- VAag.curr[,c(4,10,12)]
colnames(VAag1) <- c("country", "year", "VAag")
head(VAag1)

fao.cntry.code <- read.csv("fao.cntry.code.csv")
head(fao.cntry.code)
VAag <- merge(VAag1, fao.cntry.code, by="country", all=T)
head(VAag)

# Transparency International
# Brenton-Rule et al. 2016 - corruption (no relationship with # invasive species)
# https://royalsocietypublishing.org/doi/full/10.1098/rspb.2016.0901
# Corruption Perception index https://www.transparency.org/en/cpi/2021
cpi <- read.csv("CPI.csv")
head(cpi)

# gov expenditure on education (% GDP)
# https://data.worldbank.org/indicator/SE.XPD.TOTL.GD.ZS
expedu <- read.csv("govexpedu.csv")
head(expedu)

# Global Health Security Index
# https://www.ghsindex.org/report-model/
ghsi <- read.csv("GHSI2022.csv")
head(ghsi)

# % agricultural land
# https://data.worldbank.org/indicator/AG.LND.AGRI.ZS
agrlnd <- read.csv("pcAgrLand.csv")
head(agrlnd)

# imports of goods and services per capita (mean of last 5 years)
# https://data.worldbank.org/indicator/NE.IMP.GNFS.CD
igs <- read.csv("importGS.csv")
head(igs)
pop21 <- read.csv("pop2021.csv") # https://data.worldbank.org/indicator/SP.POP.TOTL
head(pop21)
igspc <- merge(igs, pop21, by="cntry.code", all=T)
head(igspc)
igspc$igspc <- igspc$igs/igspc$popN
head(igspc)

# scientific & technical journal articles (per capita)
# https://data.worldbank.org/indicator/IP.JRN.ARTC.SC
stja <- read.csv("stjarticles.csv")
head(stja)
stjapc <- merge(stja, pop21, by="cntry.code", all=T)
head(stjapc)
stjapc$stjapc <- stjapc$stja/stjapc$popN
head(stjapc)


ratio.gdpcpi <- merge(ratio.gdp, cpi, by="cntry.code", all=T)
head(ratio.gdpTcpi)

ratio.full1 <- merge(ratio.gdpcpi, expedu, by="cntry.code", all=T)
head(ratio.full1)

ratio.full2 <- merge(ratio.full1, VAag, by="cntry.code", all=T)
head(ratio.full2)

ratio.full3 <- merge(ratio.full2, ghsi, by="cntry.code", all=T)
head(ratio.full3)

ratio.full4 <- merge(ratio.full3, agrlnd, by="cntry.code", all=T)
head(ratio.full4)

ratio.full5 <- merge(ratio.full4, igspc, by="cntry.code", all=T)
head(ratio.full5)

ratio.full6 <- merge(ratio.full5, stjapc, by="cntry.code", all=T)
head(ratio.full6)

ratio.full <- merge(ratio.full6, pa, by="cntry.code", all=T)
head(ratio.full)

dim(ratio.full)
head(ratio.full)
length(table(ratio.full$cntry.code))
ratio.full <- ratio.full[is.na(ratio.full$cntry.code)==F,]
dim(ratio.full)
head(ratio.full)
length(table(ratio.full$cntry.code))


plot((ratio.full$CPI), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="corruption index")
fit <- lm(log10(ratio.full$dam.mgm) ~ (ratio.full$CPI))
abline(fit,lty=2, col="red")

plot((ratio.full$govexpedu), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="gov education expend %GDP")
fit <- lm(log10(ratio.full$dam.mgm) ~ (ratio.full$govexpedu))
abline(fit,lty=2, col="red")

plot((ratio.full$VAag), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="value added agriculture/fisheries/forestry %GDP")
fit <- lm(log10(ratio.full$dam.mgm) ~ (ratio.full$protareapc))
abline(fit,lty=2, col="red")

plot((ratio.full$GHSI), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="global health security index")
fit <- lm(log10(ratio.full$dam.mgm) ~ (ratio.full$GHSI))
abline(fit,lty=2, col="red")

plot((ratio.full$pcAgrLand), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="% agricultural land")
fit <- lm(log10(ratio.full$dam.mgm) ~ (ratio.full$pcAgrLand))
abline(fit,lty=2, col="red")

plot(log10(ratio.full$igspc), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 pc import goods & services")
fit <- lm(log10(ratio.full$dam.mgm) ~ log10(ratio.full$igspc))
abline(fit,lty=2, col="red")

plot(log10(ratio.full$stjapc), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 pc sci/tech articles")
fit <- lm(log10(ratio.full$dam.mgm) ~ log10(ratio.full$stjapc))
abline(fit,lty=2, col="red")


plot(log10(ratio.full$GDPpcRecent), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="log10 GDP per capita")
fit <- lm(logit(ratio.full$pdam) ~ log10(ratio.full$GDPpcRecent))
summary(fit)
abline(fit,lty=2, col="red")
plot(log10(ratio.full$GDPpcRecent), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="log10 GDP per capita")
fit <- lm(logit(ratio.full$pmgm) ~ log10(ratio.full$GDPpcRecent))
abline(fit,lty=2, col="red")

plot((ratio.full$CPI), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="corruption index")
fit <- lm(logit(ratio.full$pdam) ~ (ratio.full$CPI))
abline(fit,lty=2, col="red")
plot((ratio.full$CPI), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="corruption index")
fit <- lm(logit(ratio.full$pmgm) ~ (ratio.full$CPI))
abline(fit,lty=2, col="red")

plot((ratio.full$govexpedu), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="gov education expend %GDP")
fit <- lm(logit(ratio.full$pdam) ~ (ratio.full$govexpedu))
abline(fit,lty=2, col="red")
plot((ratio.full$govexpedu), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="gov education expend %GDP")
fit <- lm(logit(ratio.full$pmgm) ~ (ratio.full$govexpedu))
abline(fit,lty=2, col="red")

plot((ratio.full$VAag/100), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="value added primary production %GDP")
fit <- lm(logit(ratio.full$pdam) ~ (ratio.full$VAag))
abline(fit,lty=2, col="red")
plot((ratio.full$VAag), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="value added primary production %GDP")
fit <- lm(logit(ratio.full$pmgm) ~ (ratio.full$VAag))
abline(fit,lty=2, col="red")

plot((ratio.full$GHSI), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="global health security index")
fit <- lm(logit(ratio.full$pdam) ~ (ratio.full$GHSI))
abline(fit,lty=2, col="red")
plot((ratio.full$GHSI), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="global health security index")
fit <- lm(logit(ratio.full$pmgm) ~ (ratio.full$GHSI))
abline(fit,lty=2, col="red")

plot(logit(ratio.full$pcAgrLand/100), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="logit % agricultural land")
fit <- lm(logit(ratio.full$pdam) ~ logit(ratio.full$pcAgrLand/100))
abline(fit,lty=2, col="red")
plot(logit(ratio.full$pcAgrLand/100), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="logit % agricultural land")
fit <- lm(logit(ratio.full$pmgm) ~ logit(ratio.full$pcAgrLand/100))
abline(fit,lty=2, col="red")

plot(log10(ratio.full$igspc), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="log10 pc import goods & services")
fit <- lm(logit(ratio.full$pdam) ~ log10(ratio.full$igspc))
summary(fit)
abline(fit,lty=2, col="red")
plot(log10(ratio.full$igspc), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="log10 pc import goods & services")
fit <- lm(logit(ratio.full$pmgm) ~ log10(ratio.full$igspc))
abline(fit,lty=2, col="red")

plot(log10(ratio.full$stjapc), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="log10 pc sci/tech articles")
fit <- lm(logit(ratio.full$pdam) ~ log10(ratio.full$stjapc))
summary(fit)
abline(fit,lty=2, col="red")
plot(log10(ratio.full$stjapc), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="log10 pc sci/tech articles")
fit <- lm(logit(ratio.full$pmgm) ~ log10(ratio.full$stjapc))
abline(fit,lty=2, col="red")


final.dat.raw <- data.frame(ratio.full$GDPpcRecent, ratio.full$CPI, ratio.full$GHSI, ratio.full$igspc, ratio.full$pcAgrLand,
                            ratio.full$VAag, ratio.full$govexpedu, ratio.full$stjapc)
colnames(final.dat.raw) <- c("gdp", "cpi", "ghsi", "igs", "agrL", "VAag", "govexpedu", "stja")
head(final.dat.raw)
tail(final.dat.raw)
dim(final.dat.raw)
dim(na.omit(final.dat.raw))

aggr_plot <- aggr(final.dat.raw, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(final.dat.raw), cex.axis=.7, gap=3, ylab=c("histogram of missing data","pattern"))
marginplot(final.dat.raw[c(1,2)])
aggr_plot$missings
missings <- data.frame(aggr_plot$missings$Variable, round(100*aggr_plot$missings$Count/dim(final.dat.raw)[1], 1))
colnames(missings) <- c("var", "pcMiss")
missings

final.dat.rawRatio <- data.frame(ratio.full$dam.mgm, final.dat.raw)
colnames(final.dat.rawRatio)[1] <- "dam.mgm"
head(final.dat.rawRatio)
final.dat.rawRatio <- final.dat.rawRatio[is.na(final.dat.rawRatio$dam.mgm)==F,]
head(final.dat.rawRatio)
dim(final.dat.rawRatio)
aggr_plot2 <- aggr(final.dat.rawRatio, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(final.dat.raw), cex.axis=.7, gap=3, ylab=c("histogram of missing data","pattern"))

# multiple imputation
md.pattern(final.dat.raw)
final.dat.imp <- mice(final.dat.raw, m=8, maxit=500, method="pmm", seed=101)
summary(final.dat.imp)

final.dat.imp$imp$gdp
final.dat.imp$imp$cpi
final.dat.imp$imp$ghsi
final.dat.imp$imp$igs
final.dat.imp$imp$agrL
final.dat.imp$imp$VAag
final.dat.imp$imp$govexpedu
final.dat.imp$imp$stja

final.dat.compl1 <- complete(final.dat.imp, 1)
final.dat.compl2 <- complete(final.dat.imp, 2)
final.dat.compl3 <- complete(final.dat.imp, 3)
final.dat.compl4 <- complete(final.dat.imp, 4)
final.dat.compl5 <- complete(final.dat.imp, 5)
final.dat.compl6 <- complete(final.dat.imp, 6)
final.dat.compl7 <- complete(final.dat.imp, 7)
final.dat.compl8 <- complete(final.dat.imp, 8)

final.dat.compl.mn <- final.dat.raw
sub1 <- which(is.na(final.dat.raw[,1]) == T)
dat1 <- data.frame(final.dat.compl1[sub1, 1], final.dat.compl2[sub1, 1], final.dat.compl3[sub1, 1], final.dat.compl4[sub1, 1], final.dat.compl5[sub1, 1], final.dat.compl6[sub1, 1], final.dat.compl7[sub1, 1], final.dat.compl8[sub1, 1])
final.dat.compl.mn[sub1, 1] <- apply(dat1, MARGIN=1, mean)

sub2 <- which(is.na(final.dat.raw[,2]) == T)
dat2 <- data.frame(final.dat.compl1[sub2, 2], final.dat.compl2[sub2, 2], final.dat.compl3[sub2, 2], final.dat.compl4[sub2, 2], final.dat.compl5[sub2, 2], final.dat.compl6[sub2, 2], final.dat.compl7[sub2, 2], final.dat.compl8[sub2, 2])
final.dat.compl.mn[sub2, 2] <- apply(dat2, MARGIN=1, mean)

sub3 <- which(is.na(final.dat.raw[,3]) == T)
dat3 <- data.frame(final.dat.compl1[sub3, 3], final.dat.compl2[sub3, 3], final.dat.compl3[sub3, 3], final.dat.compl4[sub3, 3], final.dat.compl5[sub3, 3], final.dat.compl6[sub3, 3], final.dat.compl7[sub3, 3], final.dat.compl8[sub3, 3])
final.dat.compl.mn[sub3, 3] <- apply(dat3, MARGIN=1, mean)

sub4 <- which(is.na(final.dat.raw[,4]) == T)
dat4 <- data.frame(final.dat.compl1[sub4, 4], final.dat.compl2[sub4, 4], final.dat.compl3[sub4, 4], final.dat.compl4[sub4, 4], final.dat.compl5[sub4, 4], final.dat.compl6[sub4, 4], final.dat.compl7[sub4, 4], final.dat.compl8[sub4, 4])
final.dat.compl.mn[sub4, 4] <- apply(dat4, MARGIN=1, mean)

sub5 <- which(is.na(final.dat.raw[,5]) == T)
dat5 <- data.frame(final.dat.compl1[sub5, 5], final.dat.compl2[sub5, 5], final.dat.compl3[sub5, 5], final.dat.compl4[sub5, 5], final.dat.compl5[sub5, 5], final.dat.compl6[sub5, 5], final.dat.compl7[sub5, 5], final.dat.compl8[sub5, 5])
final.dat.compl.mn[sub5, 5] <- apply(dat5, MARGIN=1, mean)

sub6 <- which(is.na(final.dat.raw[,6]) == T)
dat6 <- data.frame(final.dat.compl1[sub6, 6], final.dat.compl2[sub6, 6], final.dat.compl3[sub6, 6], final.dat.compl4[sub6, 6], final.dat.compl5[sub6, 6], final.dat.compl6[sub6, 6], final.dat.compl7[sub6, 6], final.dat.compl8[sub6, 6])
final.dat.compl.mn[sub6, 6] <- apply(dat6, MARGIN=1, mean)

sub7 <- which(is.na(final.dat.raw[,7]) == T)
dat7 <- data.frame(final.dat.compl1[sub7, 7], final.dat.compl2[sub7, 7], final.dat.compl3[sub7, 7], final.dat.compl4[sub7, 7], final.dat.compl5[sub7, 7], final.dat.compl6[sub7, 7], final.dat.compl7[sub7, 7], final.dat.compl8[sub7, 7])
final.dat.compl.mn[sub7, 7] <- apply(dat7, MARGIN=1, mean)

sub8 <- which(is.na(final.dat.raw[,8]) == T)
dat8 <- data.frame(final.dat.compl1[sub8, 8], final.dat.compl2[sub8, 8], final.dat.compl3[sub8, 8], final.dat.compl4[sub8, 8], final.dat.compl5[sub8, 8], final.dat.compl6[sub8, 8], final.dat.compl7[sub8, 8], final.dat.compl8[sub8, 8])
final.dat.compl.mn[sub8, 8] <- apply(dat8, MARGIN=1, mean)

gdpsc <- log10(scale((final.dat.compl.mn$gdp), scale=T, center=F))
cpisc <- scale(final.dat.compl.mn$cpi, scale=T, center=F)
ghsisc <- scale((final.dat.compl.mn$ghsi), scale=T, center=F)
igssc <- log10(scale(final.dat.compl.mn$igs, scale=T, center=F))
agrLsc <- scale(logit(final.dat.compl.mn$agrL/100), scale=T, center=F)
VAagsc <- scale(logit(final.dat.compl.mn$VAag/100), scale=T, center=F)
govexpedusc <- scale(logit(final.dat.compl.mn$govexpedu/100), scale=T, center=F)
stjasc <- log10(scale((final.dat.compl.mn$stja), scale=T, center=F))
ratiosc <- log10(scale(ratio.full$dam.mgm, center=F, scale=T))
rsc <- scale(ratio.full$ratio.r, center=F, scale=T)
pdamsc <- scale(logit(ratio.full$pdam), center=F, scale=T)

head(ratio.full)
dim(ratio.full)

final.dat.imputed <- na.omit(data.frame(ratio.full$cntry.code, ratio.full$cont, ratio.full$region, ratiosc,
                                gdpsc, cpisc, ghsisc, igssc, agrLsc, VAagsc, govexpedusc, stjasc))
colnames(final.dat.imputed) <- c("cntry.code", "cont", "reg", "ratio", "gdp", "cpi", "ghsi", "igs", "agrL", "VAag", "govexpedu", "stja")
head(final.dat.imputed)
dim(final.dat.imputed)

final.datPdam.imputed <- na.omit(data.frame(ratio.full$cntry.code, ratio.full$cont, ratio.full$region, pdamsc,
                                        gdpsc, cpisc, ghsisc, igssc, agrLsc, VAagsc, govexpedusc, stjasc))
colnames(final.datPdam.imputed) <- c("cntry.code", "cont", "reg", "pdam", "gdp", "cpi", "ghsi", "igs", "agrL", "VAag", "govexpedu", "stja")
head(final.datPdam.imputed)
dim(final.datPdam.imputed)

final.datr.imputed <- na.omit(data.frame(ratio.full$cntry.code, ratio.full$cont, ratio.full$region, rsc,
                                        gdpsc, cpisc, ghsisc, igssc, agrLsc, VAagsc, govexpedusc, stjasc))
colnames(final.datr.imputed) <- c("cntry.code", "cont", "reg", "r", "gdp", "cpi", "ghsi", "igs", "agrL", "VAag", "govexpedu", "stja")
head(final.datr.imputed)
dim(final.datr.imputed)

## correlation matrix
cor.dat.ratio <- data.frame(final.dat.imputed$gdp, final.dat.imputed$cpi, final.dat.imputed$ghsi, final.dat.imputed$igs, final.dat.imputed$agrL,
                            final.dat.imputed$VAag, final.dat.imputed$govexpedu, final.dat.imputed$stja)
colnames(cor.dat.ratio) <- c("GDP","CPI","GHSI","IGS","AGRL","VAPP","EDU","STJA")
cormat.ratio <- cor(na.omit(cor.dat.ratio), method="kendall")
cormat.ratio[lower.tri(cormat.ratio)] <- NA
cormat.ratio
max(abs(cormat.ratio[cormat.ratio < 1]), na.rm=T)
median(abs(cormat.ratio[cormat.ratio < 1]), na.rm=T)
hist(abs(cormat.ratio[cormat.ratio < 1]), main="")

colNam.vec <- colnames(cormat.ratio) # column names
rowNam.vec <- rownames(cormat.ratio) # row names
samp.exp <- expand.grid(X=colNam.vec, Y=rowNam.vec) # expand to all possible row x column pairs
samp.exp$cor <- 0 # set correlation column

# loop through expanded grid to fill with correlation matrix values
for (i in 1:dim(samp.exp)[1]) {
  x.var <- as.character(samp.exp[i,1])
  y.var <- as.character(samp.exp[i,2])
  samp.exp$cor[i] <- cormat.ratio[which(colnames(cormat.ratio) == x.var), which(rownames(cormat.ratio) == y.var)]
}

# plot heatmap
ggplot(samp.exp, aes(X, Y, fill = cor)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", limits=c(-1,1)) +
  labs(x = "", y = "")



# resampled BRT loop
biter <- 100
eq.sp.points <- 100

# create storage arrays
val.arr <- pred.arr <- array(data = 0, dim = c(eq.sp.points, 5, biter), dimnames=list(paste("x",1:eq.sp.points,sep=""), attr(final.dat.imputed, "names")[c(5:9)], paste("b",1:biter,sep="")))

# create storage vectors
D2.vec <- CV.cor.vec <- CV.cor.se.vec <- GDP.ri <- govExpEdu.ri <- CPI.ri <- PA.ri <- VAag.ri <- rep(0,biter)

for (b in 1:biter) {
  # resample data among countries
  resamp.sub <- sort(sample(x = 1:dim(final.dat.imputed)[1], size = dim(final.dat.imputed)[1], replace=TRUE))
  dat.resamp <- final.dat.imputed[resamp.sub,]
  
  # boosted regression tree
  brt.fit <- gbm.step(dat.resamp, gbm.x = attr(dat.resamp, "names")[c(5:9)], gbm.y = attr(dat.resamp, "names")[4], family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.001, bag.fraction=0.75, tree.complexity = 2, silent=T, tolerance.method = "auto")
  summ.fit <- summary(brt.fit)

  # variable relative importance
  GDP.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(5:9)][1])]
  govExpEdu.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(5:9)][2])]
  CPI.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(5:9)][3])]
  PA.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(5:9)][4])]
  VAag.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(5:9)][5])]
  
  D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
  #D2 # % deviance explained
  D2.vec[b] <- D2
  CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
  CV.cor.vec[b] <- CV.cor
  CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
  CV.cor.se.vec[b] <- CV.cor.se

  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=5)
  ## output average predictions
  for (p in 1:5) {
    RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
  }
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.fit$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.fit$var.names
  
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  print(b)
  
} # end b

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 5
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] > (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  }
  print(k)
}

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)

val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median)

par(mfrow=c(3,2)) 
plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,1]),max(pred.up[,1])), lwd=2, ylab="(←lower) dam:mgm (higher→)", xlab="(←poorer) GDPpc (richer→)" )
lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2, ylab="(←lower) dam:mgm (higher→)", xlab="(←lower) education investment (higher→)" )
lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")

plot(val.med[,3],pred.med[,3],type="l",ylim=c(min(pred.lo[,3]),max(pred.up[,3])), lwd=2, ylab="(←lower) dam:mgm (higher→)", xlab="(←higher) corruption (lower→)" )
lines(val.med[,3], pred.lo[,3], type="l", lty=2, col="red")
lines(val.med[,3], pred.up[,3], type="l", lty=2, col="red")

plot(val.med[,4],pred.med[,4],type="l",ylim=c(min(pred.lo[,4]),max(pred.up[,4])), lwd=2, ylab="(←lower) dam:mgm (higher→)", xlab="(←less) protected area (more→)" )
lines(val.med[,4], pred.lo[,4], type="l", lty=2, col="red")
lines(val.med[,4], pred.up[,4], type="l", lty=2, col="red")

plot(val.med[,5],pred.med[,5],type="l",ylim=c(min(pred.lo[,5]),max(pred.up[,5])), lwd=2, ylab="(←lower) dam:mgm (higher→)", xlab="(←less) agriculture (more→)" )
lines(val.med[,5], pred.lo[,5], type="l", lty=2, col="red")
lines(val.med[,5], pred.up[,5], type="l", lty=2, col="red")
par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update <- D2.vec[1:biter]
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]
GDP.ri.update <- GDP.ri[1:biter]
govExpEdu.ri.update <- govExpEdu.ri[1:biter]
CPI.ri.update <- CPI.ri[1:biter]
PA.ri.update <- PA.ri[1:biter]
VAag.ri.update <- VAag.ri[1:biter]

for (k in 1:kappa.n) {
  D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  GDP.mean <- mean(GDP.ri.update, na.rm=T); GDP.sd <- sd(GDP.ri.update, na.rm=T)
  govExpEdu.mean <- mean(govExpEdu.ri.update, na.rm=T); govExpEdu.sd <- sd(govExpEdu.ri.update, na.rm=T)
  CPI.mean <- mean(CPI.ri.update, na.rm=T); CPI.sd <- sd(CPI.ri.update, na.rm=T)
  PA.mean <- mean(PA.ri.update, na.rm=T); PA.sd <- sd(PA.ri.update, na.rm=T)
  VAag.mean <- mean(VAag.ri.update, na.rm=T); VAag.sd <- sd(VAag.ri.update, na.rm=T)

  for (u in 1:biter) {
    D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    
    GDP.ri.update[u] <- ifelse((GDP.ri.update[u] < (GDP.mean-kappa*GDP.sd) | GDP.ri.update[u] > (GDP.mean+kappa*GDP.sd)), NA, GDP.ri.update[u])
    govExpEdu.ri.update[u] <- ifelse((govExpEdu.ri.update[u] < (govExpEdu.mean-kappa*govExpEdu.sd) | govExpEdu.ri.update[u] > (govExpEdu.mean+kappa*govExpEdu.sd)), NA, govExpEdu.ri.update[u])
    CPI.ri.update[u] <- ifelse((CPI.ri.update[u] < (CPI.mean-kappa*CPI.sd) | CPI.ri.update[u] > (CPI.mean+kappa*CPI.sd)), NA, CPI.ri.update[u])
    PA.ri.update[u] <- ifelse((PA.ri.update[u] < (PA.mean-kappa*PA.sd) | PA.ri.update[u] > (PA.mean+kappa*PA.sd)), NA, PA.ri.update[u])
    VAag.ri.update[u] <- ifelse((VAag.ri.update[u] < (VAag.mean-kappa*VAag.sd) | VAag.ri.update[u] > (VAag.mean+kappa*VAag.sd)), NA, VAag.ri.update[u])
  }
  print(k)
}

D2.med <- median(D2.update, na.rm=TRUE)
D2.lo <- quantile(D2.update, probs=0.025, na.rm=TRUE)
D2.up <- quantile(D2.update, probs=0.975, na.rm=TRUE)
print(c(D2.lo,D2.med,D2.up))

CV.cor.med <- median(CV.cor.update, na.rm=TRUE)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=TRUE)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(CV.cor.lo,CV.cor.med,CV.cor.up))

GDP.ri.lo <- quantile(GDP.ri.update, probs=0.025, na.rm=TRUE)
GDP.ri.med <- median(GDP.ri.update, na.rm=TRUE)
GDP.ri.up <- quantile(GDP.ri.update, probs=0.975, na.rm=TRUE)

govExpEdu.ri.lo <- quantile(govExpEdu.ri.update, probs=0.025, na.rm=TRUE)
govExpEdu.ri.med <- median(govExpEdu.ri.update, na.rm=TRUE)
govExpEdu.ri.up <- quantile(govExpEdu.ri.update, probs=0.975, na.rm=TRUE)

CPI.ri.lo <- quantile(CPI.ri.update, probs=0.025, na.rm=TRUE)
CPI.ri.med <- median(CPI.ri.update, na.rm=TRUE)
CPI.ri.up <- quantile(CPI.ri.update, probs=0.975, na.rm=TRUE)

PA.ri.lo <- quantile(PA.ri.update, probs=0.025, na.rm=TRUE)
PA.ri.med <- median(PA.ri.update, na.rm=TRUE)
PA.ri.up <- quantile(PA.ri.update, probs=0.975, na.rm=TRUE)

VAag.ri.lo <- quantile(VAag.ri.update, probs=0.025, na.rm=TRUE)
VAag.ri.med <- median(VAag.ri.update, na.rm=TRUE)
VAag.ri.up <- quantile(VAag.ri.update, probs=0.975, na.rm=TRUE)

ri.lo <- c(GDP.ri.lo,govExpEdu.ri.lo,CPI.ri.lo,PA.ri.lo,VAag.ri.lo)
ri.med <- c(GDP.ri.med,govExpEdu.ri.med,CPI.ri.med,PA.ri.med,VAag.ri.med)
ri.up <- c(GDP.ri.up,govExpEdu.ri.up,CPI.ri.up,PA.ri.up,VAag.ri.up)

ri.out <- as.data.frame(cbind(ri.lo,ri.med,ri.up))
colnames(ri.out) <- c("ri.lo","ri.med","ri.up")
rownames(ri.out) <- attr(final.dat.imputed, "names")[c(5:9)]
ri.sort <- ri.out[order(ri.out[,2],decreasing=T),1:3]
ri.sort

write.table(pred.med,file="BRT.boot.pred.med.csv",sep=",", row.names = T, col.names = T)
write.table(pred.lo,file="BRT.boot.pred.lo.csv",sep=",", row.names = T, col.names = T)
write.table(pred.up,file="BRT.boot.pred.up.csv",sep=",", row.names = T, col.names = T)
write.table(val.med,file="BRT.boot.val.med.csv",sep=",", row.names = T, col.names = T)
write.table(ri.sort,file="BRT.boot.ri.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.val.dat,file="BRT.val.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.pred.dat,file="BRT.pred.csv",sep=",", row.names = T, col.names = T)



##########################
## phase-based approach ##
##########################
## response: dam:mgm ratio

head(final.dat.imputed)

## PHASE 1 - wealth/capacity
# gdp, cpi, ghsi, igs

brt.fit1 <- gbm.step(final.dat.imputed, gbm.x = attr(final.dat.imputed, "names")[c(5:8)], gbm.y = attr(final.dat.imputed, "names")[4], family="gaussian", max.trees=100000, tolerance = 0.00001, learning.rate = 0.00003, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit1)
D2 <- 100 * (brt.fit1$cv.statistics$deviance.mean - brt.fit1$self.statistics$mean.resid) / brt.fit1$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit1)
gbm.plot.fits(brt.fit1)

brt1.CV.cor <- 100 * brt.fit1$cv.statistics$correlation.mean
brt1.CV.cor
brt1.CV.cor.se <- 100 * brt.fit1$cv.statistics$correlation.se
brt1.CV.cor.se
print(c(brt1.CV.cor, brt1.CV.cor.se))

# keep gdp + igs + cpi


# GLMM
# recode continent to greater region to increase per-category sample sizes
table(final.dat.imputed$cont)
final.dat.imputed$reg2 <- ifelse(final.dat.imputed$cont == "NAM" | final.dat.imputed$cont == "CAR", "NAMCAR", final.dat.imputed$cont)
final.dat.imputed$reg2 <- ifelse(final.dat.imputed$cont == "EUR" | final.dat.imputed$cont == "ME", "EURME", final.dat.imputed$reg2)
final.dat.imputed$reg2 <- ifelse(final.dat.imputed$cont == "ASIA" | final.dat.imputed$cont == "OC", "ASIAOC", final.dat.imputed$reg2)
table(final.dat.imputed$reg2)

# model set
m1 <- "ratio ~ gdp + cpi + ghsi + igs + (1|reg2)"

m2 <- "ratio ~ gdp + cpi + ghsi + (1|reg2)"
m3 <- "ratio ~ gdp + cpi + igs + (1|reg2)"
m4 <- "ratio ~ gdp + ghsi + igs + (1|reg2)"
m5 <- "ratio ~ cpi + ghsi + igs + (1|reg2)"

m6 <- "ratio ~ gdp + cpi + (1|reg2)"
m7 <- "ratio ~ gdp + ghsi + (1|reg2)"
m8 <- "ratio ~ gdp + igs + (1|reg2)"
m9 <- "ratio ~ cpi + ghsi + (1|reg2)"
m10 <- "ratio ~ cpi + igs + (1|reg2)"
m11 <- "ratio ~ ghsi + igs + (1|reg2)"

m12 <- "ratio ~ gdp + (1|reg2)"
m13 <- "ratio ~ cpi + (1|reg2)"
m14 <- "ratio ~ ghsi + (1|reg2)"
m15 <- "ratio ~ igs + (1|reg2)"

m16 <- "ratio ~ 1 + (1|reg2)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16)
length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.dat.imputed, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable1 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable1) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable1) <- as.character(mod.vec)
summary.table1 <- sumtable[order(sumtable1[,5],decreasing=F),]
summary.table1

# keep gdp + igs + cpi


## PHASE 2 - primary production
# agrL, VAag
head(final.dat.imputed)

brt.fit2 <- gbm.step(final.dat.imputed, gbm.x = attr(final.dat.imputed, "names")[c(9:10)], gbm.y = attr(final.dat.imputed, "names")[4], family="gaussian", max.trees=100000, tolerance = 0.000001, learning.rate = 0.00003, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit2)
D2 <- 100 * (brt.fit2$cv.statistics$deviance.mean - brt.fit2$self.statistics$mean.resid) / brt.fit2$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit2)
gbm.plot.fits(brt.fit2)

brt2.CV.cor <- 100 * brt.fit2$cv.statistics$correlation.mean
brt2.CV.cor
brt2.CV.cor.se <- 100 * brt.fit2$cv.statistics$correlation.se
brt2.CV.cor.se
print(c(brt2.CV.cor, brt2.CV.cor.se))

# keep both?


# GLMM
# model set
m1 <- "ratio ~ agrL + VAag + (1|reg2)"

m2 <- "ratio ~ agrL + (1|reg2)"
m3 <- "ratio ~ VAag + (1|reg2)"

m4 <- "ratio ~ 1 + (1|reg2)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)
length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.dat.imputed, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable2 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable2) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable2) <- as.character(mod.vec)
summary.table2 <- sumtable2[order(sumtable2[,5],decreasing=F),]
summary.table2

# keep agrL + VAag


## PHASE 3 - education/research
# govexpedu, stja
head(final.dat.imputed)

brt.fit3 <- gbm.step(final.dat.imputed, gbm.x = attr(final.dat.imputed, "names")[c(11:12)], gbm.y = attr(final.dat.imputed, "names")[4], family="gaussian", max.trees=100000, tolerance = 0.00001, learning.rate = 0.00003, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit3)
D2 <- 100 * (brt.fit3$cv.statistics$deviance.mean - brt.fit3$self.statistics$mean.resid) / brt.fit3$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit3)
gbm.plot.fits(brt.fit3)

brt3.CV.cor <- 100 * brt.fit3$cv.statistics$correlation.mean
brt3.CV.cor
brt3.CV.cor.se <- 100 * brt.fit3$cv.statistics$correlation.se
brt3.CV.cor.se
print(c(brt3.CV.cor, brt3.CV.cor.se))

# keep both?


# GLMM
# model set
m1 <- "ratio ~ govexpedu + stja + (1|reg2)"

m2 <- "ratio ~ govexpedu + (1|reg2)"
m3 <- "ratio ~ stja + (1|reg2)"

m4 <- "ratio ~ 1 + (1|reg2)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)
length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.dat.imputed, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable3 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable3) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable3) <- as.character(mod.vec)
summary.table3 <- sumtable3[order(sumtable3[,5],decreasing=F),]
summary.table3

# keep neither


## COMBINE PHASES
# gdp, igs, cpi, VAag, agrL

# GLMM
# model set
m1 <- "ratio ~ gdp + igs + cpi + VAag + agrL + (1|reg2)"

m2 <- "ratio ~ gdp + igs + cpi + VAag + (1|reg2)"
m3 <- "ratio ~ gdp + igs + cpi + agrL + (1|reg2)"
m4 <- "ratio ~ gdp + igs + VAag + agrL + (1|reg2)"
m5 <- "ratio ~ gdp + cpi + VAag + agrL + (1|reg2)"
m6 <- "ratio ~ igs + cpi + VAag + agrL + (1|reg2)"

m7 <- "ratio ~ gdp + igs + cpi + (1|reg2)"
m8 <- "ratio ~ gdp + igs + agrL + (1|reg2)"
m9 <- "ratio ~ gdp + VAag + agrL + (1|reg2)"
m10 <- "ratio ~ gdp + cpi + VAag + (1|reg2)"
m11 <- "ratio ~ gdp + igs + VAag + (1|reg2)"
m12 <- "ratio ~ igs + cpi + VAag + (1|reg2)"
m13 <- "ratio ~ igs + VAag + agrL + (1|reg2)"
m14 <- "ratio ~ igs + cpi + agrL + (1|reg2)"
m15 <- "ratio ~ cpi + VAag + agrL + (1|reg2)"

m16 <- "ratio ~ gdp + igs + (1|reg2)"
m17 <- "ratio ~ gdp + cpi + (1|reg2)"
m18 <- "ratio ~ gdp + VAag + (1|reg2)"
m19 <- "ratio ~ gdp + agrL + (1|reg2)"
m20 <- "ratio ~ igs + cpi + (1|reg2)"
m21 <- "ratio ~ igs + VAag + (1|reg2)"
m22 <- "ratio ~ igs + agrL + (1|reg2)"
m23 <- "ratio ~ cpi + VAag + (1|reg2)"
m24 <- "ratio ~ cpi + agrL + (1|reg2)"

m25 <- "ratio ~ gdp + (1|reg2)"
m26 <- "ratio ~ igs + (1|reg2)"
m27 <- "ratio ~ cpi + (1|reg2)"
m28 <- "ratio ~ VAag + (1|reg2)"
m29 <- "ratio ~ agrL + (1|reg2)"

m30 <- "ratio ~ 1 + (1|reg2)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30)
length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.dat.imputed, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtableF <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtableF) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtableF) <- as.character(mod.vec)
summary.tableF <- sumtableF[order(sumtableF[,5],decreasing=F),]
summary.tableF


## BRT
# gdp, igs, cpi, VAag, agrL
head(final.dat.imputed)

brt.fitF <- gbm.step(final.dat.imputed, gbm.x = attr(final.dat.imputed, "names")[c(5,6,8,9,10)], gbm.y = attr(final.dat.imputed, "names")[4], family="gaussian", max.trees=100000, tolerance = 0.00001, learning.rate = 0.00003, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fitF)
D2 <- 100 * (brt.fitF$cv.statistics$deviance.mean - brt.fitF$self.statistics$mean.resid) / brt.fitF$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fitF)
gbm.plot.fits(brt.fitF)

brtF.CV.cor <- 100 * brt.fitF$cv.statistics$correlation.mean
brtF.CV.cor
brtF.CV.cor.se <- 100 * brt.fitF$cv.statistics$correlation.se
brtF.CV.cor.se
print(c(brtF.CV.cor, brtF.CV.cor.se))


# resampled BRT loop
# gdp, igs, cpi, VAag, agrL
biter <- 100
eq.sp.points <- 100

# create storage arrays
val.arr <- pred.arr <- array(data = NA, dim = c(eq.sp.points, 5, biter), dimnames=list(paste("x",1:eq.sp.points,sep=""), attr(final.dat.imputed, "names")[c(5,6,8,9,10)], paste("b",1:biter,sep="")))

# create storage vectors
D2.vec <- CV.cor.vec <- CV.cor.se.vec <- gdp.ri <- igs.ri <- cpi.ri <- VAag.ri <- agrL.ri <- rep(NA,biter)

for (b in 1:biter) {
  # resample data among countries
  resamp.sub <- sort(sample(x = 1:dim(final.dat.imputed)[1], size = dim(final.dat.imputed)[1], replace=TRUE))
  dat.resamp <- final.dat.imputed[resamp.sub,]
  
  # boosted regression tree
  brt.fit <- gbm.step(dat.resamp, gbm.x = attr(dat.resamp, "names")[c(5,6,8,9,10)], gbm.y = attr(dat.resamp, "names")[4], family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.0001, bag.fraction=0.75, tree.complexity = 2, silent=T, tolerance.method = "auto")
  summ.fit <- summary(brt.fit)

  length(summ.fit[[1]])
  
  if (length(summ.fit[[1]]) == 5) {
    # variable relative importance
    gdp.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(5,6,8,9,10)][1])]
    cpi.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(5,6,8,9,10)][2])]
    igs.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(5,6,8,9,10)][3])]
    agrL.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(5,6,8,9,10)][4])]
    VAag.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(5,6,8,9,10)][5])]
    
    D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
    D2.vec[b] <- D2
    CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
    CV.cor.vec[b] <- CV.cor
    CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
    CV.cor.se.vec[b] <- CV.cor.se
    
    RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=5)
    ## output average predictions
    for (p in 1:5) {
      RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
      RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
    }
    RESP.val.dat <- as.data.frame(RESP.val)
    colnames(RESP.val.dat) <- brt.fit$var.names
    RESP.pred.dat <- as.data.frame(RESP.pred)
    colnames(RESP.pred.dat) <- brt.fit$var.names
    
    val.arr[, , b] <- as.matrix(RESP.val.dat)
    pred.arr[, , b] <- as.matrix(RESP.pred.dat)
    
    print(b)
  }
  
  if (length(summ.fit[[1]]) != 5) {
    b <- b+1
    print(b)
  }

} # end b

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 5
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] > (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  }
  print(k)
}

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)

val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

par(mfrow=c(3,2)) 
plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,1]),max(pred.up[,1])), lwd=2, ylab="(←lower) dam:mgm (higher→)", xlab="(←poorer) per capita GDP (richer→)")
lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2, ylab="(←lower) dam:mgm (higher→)", xlab="(←higher) corruption (lower→)")
lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")

plot(val.med[,3],pred.med[,3],type="l",ylim=c(min(pred.lo[,3]),max(pred.up[,3])), lwd=2, ylab="(←lower) dam:mgm (higher→)",  xlab="(←lower) per capita imports goods & services (higher→)" )
lines(val.med[,3], pred.lo[,3], type="l", lty=2, col="red")
lines(val.med[,3], pred.up[,3], type="l", lty=2, col="red")

plot(val.med[,4],pred.med[,4],type="l",ylim=c(min(pred.lo[,4]),max(pred.up[,4])), lwd=2, ylab="(←lower) dam:mgm (higher→)", xlab="(←less) agricultural land (more→)")
lines(val.med[,4], pred.lo[,4], type="l", lty=2, col="red")
lines(val.med[,4], pred.up[,4], type="l", lty=2, col="red")

plot(val.med[,5],pred.med[,5],type="l",ylim=c(min(pred.lo[,5]),max(pred.up[,5])), lwd=2, ylab="(←lower) dam:mgm (higher→)", xlab="(←less) value-added agric (more→)")
lines(val.med[,5], pred.lo[,5], type="l", lty=2, col="red")
lines(val.med[,5], pred.up[,5], type="l", lty=2, col="red")
par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update <- D2.vec[1:biter]
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]
gdp.ri.update <- gdp.ri[1:biter]
cpi.ri.update <- cpi.ri[1:biter]
igs.ri.update <- igs.ri[1:biter]
agrL.ri.update <- agrL.ri[1:biter]
VAag.ri.update <- VAag.ri[1:biter]

for (k in 1:kappa.n) {
  D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  gdp.mean <- mean(gdp.ri.update, na.rm=T); gdp.sd <- sd(gdp.ri.update, na.rm=T)
  cpi.mean <- mean(cpi.ri.update, na.rm=T); cpi.sd <- sd(cpi.ri.update, na.rm=T)
  igs.mean <- mean(igs.ri.update, na.rm=T); igs.sd <- sd(igs.ri.update, na.rm=T)
  agrL.mean <- mean(agrL.ri.update, na.rm=T); agrL.sd <- sd(agrL.ri.update, na.rm=T)
  VAag.mean <- mean(VAag.ri.update, na.rm=T); VAag.sd <- sd(VAag.ri.update, na.rm=T)
  
  for (u in 1:biter) {
    D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    
    gdp.ri.update[u] <- ifelse((gdp.ri.update[u] < (gdp.mean-kappa*gdp.sd) | gdp.ri.update[u] > (gdp.mean+kappa*gdp.sd)), NA, gdp.ri.update[u])
    cpi.ri.update[u] <- ifelse((cpi.ri.update[u] < (cpi.mean-kappa*cpi.sd) | cpi.ri.update[u] > (cpi.mean+kappa*cpi.sd)), NA, cpi.ri.update[u])
    igs.ri.update[u] <- ifelse((igs.ri.update[u] < (igs.mean-kappa*igs.sd) | igs.ri.update[u] > (igs.mean+kappa*igs.sd)), NA, igs.ri.update[u])
    agrL.ri.update[u] <- ifelse((agrL.ri.update[u] < (agrL.mean-kappa*agrL.sd) | agrL.ri.update[u] > (agrL.mean+kappa*agrL.sd)), NA, agrL.ri.update[u])
    VAag.ri.update[u] <- ifelse((VAag.ri.update[u] < (VAag.mean-kappa*VAag.sd) | VAag.ri.update[u] > (VAag.mean+kappa*VAag.sd)), NA, VAag.ri.update[u])
  }
  
  print(k)
}

D2.med <- median(D2.update, na.rm=TRUE)
D2.lo <- quantile(D2.update, probs=0.025, na.rm=TRUE)
D2.up <- quantile(D2.update, probs=0.975, na.rm=TRUE)
print(c(D2.lo,D2.med,D2.up))

CV.cor.med <- median(CV.cor.update, na.rm=TRUE)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=TRUE)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(CV.cor.lo,CV.cor.med,CV.cor.up))

gdp.ri.lo <- quantile(gdp.ri.update, probs=0.025, na.rm=TRUE)
gdp.ri.med <- median(gdp.ri.update, na.rm=TRUE)
gdp.ri.up <- quantile(gdp.ri.update, probs=0.975, na.rm=TRUE)

cpi.ri.lo <- quantile(cpi.ri.update, probs=0.025, na.rm=TRUE)
cpi.ri.med <- median(cpi.ri.update, na.rm=TRUE)
cpi.ri.up <- quantile(cpi.ri.update, probs=0.975, na.rm=TRUE)

igs.ri.lo <- quantile(igs.ri.update, probs=0.025, na.rm=TRUE)
igs.ri.med <- median(igs.ri.update, na.rm=TRUE)
igs.ri.up <- quantile(igs.ri.update, probs=0.975, na.rm=TRUE)

agrL.ri.lo <- quantile(agrL.ri.update, probs=0.025, na.rm=TRUE)
agrL.ri.med <- median(agrL.ri.update, na.rm=TRUE)
agrL.ri.up <- quantile(agrL.ri.update, probs=0.975, na.rm=TRUE)

VAag.ri.lo <- quantile(VAag.ri.update, probs=0.025, na.rm=TRUE)
VAag.ri.med <- median(VAag.ri.update, na.rm=TRUE)
VAag.ri.up <- quantile(VAag.ri.update, probs=0.975, na.rm=TRUE)

ri.lo <- c(gdp.ri.lo,cpi.ri.lo,igs.ri.lo,agrL.ri.lo,VAag.ri.lo)
ri.med <- c(gdp.ri.med,cpi.ri.med,igs.ri.med,agrL.ri.med,VAag.ri.med)
ri.up <- c(gdp.ri.up,cpi.ri.up,igs.ri.up,agrL.ri.up,VAag.ri.up)

ri.out <- as.data.frame(cbind(ri.lo,ri.med,ri.up))
colnames(ri.out) <- c("ri.lo","ri.med","ri.up")
rownames(ri.out) <- attr(final.dat.imputed, "names")[c(5,6,8,9,10)]
ri.sort <- ri.out[order(ri.out[,2],decreasing=T),1:3]
ri.sort


write.table(pred.med,file="BRT.boot.pred.med.csv",sep=",", row.names = T, col.names = T)
write.table(pred.lo,file="BRT.boot.pred.lo.csv",sep=",", row.names = T, col.names = T)
write.table(pred.up,file="BRT.boot.pred.up.csv",sep=",", row.names = T, col.names = T)
write.table(val.med,file="BRT.boot.val.med.csv",sep=",", row.names = T, col.names = T)
write.table(ri.sort,file="BRT.boot.ri.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.val.dat,file="BRT.val.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.pred.dat,file="BRT.pred.csv",sep=",", row.names = T, col.names = T)



## response: pdam

head(final.datPdam.imputed)

## PHASE 1 - wealth/capacity
# gdp, cpi, ghsi, igs

brt.fit4 <- gbm.step(final.datPdam.imputed, gbm.x = attr(final.datPdam.imputed, "names")[c(5:8)], gbm.y = attr(final.datPdam.imputed, "names")[4], family="gaussian", max.trees=100000, tolerance = 0.00003, learning.rate = 0.00001, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit4)
D2 <- 100 * (brt.fit4$cv.statistics$deviance.mean - brt.fit4$self.statistics$mean.resid) / brt.fit4$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit4)
gbm.plot.fits(brt.fit4)

brt4.CV.cor <- 100 * brt.fit4$cv.statistics$correlation.mean
brt4.CV.cor
brt4.CV.cor.se <- 100 * brt.fit4$cv.statistics$correlation.se
brt4.CV.cor.se
print(c(brt4.CV.cor, brt4.CV.cor.se))

# keep gdp

# GLMM
# recode continent to greater region to increase per-category sample sizes
table(final.datPdam.imputed$cont)
final.datPdam.imputed$reg2 <- ifelse(final.datPdam.imputed$cont == "NAM" | final.datPdam.imputed$cont == "CAR", "NAMCAR", final.datPdam.imputed$cont)
final.datPdam.imputed$reg2 <- ifelse(final.datPdam.imputed$cont == "EUR" | final.datPdam.imputed$cont == "ME", "EURME", final.datPdam.imputed$reg2)
final.datPdam.imputed$reg2 <- ifelse(final.datPdam.imputed$cont == "ASIA" | final.datPdam.imputed$cont == "OC", "ASIAOC", final.datPdam.imputed$reg2)
table(final.datPdam.imputed$reg2)

# model set
m1 <- "pdam ~ gdp + cpi + ghsi + igs + (1|reg2)"

m2 <- "pdam ~ gdp + cpi + ghsi + (1|reg2)"
m3 <- "pdam ~ gdp + cpi + igs + (1|reg2)"
m4 <- "pdam ~ gdp + ghsi + igs + (1|reg2)"
m5 <- "pdam ~ cpi + ghsi + igs + (1|reg2)"

m6 <- "pdam ~ gdp + cpi + (1|reg2)"
m7 <- "pdam ~ gdp + ghsi + (1|reg2)"
m8 <- "pdam ~ gdp + igs + (1|reg2)"
m9 <- "pdam ~ cpi + ghsi + (1|reg2)"
m10 <- "pdam ~ cpi + igs + (1|reg2)"
m11 <- "pdam ~ ghsi + igs + (1|reg2)"

m12 <- "pdam ~ gdp + (1|reg2)"
m13 <- "pdam ~ cpi + (1|reg2)"
m14 <- "pdam ~ ghsi + (1|reg2)"
m15 <- "pdam ~ igs + (1|reg2)"

m16 <- "pdam ~ 1 + (1|reg2)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16)
length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.datPdam.imputed, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable4 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable4) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable4) <- as.character(mod.vec)
summary.table4 <- sumtable4[order(sumtable4[,5],decreasing=F),]
summary.table4

# keep gdp


## PHASE 2 - primary production
# agrL, VAag
head(final.datPdam.imputed)
hist(final.datPdam.imputed$pdam)

brt.fit2 <- gbm.step(final.datPdam.imputed, gbm.x = attr(final.datPdam.imputed, "names")[c(9:10)], gbm.y = attr(final.datPdam.imputed, "names")[4], family="gaussian", max.trees=100000, tolerance = 0.00001, learning.rate = 0.000005, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit2)
D2 <- 100 * (brt.fit2$cv.statistics$deviance.mean - brt.fit2$self.statistics$mean.resid) / brt.fit2$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit2)
gbm.plot.fits(brt.fit2)

brt2.CV.cor <- 100 * brt.fit2$cv.statistics$correlation.mean
brt2.CV.cor
brt2.CV.cor.se <- 100 * brt.fit2$cv.statistics$correlation.se
brt2.CV.cor.se
print(c(brt2.CV.cor, brt2.CV.cor.se))

# keep agrL


# GLMM
# model set
m1 <- "pdam ~ agrL + VAag + (1|reg2)"

m2 <- "pdam ~ agrL + (1|reg2)"
m3 <- "pdam ~ VAag + (1|reg2)"

m4 <- "pdam ~ 1 + (1|reg2)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)
length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.datPdam.imputed, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable) <- as.character(mod.vec)
summary.table <- sumtable[order(sumtable[,5],decreasing=F),]
summary.table

# keep agrL


## PHASE 3 - education/research
# govexpedu, stja
head(final.datPdam.imputed)

brt.fit3 <- gbm.step(final.datPdam.imputed, gbm.x = attr(final.datPdam.imputed, "names")[c(11:12)], gbm.y = attr(final.datPdam.imputed, "names")[4], family="gaussian", max.trees=100000, tolerance = 0.00001, learning.rate = 0.00003, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit3)
D2 <- 100 * (brt.fit3$cv.statistics$deviance.mean - brt.fit3$self.statistics$mean.resid) / brt.fit3$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit3)
gbm.plot.fits(brt.fit3)

brt3.CV.cor <- 100 * brt.fit3$cv.statistics$correlation.mean
brt3.CV.cor
brt3.CV.cor.se <- 100 * brt.fit3$cv.statistics$correlation.se
brt3.CV.cor.se
print(c(brt3.CV.cor, brt3.CV.cor.se))

# keep both?


# GLMM
# model set
m1 <- "pdam ~ govexpedu + stja + (1|reg2)"

m2 <- "pdam ~ govexpedu + (1|reg2)"
m3 <- "pdam ~ stja + (1|reg2)"

m4 <- "pdam ~ 1 + (1|reg2)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)
length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.datPdam.imputed, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable) <- as.character(mod.vec)
summary.table <- sumtable[order(sumtable[,5],decreasing=F),]
summary.table

# keep govexpedu + stja


## COMBINE PHASES
# gdp, agrL, govexpedu, stja

head(final.datPdam.imputed)

brt.fitF2 <- gbm.step(final.datPdam.imputed, gbm.x = attr(final.datPdam.imputed, "names")[c(5,9,11,12)], gbm.y = attr(final.datPdam.imputed, "names")[4], family="gaussian", max.trees=100000, tolerance = 0.00001, learning.rate = 0.00003, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fitF2)
D2 <- 100 * (brt.fitF2$cv.statistics$deviance.mean - brt.fitF2$self.statistics$mean.resid) / brt.fitF2$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fitF2)
gbm.plot.fits(brt.fitF2)

brtF2.CV.cor <- 100 * brt.fitF2$cv.statistics$correlation.mean
brtF2.CV.cor
brtF2.CV.cor.se <- 100 * brt.fitF2$cv.statistics$correlation.se
brtF2.CV.cor.se
print(c(brtF2.CV.cor, brtF2.CV.cor.se))


# GLMM
# model set
m1 <- "pdam ~ gdp + agrL + govexpedu + stja + (1|reg2)"

m2 <- "pdam ~ gdp + agrL + govexpedu + (1|reg2)"
m3 <- "pdam ~ gdp + agrL + stja + (1|reg2)"
m4 <- "pdam ~ gdp + govexpedu + stja + (1|reg2)"
m5 <- "pdam ~ agrL + govexpedu + stja + (1|reg2)"

m6 <- "pdam ~ gdp + agrL + (1|reg2)"
m7 <- "pdam ~ gdp + govexpedu + (1|reg2)"
m8 <- "pdam ~ gdp + stja + (1|reg2)"
m9 <- "pdam ~ agrL + govexpedu + (1|reg2)"
m10 <- "pdam ~ agrL + stja + (1|reg2)"
m11 <- "pdam ~ govexpedu + stja + (1|reg2)"

m12 <- "pdam ~ gdp + (1|reg2)"
m13 <- "pdam ~ agrL + (1|reg2)"
m14 <- "pdam ~ govexpedu + (1|reg2)"
m15 <- "pdam ~ stja + (1|reg2)"

m16 <- "pdam ~ 1 + (1|reg2)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16)
length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.datPdam.imputed, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtableF2 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtableF2) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtableF2) <- as.character(mod.vec)
summary.tableF2 <- sumtableF2[order(sumtableF2[,5],decreasing=F),]
summary.tableF2


# resampled BRT loop
# gdp, agrL, govexpedu, stja
attr(final.datPdam.imputed, "names")[c(5,9,11,12)]

biter <- 100
eq.sp.points <- 100

# create storage arrays
val.arr2 <- pred.arr2 <- array(data = NA, dim = c(eq.sp.points, 4, biter), dimnames=list(paste("x",1:eq.sp.points,sep=""), attr(final.datPdam.imputed, "names")[c(5,9,11,12)], paste("b",1:biter,sep="")))

# create storage vectors
D2.vec2 <- CV.cor.vec2 <- CV.cor.se.vec2 <- gdp.ri2 <- agrL.ri2 <- govexpedu.ri2 <- stja.ri2 <- rep(NA,biter)

for (b in 1:biter) {
  # resample data among countries
  resamp.sub <- sort(sample(x = 1:dim(final.datPdam.imputed)[1], size = dim(final.datPdam.imputed)[1], replace=TRUE))
  dat.resamp2 <- final.datPdam.imputed[resamp.sub,]
  
  # boosted regression tree
  brt.fit2 <- gbm.step(dat.resamp2, gbm.x = attr(dat.resamp2, "names")[c(5,9,11,12)], gbm.y = attr(dat.resamp2, "names")[4], family="gaussian", max.trees=100000, tolerance = 0.00001, learning.rate = 0.00003, bag.fraction=0.75, tree.complexity = 2, silent=T, tolerance.method = "auto")
  summ.fit2 <- summary(brt.fit2)
  
  
  if (length(summ.fit2[[1]]) == 4) {
    # variable relative importance
    gdp.ri2[b] <- summ.fit2$rel.inf[which(summ.fit2$var == attr(dat.resamp2, "names")[c(5,9,11,12)][1])]
    agrL.ri2[b] <- summ.fit2$rel.inf[which(summ.fit2$var == attr(dat.resamp2, "names")[c(5,9,11,12)][2])]
    govexpedu.ri2[b] <- summ.fit2$rel.inf[which(summ.fit2$var == attr(dat.resamp2, "names")[c(5,9,11,12)][3])]
    stja.ri2[b] <- summ.fit2$rel.inf[which(summ.fit2$var == attr(dat.resamp2, "names")[c(5,9,11,12)][4])]

    D22 <- 100 * (brt.fit2$cv.statistics$deviance.mean - brt.fit2$self.statistics$mean.resid) / brt.fit2$cv.statistics$deviance.mean
    D2.vec2[b] <- D22
    CV.cor2 <- 100 * brt.fit2$cv.statistics$correlation.mean
    CV.cor.vec2[b] <- CV.cor2
    CV.cor.se2 <- 100 *brt.fit2$cv.statistics$correlation.se
    CV.cor.se.vec2[b] <- CV.cor.se2
    
    RESP.val2 <- RESP.pred2 <- matrix(data=NA, nrow=eq.sp.points, ncol=4)
    ## output average predictions
    for (p in 1:4) {
      RESP.val2[,p] <- plot.gbm(brt.fit2, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
      RESP.pred2[,p] <- plot.gbm(brt.fit2, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
    }
    RESP.val.dat2 <- as.data.frame(RESP.val2)
    colnames(RESP.val.dat2) <- brt.fit2$var.names
    RESP.pred.dat2 <- as.data.frame(RESP.pred2)
    colnames(RESP.pred.dat2) <- brt.fit2$var.names
    
    val.arr2[, , b] <- as.matrix(RESP.val.dat2)
    pred.arr2[, , b] <- as.matrix(RESP.pred.dat2)
    
    print(b)
    
  }
  
  if (length(summ.fit2[[1]]) != 4) {
    b <- b+1
    print(b)
  }
  
} # end b

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 5
pred.update2 <- pred.arr2[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean2 <- apply(pred.update2, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd2 <- apply(pred.update2, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update2[,,z] <- ifelse((pred.update2[,,z] < (boot.mean2-kappa*boot.sd2) | pred.update2[,,z] > (boot.mean2+kappa*boot.sd2)), NA, pred.update2[,,z])
  }
  print(k)
}

pred.med2 <- apply(pred.update2, MARGIN=c(1,2), median, na.rm=T)
pred.lo2 <- apply(pred.update2, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up2 <- apply(pred.update2, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)

val.med2 <- apply(val.arr2[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

par(mfrow=c(2,2)) 
plot(val.med2[,1],pred.med2[,1],type="l",ylim=c(min(pred.lo2[,1]),max(pred.up2[,1])), lwd=2, ylab="(←lower) pdam (higher→)", xlab="(←poorer) per capita GDP (richer→)")
lines(val.med2[,1], pred.lo2[,1], type="l", lty=2, col="red")
lines(val.med2[,1], pred.up2[,1], type="l", lty=2, col="red")

plot(val.med2[,2],pred.med2[,2],type="l",ylim=c(min(pred.lo2[,2]),max(pred.up2[,2])), lwd=2, ylab="(←lower) pdam (higher→)", xlab="(←lower) agricultural land (higher→)")
lines(val.med2[,2], pred.lo2[,2], type="l", lty=2, col="red")
lines(val.med2[,2], pred.up2[,2], type="l", lty=2, col="red")

plot(val.med2[,3],pred.med2[,3],type="l",ylim=c(min(pred.lo2[,3]),max(pred.up2[,3])), lwd=2, ylab="(←lower) pdam (higher→)",  xlab="(←less) gov expenditure education (more→)" )
lines(val.med2[,3], pred.lo2[,3], type="l", lty=2, col="red")
lines(val.med2[,3], pred.up2[,3], type="l", lty=2, col="red")

plot(val.med2[,4],pred.med2[,4],type="l",ylim=c(min(pred.lo2[,4]),max(pred.up2[,4])), lwd=2, ylab="(←lower) pdam (higher→)", xlab="(←less) sci/tech journal articles (more→)")
lines(val.med2[,4], pred.lo2[,4], type="l", lty=2, col="red")
lines(val.med2[,4], pred.up2[,4], type="l", lty=2, col="red")

par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update2 <- D2.vec2[1:biter]
CV.cor.update2 <- CV.cor.vec2[1:biter]
CV.cor.se.update2 <- CV.cor.se.vec2[1:biter]
gdp.ri.update2 <- gdp.ri2[1:biter]
agrL.ri.update2 <- agrL.ri2[1:biter]
govexpedu.ri.update2 <- govexpedu.ri2[1:biter]
stja.ri.update2 <- stja.ri2[1:biter]

for (k in 1:kappa.n) {
  D2.mean2 <- mean(D2.update2, na.rm=T); D2.sd2 <- sd(D2.update2, na.rm=T)
  CV.cor.mean2 <- mean(CV.cor.update2, na.rm=T); CV.cor.sd2 <- sd(CV.cor.update2, na.rm=T)
  CV.cor.se.mean2 <- mean(CV.cor.se.update2, na.rm=T); CV.cor.se.sd2 <- sd(CV.cor.se.update2, na.rm=T)
  
  gdp.mean2 <- mean(gdp.ri.update2, na.rm=T); gdp.sd2 <- sd(gdp.ri.update2, na.rm=T)
  agrL.mean2 <- mean(agrL.ri.update2, na.rm=T); agrL.sd2 <- sd(agrL.ri.update2, na.rm=T)
  govexpedu.mean2 <- mean(govexpedu.ri.update2, na.rm=T); govexpedu.sd2 <- sd(govexpedu.ri.update2, na.rm=T)
  stja.mean2 <- mean(stja.ri.update2, na.rm=T); stja.sd2 <- sd(stja.ri.update2, na.rm=T)

  for (u in 1:biter) {
    D2.update2[u] <- ifelse((D2.update2[u] < (D2.mean2-kappa*D2.sd2) | D2.update2[u] > (D2.mean2+kappa*D2.sd2)), NA, D2.update2[u])
    CV.cor.update2[u] <- ifelse((CV.cor.update2[u] < (CV.cor.mean2-kappa*CV.cor.sd2) | CV.cor.update2[u] > (CV.cor.mean2+kappa*CV.cor.sd2)), NA, CV.cor.update2[u])
    CV.cor.se.update2[u] <- ifelse((CV.cor.se.update2[u] < (CV.cor.se.mean2-kappa*CV.cor.se.sd2) | CV.cor.se.update2[u] > (CV.cor.se.mean2+kappa*CV.cor.se.sd2)), NA, CV.cor.se.update2[u])
    
    gdp.ri.update2[u] <- ifelse((gdp.ri.update2[u] < (gdp.mean2-kappa*gdp.sd2) | gdp.ri.update2[u] > (gdp.mean2+kappa*gdp.sd2)), NA, gdp.ri.update2[u])
    agrL.ri.update2[u] <- ifelse((agrL.ri.update2[u] < (agrL.mean2-kappa*agrL.sd2) | agrL.ri.update2[u] > (agrL.mean2+kappa*agrL.sd2)), NA, agrL.ri.update2[u])
    govexpedu.ri.update2[u] <- ifelse((govexpedu.ri.update2[u] < (govexpedu.mean2-kappa*govexpedu.sd2) | govexpedu.ri.update2[u] > (govexpedu.mean2+kappa*govexpedu.sd2)), NA, govexpedu.ri.update2[u])
    stja.ri.update2[u] <- ifelse((stja.ri.update2[u] < (stja.mean2-kappa*stja.sd2) | stja.ri.update2[u] > (stja.mean2+kappa*agrL.sd2)), NA, stja.ri.update2[u])
  }
  print(k)
}

D2.med2 <- median(D2.update2, na.rm=TRUE)
D2.lo2 <- quantile(D2.update2, probs=0.025, na.rm=TRUE)
D2.up2 <- quantile(D2.update2, probs=0.975, na.rm=TRUE)
print(c(D2.lo2,D2.med2,D2.up2))

CV.cor.med2 <- median(CV.cor.update2, na.rm=TRUE)
CV.cor.lo2 <- quantile(CV.cor.update2, probs=0.025, na.rm=TRUE)
CV.cor.up2 <- quantile(CV.cor.update2, probs=0.975, na.rm=TRUE)
print(c(CV.cor.lo2,CV.cor.med2,CV.cor.up2))

gdp.ri.lo <- quantile(gdp.ri.update2, probs=0.025, na.rm=TRUE)
gdp.ri.med <- median(gdp.ri.update2, na.rm=TRUE)
gdp.ri.up <- quantile(gdp.ri.update2, probs=0.975, na.rm=TRUE)

agrL.ri.lo <- quantile(agrL.ri.update2, probs=0.025, na.rm=TRUE)
agrL.ri.med <- median(agrL.ri.update2, na.rm=TRUE)
agrL.ri.up <- quantile(agrL.ri.update2, probs=0.975, na.rm=TRUE)

govexpedu.ri.lo <- quantile(govexpedu.ri.update2, probs=0.025, na.rm=TRUE)
govexpedu.ri.med <- median(govexpedu.ri.update2, na.rm=TRUE)
govexpedu.ri.up <- quantile(govexpedu.ri.update2, probs=0.975, na.rm=TRUE)

stja.ri.lo <- quantile(stja.ri.update2, probs=0.025, na.rm=TRUE)
stja.ri.med <- median(stja.ri.update2, na.rm=TRUE)
stja.ri.up <- quantile(stja.ri.update2, probs=0.975, na.rm=TRUE)

ri.lo2 <- c(gdp.ri.lo,agrL.ri.lo,govexpedu.ri.lo,stja.ri.lo)
ri.med2 <- c(gdp.ri.med,agrL.ri.med,govexpedu.ri.med,stja.ri.med)
ri.up2 <- c(gdp.ri.up,agrL.ri.up,govexpedu.ri.up,stja.ri.up)

ri.out2 <- as.data.frame(cbind(ri.lo2,ri.med2,ri.up2))
colnames(ri.out2) <- c("ri.lo","ri.med","ri.up")
rownames(ri.out2) <- attr(final.datPdam.imputed, "names")[c(5,9,11,12)]
ri.sort2 <- ri.out2[order(ri.out2[,2],decreasing=T),1:3]
ri.sort2

write.table(pred.med2,file="BRT.boot.pred.med2.csv",sep=",", row.names = T, col.names = T)
write.table(pred.lo2,file="BRT.boot.pred.lo2.csv",sep=",", row.names = T, col.names = T)
write.table(pred.up2,file="BRT.boot.pred.up2.csv",sep=",", row.names = T, col.names = T)
write.table(val.med2,file="BRT.boot.val.med2.csv",sep=",", row.names = T, col.names = T)
write.table(ri.sort2,file="BRT.boot.ri2.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.val.dat2,file="BRT.val2.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.pred.dat2,file="BRT.pred2.csv",sep=",", row.names = T, col.names = T)



##################################
## general least-squares models ##
##################################

# get world map
wmap <- getMap(resolution="high")

# get centroids
centroids <- gCentroid(wmap, byid=TRUE)

# get a data.frame with centroids
centroids.df <- as.data.frame(centroids)
head(centroids.df)
centroids.df$country <- (row.names(centroids.df))
colnames(centroids.df) <- c("lon","lat","country")
head(centroids.df)

# response = ratio
final.dat.imputed$cntry.code
cc.lab <- data.frame(cont.cntry$country,cont.cntry$cntry.code)
colnames(cc.lab) <- c("country","cntry.code")
head(cc.lab)
final.dat.imputed2 <- merge(final.dat.imputed, cc.lab, by="cntry.code")
final.dat.imputed3 <- final.dat.imputed2[!duplicated(final.dat.imputed2$cntry.code), ]
head(final.dat.imputed3)
dim(final.dat.imputed3)

final.dat.imputed4 <- merge(final.dat.imputed3, centroids.df, by="country")
head(final.dat.imputed4)
tail(final.dat.imputed4)

final.dat.imputed4$x <- latlong2grid(final.dat.imputed4[,c(15:16)])$x # equidistant coordinates
final.dat.imputed4$y <- latlong2grid(final.dat.imputed4[,c(15:16)])$y # equidistant coordinates
plot(final.dat.imputed4$x,final.dat.imputed4$y,pch=19)

## determine best correlation structure
m1 <- gls(ratio ~ gdp + igs + cpi + VAag + agrL, data = final.dat.imputed4)
vario1 <- Variogram(m1, form = ~ lon + lat, resType = "pearson")
plot(vario1, smooth = TRUE)

m2 <- gls(ratio ~ gdp + igs + cpi + VAag + agrL, correlation = corExp(form = ~ lon + lat, nugget = T), data = final.dat.imputed4)
m3 <- gls(ratio ~ gdp + igs + cpi + VAag + agrL, correlation = corGaus(form = ~ lon + lat, nugget = T), data = final.dat.imputed4)
m4 <- gls(ratio ~ gdp + igs + cpi + VAag + agrL, correlation = corSpher(form = ~ lon + lat, nugget = T), data = final.dat.imputed4)
m5 <- gls(ratio ~ gdp + igs + cpi + VAag + agrL, correlation = corRatio(form = ~ lon + lat, nugget = T), data = final.dat.imputed4)

mod.lab <- c("noCor","Exp","Gaus","Spher","Ratio")
AIC.vec <- c(AICc(m1), AICc(m2), AICc(m3), AICc(m4), AIC(m5))
dAIC.vec <- delta.IC(AIC.vec)
wAIC.vec <- weight.IC(dAIC.vec)
psR2.mcf <- c(nagelkerke(m1)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m2)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m3)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m4)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m5)$Pseudo.R.squared.for.model.vs.null[1])
psR2.cs <- c(nagelkerke(m1)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m2)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m3)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m4)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m5)$Pseudo.R.squared.for.model.vs.null[2])
psR2.cu <- c(nagelkerke(m1)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m2)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m3)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m4)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m5)$Pseudo.R.squared.for.model.vs.null[3])
results.out <- data.frame(mod.lab,AIC.vec,dAIC.vec,wAIC.vec,psR2.mcf,psR2.cs,psR2.cu)
colnames(results.out) <- c("mod","AICc","dAICc","wAICc","psR2mcf","psR2cs","psR2cu")
results.sort <- results.out[order(results.out[,4],decreasing=T),1:7]
results.sort

# percentage of variance explained by geographic coordinates
100*(1 - results.sort[3,5]/(results.sort[1,5])) # pseudo R2 - McFadden
100*(1 - results.sort[3,6]/(results.sort[1,6])) # pseudo R2 - Cox & Snell
100*(1 - results.sort[3,7]/(results.sort[1,7])) # pseudo R2 - Craig & Uhler

vario3 <- Variogram(m3, form = ~ lon + lat, resType = "pearson")
plot(vario3, smooth = TRUE)
vario3.nr <- Variogram(m3, form = ~ lon + lat, resType = "normalized")
plot(vario3.nr, smooth = TRUE)

# run with Gaussian spatial autocorrelation
# model set
m1 <- "ratio ~ gdp + igs + cpi + VAag + agrL"

m2 <- "ratio ~ gdp + igs + cpi + VAag"
m3 <- "ratio ~ gdp + igs + cpi + agrL"
m4 <- "ratio ~ gdp + igs + VAag + agrL"
m5 <- "ratio ~ gdp + cpi + VAag + agrL"
m6 <- "ratio ~ igs + cpi + VAag + agrL"

m7 <- "ratio ~ gdp + igs + cpi"
m8 <- "ratio ~ gdp + igs + agrL"
m9 <- "ratio ~ gdp + VAag + agrL"
m10 <- "ratio ~ gdp + cpi + VAag"
m11 <- "ratio ~ gdp + igs + VAag"
m12 <- "ratio ~ igs + cpi + VAag"
m13 <- "ratio ~ igs + VAag + agrL"
m14 <- "ratio ~ igs + cpi + agrL"
m15 <- "ratio ~ cpi + VAag + agrL"

m16 <- "ratio ~ gdp + igs"
m17 <- "ratio ~ gdp + cpi"
m18 <- "ratio ~ gdp + VAag"
m19 <- "ratio ~ gdp + agrL"
m20 <- "ratio ~ igs + cpi"
m21 <- "ratio ~ igs + VAag"
m22 <- "ratio ~ igs + agrL"
m23 <- "ratio ~ cpi + VAag"
m24 <- "ratio ~ cpi + agrL"

m25 <- "ratio ~ gdp"
m26 <- "ratio ~ igs"
m27 <- "ratio ~ cpi"
m28 <- "ratio ~ VAag"
m29 <- "ratio ~ agrL"

m30 <- "ratio ~ 1"

## make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30)

## define n.mod
n.mod <- length(mod.vec)

# model fitting and logLik output loop
Modnum <- length(mod.vec)
SaveCount <- BIC.vec <- AICc.vec <- LL.vec <- k.vec <- psR2mcf <- psR2cs <- psR2cu <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- gls(as.formula(mod.vec[i]), correlation = corGaus(form = ~ lon + lat, nugget = T), method="ML", data = final.dat.imputed4)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- logLik(fit)
  k.vec[i] <-  (AIC(fit) - -2*as.numeric(logLik(fit)))/2
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  psR2mcf[i] <- nagelkerke(fit)$Pseudo.R.squared.for.model.vs.null[1] # McFadden pseudo-R2
  psR2cs[i] <- nagelkerke(fit)$Pseudo.R.squared.for.model.vs.null[2] # Cox & Snell pseudo-R2
  psR2cu[i] <- nagelkerke(fit)$Pseudo.R.squared.for.model.vs.null[3] # Craig & Uhler pseudo-R2
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

# re-interpret pseudo-R2s as % deviance from maximum pseudo-R2
psR2mcf.rel <- 100 - 100*(max(psR2mcf) - psR2mcf)/max(psR2mcf)
psR2cs.rel <- 100 - 100*(max(psR2cs) - psR2cs)/max(psR2cs)
psR2cu.rel <- 100 - 100*(max(psR2cu) - psR2cu)/max(psR2cu)

sumtableGLS1 <- data.frame(mod.num,k.vec,round(LL.vec,3),round(AICc.vec,3),round(dAICc,3),round(wAICc,3),BIC.vec,round(dBIC,3),round(wBIC,3),round(psR2mcf,3),round(psR2mcf.rel,1),round(psR2cs,3),round(psR2cs.rel,1),round(psR2cu,3),round(psR2cu.rel,1))
colnames(sumtableGLS1) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","psR2mcf","mcfRel","psR2cs","csRel","psR2cu","cuRel")
row.names(sumtableGLS1) <- as.character(mod.vec)
summary.tableGLS1 <- sumtableGLS1[order(sumtableGLS1[,6],decreasing=T),1:15] # order by wBIC
summary.tableGLS1

top.mod.subGLS1 <- summary.tableGLS1[1,1]
summary(mod.list[[top.mod.subGLS1]])
variotopGLS1 <- Variogram(mod.list[[top.mod.subGLS1]], form = ~ lon + lat, resType = "pearson")
plot(variotopGLS1, smooth = F, pch=19, grid=T)
variotop.nrGLS1 <- Variogram(mod.list[[top.mod.subGLS1]], form = ~ lon + lat, resType = "normalized")
plot(variotop.nrGLS1, smooth = F, pch=19, grid=T)

plot(mod.list[[top.mod.subGLS1]],pch=19)
ggqqplot(mod.list[[top.mod.subGLS1]]$residual)
plot(variotop.nrGLS1, smooth = F, pch=19, grid=T)

plot(mod.list[[1]],pch=19)
ggqqplot(mod.list[[1]]$residual)
variosat.nrGLS1 <- Variogram(mod.list[[1]], form = ~ lon + lat, resType = "normalized")
plot(variosat.nrGLS1, smooth = F, pch=19, grid=T)



# response = pdam
final.datPdam.imputed
final.datPdam.imputed$cntry.code
final.datPdam.imputed2 <- merge(final.datPdam.imputed, cc.lab, by="cntry.code")
final.datPdam.imputed3 <- final.datPdam.imputed2[!duplicated(final.datPdam.imputed2$cntry.code), ]
head(final.datPdam.imputed3)
dim(final.datPdam.imputed3)

final.datPdam.imputed4 <- merge(final.datPdam.imputed3, centroids.df, by="country")
head(final.datPdam.imputed4)
tail(final.datPdam.imputed4)

final.datPdam.imputed4$x <- latlong2grid(final.datPdam.imputed4[,c(15:16)])$x # equidistant coordinates
final.datPdam.imputed4$y <- latlong2grid(final.datPdam.imputed4[,c(15:16)])$y # equidistant coordinates
plot(final.datPdam.imputed4$x,final.datPdam.imputed4$y,pch=19)

## determine best correlation structure
m1 <- gls(pdam ~ gdp + agrL + govexpedu + stja, data = final.datPdam.imputed4)
vario1 <- Variogram(m1, form = ~ lon + lat, resType = "pearson")
plot(vario1, smooth = TRUE)

m2 <- gls(pdam ~ gdp + agrL + govexpedu + stja, correlation = corExp(form = ~ lon + lat, nugget = T), data = final.datPdam.imputed4)
m3 <- gls(pdam ~ gdp + agrL + govexpedu + stja, correlation = corGaus(form = ~ lon + lat, nugget = T), data = final.datPdam.imputed4)
m4 <- gls(pdam ~ gdp + agrL + govexpedu + stja, correlation = corSpher(form = ~ lon + lat, nugget = T), data = final.datPdam.imputed4)
m5 <- gls(pdam ~ gdp + agrL + govexpedu + stja, correlation = corRatio(form = ~ lon + lat, nugget = T), data = final.datPdam.imputed4)

mod.lab <- c("noCor","Exp","Gaus","Spher","Ratio")
AIC.vec <- c(AICc(m1), AICc(m2), AICc(m3), AICc(m4), AIC(m5))
dAIC.vec <- delta.IC(AIC.vec)
wAIC.vec <- weight.IC(dAIC.vec)
psR2.mcf <- c(nagelkerke(m1)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m2)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m3)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m4)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m5)$Pseudo.R.squared.for.model.vs.null[1])
psR2.cs <- c(nagelkerke(m1)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m2)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m3)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m4)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m5)$Pseudo.R.squared.for.model.vs.null[2])
psR2.cu <- c(nagelkerke(m1)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m2)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m3)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m4)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m5)$Pseudo.R.squared.for.model.vs.null[3])
results.out <- data.frame(mod.lab,AIC.vec,dAIC.vec,wAIC.vec,psR2.mcf,psR2.cs,psR2.cu)
colnames(results.out) <- c("mod","AICc","dAICc","wAICc","psR2mcf","psR2cs","psR2cu")
results.sort <- results.out[order(results.out[,4],decreasing=T),1:7]
results.sort

# percentage of variance explained by geographic coordinates
100*(1 - results.sort[2,5]/(results.sort[1,5])) # pseudo R2 - McFadden
100*(1 - results.sort[2,6]/(results.sort[1,6])) # pseudo R2 - Cox & Snell
100*(1 - results.sort[2,7]/(results.sort[1,7])) # pseudo R2 - Craig & Uhler

vario3 <- Variogram(m3, form = ~ lon + lat, resType = "pearson")
plot(vario3, smooth = TRUE)
vario3.nr <- Variogram(m3, form = ~ lon + lat, resType = "normalized")
plot(vario3.nr, smooth = TRUE)

# run with Gaussian spatial autocorrelation
# model set
m1 <- "pdam ~ gdp + agrL + govexpedu + stja"

m2 <- "pdam ~ gdp + agrL + govexpedu"
m3 <- "pdam ~ gdp + agrL + stja"
m4 <- "pdam ~ gdp + govexpedu + stja"
m5 <- "pdam ~ agrL + govexpedu + stja"

m6 <- "pdam ~ gdp + agrL"
m7 <- "pdam ~ gdp + govexpedu"
m8 <- "pdam ~ gdp + stja"
m9 <- "pdam ~ agrL + govexpedu"
m10 <- "pdam ~ agrL + stja"
m11 <- "pdam ~ govexpedu + stja"

m12 <- "pdam ~ gdp"
m13 <- "pdam ~ agrL"
m14 <- "pdam ~ govexpedu"
m15 <- "pdam ~ stja"

m16 <- "pdam ~ 1"

## make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16)

## define n.mod
n.mod <- length(mod.vec)

# model fitting and logLik output loop
Modnum <- length(mod.vec)
SaveCount <- BIC.vec <- AICc.vec <- LL.vec <- k.vec <- psR2mcf <- psR2cs <- psR2cu <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- gls(as.formula(mod.vec[i]), correlation = corExp(form = ~ lon + lat, nugget = T), method="ML", data = final.datPdam.imputed4)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- logLik(fit)
  k.vec[i] <-  (AIC(fit) - -2*as.numeric(logLik(fit)))/2
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  psR2mcf[i] <- nagelkerke(fit)$Pseudo.R.squared.for.model.vs.null[1] # McFadden pseudo-R2
  psR2cs[i] <- nagelkerke(fit)$Pseudo.R.squared.for.model.vs.null[2] # Cox & Snell pseudo-R2
  psR2cu[i] <- nagelkerke(fit)$Pseudo.R.squared.for.model.vs.null[3] # Craig & Uhler pseudo-R2
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

# re-interpret pseudo-R2s as % deviance from maximum pseudo-R2
psR2mcf.rel <- 100 - 100*(max(psR2mcf) - psR2mcf)/max(psR2mcf)
psR2cs.rel <- 100 - 100*(max(psR2cs) - psR2cs)/max(psR2cs)
psR2cu.rel <- 100 - 100*(max(psR2cu) - psR2cu)/max(psR2cu)

sumtableGLS2 <- data.frame(mod.num,k.vec,round(LL.vec,3),round(AICc.vec,3),round(dAICc,3),round(wAICc,3),BIC.vec,round(dBIC,3),round(wBIC,3),round(psR2mcf,3),round(psR2mcf.rel,1),round(psR2cs,3),round(psR2cs.rel,1),round(psR2cu,3),round(psR2cu.rel,1))
colnames(sumtableGLS2) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","psR2mcf","mcfRel","psR2cs","csRel","psR2cu","cuRel")
row.names(sumtableGLS2) <- as.character(mod.vec)
summary.tableGLS2 <- sumtableGLS2[order(sumtableGLS2[,6],decreasing=T),1:15] # order by wBIC
summary.tableGLS2

top.mod.subGLS2 <- summary.tableGLS2[1,1]
summary(mod.list[[top.mod.subGLS2]])
variotopGLS2 <- Variogram(mod.list[[top.mod.subGLS2]], form = ~ lon + lat, resType = "pearson")
plot(variotopGLS2, smooth = F, pch=19, grid=T)
variotop.nrGLS2 <- Variogram(mod.list[[top.mod.subGLS2]], form = ~ lon + lat, resType = "normalized")
plot(variotop.nrGLS2, smooth = F, pch=19, grid=T)

plot(mod.list[[top.mod.subGLS2]],pch=19)
ggqqplot(mod.list[[top.mod.subGLS2]]$residual)
plot(variotop.nrGLS2, smooth = F, pch=19, grid=T)

plot(mod.list[[1]],pch=19)
ggqqplot(mod.list[[1]]$residual)
variosat.nrGLS2 <- Variogram(mod.list[[1]], form = ~ lon + lat, resType = "normalized")
plot(variosat.nrGLS2, smooth = F, pch=19, grid=T)

