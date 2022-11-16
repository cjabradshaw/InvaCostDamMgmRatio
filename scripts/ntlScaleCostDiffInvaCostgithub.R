##########################################################################################################
## national-level prediction & extrapolation of costs based on socio-economic traits
## InvaCost Workshop 2022
## November 2022
## CJA Bradshaw + Emma Hudgins + Phil Hulme + Brian Leung + Franck Courchamp
##########################################################################################################

library(invacost)
library(lme4)
library(dismo)
library(gbm)
library(boot)
library(VIM)
library(mice)

# source files
source("scripts/new_lmer_AIC_tables3.R")
source("scripts/r.squared.R")

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
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
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
  
  cntry.dam.md.boot <- cntry.mgm.md.boot <- cntry.dam.mgm.ratio.md.boot <- cntry.dam.mgm.raw.md.boot <- rep(NA,lcntry)

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
      
    
    cntry.mgm.per.md <- apply(cntry.mgm.per.md.mat, MARGIN=2, median, na.rm=T)
    cntry.mgm.per.lo <- apply(cntry.mgm.per.md.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
    cntry.mgm.per.up <- apply(cntry.mgm.per.md.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
     
    cntry.mgmPgen.per.md <- apply(cntry.mgmPgen.mat, MARGIN=2, median, na.rm=T)
    cntry.mgmPgen.per.lo <- apply(cntry.mgmPgen.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
    cntry.mgmPgen.per.up <- apply(cntry.mgmPgen.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
     
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

    
    cntry.dam <- as.vector(na.omit(c(cntry.dam.per.md,cntry.dam.per.lo,cntry.dam.per.up)))
    if (length(cntry.dam) > 0) {
      cntry.dam.md.boot[c] <- mean(sample(cntry.dam, 10000, replace=T))
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
    
    print("###################")
    print(paste(cntry.vec[c], " (", round(c/lcntry*100, 0), "% complete)", sep=""))
    print("###################")
  } # end c
  
ratio.out <- data.frame(cntry.vec,cntry.dam.mgm.ratio.md.boot,cntry.dam.mgm.raw.md.boot)
colnames(ratio.out) <- c("country","damPgen.mgmPgen","dam.mgm")
hist(log10(ratio.out$dam.mgm),main="",xlab="log10 damage:management")

# dam by mgm relationship (across countries)
dam.mgm.out <- data.frame(cntry.vec, cntry.dam.md.boot, cntry.mgm.md.boot)
colnames(dam.mgm.out) <- c("country","dam","mgm")
plot(log10(dam.mgm.out$mgm), log10(dam.mgm.out$dam), pch=19)
fit.dam.mgm <- lm(log10(dam.mgm.out$dam) ~ log10(dam.mgm.out$mgm))
summary(fit.dam.mgm)
abline(fit.dam.mgm, lty=2, col="red")


# import cntry.code
cont.cntry <- read.csv("data/continent.countryINVACOST.csv", header=T)

ratio.cntry <- merge(ratio.out, cont.cntry, by="country")
head(ratio.cntry)


# GDP per capita
gdp <- read.csv("data/GDPpc.csv")
head(gdp)

ratio.gdp <- merge(ratio.cntry, gdp, by="cntry.code", all = T)
ratio.gdp$GDPpcRecent <- as.numeric(ratio.gdp$GDPpcRecent)
head(ratio.gdp)

plot(log10(ratio.gdp$GDPpcRecent), log10(ratio.gdp$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 GDP per capita")
fit <- lm(log10(ratio.gdp$dam.mgm) ~ log10(ratio.gdp$GDPpcRecent))
summary(fit)
abline(fit,lty=2, col="red")
  
## FAO proportion of GDP devoted to agriculture
## Agriculture, forestry, and fishing, value added (% of GDP)
# Agriculture, forestry, and fishing corresponds to ISIC divisions 1-3 and includes forestry, hunting, and fishing, as well as
# cultivation of crops and livestock production. Value added is the net output of a sector after adding up all outputs and subtracting
# intermediate inputs. It is calculated without making deductions for depreciation of fabricated assets or depletion and degradation
# of natural resources. The origin of value added is determined by the International Standard Industrial Classification (ISIC),
# revision 4. Note: For VAB countries, gross value added at factor cost is used as the denominator.
FAOag <- read.csv("data/faoag.csv")
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

fao.cntry.code <- read.csv("data/fao.cntry.code.csv")
head(fao.cntry.code)
VAag <- merge(VAag1, fao.cntry.code, by="country", all=T)
head(VAag)

## proportion of land protected (% world bank) 
# https://data.worldbank.org/indicator/ER.LND.PTLD.ZS
pa <- read.csv("data/protarea.csv")
pa <- pa[,-1]
head(pa)

# Transparency International
# Brenton-Rule et al. 2016 - corruption (no relationship with # invasive species)
# https://royalsocietypublishing.org/doi/full/10.1098/rspb.2016.0901
# Corruption Perception index https://www.transparency.org/en/cpi/2021
cpi <- read.csv("data/CPI.csv")
head(cpi)

# gov expenditure on education (% GDP)
# https://data.worldbank.org/indicator/SE.XPD.TOTL.GD.ZS
expedu <- read.csv("data/govexpedu.csv")
head(expedu)

ratio.gdpT <- merge(ratio.gdp, tmpr.ann, by="cntry.code", all=T)
head(ratio.gdpT)

ratio.gdpTcpi <- merge(ratio.gdpT, cpi, by="cntry.code", all=T)
head(ratio.gdpTcpi)

ratio.full1 <- merge(ratio.gdpTcpi, expedu, by="cntry.code", all=T)
head(ratio.full1)

ratio.full2 <- merge(ratio.full1, VAag, by="cntry.code", all=T)
head(ratio.full2)

ratio.full <- merge(ratio.full2, pa, by="cntry.code", all=T)
head(ratio.full)


plot((ratio.full$Annual_temp), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="mean annual temperature")
fit <- lm(log10(ratio.full$dam.mgm) ~ (ratio.full$Annual_temp))
abline(fit,lty=2, col="red")

plot((ratio.full$CPI), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="corruption index")
fit <- lm(log10(ratio.full$dam.mgm) ~ (ratio.full$CPI))
abline(fit,lty=2, col="red")

plot((ratio.full$govexpedu), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="gov education expend %GDP")
fit <- lm(log10(ratio.full$dam.mgm) ~ (ratio.full$govexpedu))
abline(fit,lty=2, col="red")

plot((ratio.full$protareapc), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="% protected area")
fit <- lm(log10(ratio.full$dam.mgm) ~ (ratio.full$protareapc))
abline(fit,lty=2, col="red")

plot((ratio.full$VAag), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="value added agriculture/fisheries/forestry %GDP")
fit <- lm(log10(ratio.full$dam.mgm) ~ (ratio.full$protareapc))
abline(fit,lty=2, col="red")


## multiple imputation using mice
# transforms
final.dat.raw <- data.frame(ratio.full$GDPpcRecent, ratio.full$CPI, ratio.full$govexpedu, ratio.full$protareapc, ratio.full$VAag)
colnames(final.dat.raw) <- c("GDPpcRecent", "CPI", "govexpedu", "protareapc", "VAag")
head(final.dat.raw)
tail(final.dat.raw)
dim(final.dat.raw)
dim(na.omit(final.dat.raw))

aggr_plot <- aggr(final.dat.raw, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(final.dat.raw), cex.axis=.7, gap=3, ylab=c("histogram of missing data","pattern"))
marginplot(final.dat.raw[c(1,2)])

md.pattern(final.dat.raw)
final.dat.imp <- mice(final.dat.raw, m=7, maxit=500, method="pmm", seed=101)
summary(final.dat.imp)

final.dat.imp$imp$GDP
final.dat.imp$imp$CPI
final.dat.imp$imp$govexpedu
final.dat.imp$imp$protareapc
final.dat.imp$imp$VAag

final.dat.compl1 <- complete(final.dat.imp, 1)
final.dat.compl2 <- complete(final.dat.imp, 2)
final.dat.compl3 <- complete(final.dat.imp, 3)
final.dat.compl4 <- complete(final.dat.imp, 4)
final.dat.compl5 <- complete(final.dat.imp, 5)

final.dat.compl.mn <- final.dat.raw
sub1 <- which(is.na(final.dat.raw[,1]) == T)
dat1 <- data.frame(final.dat.compl1[sub1, 1], final.dat.compl2[sub1, 1], final.dat.compl3[sub1, 1], final.dat.compl4[sub1, 1], final.dat.compl5[sub1, 1])
final.dat.compl.mn[sub1, 1] <- apply(dat1, MARGIN=1, mean)

sub2 <- which(is.na(final.dat.raw[,2]) == T)
dat2 <- data.frame(final.dat.compl1[sub2, 2], final.dat.compl2[sub2, 2], final.dat.compl3[sub2, 2], final.dat.compl4[sub2, 2], final.dat.compl5[sub2, 2])
final.dat.compl.mn[sub2, 2] <- apply(dat2, MARGIN=1, mean)

sub3 <- which(is.na(final.dat.raw[,3]) == T)
dat3 <- data.frame(final.dat.compl1[sub3, 3], final.dat.compl2[sub3, 3], final.dat.compl3[sub3, 3], final.dat.compl4[sub3, 3], final.dat.compl5[sub3, 3])
final.dat.compl.mn[sub3, 3] <- apply(dat3, MARGIN=1, mean)

sub4 <- which(is.na(final.dat.raw[,4]) == T)
dat4 <- data.frame(final.dat.compl1[sub4, 4], final.dat.compl2[sub4, 4], final.dat.compl3[sub4, 4], final.dat.compl4[sub4, 4], final.dat.compl5[sub4, 4])
final.dat.compl.mn[sub4, 4] <- apply(dat4, MARGIN=1, mean)

sub5 <- which(is.na(final.dat.raw[,5]) == T)
dat5 <- data.frame(final.dat.compl1[sub5, 5], final.dat.compl2[sub5, 5], final.dat.compl3[sub5, 5], final.dat.compl4[sub5, 5], final.dat.compl5[sub5, 5])
final.dat.compl.mn[sub5, 5] <- apply(dat5, MARGIN=1, mean)

GDPsc <- log10(scale((final.dat.compl.mn$GDPpcRecent), scale=T, center=F))
CPIsc <- scale(final.dat.compl.mn$CPI, scale=T, center=F)
govexpedusc <- scale(logit(final.dat.compl.mn$govexpedu/100), scale=T, center=F)
pasc <- scale(logit(final.dat.compl.mn$protareapc/100), scale=T, center=F)
VAagsc <- scale(logit(final.dat.compl.mn$VAag/100), scale=T, center=F)
ratiosc <- log10(scale(ratio.full$dam.mgm, center=F, scale=T))

head(ratio.full)
dim(ratio.full)

final.dat.imputed <- na.omit(data.frame(ratio.full$cntry.code, ratio.full$cont, ratio.full$region, ratiosc,
                                GDPsc, CPIsc, govexpedusc, pasc, VAagsc))
colnames(final.dat.imputed) <- c("cntry.code", "cont", "reg", "ratio", "GDP","CPI","govExpEdu","PA","VAag")
head(final.dat.imputed)
dim(final.dat.imputed)


brt.fit <- gbm.step(final.dat.imputed, gbm.x = attr(final.dat.imputed, "names")[c(5:9)], gbm.y = attr(final.dat.imputed, "names")[4], family="gaussian", max.trees=100000, tolerance = 0.00001, learning.rate = 0.00003, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit)
D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit)
gbm.plot.fits(brt.fit)

brt.CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
brt.CV.cor
brt.CV.cor.se <- 100 * brt.fit$cv.statistics$correlation.se
brt.CV.cor.se
print(c(brt.CV.cor, brt.CV.cor.se))


# GLMM
# recode continent to greater region to increase per-category sample sizes
table(final.dat.imputed$cont)
final.dat.imputed$reg2 <- ifelse(final.dat.imputed$cont == "NAM" | final.dat.imputed$cont == "CAR", "NAMCAR", final.dat.imputed$cont)
final.dat.imputed$reg2 <- ifelse(final.dat.imputed$cont == "EUR" | final.dat.imputed$cont == "ME", "EURME", final.dat.imputed$reg2)
final.dat.imputed$reg2 <- ifelse(final.dat.imputed$cont == "ASIA" | final.dat.imputed$cont == "OC", "ASIAOC", final.dat.imputed$reg2)
table(final.dat.imputed$reg2)


# model set
m1 <- "ratio ~ GDP + PA + CPI + govExpEdu + VAag + (1|reg2)"

m2 <- "ratio ~ GDP + PA + CPI + govExpEdu + (1|reg2)"
m3 <- "ratio ~ GDP + PA + CPI + VAag + (1|reg2)"
m4 <- "ratio ~ GDP + PA + govExpEdu + VAag + (1|reg2)"
m5 <- "ratio ~ GDP + CPI + govExpEdu + VAag + (1|reg2)"
m6 <- "ratio ~ PA + CPI + govExpEdu + VAag + (1|reg2)"

m7 <- "ratio ~ GDP + PA + CPI + (1|reg2)"
m8 <- "ratio ~ GDP + PA + govExpEdu + (1|reg2)"
m9 <- "ratio ~ GDP + PA + VAag + (1|reg2)"
m10 <- "ratio ~ GDP + CPI + govExpEdu + (1|reg2)"
m11 <- "ratio ~ GDP + CPI + VAag + (1|reg2)"
m12 <- "ratio ~ GDP + govExpEdu + VAag + (1|reg2)"
m13 <- "ratio ~ PA + CPI + govExpEdu + (1|reg2)"
m14 <- "ratio ~ PA + CPI + VAag + (1|reg2)"
m15 <- "ratio ~ CPI + govExpEdu + VAag + (1|reg2)"

m16 <- "ratio ~ GDP + CPI + (1|reg2)"
m17 <- "ratio ~ PA + CPI + (1|reg2)"
m18 <- "ratio ~ GDP + PA + (1|reg2)"
m19 <- "ratio ~ GDP + govExpEdu + (1|reg2)"
m20 <- "ratio ~ PA + govExpEdu + (1|reg2)"
m21 <- "ratio ~ CPI + govExpEdu + (1|reg2)"
m22 <- "ratio ~ GDP + VAag + (1|reg2)"
m23 <- "ratio ~ VAag + CPI + (1|reg2)"
m24 <- "ratio ~ VAag + PA + (1|reg2)"
m25 <- "ratio ~ VAag + govExpEdu + (1|reg2)"
m26 <- "ratio ~ CPI + VAag + (1|reg2)"

m27 <- "ratio ~ GDP + (1|reg2)"
m28 <- "ratio ~ PA + (1|reg2)"
m29 <- "ratio ~ CPI + (1|reg2)"
m30 <- "ratio ~ govExpEdu + (1|reg2)"
m31 <- "ratio ~ VAag + (1|reg2)"

m32 <- "ratio ~ 1 + (1|reg2)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32)
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

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable) <- as.character(mod.vec)
summary.table <- sumtable[order(sumtable[,5],decreasing=F),]
summary.table


# resampled BRT loop
biter <- 100
eq.sp.points <- 100

# create storage arrays
val.arr <- pred.arr <- array(data = 0, dim = c(eq.sp.points, 5, iter), dimnames=list(paste("x",1:eq.sp.points,sep=""), attr(final.dat.imputed, "names")[c(5:9)], paste("b",1:iter,sep="")))

# create storage vectors
D2.vec <- CV.cor.vec <- CV.cor.se.vec <- GDP.ri <- govExpEdu.ri <- CPI.ri <- PA.ri <- VAag.ri <- rep(0,iter)

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
pred.update <- pred.arr

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:iter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] > (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  }
  
  print(k)
}

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)

val.med <- apply(val.arr, MARGIN=c(1,2), median)

par(mfrow=c(3,2)) 
plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,1]),max(pred.up[,1])), lwd=2, ylab="(←higher) dam:mgm (lower→)", xlab="(←poorer) GDPpc (richer→)" )
lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2, ylab="(←higher) dam:mgm (lower→)", xlab="(←lower) education investment (higher→)" )
lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")

plot(val.med[,3],pred.med[,3],type="l",ylim=c(min(pred.lo[,3]),max(pred.up[,3])), lwd=2, ylab="(←higher) dam:mgm (lower→)", xlab="(←higher) corruption (lower→" )
lines(val.med[,3], pred.lo[,3], type="l", lty=2, col="red")
lines(val.med[,3], pred.up[,3], type="l", lty=2, col="red")

plot(val.med[,4],pred.med[,4],type="l",ylim=c(min(pred.lo[,4]),max(pred.up[,4])), lwd=2, ylab="(←higher) dam:mgm (lower→)", xlab="(←less) protected area (more→)" )
lines(val.med[,4], pred.lo[,4], type="l", lty=2, col="red")
lines(val.med[,4], pred.up[,4], type="l", lty=2, col="red")

plot(val.med[,5],pred.med[,5],type="l",ylim=c(min(pred.lo[,5]),max(pred.up[,5])), lwd=2, ylab="(←higher) dam:mgm (lower→)", xlab="(←less) agriculture (more→)" )
lines(val.med[,5], pred.lo[,5], type="l", lty=2, col="red")
lines(val.med[,5], pred.up[,5], type="l", lty=2, col="red")
par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update <- D2.vec
CV.cor.update <- CV.cor.vec
CV.cor.se.update <- CV.cor.se.vec
GDP.ri.update <- GDP.ri
govExpEdu.ri.update <- govExpEdu.ri
CPI.ri.update <- CPI.ri
PA.ri.update <- PA.ri
VAag.ri.update <- VAag.ri

for (k in 1:kappa.n) {
  D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  GDP.mean <- mean(GDP.ri.update, na.rm=T); GDP.sd <- sd(GDP.ri.update, na.rm=T)
  govExpEdu.mean <- mean(govExpEdu.ri.update, na.rm=T); govExpEdu.sd <- sd(govExpEdu.ri.update, na.rm=T)
  CPI.mean <- mean(CPI.ri.update, na.rm=T); CPI.sd <- sd(CPI.ri.update, na.rm=T)
  PA.mean <- mean(PA.ri.update, na.rm=T); PA.sd <- sd(PA.ri.update, na.rm=T)
  VAag.mean <- mean(VAag.ri.update, na.rm=T); VAag.sd <- sd(VAag.ri.update, na.rm=T)

  for (u in 1:iter) {
    D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    
    GDP.ri.update[u] <- ifelse((GDP.ri.update[u] < (GDP.mean-kappa*GDP.sd) | GDP.ri.update[u] > (GDP.mean+kappa*GDP.sd)), NA, GDP.ri.update[u])
    govExpEdu.ri.update[u] <- ifelse((govExpEdu.ri.update[u] < (govExpEdu.mean-kappa*govExpEdu.sd) | govExpEdu.ri.update[u] > (govExpEdu.mean+kappa*govExpEdu.sd)), NA, govExpEdu.ri.update[u])
    CPI.ri.update[u] <- ifelse((CPI.ri.update[u] < (CPI.mean-kappa*CIP.sd) | CPI.ri.update[u] > (CPI.mean+kappa*CPI.sd)), NA, CPI.ri.update[u])
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

write.table(pred.med,file="BRT.boot.pred.med.CHLTH.HS.csv",sep=",", row.names = T, col.names = T)
write.table(pred.lo,file="BRT.boot.pred.lo.CHLTH.HS.csv",sep=",", row.names = T, col.names = T)
write.table(pred.up,file="BRT.boot.pred.up.CHLTH.HS.csv",sep=",", row.names = T, col.names = T)
write.table(val.med,file="BRT.boot.val.med.CHLTH.HS.csv",sep=",", row.names = T, col.names = T)
write.table(ri.sort,file="BRT.boot.ri.CHLTH.HS.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.val.dat,file="BRT.val.CHLTH.HS.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.pred.dat,file="BRT.pred.CHLTH.HS.csv",sep=",", row.names = T, col.names = T)

save.image("ratioOut.RData")
