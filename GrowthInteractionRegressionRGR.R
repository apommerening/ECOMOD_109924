# Functions for shot-noise regression. Based on Haebel et al. (2019). Updated on 04.08.2021.

rm(list = ls())
library(Rcpp)
options(digits = 14, width = 100)

dataFile <- "/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/MostCompetitiveNeighbour/Data/"
codeFile <- "/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/MostCompetitiveNeighbour/Code/"
outputFile <- "/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/MostCompetitiveNeighbour/Results/"


calcBasalArea <- function(dbh) 
  return(pi * (dbh / 200)^2)

#-- Read in data for regressions --#
#----------------------------------#
timeSeries <- "giRegresEucalypt4000.txt"
TreeList  <- read.table(paste(dataFile, timeSeries, sep = ""), header = T)
names(TreeList)
str(TreeList)
range(TreeList$dbh)
range(TreeList$year)
TreeList[200 : 300, c(2 : 3, 6 : 7, 11 : 14)]
table(TreeList$Species)

# Checking whether dbh RGR is really half of basal-area RGR
par(mar = c(2, 3.2, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeList$RGR, TreeList$RGR.g / 2, ylab = "", xlab = "", 
     cex = .9, col = "black", pch = 16, axes = FALSE, xlim = c(0, 0.5), ylim = c(0, 0.5))
abline(0, 1, col = "red")
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

par(mar = c(2, 5.0, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeList$dbh, TreeList$RGR - TreeList$RGR.g / 2, cex = .9, ylab = "", xlab = "", col = "black", pch = 16, axes = FALSE)
abline(h = 0, col = "red")
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)


# Austria (Litschau)
if(substr(timeSeries, 9, 12) == "FLNR") {
  table(TreeList$Species)
  TreeList$Species[TreeList$Species == 7] <- 4
  spec <- table(TreeList$Species)
  spec <- sort(spec, decreasing = T)
  speciesCodes <- as.numeric(names(spec))
} else {
  spec <- table(TreeList$Species)
  spec <- sort(spec, decreasing = T)
  if((substr(timeSeries, 9, 10) == "DF") | (substr(timeSeries, 9, 10) == "NS")) 
    TreeList$Species <- rep(1, length(TreeList$Species)) 
  else 
    TreeList$Species <- rep(as.numeric(names(spec)[1]), length(TreeList$Species)) 
  spec <- table(TreeList$Species)
  spec <- sort(spec, decreasing = T)
  speciesCodes <- as.numeric(names(spec))
}
spec
speciesCodes

#-- Regression mark potential --#
#-------------------------------#
TreeLista <- TreeList[TreeList$AGR > 0, ]
TreeLista <- TreeLista[!is.na(TreeLista$AGR) & !is.na(TreeLista$dbh),]
TreeLista <- TreeLista[!(TreeLista$AGR == Inf),]
TreeLista <- TreeLista[!(TreeLista$dbh == Inf),]
range(TreeLista$AGR)
range(TreeLista$dbh)
table(TreeList$Species)
TreeList$Species <- as.integer(TreeList$Species)
typeof(TreeList$Species)
table(TreeList$year)
# TreeLista <- TreeLista[TreeLista$year > 2002,] # For Hirschlacke only


# Quantile regression for potential AGR
#-------------------------------------------------------#
# install.packages("quantreg", dep = T)
library(quantreg)

# AGR
dpot <- function(dbh, xA, xk, xp) {				
  return(xA * dbh^xk * exp(-xp *  dbh)) # Modified from Zeide (1993, p. 609, Eq. 7))
}

# Mono-specific case with quantreg (Zeide)
nlsout <- nlrq(AGR.g ~ dpot(dbh, A, k, p), data = TreeLista, start = list(A = 0.18, k = 1.54, p = 0.01), # NS, BE, DF, FLNR, HK, PP, SP (tau = 0.975)
              tau = 0.975, trace = TRUE)

summary(nlsout)
param <- c(summary(nlsout)$coefficients[1], summary(nlsout)$coefficients[2], summary(nlsout)$coefficients[3])
rm(nlsout)

maxG <- ceiling(max(TreeList$AGR.g, na.rm = T) * 10000) / 10000
# pdf(file = paste(outputFile, substr(timeSeries, 1, (nchar(timeSeries) - 4)), "_AGR_Potential.pdf", sep = ""))
par(mar = c(2, 5, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeLista$dbh, TreeLista$AGR.g, las = 1, ylab = "", xlab = "", 
     cex = .9, col = "black", pch = 16, axes = FALSE, xlim = c(0, ceiling(max(TreeList$dbh, na.rm = T))), ylim = c(0, 0.030))#, ylim = c(0, 1)
# points(x, y, pch = 16, col = "red")
curve(dpot(x, param[1], param[2], param[3]), from = min(TreeLista$dbh), to = max(TreeLista$dbh), lwd = 4, lty = 1, col = "red", add = TRUE)
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)
# dev.off()
rm(maxG)

dbhList <- seq(min(TreeLista$dbh), max(TreeLista$dbh), by = 0.5)
ik <- dpot(dbhList, param[1], param[2], param[3])
# pkplus <- ik / dbhList # When using diameter RGR
pkplus <- ik / calcBasalArea(dbhList) # When using basal-area RGR
pk <- log(1 + pkplus)

# pdf(file = paste(outputFile, substr(timeSeries, 1, (nchar(timeSeries) - 4)), "_RGR_Potential.pdf", sep = ""))
par(mar = c(2, 4.2, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeLista$dbh, TreeLista$RGR.g, las = 1, ylab = "", xlab = "", 
     cex = .9, col = "black", pch = 16, axes = FALSE, xlim = c(0, ceiling(max(TreeList$dbh, na.rm = T))), ylim = c(0, 0.25)) # , ylim = c(0, 1)
lines(dbhList, pk, lwd = 4, col = "red")
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)
# dev.off()


#-- Growth interaction regressions --#
#------------------------------------#

# Loss function for estimating interaction parameters abdn
#
# @param abdn A parameter vector containing interaction parameters (standardization (=first parameter) + kernel)
# @param xTreeList A suitable data frame?
# @param xparam The potential mark growth parameters k, p, q

loss.L2 <- function(abdn, xTreeList, xparam, specCodes, fact, denom, varr, relascope) {
  gi <- estimateRelativeGrowthRate(xData = xTreeList,  param = xparam, abdn = abdn, species = xTreeList$Species, codes = specCodes, 
                                   rep(pi, length(xTreeList$dbh)), relascope)
  dev <- 1 / varr * (xTreeList$RGR.g - gi)^2
  # is.na(dev) <- sapply(dev, is.infinite)
  if(denom == 0)
     ret <- sum(dev, na.rm = TRUE) + fact * any(abdn[1 : 2] < 0) else ret <- sum(dev, na.rm = T) + fact * any(abdn[1 : 2] < 0) + sum(abs(abdn[1 : 2]) / denom)
  # ret <- sum(dev, na.rm = T)
  return(ret)
}

loss.ML <- function(abdn, xTreeList, xparam, specCodes, fact, denom, relascope) {
  gi <- estimateRelativeGrowthRate(xData = xTreeList,  param = xparam, abdn = abdn[1 : 2], species = xTreeList$Species, 
                                   codes = specCodes, rep(pi, length(xTreeList$dbh)), relascope)
  dev <- dnorm(xTreeList$RGR.g, mean = gi, sd = exp(abdn[length(abdn)]), log = TRUE)
  # is.na(dev) <- sapply(dev, is.infinite)
  if(denom == 0)
      ret <- sum(dev, na.rm = TRUE) else ret <- sum(dev, na.rm = TRUE) - 2 * sum(abdn[1 : 2] / denom) - fact * any(abdn[1 : 2] < 0)
  # ret <- sum(dev, na.rm = T)
  return(ret)
}

#-- Least squares
#----------------

TreeList <- TreeList[TreeList$dbh > 0, ] 
TreeList$RGR[is.infinite(TreeList$AGR)] <- NA
TreeList$RGR[TreeList$AGR < 0] <- NA
# TreeList <- TreeList[TreeList$year > 2002,] # Hirschlacke only

sourceCpp(paste(codeFile, "GrowthInteraction.cpp", sep = ""))
#-----------------------------#
# SNx: nu, alpha, beta, delta #
#-----------------------------#
# Starting weights (Run before first LS regression)
w <- rep(1, length(TreeList$dbh))
# Weights for weighted regression (Run before second LS regression and modify number in next line if necessary)
group <- as.numeric(cut(TreeList$dbh, 50))
y <- as.numeric(tapply(abs(error.hats), group, var, na.rm = T))
x <- as.numeric(tapply(TreeList$dbh, group, median, na.rm = T))
w <- group
w <- y[w]
  
par(mar = c(2, 3.3, 0.5, 0.5), mfrow = c(1, 1))
plot(x, y, las = 1, ylab = "", xlab = "", 
     cex = .9, col = "black", pch = 16, axes = FALSE)
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

relas <- F
abdn0 <- c(0.0547824493386750, 0.0605973818356007)
if(!relas)
  abdn0 <- c(0.0547824493386750, 0.0605973818356007, 0.074714727055536, 0.381539387191665)
  # abdn0 <- c(0.20, 1.979, 21.929, 3.357) # SP (Haebel et al., 2019)

#-- LS
#-----
# Optimization
myFact <- 500 
myDenom <- 50 # 50

abdn.L2 <- optim(abdn0, loss.L2, xTreeList = TreeList, xparam = as.matrix(param), specCodes = speciesCodes, fact = myFact, 
                 denom = myDenom, varr = w, relascope = relas, control = list(maxit = 900000, reltol = .Machine$double.eps, 
                 temp = 80000, trace = 6, REPORT = 5000)) #, method = "L-BFGS-B", lower = rep(0, 4), upper = rep(Inf, 4))
# Parameter values
abdn.L2$par

ci <- estTreeInteraction(abdn.L2$par, TreeList, relas)
hist(ci)
range(ci, na.rm = TRUE)

yhat <- estimateRelativeGrowthRate(xData = TreeList,  param = as.matrix(param), abdn = abdn.L2$par, 
                 species = TreeList$Species, codes = speciesCodes, pi = rep(pi, length(TreeList$dbh)), choice = relas) / 2
(varres <- var(TreeList$RGR - yhat, na.rm = TRUE)) 
(bias <- mean(TreeList$RGR - yhat, na.rm = TRUE)) # Should be: mean(yhat - TreeList$RGR, na.rm = TRUE)
bias / mean(TreeList$RGR, na.rm = TRUE)
(rmse <- sqrt(varres + bias^2))
rmse / mean(TreeList$RGR, na.rm = TRUE)
error.hats <- TreeList$RGR - yhat 
(Eff <- 1 - (sum((TreeList$RGR - yhat)^2, na.rm = T) / sum((TreeList$RGR - mean(TreeList$RGR, na.rm = T))^2, na.rm = T))) # Efficiency

# AIC
length(TreeList$RGR) * log(sum((TreeList$RGR - yhat)^2, na.rm = T) / length(TreeList$RGR)) + 2 * length(abdn.L2$par) 

par(mar = c(2, 4.0, 0.5, 0.5), mfrow = c(1, 1))
plot(yhat, TreeList$RGR, ylab = "", xlab = "", cex = .9, col = "black", pch = 16, axes = FALSE, ylim = c(0, 0.3), 
     xlim = c(0, 0.3))
axis(1, cex.axis = 1.8)
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
abline(0, 1)
box(lwd = 2)

par(mar = c(2, 4.0, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeList$dbh, error.hats, ylab = "", xlab = "", cex = .9, col = "black", pch = 16, axes = FALSE) #, ylim = c(0, 0.3), xlim = c(0, 0.3))
axis(1, cex.axis = 1.8)
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
abline(h = 0, col = "red")
box(lwd = 2)

qqnorm(error.hats, las = 1, cex.axis = 0.9)
qqline(error.hats)
library(car)
qqPlot(error.hats)

# pdf(paste(outputFile, "Signal.pdf", sep = ""))
maxS <- ceiling((60 / (1 + 2 * abs(0)))^abdn.L2$par[1])
my.labels <- seq(-15, 15, 5)
par(lab = c(length(my.labels), 5, 7), mar = c(2, 4, 0.5, 0.5), mfrow = c(1, 1))
curve((60 / (1 + 2 * abs(x)))^abdn.L2$par[1], from = -max(my.labels), to = max(my.labels), xlab = "", ylab = "", col = "black", xaxt = "n", axes = FALSE, lwd = 2, xlim = c(-max(my.labels), max(my.labels)), ylim = c(0, maxS))
curve((20 / (1 + 2 * abs(x)))^abdn.L2$par[1], from = -max(my.labels), to = max(my.labels), col = "red", lwd = 2, add = T) 
curve((10 / (1 + 2 * abs(x)))^abdn.L2$par[1], from = -max(my.labels), to = max(my.labels), col = "blue", lwd = 2, add = T) 
axis(side = 1, at = my.labels, labels = abs(my.labels), cex.axis = 1.8)
axis(side = 2, las = 1, cex.axis = 1.8)
box(lwd = 2)
# dev.off()

# (pow(mark[j], abd[0]) * exp(-abd[2] * pow(distance, 1) / pow(mark[j], abd[1])))
# pdf(paste(outputFile, "Signal.pdf", sep = ""))
maxS <- ceiling(60^abdn.L2$par[1] * exp(-abdn.L2$par[3] * abs(0) / 60^abdn.L2$par[2]))
my.labels <- seq(-15, 15, 5)
par(lab = c(length(my.labels), 5, 7), mar = c(2, 3.5, 0.5, 0.5), mfrow = c(1, 1))
curve(60^abdn.L2$par[1] * exp(-abdn.L2$par[3] * abs(x) / 60^abdn.L2$par[2]), from = -max(my.labels), to = max(my.labels), xlab = "", ylab = "", col = "black", xaxt = "n", axes = FALSE, lwd = 2, xlim = c(-max(my.labels), max(my.labels)), ylim = c(0, maxS))
curve(20^abdn.L2$par[1] * exp(-abdn.L2$par[3] * abs(x) / 20^abdn.L2$par[2]), from = -max(my.labels), to = max(my.labels), col = "red", lwd = 2, add = T) 
curve(10^abdn.L2$par[1] * exp(-abdn.L2$par[3] * abs(x) / 10^abdn.L2$par[2]), from = -max(my.labels), to = max(my.labels), col = "blue", lwd = 2, add = T) 
axis(side = 1, at = my.labels, labels = abs(my.labels), cex.axis = 1.8)
axis(side = 2, las = 1, cex.axis = 1.8)
box(lwd = 2)
# dev.off()



#-- Maximum Likelihood
#---------------------
abdn0 <- c(abdn0, 5)
# Optimization
abdn.ML <- optim(abdn0, loss.ML, xTreeList = TreeList, xparam = as.matrix(param), specCodes = speciesCodes, fact = myFact,
                 denom = myDenom, relascope = relas, hessian = TRUE, control = list(maxit = 900000, reltol = .Machine$double.eps, 
                 temp = 80000, trace = 6, fnscale = -1, REPORT = 5000))

# Comparison of LS and ML
abdn.ML$par
abdn.L2$par
# Also, if you add the hessian = TRUE then you can compute this to get asymptotic estimates of the standard errors of the parameter estimates.
sqrt(diag(solve(-abdn.ML$hessian))) 
# Evaluate the L2 loss function at the ML location (ignoring variance)
loss.L2(abdn.ML$par[-length(abdn.ML$par)], xTreeList = TreeList, xparam = as.matrix(param), specCodes = speciesCodes, fact = myFact, denom = myDenom, varr = w, relascope = relas) 
abdn.L2$value

# AIC
-2 * abdn.ML$value + 2 * length(abdn.ML$par) 

ci <- estTreeInteraction(abdn.ML$par, TreeList, relas)
hist(ci)
range(ci, na.rm = TRUE)

yhat <- estimateRelativeGrowthRate(xData = TreeList,  param = as.matrix(param), abdn = abdn.ML$par, 
                species = TreeList$Species, codes = speciesCodes, pi = rep(pi, length(TreeList$dbh)), choice = relas) / 2
(varres <- var(TreeList$RGR - yhat, na.rm = TRUE)) 
(bias <- mean(TreeList$RGR - yhat, na.rm = TRUE)) # Should be: mean(yhat - TreeList$RGR, na.rm = TRUE)
bias / mean(TreeList$RGR, na.rm = TRUE)
(rmse <- sqrt(varres + bias^2))
rmse / mean(TreeList$RGR, na.rm = TRUE)
error.hats <- TreeList$RGR - yhat 
(Eff <- 1 - (sum((TreeList$RGR - yhat)^2, na.rm = T) / sum((TreeList$RGR - mean(TreeList$RGR, na.rm = T))^2, na.rm = T))) # Efficiency

par(mar = c(2, 4.0, 0.5, 0.5), mfrow = c(1, 1))
plot(yhat, TreeList$RGR, ylab = "", xlab = "", cex = .9, col = "black", pch = 16, axes = FALSE, ylim = c(0, 0.3),
     xlim = c(0, 0.3))
axis(1, cex.axis = 1.8)
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
abline(0, 1)
box(lwd = 2)

par(mar = c(2, 4.0, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeList$dbh, error.hats, ylab = "", xlab = "", cex = .9, col = "black", pch = 16, axes = FALSE) #, ylim = c(0, 0.3), xlim = c(0, 0.3))
axis(1, cex.axis = 1.8)
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
abline(h = 0, col = "red")
box(lwd = 2)

# pdf(paste(dataFile, "Paper/Signal.pdf", sep = ""))
maxS <- ceiling((60 / (1 + 2 * abs(0)))^abdn.ML$par[1])
my.labels <- seq(-15, 15, 5)
par(lab = c(length(my.labels), 5, 7), mar = c(2, 4, 0.5, 0.5), mfrow = c(1, 1))
curve((60 / (1 + 2 * abs(x)))^abdn.ML$par[1], from = -max(my.labels), to = max(my.labels), xlab = "", ylab = "", col = "black", xaxt = "n", axes = FALSE, lwd = 2, ylim = c(0, maxS), xlim = c(-max(my.labels), max(my.labels))) # ylim = c(0, maxS),
curve((20 / (1 + 2 * abs(x)))^abdn.ML$par[1], from = -max(my.labels), to = max(my.labels), col = "red", lwd = 2, add = T) 
curve((10 / (1 + 2 * abs(x)))^abdn.ML$par[1], from = -max(my.labels), to = max(my.labels), col = "blue", lwd = 2, add = T) 
axis(side = 1, at = my.labels, labels = abs(my.labels), cex.axis = 1.8)
axis(side = 2, las = 1, cex.axis = 1.8)
box(lwd = 2)
# dev.off()

# (pow(mark[j], abd[0]) * exp(-abd[2] * pow(distance, 1) / pow(mark[j], abd[1])))
# pdf(paste(outputFile, "Signal.pdf", sep = ""))
maxS <- ceiling(60^abdn.ML$par[1] * exp(-abdn.ML$par[3] * abs(0) / 60^abdn.ML$par[2]))
my.labels <- seq(-15, 15, 5)
par(lab = c(length(my.labels), 5, 7), mar = c(2, 3.5, 0.5, 0.5), mfrow = c(1, 1))
curve(60^abdn.ML$par[1] * exp(-abdn.ML$par[3] * abs(x) / 60^abdn.ML$par[2]), from = -max(my.labels), to = max(my.labels), xlab = "", ylab = "", col = "black", xaxt = "n", axes = FALSE, lwd = 2, xlim = c(-max(my.labels), max(my.labels)), ylim = c(0, maxS))
curve(20^abdn.ML$par[1] * exp(-abdn.ML$par[3] * abs(x) / 20^abdn.ML$par[2]), from = -max(my.labels), to = max(my.labels), col = "red", lwd = 2, add = T) 
curve(10^abdn.ML$par[1] * exp(-abdn.ML$par[3] * abs(x) / 10^abdn.ML$par[2]), from = -max(my.labels), to = max(my.labels), col = "blue", lwd = 2, add = T) 
axis(side = 1, at = my.labels, labels = abs(my.labels), cex.axis = 1.8)
axis(side = 2, las = 1, cex.axis = 1.8)
box(lwd = 2)
# dev.off()


#-- Maximum likelihood with LS-starting values
#--------------------------------------------
# Fit the ML function using the LS parameter estimates as starting points and an arbitrary variance
abdn.ML2 <- optim(c(abdn.L2$par, 5), loss.ML, xTreeList = TreeList, xparam = as.matrix(param), specCodes = speciesCodes, 
                  fact = myFact, denom = myDenom, relascope = relas,
                  control = list(maxit = 900000, reltol = .Machine$double.eps, temp = 80000, trace = 6, fnscale = -1, 
                  REPORT = 5000))

# Compare the three parameter estimates
abdn.ML2$par
abdn.L2$par
abdn.ML$par

# Compare the L2 loss function evaluated at each set
loss.L2(abdn.ML2$par[-length(abdn.ML2$par)], xTreeList = TreeList, xparam = param, specCodes = speciesCodes, fact = myFact, denom = myDenom, varr = w, relascope = relas)
loss.L2(abdn.ML$par[-length(abdn.ML$par)], xTreeList = TreeList, xparam = param, specCodes = speciesCodes, fact = myFact, denom = myDenom, varr = w, relascope = relas)
abdn.L2$value

# AIC
-2 * abdn.ML2$value + 2 * length(abdn.ML2$par) 

ci <- estTreeInteraction(abdn.ML2$par, TreeList, relas)
hist(ci)
range(ci, na.rm = TRUE)

yhat <- estimateRelativeGrowthRate(xData = TreeList,  param = as.matrix(param), abdn = abdn.ML2$par, 
                species = TreeList$Species, codes = speciesCodes, pi = rep(pi, length(TreeList$dbh)), choice = relas) / 2
(varres <- var(TreeList$RGR - yhat, na.rm = TRUE)) 
(bias <- mean(TreeList$RGR - yhat, na.rm = TRUE)) # Should be: mean(yhat - TreeList$RGR, na.rm = TRUE)
bias / mean(TreeList$RGR, na.rm = TRUE)
(rmse <- sqrt(varres + bias^2))
rmse / mean(TreeList$RGR, na.rm = TRUE)
error.hats <- TreeList$RGR - yhat 
(Eff <- 1 - (sum((TreeList$RGR - yhat)^2, na.rm = T) / sum((TreeList$RGR - mean(TreeList$RGR, na.rm = T))^2, na.rm = T))) # Efficiency

par(mar = c(2, 4.0, 0.5, 0.5), mfrow = c(1, 1))
plot(yhat, TreeList$RGR, ylab = "", xlab = "", cex = .9, col = "black", pch = 16, axes = FALSE, ylim = c(0, 0.3), 
     xlim = c(0, 0.3))
axis(1, cex.axis = 1.8)
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
abline(0, 1)
box(lwd = 2)

par(mar = c(2, 4.0, 0.5, 0.5), mfrow = c(1, 1))
plot(TreeList$dbh, error.hats, ylab = "", xlab = "", cex = .9, col = "black", pch = 16, axes = FALSE) #, ylim = c(0, 0.3), xlim = c(0, 0.3))
axis(1, cex.axis = 1.8)
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
abline(h = 0, col = "red")
box(lwd = 2)

# pdf(paste(dataFile, "Paper/Signal.pdf", sep = ""))
maxS <- ceiling((60 / (1 + 2 * abs(0)))^abdn.ML2$par[1])
my.labels <- seq(-15, 15, 5)
par(lab = c(length(my.labels), 5, 7), mar = c(2, 3.5, 0.5, 0.5), mfrow = c(1, 1))
curve((60 / (1 + 2 * abs(x)))^abdn.ML2$par[1], from = -max(my.labels), to = max(my.labels), xlab = "", ylab = "", col = "black", xaxt = "n", axes = FALSE, lwd = 2, ylim = c(0, maxS), xlim = c(-max(my.labels), max(my.labels))) # ylim = c(0, maxS),
curve((20 / (1 + 2 * abs(x)))^abdn.ML2$par[1], from = -max(my.labels), to = max(my.labels), col = "red", lwd = 2, add = T) 
curve((10 / (1 + 2 * abs(x)))^abdn.ML2$par[1], from = -max(my.labels), to = max(my.labels), col = "blue", lwd = 2, add = T) 
axis(side = 1, at = my.labels, labels = abs(my.labels), cex.axis = 1.8)
axis(side = 2, las = 1, cex.axis = 1.8)
box(lwd = 2)
# dev.off()

# (pow(mark[j], abd[0]) * exp(-abd[2] * pow(distance, 1) / pow(mark[j], abd[1])))
# pdf(paste(outputFile, "Signal.pdf", sep = ""))
maxS <- ceiling(60^abdn.ML2$par[1] * exp(-abdn.ML2$par[3] * abs(0) / 60^abdn.ML2$par[2]))
my.labels <- seq(-15, 15, 5)
par(lab = c(length(my.labels), 5, 7), mar = c(2, 3.5, 0.5, 0.5), mfrow = c(1, 1))
curve(60^abdn.ML2$par[1] * exp(-abdn.ML2$par[3] * abs(x) / 60^abdn.ML2$par[2]), from = -max(my.labels), to = max(my.labels), xlab = "", ylab = "", col = "black", xaxt = "n", axes = FALSE, lwd = 2, xlim = c(-max(my.labels), max(my.labels)), ylim = c(0, maxS))
curve(20^abdn.ML2$par[1] * exp(-abdn.ML2$par[3] * abs(x) / 20^abdn.ML2$par[2]), from = -max(my.labels), to = max(my.labels), col = "red", lwd = 2, add = T) 
curve(10^abdn.ML2$par[1] * exp(-abdn.ML2$par[3] * abs(x) / 10^abdn.ML2$par[2]), from = -max(my.labels), to = max(my.labels), col = "blue", lwd = 2, add = T) 
axis(side = 1, at = my.labels, labels = abs(my.labels), cex.axis = 1.8)
axis(side = 2, las = 1, cex.axis = 1.8)
box(lwd = 2)
# dev.off()

SNx.parameters <- list(L2 = abdn.L2, ML = abdn.ML, ML2 = abdn.ML2)



