# Model with weather/climate data----
library(rjags)
library(coda)
library(readxl)
library(lubridate)
library(ggplot2)
library(tidyverse)
#########################################
# pulling inspiration from Exercise 5b from the EFI activities.

# Borrowing code from the EFI course. Using hind-casted data, will need to step through time
google.path <- file.path("G:/Shared drives/Urban Ecological Drought")

# loading in NDVI data

low.urban.ndvi <- read.csv(file.path(google.path, "Neighborhood remote sensing analysis/Landsat NDVI/Chi-NDVI-UrbanLow-Landsat_5-9.csv"), header=T)
head(low.urban.ndvi)

low.urban.ndvi$time <- as_date(low.urban.ndvi$time)
summary(low.urban.ndvi)

low.urban.ndvi$month <- month(low.urban.ndvi$time)

ggplot(data=low.urban.ndvi) +
  geom_line(aes(x=time, y=NDVI))

ggplot(data=low.urban.ndvi[is.na(low.urban.ndvi$NDVI),])+
  geom_histogram(aes(x=month), bins=12, col="green") # how can we handle NA's in the data frame?


# looking to see where NDVI dates are falling


mort.clim.dat <- read.csv(file.path(google.path,"data/r_files/input_files/data_explore/morton_stacked_drought_indices.csv"), header=T)
head(mort.clim.dat)


mort.vpd.index <- mort.clim.dat[mort.clim.dat$index=="vpd.index.value",]
# creating a date variable
mort.vpd.index$date <- as.Date(paste(mort.vpd.index$year, mort.vpd.index$month, "01", sep="-"))
summary(mort.vpd.index)

mort.vpd.index2 <- mort.vpd.index[mort.vpd.index$date %in% thorn.dat.month$date,]
head(mort.vpd.index2)

# creating ndvi index
mean.ndvi <- aggregate(meanNDVI ~ month, FUN=mean, data=thorn.dat.month)
head(mean.ndvi)

thorn.dat.month$mean.ndvi.agg <- mean.ndvi$meanNDVI[match(thorn.dat.month$month, mean.ndvi$month)]
head(thorn.dat.month)

thorn.dat.month$ndvi.index <- thorn.dat.month$meanNDVI/thorn.dat.month$mean.ndvi.agg

# making some of the data blank, because this is what we will be modeling
#thorn.dat.month$model.ndvi.index <- ifelse(thorn.dat.month$date >= "2020-01-01", NA, thorn.dat.month$ndvi.index)



library(ggplot2)
ggplot(data=thorn.dat.month) +
  geom_line(aes(x=date, y=meanNDVI ))

acf(thorn.dat.month$meanNDVI) # not much autocorrelation, but the detrending is pretty harsh here.

# need to combine ndvi and vpd index into a single list file
data <- list(ndvi = thorn.dat.month$model.ndvi[order(thorn.dat.month$date, decreasing = F) & thorn.dat.month$date < "2020-01-01"],
             vpd = mort.vpd.index2$mean.month[order(thorn.dat.month$date, decreasing = F)& thorn.dat.month$date < "2020-01-01"],
             #time = thorn.dat.month$date[order(thorn.dat.month$date, decreasing = F)],
             nt = length(thorn.dat.month$date[thorn.dat.month$date < "2020-01-01"]))

# quick lm to see what the relationship might look like 
test <- lm(thorn.dat.month$model.ndvi[order(thorn.dat.month$date, decreasing = F)] ~ mort.vpd.index2$mean.month[order(thorn.dat.month$date, decreasing = F)]-1)
summary(test)

# model from EFI course with a climate driver
# Sets the model structure
# The s loops won't apply here, because we aren't looking across sites or land cover type yet.
univariate_regression <- "
model{

  beta ~ dmnorm(b0,Vb)  	## multivariate Normal prior on vector of regression params
  prec ~ dgamma(s1,s2)    ## prior precision
  meow[1] <- 0.05356667 ## setting the seed value for NDVI
  
  for(i in 2:nt){
	  meow[i] <- beta[1] + meow[i-1] + beta[2]*vpd[i]   	## process model
	  ndvi[i]  ~ dnorm(meow[i],prec)		        ## data model
  }
}
"

## specify priors- uninformative

data$b0 <- as.vector(c(0,0))      ## regression beta means
data$Vb <- solve(diag(10000,2))   ## regression beta precisions
data$s1 <- 0.1                    ## error prior n/2
data$s2 <- 0.1                    ## error prior SS/2

## initial conditions
nchain = 3
inits <- list()
for(i in 1:nchain){
  inits[[i]] <- list(beta = rnorm(2,0,5), prec = runif(1,1/100,1/20))
}

j.model   <- jags.model(file = textConnection(univariate_regression),
                        data = data,
                        inits = inits,
                        n.chains = nchain)

var.out   <- coda.samples (model = j.model,
                           variable.names = c("beta","prec"),
                           n.iter = 10000)
## remember to assess convergence and remove burn-in before doing other diagnostics
GBR <- gelman.plot(var.out)
burnin <- 8000
var.burn <- window(var.out,start=burnin)
## convert to matrix
var.mat      <- as.matrix(var.burn)

## Pairwise scatter plots & correlation
pairs(var.mat)	## pairs plot to evaluate parameter correlation
cor(var.mat)    ## correlation matrix among model parameters


xpred <- mort.vpd.index2$mean.month[mort.vpd.index2$date >= "2020-01-01"]
plot(data$vpd, data$ndvi)
for(i in 1:10){
  lines(xpred, var.mat[i,"beta[1]"] + var.mat[i,"beta[2]"]*xpred)
}


nsamp <- 5000
samp <- sample.int(nrow(var.mat),nsamp)
xpred <- mort.vpd.index2$mean.month[mort.vpd.index2$date >= "2020-01-01"]					## sequence of x values we're going to
npred <- length(xpred)				##      make predictions for
ypred <- matrix(0.0,nrow=nsamp,ncol=npred)	## storage for predictive interval
ycred <- matrix(0.0,nrow=nsamp,ncol=npred)	## storage for credible interval

for(g in seq_len(nsamp)){
  theta = var.mat[samp[g],]
  ycred[g,] <- theta["beta[1]"] + theta["beta[2]"]*xpred
  ypred[g,] <- rnorm(npred,ycred[g,],1/sqrt(theta["prec"]))
}

ci <- apply(ycred,2,quantile,c(0.025,0.5,0.975))  ## credible interval and median
pi <- apply(ypred,2,quantile,c(0.025,0.975))		## prediction interval

plot(data$vpd,data$ndvi,cex=0.5,xlim=c(0,max(data$vpd)),ylim=c(0,max(data$ndvi, na.rm=T)))
lines(xpred,ci[1,],col=3,lty=2)	## lower CI
lines(xpred,ci[2,],col=3,lwd=3)	## median
lines(xpred,ci[3,],col=3,lty=2)	## upper CI
lines(xpred,pi[1,],col=4,lty=2)	## lower PI
lines(xpred,pi[2,],col=4,lty=2)	## upper PI
abline(b0,b1)				## true model
