# Initial code for random walk of NDVI model
library(rjags)
library(coda)
library(readxl)
library(lubridate)
library(ggplot2)
library(ecoforecastR)
# going to work with Thornhill data for now to test out the jags code

thorn.dat <- read.csv("input_data/NDVI_Thornhill.csv", header=T)
summary(thorn.dat)

#creating a date variable
thorn.dat$date <- as.Date(thorn.dat$system.index, format="%Y_%m_%d")
head(thorn.dat$date)

# creating new object that's parsed down a bit
thorn.dat2 <- thorn.dat[,c("date", "meanNDVI", "month", "year")]
summary(thorn.dat2)

thorn.dat2$doy <- yday(thorn.dat2$date) 

# creating a monthly mean value for NDVI
thorn.dat.month <- aggregate(meanNDVI~year + month, data=thorn.dat2, FUN="mean")
thorn.dat.month$meanNDVI.sd <- aggregate(meanNDVI~year + month, data=thorn.dat2, FUN="sd")$meanNDVI

# adding a date variable
thorn.dat.month$date <- paste(thorn.dat.month$year, thorn.dat.month$month, "01", sep="-")
thorn.dat.month$date <- as.Date(thorn.dat.month$date, format="%Y-%m-%d")
# making negative numbers equal to zero
thorn.dat2$meanNDVI[thorn.dat2$meanNDVI < 0] <- 0

saveRDS(thorn.dat2, "processed_data/thornhill_ndvi.rds")
saveRDS(thorn.dat.month, "processed_data/thornhill_ndvi_monthly.rds")

ggplot(data=thorn.dat.month) +
  geom_density(aes(x=meanNDVI))
ggplot(data=thorn.dat.month) +
  geom_line(aes(x=date, y=meanNDVI))

# Random Walk Code Example from mike's stateSpace script in EFI Activities

# Defining data for JAGS to reach out to
# ordering the df so it plots correctly
thorn.dat.month <- thorn.dat.month[order(thorn.dat.month$date, decreasing=F),]

time = as.Date(thorn.dat.month$date)
y = thorn.dat.month$meanNDVI
plot(time,y,type='l',ylab="Mean NDVI",lwd=2)


# Define the process model
RandomWalk = "
model{
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs)
  }
  
  #### Process Model
  for(t in 2:n){
    x[t]~dnorm(x[t-1],tau_add)
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
}
"

#Define the Data model
data <- list(y=y,n=length(y),           ## data
             x_ic=0.5,tau_ic=1,      ## initial condition prior # ROSS QUESTION: What sets the initial conditions?
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)

# Define the initial conditions
nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y,length(y),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(y.samp)),  ## initial guess on process precision
                    tau_obs=5/var(y.samp))        ## initial guess on obs precision
}

# Sending model info to JAGS
j.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                         inits = init,
                         n.chains = 3)

# Assessing model convergence
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","tau_obs"),
                            n.iter = 1000)
plot(jags.out)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 10000)

# Calculating credible interval and plotting output
time.rng = c(1,length(time))       ## adjust to zoom in and out
out <- as.matrix(jags.out)         ## convert from coda to matrix  
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975))

plot(time,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="Mean NDVI",xlim=time[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(time,y,pch="+",cex=0.5)


# Looking at posterior distributions - translating back into standard deviations from precisions.
hist(1/sqrt(out[,1]),main=colnames(out)[1])
hist(1/sqrt(out[,2]),main=colnames(out)[2])

# Looking at joint distribution of parameters and assessing covariance between the two
plot(out[,1],out[,2],pch=".",xlab=colnames(out)[1],ylab=colnames(out)[2])
cor(out[,1:2])


# Model with weather/climate data----
#########################################
# Using Mike's Dynamic linear model framework from Exercise 6

## grab weather data--Using VDP index that was calculated for the Morton
mort.clim.dat <- read.csv("../data_explore/input_data/morton_stacked_drought_indices.csv", header=T)

mort.vpd.index <- mort.clim.dat[mort.clim.dat$index=="vpd.index.value",]
# creating a date variable
mort.vpd.index$date <- as.Date(paste(mort.vpd.index$year, mort.vpd.index$month, "01", sep="-"))
summary(mort.vpd.index)


# refreshing on the format of object 'data'
head(data)

data$vpd = -1*mort.vpd.index$mean.month[match(time,mort.vpd.index$date)] # we have a monthly timestep in both the NDVI and the VPD data

## fit the model
ef.out <- ecoforecastR::fit_dlm(model=list(obs="y",fixed="~ 1 + X + vpd"),data)
names(ef.out)


## parameter diagnostics
params <- window(ef.out$params,start=1000) ## remove burn-in
plot(params)
summary(params)
cor(as.matrix(params))
pairs(as.matrix(params))

## confidence interval
out <- as.matrix(ef.out$predict)
ci <- apply(out,2,quantile,c(0.025,0.5,0.975))
plot(time,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="NDVI",xlim=time[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(time,y,pch="+",cex=0.5)
