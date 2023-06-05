# Model with weather/climate data----
library(rjags)
library(coda)
library(readxl)
library(lubridate)
library(ggplot2)
library(tidyverse)
#########################################
# Borrowing code from teh EFI course. Using hind-casted data, will need to step through time

# loading in NDVI data
thorn.dat.month <- readRDS("processed_data/thornhill_ndvi_monthly.rds")
summary(thorn.dat.month)
mort.clim.dat <- read.csv("../data_explore/input_data/morton_stacked_drought_indices.csv", header=T)

mort.vpd.index <- mort.clim.dat[mort.clim.dat$index=="vpd.index.value",]
# creating a date variable
mort.vpd.index$date <- as.Date(paste(mort.vpd.index$year, mort.vpd.index$month, "01", sep="-"))
summary(mort.vpd.index)

mort.vpd.index2 <- mort.vpd.index[mort.vpd.index$date %in% thorn.dat.month$date,]


# model from EFI course with a climate driver
# Sets the model structure
model_RandomWalk <- function() {
  "
model{
  
  for (s in 1:ns) {
    #### Data Model
    for(t in 1:nt){
      ndvi[t, s] ~ dnorm(x[t, s],tau_obs_ndvi)
      vpd[t, s] ~ dnorm(x[t, s],tau_obs_vpd)
    }
    
    #### Process Model
    for(t in 2:nt){
      x[t, s]~dnorm(x[t-1, s],tau_add)
    }
  }
  
  
  #### Priors
  for (s in 1:ns) {
     x[1, s] ~ dnorm(x_ic[s],tau_ic[s])
  }
  tau_obs_ndvi ~ dgamma(a_obs_ndvi,r_obs_ndvi)
  tau_obs_vpd ~ dgamma(a_obs_vpd,r_obs_vpd)
  tau_add ~ dgamma(a_add,r_add)
}
"
}

# creating a list of dates to simulate through
# date_list <- unique(thorn.dat.month$date)

all_dat <- thorn.dat.month
all_dat$vpd <- mort.vpd.index2$mean.month[match(mort.vpd.index2$date,thorn.dat.month$date)]
outdir <- "processed_data/forecast_test/"

# starting point
mindate = min(all_dat$date)
batch = 1
date_list = seq(mindate + months(batch), max(thorn.dat.month$date), by = paste0(batch, " month"))
RandomWalk = model_RandomWalk()
# forecasts = run_model(date_list, RandomWalk, all_dat, batch, mindate)
# tar_render(EDA, "docs/EDA.Rmd", output_format = "all"),
# tar_render(README, "README.Rmd", output_format = "all")

# Runs the model through time----
run_model <- function(date_list, RandomWalk, all_dat, batch, mindate, outdir = "./forecasts/") {
  for (d in 1:length(date_list)) {
    today <- date_list[d]
    
    # create dir
    dir.create(paste0(outdir, today), recursive = T)
    
    # subset data
    dat_new <- all_dat %>%
      filter(date >= today-months(batch)) %>%
      mutate( # set data after today to NA
        ndvi = case_when(date < today ~ meanNDVI),
        #ndvi_sd = case_when(date < today ~ ndvi_sd),
        vpd2 = case_when(date < today ~ vpd),
        #vpd_sd = case_when(date < today ~ vpd_sd)
      )
    
    # make matrices for observations
    ndvi <- dat_new %>%
      dplyr::select(date, ndvi) %>%
      #spread(key = "site", value = "ndvi_90") %>%
      dplyr::select(-date) %>%
      as.matrix()
    vpd <- dat_new %>%
      dplyr::select(date, vpd2) %>%
      #spread(key = "site", value = "vpd") %>%
      dplyr::select(-date) %>%
      as.matrix()
    
    # load prvpdous data
    if (d ==1) {
      
      dat_old <- all_dat %>%
        filter(date>=mindate) %>% 
        head(0)
    } else {
      prev_day <- date_list[d - 1]
      dat_old<-read_rds(paste0(outdir, prev_day, "/data.rds")) %>% 
        filter(date<prev_day)
    }
    
    # load priors
    if (d == 1) {
      # initialize
      priors <- list(
        x_ic = rep(0, ncol(ndvi)), tau_ic = rep(10, ncol(ndvi)), ## initial condition prior
        a_obs_ndvi = 1, r_obs_ndvi = 1, ## obs error prior
        a_obs_vpd = 1, r_obs_vpd = 1, ## obs error prior
        a_add = 1, r_add = 1 ## process error prior
      )
    } else {
      # load prvpdous model
      prev_day <- date_list[d - 1]
      priors <- read_rds(paste0(outdir, prev_day, "/posterior.rds"))
    }
    
    
    # Setting data and prior for the model
    data <- c(
      list(
        ndvi = ndvi,
        vpd = vpd,
        nt = nrow(ndvi),
        ns = ncol(ndvi) ## data
      ),
      priors
    )
    
    
    ####
    # Defining the initial state of each of the chains we are going to run
    
    nchain <- 3
    init <- list()
    for (i in 1:nchain) {
      # y.samp = sample(y,length(y),replace=TRUE)
      # init[[i]] <- list(tau_add=1/var(diff(y.samp)),  ## initial guess on process precision
      #                   tau_obs=5/var(y.samp))        ## initial guess on obs precision
      init[[i]] <- list(
        tau_add = 500, ## initial guess on process precision
        tau_obs_ndvi = 500,
        tau_obs_vpd = 500
      ) ## initial guess on obs precision
    }
    
    ####
    # Presenting basic model to JAGS
    j.model <- jags.model(
      file = textConnection(RandomWalk),
      data = data,
      inits = init,
      n.chains = 3
    )
    
    # run MCMC
    jags.out <- coda.samples(
      model = j.model,
      variable.names = c("x", "tau_add", "tau_obs_ndvi", "tau_obs_vpd"),
      n.iter = 1000
    )
    out <- as.matrix(jags.out) ## convert from coda to matrix
    
    # Save posteriors
    x_ic <- tau_ic <- rep(NA, ncol(ndvi))
    for (s in 1:ncol(ndvi)) {
      x.col <- paste0("x[", batch+1, ",", s, "]")
      x_ic[s] <- mean(out[, x.col])
      tau_ic[s] <- 1 / sd(out[, x.col])^2
    }
    mu_obs_ndvi<-mean(out[, "tau_obs_ndvi"])
    sd_obs_ndvi<-sd(out[, "tau_obs_ndvi"])
    mu_obs_vpd<-mean(out[, "tau_obs_vpd"])
    sd_obs_vpd<-sd(out[, "tau_obs_vpd"])
    mu_add<-mean(out[, "tau_add"])
    sd_add<-sd(out[, "tau_add"])
    
    posteriors <- list(
      x_ic = x_ic, tau_ic = tau_ic, ## initial condition prior
      a_obs_ndvi = (mu_obs_ndvi/sd_obs_ndvi)^2,
      r_obs_ndvi =mu_obs_ndvi/sd_obs_ndvi^2, ## obs error prior
      a_obs_vpd = (mu_obs_vpd/sd_obs_vpd)^2,
      r_obs_vpd = mu_obs_vpd/sd_obs_vpd^2, ## obs error prior
      a_add = (mu_add/sd_add)^2,
      r_add = mu_add/sd_add^2
    )
    write_rds(posteriors, paste0(outdir, today, "/posterior.rds"))
    
    # plot
    x.cols <- grep("^x", colnames(out)) ## grab all columns that start with the letter x
    ci <- apply(out[, x.cols], 2, quantile, c(0.025, 0.5, 0.975)) ## model was fit on log scale
    
    dat_fitted <- bind_rows(dat_old,
                            dat_new %>%
                              cbind(t(ci))
    )
    write_rds(dat_fitted, paste0(outdir, today, "/data.rds"))
    
    p <- ggplot(dat_fitted) +
      geom_line(aes(x = date, y = ndvi), col = "dark green") +
      geom_point(aes(x = date, y = vpd), col = "light green") +
      geom_ribbon(aes(x = date, ymin = `2.5%`, ymax = `97.5%`), fill = "blue", alpha = 0.5) +
      geom_line(aes(x = date, y = `50%`), col = "blue") +
      scale_y_continuous(limits = c(-80, 80))
      #facet_wrap(. ~ site) +
      theme_classic()
    ggsave(paste0(outdir, today, "/plot.pdf"), p)
    ggsave(paste0("figures/test_forecast/", today, "_plot.png"), p)
    print (paste0(today, " completed."))
  }
}

