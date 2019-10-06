##### Code to produce simulations for models A-F
##### requires inputs (inputs.RData) obtained from trial data

# clear background objects
rm(list=ls())

# libraries
library(survival)
library(dplyr)
library(magrittr)
library(ggplot2)
library(reshape2)


# set directory path
PATH = "V:/STATISTICS/STUDY FOLDER/Vasculitis/Ritazarem/Statistics/3 Analysis/Time Dependant Simulations"
setwd(PATH)


# load point estimate parameters, obtained from trial data via Weibull parametrisation
load("inputs.RData")


# simulate according to the two distributions using the point estimates loaded above

# Function to run simulations
# nSim = number of simulations
# N = sample size
# r = parametrisation of true HR (1=null hypothesis): HR = r^(-1/model$scale)
# tau = true change point (i.e. when HR changes)
# tend = analysis change point (i.e. assumed end of treatment, e.g. month 20)
run_sims <- function(nSim=1000, N=167, r=1, tau=20, tend=20) {
  
  #set up an array to store outputs
  cox.array <- array(0, dim=c(nSim,7))
  cox.t.array <- array(0, dim=c(nSim,14))
  
  # colnames
  colnames(cox.array) <- c("Aza", "Rit", "Total", "HR","lower CI","higher CI","p-value")
  colnames(cox.t.array) <- c("Aza.early", "Rit.early",
                             "Aza.late", "Rit.late",
                             "Total", "HR.early", "HR.late",
                             "lower CI.early", "lower CI.late",
                             "higher CI.early", "higher CI.late",
                             "p-value.early", "p-value.late",
                             "convergence")
  
  # simulated data
  for (i in 1:nSim) {
    
    # Now simulate according to the two distributions using the point estimates
    #sim_data <- data.frame(rx=primary0$rx)
    sim_data <- data.frame(rx=sample(c("Azathioprine", "Rituximab"), size = N, replace = TRUE, prob=c(1/2,1/2)))
    n <- nrow(sim_data)
    alpha = 1/m_relapse$scale
    lambda.e = rep(exp(m_relapse$coefficients), n)
    lambda.e = lambda.e*ifelse(sim_data$rx=="Rituximab", r[1], 1)
    relapse <- rweibull( n,shape=1/m_relapse$scale, scale=lambda.e )
    
    if(length(r)!=1){
      lambda.l = r[2]*exp(m_relapse$coefficients)
      U <- runif(n, 0, 1)
      relapse <- ifelse(sim_data$rx=="Rituximab",
                        ifelse(relapse>tau,
                               # conditional weibull distribution
                               lambda.l*( -log(U)+ (tau/lambda.l)^alpha)^(1/alpha),
                               relapse),
                        relapse)
    }
    
    censor  <- rweibull( n,shape=1/m_censor$scale,  scale=exp(m_censor$coefficients) )
    sim_data$subjid <- seq(1, nrow(sim_data))
    sim_data$time   <- pmin(relapse, censor)
    sim_data$status <- 1*(relapse<censor)
    
    
    ##### fit cox model (univariate) #####
    # fit cox proportional hazard model
    cox <- coxph(Surv(time, status)~rx, data=sim_data, timefix = FALSE)
    
    # store parameters
    obj_summary <- summary(cox)
    
    # events per arm
    cox.array[i,1:2] <- xtabs(~status+rx, sim_data)["1",]
    # total number of events
    cox.array[i,3] <- obj_summary$nevent
    # HR
    cox.array[i,4] <- obj_summary$coefficients[,2]
    # lower ci
    cox.array[i,5] <- obj_summary$conf.int[,3]
    # upper ci
    cox.array[i,6] <- obj_summary$conf.int[,4]
    # p-value
    cox.array[i,7] <- obj_summary$coefficients[,5]
    
    
    ##### fit cox model with time-varying covariate #####
    # set treatment endtime 
    sim_data$treatend <- tend
    
    # build dataframe with time-varing covariate (period)
    sim_data.t <-
      tmerge(sim_data,
             sim_data %>% mutate(period=treatend),
             id = subjid,
             status = event(time, status),
             period = tdc(period)) %>%
      mutate(
        period = factor(period, levels=0:1, labels=c("early","late"))
      )
    
    # use tryCatch to catch when the cox does not converge (warnings) 
    # & when there are errors (simulate time interval of size 0)
    tryCatch(
      expr = {
        # fit cox proportional hazard model with time-dependent covariate
        cox.t <- coxph(Surv(tstart, tstop, status) ~ period / rx, data = sim_data.t, timefix = FALSE)
        # store parameters
        obj_summary <- summary(cox.t)
        # events per arm
        cox.t.array[i, 1:2]   <- xtabs( ~ period + rx + status, sim_data.t)["early", , "1"]
        cox.t.array[i, 3:4]   <- xtabs( ~ period + rx + status, sim_data.t)["late", , "1"]
        # total number of events
        cox.t.array[i, 5]     <- obj_summary$nevent
        # hr
        cox.t.array[i, 6:7]   <- obj_summary$coefficients[2:3, 2]
        # lower ci
        cox.t.array[i, 8:9]   <- obj_summary$conf.int[2:3, 3]
        # upper ci
        cox.t.array[i, 10:11] <- obj_summary$conf.int[2:3, 4]
        # p-value
        cox.t.array[i, 12:13] <- obj_summary$coefficients[2:3, 5]
        # convergence
        cox.t.array[i, 14]    <- 1
      },
      warning = function(e) {
        cox.t.array[i, 1:2]   <- xtabs( ~ period + rx + status, sim_data.t)["early", , "1"]
        cox.t.array[i, 3:4]   <- xtabs( ~ period + rx + status, sim_data.t)["late", , "1"]
        cox.t.array[i, 5]     <- addmargins(xtabs( ~ period + status, sim_data.t))["Sum", "1"]
        cox.t.array[i, 6:13]  <- NA
        cox.t.array[i, 14]    <- 0
      },
      error = function(e) {
        cox.t.array[i, 1:2]  <- xtabs( ~ period + rx + status, sim_data.t)["early", , "1"]
        cox.t.array[i, 3:4]  <- xtabs( ~ period + rx + status, sim_data.t)["late", , "1"]
        cox.t.array[i, 5]    <- addmargins(xtabs( ~ period + status, sim_data.t))["Sum", "1"]
        cox.t.array[i, 6:13] <- NA
        cox.t.array[i, 14]   <- 0
      }
    )
    
    
    # change array to dataframe
    cox.array <- as.data.frame(cox.array)
    cox.t.array <- as.data.frame(cox.t.array)
    
  }
  
  list(nSim=nSim, N=N, tend=tend, r=r, tau=tau, cox=cox.array, cox.t=cox.t.array)
  
}


# examples (null hypothesis), sims.1 should yield same results as sim.2
set.seed(1549)
sims.1 <- run_sims(nSim=1000, N=167, r=1, tau=20, tend=20)
set.seed(1549)
sims.2 <- run_sims(nSim=1000, N=167, r=c(1,1), tend=20, tau=20)


# paremeters for simulations
nSim <- c(100, 1000, 10000)
N    <- c(100, 167, 300, 500)
tend <- seq(10, 30, by= 1)
r    <- rev(c(seq(0.3, 1, by=0.2),1)^(-m_relapse$scale))
tau  <- seq(10, 30, by= 1)


####### Model 1: simulate the null hypothesis #######
# explores the effect of sample size

# set seed
set.seed(1549)

# simulate using different values of N and tend
sims <- vector("list", 0)
for (i in 1:length(tend)) {
  tmp <- vector("list", length(N))
  for (j in 1:length(N)) {
    tmp[[j]] <- run_sims(nSim=1000, N=N[j], r=1, tau=20, tend=tend[i])
  }
  sims <- c(sims, tmp)
  print(tend[i])
}
rm(tmp)

# save simulations
save(sims, file="sims_1.RData")



####### Model 2: constant overall treatment effect #######
# explores the effect of sample size on different treatment effects

# set seed
set.seed(1549)

# simulate using different values of N, tend and r
sims <- vector("list", 0)
for (k in 1:length(r)) {
  for (i in 1:length(tend)) {
    tmp <- vector("list", length(N))
    for (j in 1:length(N)) {
      tmp[[j]] <- run_sims(nSim=1000, N=N[j], r=r[k], tau=20, tend=tend[i])
    }
    sims <- c(sims, tmp)
  }
  rm(tmp)
  print(r[k])
}
rm(tmp)

# save simulations
save(sims, file="sims_2.RData")



####### Model 3: late efficacy (HRearly = 1), treatment effect only after month 20) #######

# set seed
set.seed(1549)

# simulate using different values of N, r and tend
sims <- vector("list", 0)
for (k in 1:length(N)) {
  for (i in 1:length(tend)) {
    tmp <- vector("list", length(r))
    for (j in 1:length(r)) {
      tmp[[j]] <- run_sims(nSim=1000, N=N[k], r=c(1,r[j]), tau=20, tend=tend[i])
    }
    sims <- c(sims, tmp)
    print(paste0(N[k],", ", tend[i]))
  }
}
rm(tmp)

# save simulations
save(sims, file="sims_3.RData")



####### Model 4: early efficacy (HRlate = 1), treatment effect up until month 20) #######

# set seed
set.seed(1549)

# simulate using different values of N, r and tend
sims <- vector("list", 0)
for (k in 1:length(N)) {
  for (i in 1:length(tend)) {
    tmp <- vector("list", length(r))
    for (j in 1:length(r)) {
      tmp[[j]] <- run_sims(nSim=1000, N=N[k], r=c(r[j],1), tau=20, tend=tend[i])
    }
    sims <- c(sims, tmp)
    print(paste0(N[k],", ", tend[i]))
  }
}
rm(tmp)

# save simulations
save(sims, file="sims_4.RData")




####### Model 5: different treatment effects before/after month 20 #######

# set seed
set.seed(1549)

# using different values of N, r and tend
sims <- vector("list", 0)
for (k in 1:length(N)) {
  for (i in 1:length(tend)) {
    tmp <- vector("list", length(r))
    for (j in 1:length(r)) {
      tmp[[j]] <- run_sims(nSim=1000, N=N[k], r=c(r[j],rev(r)[j]), tau=20, tend=tend[i])
    }
    sims <- c(sims, tmp)
    print(paste0(N[k],", ", tend[i]))
  }
}
rm(tmp)

# save simulations
save(sims, file="sims_5.RData")




####### Model 6: different treatment effects, wrong change point tau #######

# set seed
set.seed(1549)

# using different values of N, r and tau
sims <- vector("list", 0)
for (k in 1:length(N)) {
  for (i in 1:length(tau)) {
    tmp <- vector("list", length(r))
    for (j in 1:length(r)) {
      tmp[[j]] <- run_sims(nSim=1000, N=N[k], r=c(r[j],rev(r)[j]), tau=tau[i], tend=20)
    }
    sims <- c(sims, tmp)
    print(paste0(N[k],", ", tau[i]))
  }
}
rm(tmp)

# save simulations
save(sims, file="sims_6.RData")

