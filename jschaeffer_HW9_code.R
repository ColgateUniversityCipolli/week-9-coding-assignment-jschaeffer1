library(tidyverse)
library(patchwork)
library(nleqslv)

################################################################################
# Precipitation in Madison County
################################################################################
dat.precip <- read_csv(file = "agacis.csv")

#####################################
# Clean Data
#####################################
dat.precip.long <- dat.precip |>    
  dplyr::select(-Annual) |>                   # Remove annual column 
  pivot_longer(cols = c(Jan, Feb, Mar, Apr,   # pivot the column data into one col
                        May, Jun, Jul, Aug, 
                        Sep, Oct, Nov, Dec), 
               values_to = "Precipitation",   # store the values in Precipitation
               names_to = "Month") |>         # store the months in Month
  mutate(Precipitation = case_when(Precipitation == "M" ~ NA_character_,
                                   TRUE                 ~ Precipitation))|>
  mutate(Precipitation = as.numeric(Precipitation))

view(dat.precip.long)

#######
#Known for Weibull
#######
weibull.alpha.hat = 2.1871
weibull.sigma.hat = 3.9683

logL_weibull = -2166.496



#################
#A: Gamma distribution
##################
###################
# MLE
###################
llgamma <- function(data, par, neg=F){
  alpha <- par[1]
  beta <- par[2]
  
  loglik <- sum(log(dgamma(x=data, shape=alpha, rate=beta)),  na.rm=T)
  
  return(ifelse(neg, -loglik, loglik))
}

(mles <- optim(par = c(1,1),
               fn = llgamma,
               data=dat.precip.long$Precipitation,
               neg=T))
alpha.hat.mle <- mles$par[1]
beta.hat.mle <- mles$par[2]


#################
#B: Log-Normal distribution
##################
###################
# MLE
###################
n=sum(!is.na(dat.precip.long$Precipitation))
mu.hat.mle    <- mean(log(dat.precip.long$Precipitation), na.rm=T)
sigma.hat.mle <- sqrt((1/n) * sum((log(dat.precip.long$Precipitation) 
                                   - mean(log(dat.precip.long$Precipitation), na.rm=T))^2, na.rm=T))

###################################
#Calculaing likelihoods
####################################

logL_gamma <- sum(dgamma(x = dat.precip.long$Precipitation, 
                         shape = alpha.hat.mle, 
                         rate = beta.hat.mle, 
                         log = TRUE), 
                  na.rm = TRUE)


logL_lognormal <- sum(dlnorm(x = dat.precip.long$Precipitation, 
                             meanlog = mu.hat.mle, 
                             sdlog = sigma.hat.mle, 
                             log = TRUE), 
                      na.rm = TRUE)


#######################
# C: Weibull and Gamma likelihood ratio
#######################
weibull_gamma_comp = exp(logL_weibull - logL_gamma)

########################
# D: Weibull and log-normal
########################
weibull_log_comp = exp(logL_weibull - logL_lognormal)

########################
# E: Gamma and log-normal
########################
gamma_log_comp = exp(logL_gamma - logL_lognormal)






#view(dat.precip.long)

