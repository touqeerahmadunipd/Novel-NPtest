rm(list=ls());gc() 
setwd("~/Documents/Touqeer-Docs/Partition D/Joint works/Dr. Abid/Githubcode")
source("Functions_tests.R")


#install.packages("survival") 
#install.packages("flexsurv")
#install.packages("ggamma")
#install.packages("PHInfiniteEstimates")
#install.packages("tpwb")

library(survival)  # for survdiff, Surv
library(ggamma)
library(flexsurv)
library(PHInfiniteEstimates)
library(tpwb)
######################################
##.    Case I------------------------
#####################################

# Parameters
n_sim <- 10000
n1 <- 200; n2 <- 200
#parameters for T1 and T2
#T1
scale_weibull_T1<- 2
shape_weibull_T1<- 2
#T2
scale_gamma_T2<- 0.557706
shape_gamma_T2<- 3.12154

#Parameters for C1 and C2
#C1
location_weibull_C1 = 0
scale_weibull_C1<- 5.795
shape_weibull_C1<- 2
#C2
location_weibull_C2 = 0
scale_weibull_C2<- 5.847
shape_weibull_C2<- 2

# ------------------------
# Check if censoring is applied
# ------------------------

apply_censoring <- TRUE # set to FALSE for no censoring

# Initialize matrix to store p-values
pvals <- matrix(NA, nrow = n_sim, ncol = 5)
colnames(pvals) <- c("Gehan","CoxMantel", "Logrank", "PetoPeto", "ModifiedWcox")
cen_per <- matrix(NA, nrow = n_sim, ncol = 3)
colnames(cen_per) <- c("% censoring X", "% censoring in Y", "overall % censoring")

for(i in 1:n_sim){
  
  # ------------------------
  # Event times
  # ------------------------
  T1 <- rweibull(n1, scale = scale_weibull_T1, shape = shape_weibull_T1) #weibull
  T2 <- rgamma(n2, scale = scale_gamma_T2, shape = shape_gamma_T2)       #gamma
  
  #Censoring application
  
  if(apply_censoring){
    # With censoring
    #C1 <- rweibull(n1, scale = scale_weibull_C1, shape = shape_weibull_C1)
    #C2 <- rweibull(n2, scale = scale_weibull_C2, shape = shape_weibull_C2)
    
    C1 <- rtpwb(n1, location = location_weibull_C1, scale = scale_weibull_C1, shape = shape_weibull_C1 )
    C2 <- rtpwb(n2, location = location_weibull_C2, scale = scale_weibull_C2, shape = shape_weibull_C2 )
    
    timeX <- pmin(T1, C1)
    statusX <- as.integer(T1 <= C1)
    timeY <- pmin(T2, C2)
    statusY <- as.integer(T2 <= C2)
    
  } else {
    # No censoring
    timeX <- T1
    statusX <- rep(1, n1)
    timeY <- T2
    statusY <- rep(1, n2)
  }
  
  # ------------------------
  # Calculate censoring percentages (optional if no censoring)
  # ------------------------
  cen_per[i, 1] <- mean(statusX == 0) * 100
  cen_per[i, 2] <- mean(statusY == 0) * 100
  cen_per[i, 3] <- mean(c(statusX, statusY) == 0) * 100
  
  # ------------------------
  # Combine data
  # ------------------------
  data <- data.frame(
    time = c(timeX, timeY),
    status = c(statusX, statusY),
    group = factor(c(rep("X", n1), rep("Y", n2)))
  )
  
  # ------------------------
  # Run survival tests
  # ------------------------
  pvals[i, "Gehan"] <- gehan.wilcoxon.test(Surv(time, status) ~ group, data = data, gehan = TRUE)$p.value
  pvals[i, "CoxMantel"] <- cox_mantel_test(timeX, timeY, statusX, statusY)$p.value
  pvals[i, "Logrank"] <- logrank_chisq(timeX, timeY, statusX, statusY)$p.value
  test_peto <- survdiff(Surv(time, status) ~ group, rho = 1, data = data)
  pvals[i, "PetoPeto"] <- 1 - pchisq(test_peto$chisq, length(test_peto$n) - 1)
  pvals[i, "ModifiedWcox"] <- modified_mw_test(timeX, timeY, statusX, statusY)$p.value
}


# Rejection rates
list(colMeans(pvals <= 0.05, na.rm = TRUE),
     colMeans(cen_per, na.rm = TRUE))





######################################
##.    Case II------------------------
######################################

# Parameters
n_sim <- 5000
n1 <- 200; n2 <- 200
#parameters for T1 and T2
#T1
scale_weibull_T1<- 2
shape_weibull_T1<- 2
#T2
mean_lognorm_T2  <- 0.4096
scale_lognorm_T2 <- 0.6179

#Parameters for C1 and C2
#C1
location_weibull_C1 = 0
scale_weibull_C1<- 5.795
shape_weibull_C1<- 2
#C2
location_weibull_C2 = 0.410
scale_weibull_C2<- 5.220
shape_weibull_C2<- 2

# ------------------------
# Check if censoring is applied
# ------------------------

apply_censoring <- TRUE # set to FALSE for no censoring

# Initialize matrix to store p-values
pvals <- matrix(NA, nrow = n_sim, ncol = 5)
colnames(pvals) <- c("Gehan","CoxMantel", "Logrank", "PetoPeto", "ModifiedWcox")
cen_per <- matrix(NA, nrow = n_sim, ncol = 3)
colnames(cen_per) <- c("% censoring X", "% censoring in Y", "overall % censoring")

for(i in 1:n_sim){
  
  # ------------------------
  # Event times
  # ------------------------
  T1 <- rweibull(n1, scale = scale_weibull_T1, shape = shape_weibull_T1) #weibull
  T2 <- rlnorm(n2, mean_lognorm_T2, scale_lognorm_T2)
  
  #Censoring
  
  if(apply_censoring){
    # With censoring
    C1 <- rtpwb(n1, location = location_weibull_C1, scale = scale_weibull_C1, shape = shape_weibull_C1 )
    C2 <- rtpwb(n2, location = location_weibull_C2, scale = scale_weibull_C2, shape = shape_weibull_C2 )
    
    timeX <- pmin(T1, C1)
    statusX <- as.integer(T1 <= C1)
    timeY <- pmin(T2, C2)
    statusY <- as.integer(T2 <= C2)
    
  } else {
    # No censoring
    timeX <- T1
    statusX <- rep(1, n1)
    timeY <- T2
    statusY <- rep(1, n2)
  }
  
  # ------------------------
  # Calculate censoring percentages (optional if no censoring)
  # ------------------------
  cen_per[i, 1] <- mean(statusX == 0) * 100
  cen_per[i, 2] <- mean(statusY == 0) * 100
  cen_per[i, 3] <- mean(c(statusX, statusY) == 0) * 100
  
  # ------------------------
  # Combine data
  # ------------------------
  data <- data.frame(
    time = c(timeX, timeY),
    status = c(statusX, statusY),
    group = factor(c(rep("X", n1), rep("Y", n2)))
  )
  
  # ------------------------
  # Run survival tests
  # ------------------------
  pvals[i, "Gehan"] <- gehan.wilcoxon.test(Surv(time, status) ~ group, data = data, gehan = TRUE)$p.value
  pvals[i, "CoxMantel"] <- cox_mantel_test(timeX, timeY, statusX, statusY)$p.value
  pvals[i, "Logrank"] <- logrank_chisq(timeX, timeY, statusX, statusY)$p.value
  test_peto <- survdiff(Surv(time, status) ~ group, rho = 1, data = data)
  pvals[i, "PetoPeto"] <- 1 - pchisq(test_peto$chisq, length(test_peto$n) - 1)
  pvals[i, "ModifiedWcox"] <- modified_mw_test(timeX, timeY, statusX, statusY)$p.value
}


# Rejection rates
list(colMeans(pvals <= 0.05, na.rm = TRUE),
     colMeans(cen_per, na.rm = TRUE))









######################################
##.    Case III------------------------
######################################

# Parameters
n_sim <- 5000
n1 <- 200; n2 <- 200
#parameters for T1 and T2
#T1
scale_weibull_T1<- 2
shape_weibull_T1<- 2
#T2
scale_weibull_T2<- 2.3
shape_weibull_T2<- 2.4

#Parameters for C1 and C2
#C1
location_weibull_C1 = 0
scale_weibull_C1<- 5.795
shape_weibull_C1<- 2
#C2
location_weibull_C2 = 0.0
scale_weibull_C2<- 3.736
shape_weibull_C2<- 6.340

# ------------------------
# Check if censoring is applied
# ------------------------

apply_censoring <- TRUE # set to FALSE for no censoring

# Initialize matrix to store p-values
pvals <- matrix(NA, nrow = n_sim, ncol = 5)
colnames(pvals) <- c("Gehan","CoxMantel", "Logrank", "PetoPeto", "ModifiedWcox")
cen_per <- matrix(NA, nrow = n_sim, ncol = 3)
colnames(cen_per) <- c("% censoring X", "% censoring in Y", "overall % censoring")

for(i in 1:n_sim){
  
  # ------------------------
  # Event times
  # ------------------------
  T1 <- rweibull(n1, scale = scale_weibull_T1, shape = shape_weibull_T1) #weibull
  T2 <- rweibull(n1, scale = scale_weibull_T2, shape = shape_weibull_T2) #weibull
  
  #Censoring
  
  if(apply_censoring){
    # With censoring
    C1 <- rtpwb(n1, location = location_weibull_C1, scale = scale_weibull_C1, shape = shape_weibull_C1 )
    C2 <- rtpwb(n2, location = location_weibull_C2, scale = scale_weibull_C2, shape = shape_weibull_C2 )
    
    timeX <- pmin(T1, C1)
    statusX <- as.integer(T1 <= C1)
    timeY <- pmin(T2, C2)
    statusY <- as.integer(T2 <= C2)
    
  } else {
    # No censoring
    timeX <- T1
    statusX <- rep(1, n1)
    timeY <- T2
    statusY <- rep(1, n2)
  }
  
  # ------------------------
  # Calculate censoring percentages (optional if no censoring)
  # ------------------------
  cen_per[i, 1] <- mean(statusX == 0) * 100
  cen_per[i, 2] <- mean(statusY == 0) * 100
  cen_per[i, 3] <- mean(c(statusX, statusY) == 0) * 100
  
  # ------------------------
  # Combine data
  # ------------------------
  data <- data.frame(
    time = c(timeX, timeY),
    status = c(statusX, statusY),
    group = factor(c(rep("X", n1), rep("Y", n2)))
  )
  
  # ------------------------
  # Run survival tests
  # ------------------------
  pvals[i, "Gehan"] <- gehan.wilcoxon.test(Surv(time, status) ~ group, data = data, gehan = TRUE)$p.value
  pvals[i, "CoxMantel"] <- cox_mantel_test(timeX, timeY, statusX, statusY)$p.value
  pvals[i, "Logrank"] <- logrank_chisq(timeX, timeY, statusX, statusY)$p.value
  test_peto <- survdiff(Surv(time, status) ~ group, rho = 1, data = data)
  pvals[i, "PetoPeto"] <- 1 - pchisq(test_peto$chisq, length(test_peto$n) - 1)
  pvals[i, "ModifiedWcox"] <- modified_mw_test(timeX, timeY, statusX, statusY)$p.value
}


# Rejection rates
list(colMeans(pvals <= 0.05, na.rm = TRUE),
     colMeans(cen_per, na.rm = TRUE))






######################################
##.    Case IV------------------------
######################################

# Parameters
n_sim <- 5000
n1 <- 200; n2 <- 200
#parameters for T1 and T2
#T1
scale_weibull_T1<- 2.1
shape_weibull_T1<- 2.1
#T2
scale_weibull_T2<- 1.75
shape_weibull_T2<- 2.1

#Parameters for C1 and C2
#C1
location_weibull_C1 = 0
scale_weibull_C1<- 6.227
shape_weibull_C1<- 2
#C2
location_weibull_C2 = 0.0
scale_weibull_C2<- 5.301
shape_weibull_C2<- 2

# ------------------------
# Check if censoring is applied
# ------------------------

apply_censoring <- TRUE # set to FALSE for no censoring

# Initialize matrix to store p-values
pvals <- matrix(NA, nrow = n_sim, ncol = 5)
colnames(pvals) <- c("Gehan","CoxMantel", "Logrank", "PetoPeto", "ModifiedWcox")
cen_per <- matrix(NA, nrow = n_sim, ncol = 3)
colnames(cen_per) <- c("% censoring X", "% censoring in Y", "overall % censoring")

for(i in 1:n_sim){
  
  # ------------------------
  # Event times
  # ------------------------
  T1 <- rweibull(n1, scale = scale_weibull_T1, shape = shape_weibull_T1) #weibull
  T2 <- rweibull(n1, scale = scale_weibull_T2, shape = shape_weibull_T2) #weibull
  
  #Censoring
  
  if(apply_censoring){
    # With censoring
    C1 <- rtpwb(n1, location = location_weibull_C1, scale = scale_weibull_C1, shape = shape_weibull_C1 )
    C2 <- rtpwb(n2, location = location_weibull_C2, scale = scale_weibull_C2, shape = shape_weibull_C2 )
    
    timeX <- pmin(T1, C1)
    statusX <- as.integer(T1 <= C1)
    timeY <- pmin(T2, C2)
    statusY <- as.integer(T2 <= C2)
    
  } else {
    # No censoring
    timeX <- T1
    statusX <- rep(1, n1)
    timeY <- T2
    statusY <- rep(1, n2)
  }
  
  # ------------------------
  # Calculate censoring percentages (optional if no censoring)
  # ------------------------
  cen_per[i, 1] <- mean(statusX == 0) * 100
  cen_per[i, 2] <- mean(statusY == 0) * 100
  cen_per[i, 3] <- mean(c(statusX, statusY) == 0) * 100
  
  # ------------------------
  # Combine data
  # ------------------------
  data <- data.frame(
    time = c(timeX, timeY),
    status = c(statusX, statusY),
    group = factor(c(rep("X", n1), rep("Y", n2)))
  )
  
  # ------------------------
  # Run survival tests
  # ------------------------
  pvals[i, "Gehan"] <- gehan.wilcoxon.test(Surv(time, status) ~ group, data = data, gehan = TRUE)$p.value
  pvals[i, "CoxMantel"] <- cox_mantel_test(timeX, timeY, statusX, statusY)$p.value
  pvals[i, "Logrank"] <- logrank_chisq(timeX, timeY, statusX, statusY)$p.value
  test_peto <- survdiff(Surv(time, status) ~ group, rho = 1, data = data)
  pvals[i, "PetoPeto"] <- 1 - pchisq(test_peto$chisq, length(test_peto$n) - 1)
  pvals[i, "ModifiedWcox"] <- modified_mw_test(timeX, timeY, statusX, statusY)$p.value
}


# Rejection rates
list(colMeans(pvals <= 0.05, na.rm = TRUE),
     colMeans(cen_per, na.rm = TRUE))






######################################
##.    Case V------------------------
######################################

# Parameters
n_sim <- 5000
n1 <- 200; n2 <- 200
#parameters for T1 and T2
#T1
scale_weibull_T1<- 1
shape_weibull_T1<- 1.1
#T2
scale_weibull_T2<- 0.7
shape_weibull_T2<- 0.9

#Parameters for C1 and C2
#C1
location_weibull_C1 = 0
scale_weibull_C1<- 3.724
shape_weibull_C1<- 2
#C2
location_weibull_C2 = 0.0
scale_weibull_C2<- 2.819
shape_weibull_C2<- 2

# ------------------------
# Check if censoring is applied
# ------------------------

apply_censoring <- TRUE # set to FALSE for no censoring

# Initialize matrix to store p-values
pvals <- matrix(NA, nrow = n_sim, ncol = 5)
colnames(pvals) <- c("Gehan","CoxMantel", "Logrank", "PetoPeto", "ModifiedWcox")
cen_per <- matrix(NA, nrow = n_sim, ncol = 3)
colnames(cen_per) <- c("% censoring X", "% censoring in Y", "overall % censoring")

for(i in 1:n_sim){
  
  # ------------------------
  # Event times
  # ------------------------
  T1 <- rweibull(n1, scale = scale_weibull_T1, shape = shape_weibull_T1) #weibull
  T2 <- rweibull(n1, scale = scale_weibull_T2, shape = shape_weibull_T2) #weibull
  
  #Censoring
  
  if(apply_censoring){
    # With censoring
    C1 <- rtpwb(n1, location = location_weibull_C1, scale = scale_weibull_C1, shape = shape_weibull_C1 )
    C2 <- rtpwb(n2, location = location_weibull_C2, scale = scale_weibull_C2, shape = shape_weibull_C2 )
    
    timeX <- pmin(T1, C1)
    statusX <- as.integer(T1 <= C1)
    timeY <- pmin(T2, C2)
    statusY <- as.integer(T2 <= C2)
    
  } else {
    # No censoring
    timeX <- T1
    statusX <- rep(1, n1)
    timeY <- T2
    statusY <- rep(1, n2)
  }
  
  # ------------------------
  # Calculate censoring percentages (optional if no censoring)
  # ------------------------
  cen_per[i, 1] <- mean(statusX == 0) * 100
  cen_per[i, 2] <- mean(statusY == 0) * 100
  cen_per[i, 3] <- mean(c(statusX, statusY) == 0) * 100
  
  # ------------------------
  # Combine data
  # ------------------------
  data <- data.frame(
    time = c(timeX, timeY),
    status = c(statusX, statusY),
    group = factor(c(rep("X", n1), rep("Y", n2)))
  )
  
  # ------------------------
  # Run survival tests
  # ------------------------
  pvals[i, "Gehan"] <- gehan.wilcoxon.test(Surv(time, status) ~ group, data = data, gehan = TRUE)$p.value
  pvals[i, "CoxMantel"] <- cox_mantel_test(timeX, timeY, statusX, statusY)$p.value
  pvals[i, "Logrank"] <- logrank_chisq(timeX, timeY, statusX, statusY)$p.value
  test_peto <- survdiff(Surv(time, status) ~ group, rho = 1, data = data)
  pvals[i, "PetoPeto"] <- 1 - pchisq(test_peto$chisq, length(test_peto$n) - 1)
  pvals[i, "ModifiedWcox"] <- modified_mw_test(timeX, timeY, statusX, statusY)$p.value
}


# Rejection rates
list(colMeans(pvals <= 0.05, na.rm = TRUE),
     colMeans(cen_per, na.rm = TRUE))








