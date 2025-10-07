rm(list=ls());gc() 
setwd("~/Documents/Touqeer-Docs/Partition D/Joint works/Dr. Abid/Githubcode")
source("Functions_tests.R")


#Figure 1-----------------
#Case I
t <- seq(0,5,length.out=1000)
p_H1 <- plot_survival(
  t,
  surv_list = list(
    1 - pweibull(t, shape=2, scale=2),
    pgamma(t, shape=3.12154, scale=0.557706, lower.tail=FALSE)
  ),
  sample_names = c("Weibull", "Gamma"),
  title = "Case I"
)
p_H1

#Case II
p_H2 <- plot_survival(
  t,
  surv_list = list(
    1 - pweibull(t, shape=2, scale=2),
    plnorm(t, meanlog=0.4096, sdlog=0.6179, lower.tail=FALSE)
  ),
  sample_names = c("Weibull", "Lognormal"),
  title = "Case II"
)
p_H2

# Case III

p_H3 <- plot_survival(
  t,
  surv_list = list(
    1- pweibull(t, 2, 2),
    1- pweibull(t, 2.4, 2.3)
  ),
  sample_names = c("Weibull", "Weibull."),
  title = "Case III"
)
p_H3

#Case IV


p_H4 <- plot_survival(
  t,
  surv_list = list(
    1- pweibull(t, 2.1, 2.1),
    1- pweibull(t, 2.1, 1.75)
  ),
  sample_names = c("Weibull", "Weibull."),
  title = "Case IV"
)
p_H4

#Case V

p_H5 <- plot_survival(
  t,
  surv_list = list(
    1- pweibull(t, 1.1, 1),
    1- pweibull(t, 0.9, 0.7)
  ),
  sample_names = c("Weibull", "Weibull."),
  title = "Case V"
)
p_H5

