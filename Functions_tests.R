#####################
##. gehan test
####################
gehan_test <- function(X, Y, statusX, statusY,
                       alternative = c("two.sided", "greater", "less")) {
  # Gehan's Generalized Wilcoxon Test (Lee & Wang 5.1.1)
  #
  # Args:
  #   X, Y     : numeric vectors of times for two groups
  #   statusX  : 1=event, 0=censored for group X
  #   statusY  : 1=event, 0=censored for group Y
  #   alternative : "two.sided" (default), "greater" (X > Y), "less" (X < Y)
  #
  # Returns: W, Var, Z, p.value, U
  
  alternative <- match.arg(alternative)
  
  n1 <- length(X)
  n2 <- length(Y)
  n <- n1 + n2
  
  # Combine
  time <- c(X, Y)
  status <- c(statusX, statusY)
  group <- c(rep("X", n1), rep("Y", n2))
  
  # --- Pairwise scoring matrix ---
  S <- matrix(0L, n, n)
  for(i in 1:n) {
    for(j in 1:n) {
      if(i == j) next
      if(time[i] > time[j] && status[j] == 1) {
        S[i,j] <- 1
      } else if(time[i] < time[j] && status[i] == 1) {
        S[i,j] <- -1
      } else if(time[i] == time[j]) {  # tie
        if(status[i] == 0 && status[j] == 1) {
          S[i,j] <- 1
        } else if(status[i] == 1 && status[j] == 0) {
          S[i,j] <- -1
        }
      }
    }
  }
  
  # --- Compute statistics ---
  U <- rowSums(S)
  W <- sum(U[group == "X"])
  VarW <- (n1 * n2) / (n * (n - 1)) * sum(U^2)
  Z <- W / sqrt(VarW)
  
  # p-value depending on alternative
  
  if (alternative == "two.sided") {
    pval <- 2 * (1 - pnorm(abs(Z)))
  } else if (alternative == "greater") {
    pval <- 1 - pnorm(Z)    # H1: X > Y
  } else {
    pval <- pnorm(Z)        # H1: X < Y
  }
  
  return(list(W = W, Var = VarW, Z = Z,
              p.value = pval, U = U))
}



#####################
##. cox_mantel test
####################



cox_mantel_test<- function(X, Y, statusX = NULL, statusY = NULL,
                           alternative = c("two.sided","greater","less"),
                           return_table = FALSE) {
  # Cox-Mantel test implemented as in Lee & Wang, Section 5.1.2 (eqs 5.1.5-5.1.6)
  # X, Y : numeric survival times for group1 and group2
  # statusX, statusY : optional vectors of 0/1 (1=event). If NULL, assume all events.
  # alternative : "two.sided" (default), "greater" (group1 > group2), or "less"
  # return_table : if TRUE, return per-failure-time table (m_i, n1_i, n2_i, r_i, A_i)
  #
  # Returns a list with U, I, C (stat), p.value (depending on alternative),
  # p.value.two.sided, and optionally the table.
  
  alternative <- match.arg(alternative)
  n1 <- length(X)
  n2 <- length(Y)
  if (is.null(statusX)) statusX <- rep(1L, n1)
  if (is.null(statusY)) statusY <- rep(1L, n2)
  if (length(statusX) != n1 || length(statusY) != n2)
    stop("statusX and statusY must match lengths of X and Y")
  
  # combine
  time <- c(X, Y)
  status <- c(as.integer(statusX), as.integer(statusY))
  group <- c(rep(1L, n1), rep(2L, n2))   # 1 = group X, 2 = group Y
  
  # distinct failure times (only where status==1), sorted ascending
  fail_times <- sort(unique(time[status == 1]))
  k <- length(fail_times)
  if (k == 0) stop("No failures (status==1) found in the pooled data.")
  
  # containers
  m <- numeric(k)
  n1_at <- numeric(k)
  n2_at <- numeric(k)
  r_at <- numeric(k)
  A <- numeric(k)
  
  for (i in seq_along(fail_times)) {
    t <- fail_times[i]
    m[i] <- sum(time == t & status == 1)         # multiplicity at t
    at_risk <- time >= t                         # risk set: time >= t
    r_at[i] <- sum(at_risk)
    n1_at[i] <- sum(at_risk & group == 1L)
    n2_at[i] <- r_at[i] - n1_at[i]
    # avoid division by zero (shouldn't happen if r_at[i] >= 1)
    A[i] <- if (r_at[i] > 0) n2_at[i] / r_at[i] else 0
  }
  
  # r2 = total failures in group 2 (sum of events in group 2)
  r2 <- sum(status == 1 & group == 2L)
  
  # U per eq (5.1.5)
  U <- r2 - sum(m * A)
  
  # I per eq (5.1.6): sum m_i * (r_i - m_i)/(r_i - 1) * A_i(1 - A_i)
  I_terms <- numeric(k)
  for (i in seq_len(k)) {
    ri <- r_at[i]
    if (ri <= 1) {
      I_terms[i] <- 0  # safeguard: contribution is effectively 0 or undefined; set 0
    } else {
      I_terms[i] <- m[i] * (ri - m[i])/(ri - 1) * (A[i] * (1 - A[i]))
    }
  }
  I <- sum(I_terms)
  
  if (I <= 0) {
    C <- NaN
  } else {
    C <- U / sqrt(I)
  }
  
  # p-values
  p_two <- NA_real_
  p_alt <- NA_real_
  if (!is.nan(C)) {
    p_two <- 2 * (1 - pnorm(abs(C)))
    if (alternative == "two.sided") p_alt <- p_two
    else if (alternative == "greater") p_alt <- 1 - pnorm(C)    # H1: group1 > group2
    else if (alternative == "less")    p_alt <- pnorm(C)        # H1: group1 < group2
  }
  
  out <- list(U = U, I = I, C = C, p.value = p_alt, p.value.two.sided = p_two)
  if (return_table) {
    tab <- data.frame(time = fail_times, m = m, n1 = n1_at, n2 = n2_at,
                      r = r_at, A = A, I_term = I_terms)
    out$table <- tab
  }
  return(out)
}



#####################
##. logrank test
####################



logrank_chisq <- function(X, Y, statusX, statusY) {
  # Chi-square version of logrank test (Lee & Wang, Eq. 5.1.10, 5.1.11)
  #
  # Args:
  #   X, Y       : survival times for groups 1 and 2
  #   statusX/Y  : 1 = event, 0 = censored
  #
  # Returns: O1, O2, E1, E2, X2, p.value
  
  # Combine data
  time   <- c(X, Y)
  status <- c(statusX, statusY)
  group  <- c(rep(1L, length(X)), rep(2L, length(Y)))
  df <- data.frame(time, status, group)
  
  # Unique event times
  event_times <- sort(unique(df$time[df$status == 1]))
  
  O1 <- sum(statusX)  # observed deaths in group 1
  O2 <- sum(statusY)  # observed deaths in group 2
  E1 <- 0
  E2 <- 0
  
  for (t in event_times) {
    # risk set at time t
    risk <- df$time >= t
    n1 <- sum(risk & df$group == 1)
    n2 <- sum(risk & df$group == 2)
    n  <- n1 + n2
    
    # deaths at t
    d  <- sum(df$time == t & df$status == 1)
    
    if (n > 0 && d > 0) {
      e1 <- d * (n1 / n)
      e2 <- d * (n2 / n)
      E1 <- E1 + e1
      E2 <- E2 + e2
    }
  }
  
  # Chi-square statistic
  chi_square <- (O1 - E1)^2 / E1 + (O2 - E2)^2 / E2
  pval <- 1 - pchisq(chi_square, df = 1)
  
  return(list(O1 = O1, O2 = O2,
              E1 = E1, E2 = E2,
              chi_square = chi_square, p.value = pval))
}






#####################
##. logrank Z test
####################



logrank_chisq_Z <- function(X, Y, statusX, statusY,
                            alternative = c("two.sided", "greater", "less")) {
  # Logrank Test (Chi-square + directional Z version)
  #
  # Args:
  #   X, Y       : survival times for groups 1 and 2
  #   statusX/Y  : 1 = event, 0 = censored
  #   alternative: "two.sided", "greater", "less"
  #
  # Returns: list with O1, O2, E1, E2, chi_square, Z, p.value
  
  alternative <- match.arg(alternative)
  
  if (length(X) != length(statusX))
    stop("X and statusX must have the same length")
  if (length(Y) != length(statusY))
    stop("Y and statusY must have the same length")
  
  # Combine data
  time   <- c(X, Y)
  status <- c(statusX, statusY)
  group  <- c(rep(1L, length(X)), rep(2L, length(Y)))
  df <- data.frame(time, status, group)
  
  # Unique event times
  event_times <- sort(unique(df$time[df$status == 1]))
  
  # Observed deaths
  O1 <- sum(statusX)
  O2 <- sum(statusY)
  E1 <- 0
  E2 <- 0
  
  # Compute expected deaths
  for (t in event_times) {
    risk <- df$time >= t
    n1 <- sum(risk & df$group == 1)
    n2 <- sum(risk & df$group == 2)
    n  <- n1 + n2
    d  <- sum(df$time == t & df$status == 1)
    
    if (n > 0 && d > 0) {
      e1 <- d * (n1 / n)
      e2 <- d * (n2 / n)
      E1 <- E1 + e1
      E2 <- E2 + e2
    }
  }
  
  # Chi-square statistic
  chi_square <- if (E1 > 0 && E2 > 0) {
    (O1 - E1)^2 / E1 + (O2 - E2)^2 / E2
  } else {
    NA
  }
  
  # Variance for Z (from 2x2 structure: O1+O2 = E1+E2 = total deaths)
  VarO1 <- if ((E1 + E2) > 0) (E1 * E2) / (E1 + E2) else NA
  Z <- if (!is.na(VarO1) && VarO1 > 0) (O1 - E1) / sqrt(VarO1) else NA
  
  # p-value depending on alternative
  pval <- NA
  if (!is.na(Z)) {
    if (alternative == "two.sided") {
      pval <- 2 * (1 - pnorm(abs(Z)))
    } else if (alternative == "greater") {
      pval <- 1 - pnorm(Z)   # group 1 worse survival
    } else if (alternative == "less") {
      pval <- pnorm(Z)       # group 1 better survival
    }
  }
  
  return(list(
    O1 = O1, E1 = E1,
    O2 = O2, E2 = E2,
    chi_square = chi_square,
    Z = Z,
    df = 1,
    p.value = pval,
    alternative = alternative
  ))
}


#####################
##. Novel Mann Whitney  test
####################


modified_mw_test <- function(X, Y, statusX, statusY,
                               alternative = c("two.sided", "greater", "less")) {
  
  alternative <- match.arg(alternative)
  
  # Counts
  n11 <- sum(statusX == 1)  # uncensored in X
  n22 <- sum(statusX == 0)  # censored in X
  m11 <- sum(statusY == 1)  # uncensored in Y
  m22 <- sum(statusY == 0)  # censored in Y
  
  # Combine data
  Z <- data.frame(
    Value  = c(X, Y),
    Status = c(statusX, statusY),
    Group  = c(rep("X", length(X)), rep("Y", length(Y)))
  )
  
  # Sort: uncensored first, then censored; within each sort by Value
  Z_sorted <- Z[order(-Z$Status, Z$Value), ]
  
  # Rank separately within uncensored/censored groups
  
  Z_sorted$Rank <- ave(Z_sorted$Value, Z_sorted$Status, FUN = rank)
  
  
  #Z_sorted$Rank <- ave(Z_sorted$Value, Z_sorted$Status,FUN = function(x) rank(x, ties.method = "first"))
  
  #print(Z_sorted)
  # Test statistic W = sum of ranks for Y
  W <- sum(Z_sorted$Rank[Z_sorted$Group == "Y"])
  
  # Mean and Variance
  Mean <- (m11 * (n11 + m11 + 1) / 2) + (m22 * (n22 + m22 + 1) / 2)
  Var  <- (n11 * m11 * (n11 + m11 + 1) / 12) + (n22 * m22 * (n22 + m22 + 1) / 12)
  
  # T statistic
  T_stat <- (W - Mean) / sqrt(Var)
  
  # p-value depending on alternative
  if (alternative == "two.sided") {
    pval <- 2 * (1 - pnorm(abs(T_stat)))
  } else if (alternative == "greater") {
    pval <- 1 - pnorm(T_stat)    # H1: Y > X
  } else {
    pval <- pnorm(T_stat)        # H1: Y < X
  }
  
  return(list(
    W = W,
    Mean = Mean,
    Variance = Var,
    T_stat = T_stat,
    p.value = pval
  ))
}






####Cases Function for data generation 


# Generate survival times from piecewise exponential hazards
rpexp <- function(n, breaks, rates) {
  # n: number of observations
  # breaks: time cut points (vector, increasing)
  # rates: hazard rates for each interval
  times <- numeric(n)
  
  for (i in 1:n) {
    u <- runif(1)
    t <- 0
    H <- 0
    for (j in seq_along(rates)) {
      if (j < length(rates)) {
        interval_len <- breaks[j+1] - breaks[j]
      } else {
        interval_len <- Inf  # last interval open-ended
      }
      
      # Cumulative hazard if we survive full interval
      deltaH <- rates[j] * interval_len
      if (-log(u) > H + deltaH) {
        H <- H + deltaH
        t <- breaks[j+1]
      } else {
        t <- breaks[j] + (-log(u) - H) / rates[j]
        break
      }
    }
    times[i] <- t
  }
  return(times)
}

# Wrapper to generate one dataset for a case
generate_dataset_cases <- function(n1, n2, case) {
  if (case == 1) {
    breaks <- c(0)       # no cutpoints, constant hazard
    h1 <- 1; h2 <- 2
    times1 <- rexp(n1, h1)
    times2 <- rexp(n2, h2)
  }
  
  if (case == 2) {
    breaks <- c(0, 0.4, 0.6)
    rates1 <- c(0.75, 3, 1)
    rates2 <- c(3, 0.75, 1)
    times1 <- rpexp(n1, breaks, rates1)
    times2 <- rpexp(n2, breaks, rates2)
  }
  
  if (case == 3) {
    breaks <- c(0, 0.5)
    rates1 <- c(2, 0.4)
    rates2 <- c(2, 4)
    times1 <- rpexp(n1, breaks, rates1)
    times2 <- rpexp(n2, breaks, rates2)
  }
  
  if (case == 4) {
    breaks <- c(0, 0.2, 0.6, 0.9)
    rates1 <- c(2, 0.75, 3, 1)
    rates2 <- c(2, 3, 0.75, 1)
    times1 <- rpexp(n1, breaks, rates1)
    times2 <- rpexp(n2, breaks, rates2)
  }
  
  list(times1, times2)
  # # Add censoring
  # censor1 <- runif(N1, 0, 2)
  # censor2 <- runif(N2, 0, 2)
  # 
  # obs_time1 <- pmin(times1, censor1)
  # obs_time2 <- pmin(times2, censor2)
  # status1 <- as.numeric(times1 <= censor1)
  # status2 <- as.numeric(times2 <= censor2)
  # 
  # data.frame(
  #   time = c(obs_time1, obs_time2),
  #   status = c(status1, status2),
  #   group = rep(c(1, 2), c(N1, N2))
  # )
}





###Survival plots Figure 1
#
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to plot survival curves with optional intersections
plot_survival <- function(t, surv_list, sample_names, title) {
  # Prepare data
  df <- data.frame(
    t = rep(t, length(surv_list)),
    S = unlist(surv_list),
    Sample = rep(sample_names, each=length(t))
  )
  
  # Wide form to compute differences
  df_wide <- df %>%
    pivot_wider(names_from = Sample, values_from = S)
  
  # Compute pairwise differences and find intersections
  # For simplicity, we only compare the first two samples
  col1 <- sample_names[1]
  col2 <- sample_names[2]
  
  df_wide <- df_wide %>%
    mutate(diff = .data[[col1]] - .data[[col2]])
  
  # Find indices where sign changes (approx intersections)
  cross_idx <- which(diff(df_wide$diff > 0) != 0)
  
  if(length(cross_idx) > 0){
    t_int <- df_wide$t[cross_idx]
    S_int <- df_wide[[col1]][cross_idx]
    df_int <- data.frame(
      t = t_int,
      S = S_int,
      label = paste0("(", round(t_int,2), ", ", round(S_int,2), ")")
    )
  } else {
    df_int <- data.frame(t=numeric(0), S=numeric(0), label=character(0))
  }
  
  # Plot
  p <- ggplot(df, aes(x=t, y=S, color=Sample, linetype=Sample)) +
    geom_line(size=1) +
    # Add intersection points if any
    # geom_point(
    #   data = df_int,
    #   aes(x=t, y=S),
    #   color="black",
    #   size=3,
    #   shape=21,
    #   fill="white",
    #   inherit.aes = FALSE
    # ) +
    # geom_text(
    #   data = df_int,
    #   aes(x=t, y=S, label=label),
    #   vjust=-1,
    #   size=3,
    #   inherit.aes = FALSE
    # ) +
    labs(title=title, x="Time", y="Survival Probabilities") +
    ylim(0, 1.05) +
    theme_minimal(base_size=14) +
    theme(legend.position="none")
  
  return(p)
}



