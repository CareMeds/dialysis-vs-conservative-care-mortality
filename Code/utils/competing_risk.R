# Extract state occupation probabilities 
Pprobtrans <- function(bh12, bh13, bh23) {
  N <- length(bh12)
  out_P11 <- rep(NA_real_, N)
  out_P12 <- rep(NA_real_, N)
  out_P13 <- rep(NA_real_, N)
  P12_1 <- 0
  P13_1 <- 0
  P11_1 <- 1
  for (i in 1:N) {
    P13_2 <- P11_1*bh13[i] + P12_1*bh23[i] + P13_1
    P11_2 <- P11_1*(1 - bh12[i] - bh13[i])
    P12_2 <- P11_1*bh12[i] + P12_1*(1 - bh23[i])
    P11_1 <- out_P11[i] <- P11_2
    P12_1 <- out_P12[i] <- P12_2
    P13_1 <- out_P13[i] <- P13_2
  }
  df <- data.frame(P11 = out_P11, P12 = out_P12, P13 = out_P13)
  return(df)
}

# Calculate state probabilities using weighted competing risk model
state_probabilities <- function(dt, horizon){
  # transition matrix
  tmat <- mstate::trans.illdeath(c("Decision","KRT","Death"))
  
  # set status and time at horizon
  dt[, KRT_event := fifelse(time2event_KRT_inf <= horizon, event_KRT_inf, 0)]
  dt[, KRT_time := fifelse(time2event_KRT_inf <= horizon, time2event_KRT_inf, horizon)]
  dt[, death_event := fifelse(time2event_death_inf <= horizon, event_death_inf, 0)]
  dt[, death_time := fifelse(time2event_death_inf <= horizon, time2event_death_inf, horizon)]
  
  # prepare data for competing risk model
  msdia <- mstate::msprep(
    time = c(NA, "KRT_time", "death_time"),
    status = c(NA, "KRT_event", "death_event"),
    data = dt,
    trans = tmat,
    keep = "sw_IPTW"
  )
  
  # fit competing risk model
  cox_dia <- survival::coxph(
    Surv(Tstart, Tstop, status) ~ strata(trans),
    data = msdia,
    weights = msdia$sw_IPTW,
    method = "breslow"
  )
  
  # Manually calculate state occupation probabilities (%)
  # extract baseline hazards for all transitions
  bh <- survival::basehaz(cox_dia, centered = F)
  bh_11 <- bh[bh$strata == "trans=1", ]
  bh_12 <- bh[bh$strata == "trans=2", ]
  bh_13 <- bh[bh$strata == "trans=3", ]
  
  # ensure that baseline hazard is extracted at uniform times
  alltimes <- 0:(horizon + 1)
  
  # fill the NAs of the cumulative hazard with the previous value
  bh_11_allt <- merge(bh_11, data.frame(time = alltimes), all = T) |>
    tidyr::fill(hazard, .direction = "down") 
  bh_12_allt <- merge(bh_12, data.frame(time = alltimes), all = T) |>
    tidyr::fill(hazard, .direction = "down") 
  bh_13_allt <- merge(bh_13, data.frame(time = alltimes), all = T) |>
    tidyr::fill(hazard, .direction = "down")
  
  # if there is a remaining NA at the iniatial times, that should be a 0
  bh_11_allt$hazard[is.na(bh_11_allt$hazard)] <- 0
  bh_12_allt$hazard[is.na(bh_12_allt$hazard)] <- 0
  bh_13_allt$hazard[is.na(bh_13_allt$hazard)] <- 0
  
  # extract hazard from cumulative hazard
  bh_11_allt$haz <- diff(c(0, bh_11_allt$hazard))
  bh_12_allt$haz <- diff(c(0, bh_12_allt$hazard))
  bh_13_allt$haz <- diff(c(0, bh_13_allt$hazard))
  
  # obtain state probabilities
  state_prob <- Pprobtrans(bh12 = bh_11_allt$haz,
                           bh13 = bh_12_allt$haz, 
                           bh23 = bh_13_allt$haz)
  # P11 = probability of remaining event-free 
  # P12 = probability of transitioning from decision to KRT
  # P13 = probability of transitioning from decision to death state
  
  # Average time spent in the decison state in the next two years
  RMST_11 <- sum(state_prob$P11 * diff(c(alltimes, horizon)))/30.5
  # Average time spent in dialysis in the next two years
  RMST_12 <- sum(state_prob$P12 * diff(c(alltimes, horizon)))/30.5
  # Average time spent alive in the next two years
  RMST_13 <- sum((1-state_prob$P13) * diff(c(alltimes, horizon)))/30.5
  
  return(list(times = alltimes, 
              state_prob = state_prob,
              RMST = list(RMST_11 = RMST_11,
                          RMST_12 = RMST_12,
                          RMST_13 = RMST_13)))  
}