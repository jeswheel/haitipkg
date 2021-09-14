#' Make covariate table for Model 1
#'
#' Generate a data frame containing weekly values for periodic basis splines for the seasonal transmission term beta.
#'
#' @param tmin first observation time
#' @param tmax last observation time
#' @param byt number of steps between each time
#' @param nbasis number of basis spline terms
#' @param degree degrees of freedom for the basis splines
#' @param settings settings for the vaccination campaign scenario
#' @export

#### construct vaccination covariate table
covars <- function(tmin, tmax, byt = 1, nbasis = 6, degree = 6, per = 52.14, data, settings) {
  haiti_dat <- data
  ndept <- settings$nd
  nweeks <- settings$nw
  coverage_2dose <- settings$c2
  coverage_1dose <- settings$c1
  first_vac_t <- nrow(haiti_dat) + 4
  ve_scen <- settings$vescen

  ## deployment scenario 1 (fast national): one departmental campaign every 10 weeks for all depts
  tbasis <- seq(from=tmin,to=tmax,by=byt)
  covar_tab <- data.frame(cbind(time = tbasis,
                                pomp::periodic.bspline.basis(
                                  x = tbasis,
                                  nbasis = nbasis,
                                  degree = degree,
                                  period = per,
                                  names = "seas%d"
                                )))

  ## no vaccinations
  if (ve_scen == "ve_s0") {
    return (covar_tab)
  }

  time_check <- c()
  for(i in 1:ndept) {
    time_check <- c(time_check, rep(0, nweeks-1), rep(i, 1))
  } ## 10 week vacc campaigns (nweeks)

  ## number of vaccines per week by department
  pop_dept <- data.frame(ocv_order = 1:10,
                         dept = c("Centre", "Artibonite", "Ouest", "Nord_Ouest",
                                  "Nord", "Sud", "Nippes", "Nord_Est", "Sud_Est",
                                  "Grand_Anse"),
                         pop = c(746236, 1727524, 4029705, 728807, 1067177,
                                 774976, 342525, 393967, 632601, 468301)) %>%
    dplyr::mutate(num_vacc = (coverage_2dose+coverage_1dose)*pop/1) ## pulse vaccinees in last 1 week of campaign
  # browser()
  ## create dataframe with number vaccinated for each campaign
  vactab <- data.frame(time = tbasis, vac_tcheck = 0)
  vactab[which(vactab$time %in% first_vac_t:(first_vac_t+length(time_check)-1)),]$vac_tcheck <- time_check
  vactab2 <- dplyr::left_join(vactab, pop_dept %>%
                                dplyr::select(ocv_order, num_vacc), by = c("vac_tcheck"="ocv_order")) %>%
    dplyr::mutate(num_vacc = ifelse(is.na(num_vacc), 0, round(num_vacc)))

  ## ve decay after x weeks
  veDecay_mo <- ve_decay
  veDecay_adult <- dplyr::bind_rows(veDecay_mo, veDecay_mo, veDecay_mo, veDecay_mo) %>%
    dplyr::arrange(month) %>%
    dplyr::select(!!ve_scen) %>%
    unlist %>% unname
  ## adjust for lower VE in U5 population (U5 VE is 0.4688*adult VE; roughly 11% of the population is 0-4 years old according to UN World Population Prospects 2017)
  ## adjust ve for age (population VE = adult VE * (1-(1-0.4688)*proportion under-5))
  veDecay <- veDecay_adult * (1-(1-0.4688)*0.11)
  ## adjust for one-dose decay after 52 weeks
  veDecay[53:length(veDecay)] <- coverage_2dose/(coverage_2dose+coverage_1dose)*veDecay[53:length(veDecay)]

  ## add vaccine immunity decay, time shifted for each campaign
  decay_times <- vactab2 %>%
    dplyr::filter(vac_tcheck > 0) %>%
    dplyr::group_by(vac_tcheck) %>%
    dplyr::filter(time == max(time)) %>%
    dplyr::select(-num_vacc) %>%
    dplyr::ungroup()
  for(i in 1:ndept) {
    decay_start <- decay_times %>%
      dplyr::filter(vac_tcheck == i) %>%
      dplyr::select(time) %>% unlist %>% unname
    decay_df <- tibble::as_tibble(data.frame(time = (seq_along(veDecay)+decay_start))) %>%
      dplyr::mutate(!!paste0("ve_d", i) := veDecay)
    vactab2 <- dplyr::left_join(vactab2, decay_df, by = c("time"))
  }

  vactab3 <- vactab2 %>%
    dplyr::mutate_at(dplyr::vars(contains("ve_")), list(~ifelse(is.na(.), 0, .)))


  covartab.all <- covar_tab %>%
    dplyr::mutate(time_endfc = ifelse(time < 232, NA, time)) %>%
    dplyr::mutate(time_fc = ifelse(time < max(haiti_dat), NA, time)) %>%
    dplyr::mutate(time_end = ifelse(time <= max(haiti_dat$week), time, NA))
  covartab.fc <- covartab.all %>%
    dplyr::filter(!is.na(time_fc)) %>%
    dplyr::select(time_fc, contains("seas"))
  covartab.all2 <- covartab.all %>%
    dplyr::select(time, contains("seas")) %>%
    dplyr::full_join(vactab3, by = c("time"))

  return(covartab.all2)
}
