#' Make covariate table for Model 1
#'
#' Generate a data frame containing weekly values for periodic basis splines for the seasonal transmission term beta and vaccination efficacy values
#'
#' NOTE: # represents the code of a different vaccination deployment strategy
#' 0 = no vaccinations,
#' 1 = fast national,
#' 2 = 2-department,
#' 3 = slow national,
#' 4 = 3-department,
#' 25 = fast national high coverage, and all others represent combinations
#'
#' @param scen_code Code for vaccination scenario
#' @return A list of codes specific to the given vaccination scenario
#' @export

vac_scen <- function(scen_code) {
  scencode <- scen_code
  fc_nd = NA; fc_nw = NA; fc_c2 = NA; fc_c1 = NA; fc_vescen = NA

  ## grouped by deployment strategy
  if (scencode %in% paste0("id", seq(2, 34, by = 4))){ ## 2 dept
    fc_nd = 2; fc_nw = 52
  } else if (scencode %in% paste0("id", seq(4, 36, by = 4))){ ## 3 dept
    fc_nd = 3; fc_nw = 33
  } else if (scencode %in% paste0("id", seq(3, 35, by = 4))){ ## slow national
    fc_nd = 10; fc_nw = 26
  } else if (scencode %in% paste0("id", seq(1, 33, by = 4))){ ## fast national
    fc_nd = 10; fc_nw = 10
  } else {
    fc_nd = 1; fc_nw = 26
  }

  ## grouped by coverage level
  if (scencode %in% paste0("id", 1:12)){
    fc_c2 = 0.7; fc_c1 = 0.1
  } else if (scencode %in% paste0("id", 13:24)){
    fc_c2 = 0.4; fc_c1 = 0.2
  } else if (scencode %in% paste0("id", 25:36)){
    fc_c2 = 0.95; fc_c1 = 0.0167
  } else {
    fc_c2 = 0; fc_c1 = 0
  }

  ## grouped by ve scenario
  if (scencode %in% paste0("id", c(1:4, 13:16, 25:28))){
    fc_vescen = "ve_s1"
  } else if (scencode %in% paste0("id", c(5:8, 17:20, 29:32))){
    fc_vescen = "ve_s2"
  } else if (scencode %in% paste0("id", c(9:12, 21:24, 33:36))){
    fc_vescen = "ve_s3"
  } else {
    fc_vescen = "ve_s0"
  }

  return(list(scode = scencode, nd = fc_nd, nw = fc_nw,
              c2 = fc_c2, c1 = fc_c1, vescen = fc_vescen))
}
