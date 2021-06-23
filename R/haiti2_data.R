#' Clean haiti cases data for Model 2
#'
#' Generate a data frame of department level reported cases of cholera
#' in Haiti.
#'
#' @importFrom magrittr %>%
#' @return Data frame
#' @examples
#' haiti2_data()

haiti2_data <- function(){
  haiti <- haiti_case_data
  haiti <- haiti %>% dplyr::select(-report)
  haiti <- reshape2::melt(haiti, id.vars="date_sat_orig")
  colnames(haiti) <- c("year","department","cases")
  haiti$year <- lubridate::decimal_date(lubridate::ymd(haiti$year))
  haiti$department <- as.character(haiti$department)

  return(haiti)
}

