#' Clean haiti cases data for Model 2
#'
#' Generate a data frame of department level reported cases of cholera
#' in Haiti.
#'
#' @importFrom magrittr %>%
#' @return Data frame
#' @examples
#' haiti2_data()
#' @export

haiti2_data <- function(){
  haiti <- haitiCholera
  haiti <- haiti %>% dplyr::select(-report)

  dates <- c("2017-04-22","2018-01-27","2018-09-22","2018-12-01")
  r <- rbind(rep(NA,10),rep(NA,10),rep(NA,10),rep(NA,10))

  final <- cbind(dates,r) %>% as.data.frame()
  colnames(final) <- colnames(haiti)

  final <- rbind(haiti, final)
  final <- dplyr::arrange(final, date_saturday)

  haiti <- reshape2::melt(final, id.vars="date_saturday")
  colnames(haiti) <- c("year","department","cases")
  haiti$year <- lubridate::decimal_date(lubridate::ymd(haiti$year))
  haiti$department <- as.character(haiti$department)
  haiti$cases <- as.numeric(haiti$cases)

  return(haiti)
}

