#' Clean haiti cases data for Model 1
#'
#' Generate a data frame of the reported cases of cholera in Haiti from October 2010 to
#' January 2019, aggregated over all departments.
#'
#' Code adapted from Model 1 Team at JHU Bloomberg School of Public Health
#'
#' @importFrom magrittr %>%
#' @return A data frame
#' @examples
#' haiti1_agg_data()
#' @export

haiti1_agg_data <- function(){
  allDat <- haitiCholera
  ## must switch back to numeric from NA because sum() will ignore other values in row if we try to use na.rm = TRUE
  allDat[allDat$date_saturday == '2016-10-01', "Artibonite"] <- 0
  allDat[allDat$date_saturday == '2017-11-11', 'Artibonite'] <- 0
  allDat[allDat$date_saturday == '2017-01-07', 'Ouest'] <- 0
  allDat <- allDat %>%
    dplyr::mutate(date_sat_orig = date_saturday)

  splitDate <- strsplit(allDat$date_sat_orig, "-")
  data.table::setattr(splitDate[[1]], 'names', c("year", "month", "day"))
  dateDf <- tibble::as_tibble(as.data.frame(do.call(rbind, splitDate))) %>%
    dplyr::mutate(month = as.character(month)) %>%
    dplyr::mutate(day = as.character(day)) %>%
    dplyr::mutate(year = as.character(year)) %>%
    dplyr::mutate(month = ifelse(nchar(month) == 1, paste0("0", month), month)) %>%
    dplyr::mutate(day = ifelse(nchar(day) == 1, paste0("0", day), day)) %>%
    dplyr::mutate(date_sat = as.Date(paste(year, month, day, sep = "-"), origin = "1900-01-01"))

  fullDateVec <- data.frame(date_sat = seq(min(dateDf$date_sat), max(dateDf$date_sat), by=7))

  cleanDat <- allDat %>%
    dplyr::mutate(date_sat = as.Date(dateDf$date_sat, origin = "1900-01-01")) %>%
    dplyr::select(-date_sat_orig) %>%
    dplyr::full_join(fullDateVec, by = c("date_sat")) %>%
    dplyr::arrange(date_sat) %>%
    dplyr::mutate(week = seq_along(date_sat)) %>%
    dplyr::rename(Grand_Anse = Grand.Anse,
                  Sud_Est = Sud.Est) %>%
    tidyr::gather(department, cases, Artibonite:Sud_Est)

  aggDat <- cleanDat %>%
    dplyr::group_by(week) %>% ## "day" or "week"
    dplyr::summarise(cases = sum(cases))

  return(aggDat)
}
