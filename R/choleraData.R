#' Haiti Cholera Data
#'
#' This dataset contains observed case numbers for each departement in Haiti
#' from October, 2010, to January, 2019. This dataset is used by each of the
#' models considered by Lee et al, and our re-analysis.
#'
#' @format The dataset contains 426 rows and 12 columns:
#' \describe{
#'    \item{date_saturday}{Date stored as year, month, day.}
#'    \item{report}{TODO: probably remove this column.}
#'    \item{Artibonite}{Reported weekly cases of Cholera in Artibonite.}
#'    ...
#' }
#'
#' This dataset was created by loading the .csv file saved as
#' inst/extdata/haiti-data-from-2010-10-to-2019-01.csv. After loading the data,
#' three datapoints were removed because it was determined that these points
#' were outliers:
#'
#' \itemize{
#'    \item 2016-10-01, Artibonite.
#'    \item 2017-11-11, Artibonite.
#'    \item 2017-01-07, Ouest.
#' }
"haitiCholera"