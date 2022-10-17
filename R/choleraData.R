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

#' Lee et al (2020) Model 1 Epidemic Parameters
#'
#' This dataset contains parameters similar to those used by Lee et al (2020)
#' to project cholera incidence in Haiti using the epidemic version of Model 1.
#'
#' @format The dataset has 152 rows and 23 columns. Each column corresponds to
#' a parameter value for the epidemic version of Model 1, except for the column
#' `parid` which identifies the parameter set.
#'
#' @source The data was created by running the source code of Lee et al (2020).
#'    Because of differences in package versions and other possible numeric
#'    considerations, the set of parameters are not the exact same as those
#'    used by Lee et al (2020), but the same code was used to generate these
#'    data, resulting in very similar sets of parameters. After generating the
#'    data using the Lee et al (2020) souce code, some parameters with low log
#'    likelihoods were removed, as well as other outlying values that were
#'    inconsistent with those shown in the Lee et al (2020) Figures S8 and S9.
#'
"h1LeeStartsEpi"


#' Lee et al (2020) Model 1 Endemic Parameters
#'
#' This dataset contains parameters similar to those used by Lee et al (2020)
#' to project cholera incidence in Haiti using the endemic version of Model 1.
#'
#' @format The dataset has 147 rows and 24 columns. Each column corresponds to
#' a parameter value for the epidemic version of Model 1, except for the column
#' `parid` which identifies the parameter set.
#'
#' @source The data was created by running the source code of Lee et al (2020).
#'    Because of differences in package versions and other possible numeric
#'    considerations, the set of parameters are not the exact same as those
#'    used by Lee et al (2020), but the same code was used to generate these
#'    data, resulting in very similar sets of parameters. After generating the
#'    data using the Lee et al (2020) souce code, some parameters with low log
#'    likelihoods were removed, as well as other outlying values that were
#'    inconsistent with those shown in the Lee et al (2020) Figures S8 and S9.
#'
"h1LeeStartsEnd"
