library(tidyverse)
library(haitipkg)

path <- system.file("extdata", "rainfall.csv", package = 'haitipkg')
haitiRainfall <- readr::read_csv(path)

usethis::use_data(
  haitiRainfall,
  overwrite = TRUE
)
