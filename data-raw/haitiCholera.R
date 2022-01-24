library(haitipkg)
library(tidyverse)

path <- system.file("extdata", "haiti-data-from-2010-10-to-2019-01.csv", package = 'haitipkg')
haitiCholera <- read.csv(path)

haitiCholera[haitiCholera$date_saturday == '2016-10-01', "Artibonite"] <- NA
haitiCholera[haitiCholera$date == '2017-11-11', 'Artibonite'] <- NA
haitiCholera[haitiCholera$date == '2017-01-07', 'Ouest'] <- NA

usethis::use_data(
  haitiCholera,
  overwrite = TRUE
)
