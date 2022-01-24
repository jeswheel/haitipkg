library(haitipkg)
library(tidyverse)

path <- system.file("extdata", "ve_decay_bymonth.csv", package = 'haitipkg')
ve_decay <- read.csv(path)

usethis::use_data(
  ve_decay,
  overwrite = TRUE
)
