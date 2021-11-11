## code to prepare `MODEL3POMPDATA` dataset goes here

library(tidyverse)

MODEL3_INPUT_PARAMETERS <- yaml::read_yaml("https://raw.githubusercontent.com/ionides/haiti/main/model3/input/haiti-data/input_parameters.yaml?token=AQ356PKKXO27YIXW4X6WB3LBS2WS6")
MODEL3_CASES <- readr::read_csv("https://raw.githubusercontent.com/ionides/haiti/main/model3/input/haiti-data/fromAzman/cases.csv?token=AQ356PJQTKMQMRK7GL5MKTTBS2WVK")
MODEL3_RAIN <- readr::read_csv("https://raw.githubusercontent.com/ionides/haiti/main/model3/input/haiti-data/fromAzman/rainfall.csv?token=AQ356PPCCFFDQMHFHA7GQFTBS2WX4")
MODEL3_PARAMS <- read.csv('https://raw.githubusercontent.com/ionides/haiti/main/model3/output/best_artibonite_params.csv?token=AQ356POMRAQP4HFCLVH7ZMDBS2W2I')

MODEL3_CASES[MODEL3_CASES$date == '2016-10-01', 'Artibonite'] <- NA
MODEL3_CASES[MODEL3_CASES$date == '2017-11-11', 'Artibonite'] <- NA
MODEL3_CASES[MODEL3_CASES$date == '2017-01-07', 'Ouest'] <- NA

usethis::use_data(
  MODEL3_INPUT_PARAMETERS,
  MODEL3_CASES,
  MODEL3_RAIN,
  MODEL3_PARAMS,
  overwrite = TRUE
)
