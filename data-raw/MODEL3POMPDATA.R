## code to prepare `MODEL3POMPDATA` dataset goes here

library(magrittr)

MODEL3_INPUT_PARAMETERS <- yaml::read_yaml("https://raw.githubusercontent.com/ionides/haiti/main/model3/input/haiti-data/input_parameters.yaml?token=AQ356PMPCSG54WKPJDZZ6I3A64K2C")
MODEL3_CASES <- readr::read_csv("https://raw.githubusercontent.com/ionides/haiti/main/model3/input/haiti-data/fromAzman/cases_corrected.csv?token=AQ356PNOJAOF4MGB2BIWCK3A64K7A")
MODEL3_RAIN <- readr::read_csv("https://raw.githubusercontent.com/ionides/haiti/main/model3/input/haiti-data/fromAzman/rainfall.csv?token=AQ356PM3DCWJXDUWFXGA5ETA64LLM")
MODEL3_PARAMS <- read.csv('https://raw.githubusercontent.com/ionides/haiti/main/model3/output/best_artibonite_params.csv?token=AQ356PMHAND3J5L5OSGYLR3A64QJA')

usethis::use_data(
  MODEL3_INPUT_PARAMETERS,
  MODEL3_CASES,
  MODEL3_RAIN,
  MODEL3_PARAMS,
  overwrite = TRUE
)



# usethis::use_data(MODEL3POMPDATA, overwrite = TRUE)
