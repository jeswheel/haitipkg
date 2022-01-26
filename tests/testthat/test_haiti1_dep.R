library(pomp)
library(haitipkg)
library(dplyr)
library(testthat)

#### Build Models ##############################################################
#### Goal: test an assortment of possible model builds/configurations

## default model: Artibonite, no vaccinations
def_mod <- haiti1_dep()

## Sud, no vaccinations
sud_novac <- haiti1_dep(dept = "Sud")

## Sud, 3 department vaccination campaign (Sud Est is NOT included)
sud_vac <- haiti1_dep(dept = "Sud", vacscen = "id4")

## Centre, 3 department vaccination campaign (Centre is included)
cen_vac <- haiti1_dep(dept = "Centre", vacscen = "id4")

## Nord, national vaccination campaign
nord_vac <- haiti1_dep(dept = "Nord", vacscen = "id5")

#### Tests #####################################################################

#### Test 1: check that the models have the appropriate data for their dept
test_that("check data", {
  ## get random week
  week <- sample(1:430, 1)
  ## get dept data for that week
  haiti_dat <- haitiCholera[week, ]

  ## default model
  expect_equal(def_mod@data[week], haiti_dat[, "Artibonite"])
  ## Sud, no vac
  expect_equal(sud_novac@data[week], haiti_dat[, "Sud"])
  ## Sud, vac
  expect_equal(sud_vac@data[week], haiti_dat[, "Sud"])
  ## Centre, vac
  expect_equal(cen_vac@data[week], haiti_dat[, "Centre"])
  ## Nord, vac
  expect_equal(nord_vac@data[week], haiti_dat[, "Nord"])
})

#### Test 2: check that models have the appropriate parameters
####         this test also will implicitly test whether the given dept
####         receives vaccinations because additional parameters are included
####         in the model if the campaign covers the dept
test_that("check parameters", {
  art_params <- MODEL1_INPUT_PARAMETERS$dep_params[["Artibonite"]] %>%
    unlist()
  sud_params <- MODEL1_INPUT_PARAMETERS$dep_params[["Sud"]] %>%
    unlist()

  cen_params <- MODEL1_INPUT_PARAMETERS$dep_params[["Centre"]] %>%
    unlist()
  par_names <- names(cen_params)
  cen_params <- c(cen_params, rep(0.0, 5))
  names(cen_params) <- c(par_names, "S1_0", "E1_0", "I1_0", "A1_0", "R1_0")
  nord_params <- MODEL1_INPUT_PARAMETERS$dep_params[["Nord"]] %>%
    unlist()
  par_names <- names(nord_params)
  nord_params <- c(nord_params, rep(0.0, 5))
  names(nord_params) <- c(par_names, "S1_0", "E1_0", "I1_0", "A1_0", "R1_0")

  ## default model
  expect_identical(def_mod@params, art_params)
  ## Sud, no vac
  expect_identical(sud_novac@params, sud_params)
  ## Sud, vac
  expect_identical(sud_novac@params, sud_params)
  ## Centre, vac
  expect_identical(cen_vac@params, cen_params)
  ## Nord, vac
  expect_identical(nord_vac@params, nord_params)
})
