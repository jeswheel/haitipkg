% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/haiti3_dep_data.R
\name{haiti3_dep_data}
\alias{haiti3_dep_data}
\title{Function used to get department data for individual departments.}
\usage{
haiti3_dep_data(departement = "Artibonite", start_time = "2014-03-01")
}
\arguments{
\item{departement}{Name of department to fit POMP model to}

\item{start_time}{Time of which to start the time series. In the
original paper, the authors used a start time of 2014-03-01 for all
departements except for Ouest, which started with a time of 2017-06-10.
This means that there wasn't a national-coupled model fit to the
data until 2017-06-01, which means that there was less than 2 years
of data to fit the model. Still, the authors fit the remaining departemental
models to data from 2014-03-01. This added argument allows us to
fit the Ouest model earlier, or fit the remaining departement models
at 2017-06-01}
}
\description{
Function used to get department data for individual departments.
}
