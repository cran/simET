% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_meteorlogy.R
\name{cal_skySolarRadiation_withas_elevation}
\alias{cal_skySolarRadiation_withas_elevation}
\title{Calculating clear sky solar radiation with elevation}
\usage{
cal_skySolarRadiation_withas_elevation(z, Ra)
}
\arguments{
\item{z}{station elevation above sea level [m].}

\item{Ra}{extraterrestrial radiation [MJ m-2 day-1].
From cal_extraterrestrialRadiation_for_daily()}
}
\value{
A vector for clear-sky solar radiation [MJ m-2 day-1].
}
\description{
The calculation of the clear-sky radiation, Rso, when n = N,
   is required for computing net longwave radiation.
}
\references{
Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
   FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
   Organization of the United Nations, 1998.
}
