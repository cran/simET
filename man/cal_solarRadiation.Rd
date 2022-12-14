% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_meteorlogy.R
\name{cal_solarRadiation}
\alias{cal_solarRadiation}
\title{Calculating Solar radiation}
\usage{
cal_solarRadiation(as = 0.25, bs = 0.5, n, N, Ra)
}
\arguments{
\item{as}{regression constant, expressing the fraction of extraterrestrial
radiation reaching the earth on overcast days (n = 0).Default is 0.25.}

\item{bs}{as+bs is fraction of extraterrestrial radiation reaching the earth
on clear days (n = N). Default is 0.50.}

\item{n}{actual duration of sunshine [hour].}

\item{N}{maximum possible duration of sunshine or daylight hours [hour].from cal_daylightHours()}

\item{Ra}{extraterrestrial radiation [MJ m-2 day-1]. From cal_extraterrestrialRadiation_for_daily()}
}
\value{
A vector for solar or shortwave radiation [MJ m-2 day-1]
}
\description{
If the solar radiation, Rs, is not measured, it can be calculated
   with the Angstrom formula, which relates solar radiation to extraterrestrial
   radiation and relative sunshine duration. This is a shortwave radiation.
}
\note{
Rs is expressed in the above equation in MJ m-2 day-1. The corresponding
   equivalent evaporation in mm day-1 is obtained by multiplying Rs by 0.408
   (Equation 20). Depending on atmospheric conditions (humidity, dust) and
   solar declination (latitude and month), the Angstrom values as and bs will
   vary. Where no actual solar radiation data are available and no calibration
   has been carried out for improved as and bs parameters, the values
   as = 0.25 and bs = 0.50 are recommended.
}
\references{
Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
   FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
   Organization of the United Nations, 1998.
}
