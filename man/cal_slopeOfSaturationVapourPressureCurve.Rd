% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_meteorlogy.R
\name{cal_slopeOfSaturationVapourPressureCurve}
\alias{cal_slopeOfSaturationVapourPressureCurve}
\title{Calculating slope of saturation vapour pressure curve}
\usage{
cal_slopeOfSaturationVapourPressureCurve(Tem)
}
\arguments{
\item{Tem}{is air temperature (degrees Celsius).}
}
\value{
A vector for slope of saturation vapour pressure curve at air temperature T
}
\description{
Calculating slope of saturation vapour pressure curve
}
\details{
In the FAO Penman-Monteith equation, where it occurs in the numerator
   and denominator, the slope of the vapour pressure curve is calculated
   using mean air temperature.
}
\references{
Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
   FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
   Organization of the United Nations, 1998.
}
