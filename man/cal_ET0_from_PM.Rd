% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_ReferenceEvapotranspiration.R
\name{cal_ET0_from_PM}
\alias{cal_ET0_from_PM}
\title{calculating reference evapotranspiration from Penman-Monteith method}
\usage{
cal_ET0_from_PM(delta, Rn, G, gamma, Tem, u2, es, ea)
}
\arguments{
\item{delta}{slope vapour pressure curve (kPa &deg;C). From cal_slopeOfSaturationVapourPressureCurve()}

\item{Rn}{net Radiation at the crop surface [MJ m-2 day-1]. From cal_netRadiation()}

\item{G}{soil heat flux density [MJ m-2 day-1].}

\item{gamma}{psychrometric constant (kPa &deg;C).}

\item{Tem}{air temperature at 2 m height [&deg;C].}

\item{u2}{wind speed at 2 m height [m s-1].}

\item{es}{saturation vapour pressure [kPa].}

\item{ea}{actual vapour pressure [kPa].}
}
\value{
A vector for reference evapotranspiration [mm day-1].
}
\description{
The FAO Penman-Monteith method is maintained as the sole standard
   method for the computation of ETo from meteorological data.
}
\note{
Ten-day or monthly time step :

 Notwithstanding the non-linearity in the Penman-Monteith equation and some weather
 parameter methods, mean ten-day or monthly weather data can be used to compute the mean
 ten-day or monthly values for the reference evapotranspiration. The value of the reference
 evapotranspiration calculated with mean monthly weather data is indeed very similar to the
 average of the daily ETo values calculated with daily average weather data for that month.

 When the soil is warming (spring) or cooling (autumn), the soil heat flux (G) for monthly
 periods may become significant relative to the mean monthly Rn. In these cases G cannot be
 ignored and its value should be determined from the mean monthly air temperatures of the
 previous and next month.

 Daily time step:

 Calculation of ETo with the Penman-Monteith equation on 24-hour time scales will generally
 provide accurate results.

 As the magnitude of daily soil heat flux (G) beneath the reference grass surface is relatively
 small, it may be ignored for 24-hour time steps.
}
\references{
Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
   FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
   Organization of the United Nations, 1998.
}
