% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_meteorlogy.R
\name{cal_soilHeatFlux_general}
\alias{cal_soilHeatFlux_general}
\title{Calculating soil heat flux (G) for general}
\usage{
cal_soilHeatFlux_general(cs, T1, T0, delta_t, delta_z)
}
\arguments{
\item{cs}{soil heat capacity [MJ m-3 degrees Celsius-1].}

\item{T1}{air temperature at time i [degrees Celsius].}

\item{T0}{air temperature at time i-1 [degrees Celsius].}

\item{delta_t}{length of time interval [day].}

\item{delta_z}{effective soil depth [m].}
}
\value{
A vector for soil heat flux [MJ m-2 day-1]
}
\description{
Complex models are available to describe soil heat flux. Because
   soil heat flux is small compared to Rn, particularly when the surface is
   covered by vegetation and calculation time steps are 24 hours or longer,
   a simple calculation procedure is presented here for long time steps,
   based on the idea that the soil temperature follows air temperature.
}
\note{
Complex models are available to describe soil heat flux. Because soil
   heat flux is small compared to Rn, particularly when the surface is
   covered by vegetation and calculation time steps are 24 hours or longer,
   a simple calculation procedure is presented here for long time steps,
   based on the idea that the soil temperature follows air temperature.
}
\references{
Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
   FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
   Organization of the United Nations, 1998.
}
