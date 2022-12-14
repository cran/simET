% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_meteorlogy.R
\name{cal_extraterrestrialRadiation_for_daily}
\alias{cal_extraterrestrialRadiation_for_daily}
\title{Calculating extraterrestrial radiation for daily periods}
\usage{
cal_extraterrestrialRadiation_for_daily(J, lat)
}
\arguments{
\item{J}{is the number of the day in the year between 1 (1 January) and 365
or 366 (31 December).}

\item{lat}{latitude (Radian), positive for the northern hemisphere and
negative for the southern hemisphere.}
}
\value{
A vector for extraterrestrial radiation for daily(MJ m-2 day-1)
}
\description{
The extraterrestrial radiation, Ra, for each day of the year and
   for different latitudes can be estimated from the solar constant, the
   solar declination and the time of the year.
}
\references{
Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
   FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
   Organization of the United Nations, 1998.
}
