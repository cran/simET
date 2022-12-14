% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Estimate.R
\name{estimate_Rs_from_airTemDiff}
\alias{estimate_Rs_from_airTemDiff}
\title{Estimating solar radiation data derived from air temperature differences}
\usage{
estimate_Rs_from_airTemDiff(Ra, Tmax, Tmin, locations)
}
\arguments{
\item{Ra}{extraterrestrial radiation [MJ m-2 d-1].}

\item{Tmax}{maximum air temperature.}

\item{Tmin}{minimum air temperature.}

\item{locations}{The adjustment coefficient kRs is empirical and differs for
interior' or 'coastal' regions.}
}
\value{
A vector for solar radiation
}
\description{
The difference between the maximum and minimum air temperature
   is related to the degree of cloud cover in a location. Clear-sky conditions
   result in high temperatures during the day (Tmax) because the atmosphere
   is transparent to the incoming solar radiation and in low temperatures
   during the night (Tmin) because less outgoing longwave radiation is absorbed
   by the atmosphere. On the other hand, in overcast conditions, Tmax is
   relatively smaller because a significant part of the incoming solar radiation
   never reaches the earth's surface and is absorbed and reflected by the clouds.
   Similarly, Tmin will be relatively higher as the cloud cover acts as a
   blanket and decreases the net outgoing longwave radiation. Therefore, the
   difference between the maximum and minimum air temperature (Tmax - Tmin)
   can be used as an indicator of the fraction of extraterrestrial radiation
   that reaches the earth's surface. This principle has been utilized by
   Hargreaves and Samani to develop estimates of ETo using only air temperature
   data.
}
\note{
The temperature difference method is recommended for locations where
   it is not appropriate to import radiation data from a regional station,
   either because homogeneous climate conditions do not occur, or because
   data for the region are lacking. For island conditions, the methodology
   of Equation 50 is not appropriate due to moderating effects of the
   surrounding water body. Caution is required when daily computations of
   ETo are needed. The advice given for Equation 49 fully applies. It is
   recommended that daily estimates of ETo that are based on estimated Rs be
   summed or averaged over a several-day period, such as a week, decade or
   month to reduce prediction error.
}
\references{
Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
   FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
   Organization of the United Nations, 1998.
}
