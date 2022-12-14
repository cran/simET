% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_soil.R
\name{cal_DeepPercolation}
\alias{cal_DeepPercolation}
\title{Calculating Deep percolation}
\usage{
cal_DeepPercolation(Wa, Wfc, a, b, t)
}
\arguments{
\item{Wa}{actual soil water storage in the root zone (mm)}

\item{Wfc}{soil water storage to maximum root depth (Zr ) at field capacity (mm)}

\item{a}{A water storage value comprised between WFc and Wa at saturation.}

\item{b}{b<-0.0173 for soils draining quickly.Otherwise b>-0.0173.}

\item{t}{time after an irrigation or rain that produced a storage above field capacity (days)}
}
\value{
A vector for deep percolation(mm/day).
}
\description{
Calculating Deep percolation
}
\references{
Liu Y, Pereira L S, Fernando R M. Fluxes through the bottom boundary of the root zone in silty soils:
 Parametric approaches to estimate groundwater contribution and percolation[J].
 Agricultural Water Management, 2006, 84(1):27-40.
}
