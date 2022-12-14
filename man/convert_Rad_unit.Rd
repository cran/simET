% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.R
\name{convert_Rad_unit}
\alias{convert_Rad_unit}
\title{Convert radiation unit}
\usage{
convert_Rad_unit(rad, type)
}
\arguments{
\item{rad}{Radiation data need to be converted from one unit to another unit.}

\item{type}{Used to specify how to convert.}
}
\value{
A vector for radiation converted unit
}
\description{
Type has the following types:
 MJ_m2_day_to_J_cm2_day;
 MJ_m2_day_to_cal_cm2_day;
 MJ_m2_day_to_W_m2;
 cal_cm2_day_to_MJ_m2_day
 cal_cm2_day_to_J_cm2_day
 cal_cm2_day_to_W_m2
 cal_cm2_day_to_mm_day
 W_m2_to_MJ_m2_day
 W_m2_to_J_cm2_day
 W_m2_to_cal_cm2_day
 W_m2_to_mm_day
 mm_day_to_MJ_m2_day
 mm_day_to_J_cm2_day
 mm_day_to_cal_cm2_day
 mm_day_to_W_m2
}
\references{
Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
   FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
   Organization of the United Nations, 1998.
}
