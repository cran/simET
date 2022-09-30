#===============================================================================
#Calculating meteorological parameters
#===============================================================================

#'@title Calculating solar declination
#'@note The solar declination actually varies throughout the day too but
#'    its variation is very small; thus, it is often ignored.
#'    Negative angles occur when the angle is below the equator plane,
#'    positive for above the equator.
#'@param td is the day of year.
#'@return A vector for solar declination (Radian)
#'@examples cal_solarDeclination(34)
#'@export
#'@references Teh CBS.Introduction to mathematical modeling of crop growth:
#'            How the equations are derived and assembled into a computer model.
#'            Brown Walker Press, 2006.

cal_solarDeclination<-function(td){
  delta<--0.4093*cos(2*pi*(td+10)/365)
  return(delta)
}

#'@title Calculating local solar time
#'@description Local solar time is different with local time.
#'@param td The day of year.
#'@param t is the local time.
#'@param gamma is the local longitude (Radian).
#'@param gamma_sm is the standard longitude (Radian).
#'@return A vector for local solar time(Hour)
#'@examples cal_localDolarTime(td=1,t=12,gamma=0.52,gamma_sm=2.09)
#'@export
#'@references Teh CBS.Introduction to mathematical modeling of crop growth:
#'            How the equations are derived and assembled into a computer model.
#'            Brown Walker Press, 2006.


cal_localDolarTime<-function(td,t,gamma,gamma_sm){
  B<-2*pi*(td-81)/364
  EoT<-9.87*sin(2*B)-7.53*cos(B)-1.5*sin(B)
  th<-t+((gamma_sm-gamma)/(pi/12))+(EoT/60)
  return(th)
}


#'@title Calculating hour angle
#'@param th is the local solar time.
#'@return A vector for hour angle (Radian)
#'@examples cal_hourAngle(12)
#'@export
#'@references Teh CBS.Introduction to mathematical modeling of crop growth:
#'            How the equations are derived and assembled into a computer model.
#'            Brown Walker Press, 2006.
cal_hourAngle<-function(th){
  hourangle<-(pi/12)*(th-12)
  return(hourangle)
}

#' @title Calculating solar inclination
#' @description A parameter used to determine the position of the sun relative
#'   to the observer (the other one is the angle from south).Conversion
#'   relationship with solar altitude angle:
#'   solar inclination=pi/2-solar altitude.
#'@param solar_declination is solar declination anger. It can be calculated from cal_solardeclination().
#'@param latitude is the latitude data (Radian).
#'@param hour_anger is hour anger. It can be calculated from cal_hourangle().
#'@return A vector for solar inclination (Radian)
#'@examples cal_solarInclination(solar_declination=-0.297,latitude=30,hour_anger=0)
#'@export
#'@references Teh CBS.Introduction to mathematical modeling of crop growth:
#'            How the equations are derived and assembled into a computer model.
#'            Brown Walker Press, 2006.

cal_solarInclination<-function(solar_declination,latitude,hour_anger){
  solarinclination<-acos(sin(solar_declination)*sin(latitude)+cos(solar_declination)*sin(latitude)*cos(hour_anger))
  return(solarinclination)
}

#'@title Calculating anger from south
#'@description A parameter used to determine the position of the sun relative
#'    to the observer (the other one is solar inclination).
#'@param latitude is the latitude data (Radian).
#'@param solar_altitude It can be calculated from pi/2-cal_solarinclination.
#'@param solar_declination is solar declination anger. It can be calculated from cal_solardeclination()
#'@return A vector for anger from south (Radian)
#'@details The minus and positive signs are taken before and after solar noon,
#'    receptively. The reason for having the positive-and-negative signs is
#'    merely an artificial convention so that we are able to distinguish between
#'    the sun lying westwards (posite angles and after solar noon) and eastwards
#'    (negative angles and before solar noon).
#'@examples cal_angerFromSouth(latitude=0.52,solar_altitude=-0.715,solar_declination=-0.2974005)
#'@export
#'@references Teh CBS.Introduction to mathematical modeling of crop growth:
#'            How the equations are derived and assembled into a computer model.
#'            Brown Walker Press, 2006.

cal_angerFromSouth<-function(latitude,solar_altitude,solar_declination){
  angerfromsouth<-acos((sin(latitude)*sin(solar_altitude)-sin(solar_declination))/(cos(latitude)*cos(solar_altitude)))
  return(angerfromsouth)
}


#'@title Calculating the local solar time for sunset/sunrise
#'@description Calculating the local solar time for sunset/sunrise.
#'@param solar_declination can be calculated by cal_solardeclination().
#'@param latitude is latitude data(Radian).
#'@return A vector for the local solar time for sunset/sunrise
#'@details Knowing the time of sunrise can calculate the time of sunset.
#'    sunrise_time=24-sunset_time. Day_length=2*(sunset_time-12).
#'@export
#'@references Teh CBS.Introduction to mathematical modeling of crop growth:
#'            How the equations are derived and assembled into a computer model.
#'            Brown Walker Press, 2006.

cal_sunsetTime<-function(solar_declination,latitude){
  tss<-12+(12/pi)*acos(-(sin(solar_declination)*sin(latitude))/(cos(solar_declination)*cos(latitude)))
  return(tss)
}



#----FAO56----

#'@title Calculating atmospheric pressure
#'@description The atmospheric pressure, P, is the pressure exerted by the
#'    weight of the earth's atmosphere.
#'@param elevation elevation above sea level (m)
#'@return A vector for atmospheric pressure (Kpa)
#'@details Assuming 20°C for a standard atmosphere. Evaporation at high altitudes
#'    is promoted due to low atmospheric pressure as expressed in thepsychrometric
#'    constant. The effect is, however, small and in the calculation procedures,
#'    the average value for a location is sufficient.
#'@examples cal_atmosphericPressure(100)
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_atmosphericPressure<-function(elevation){
  P<-101.3*((293-0.0065*elevation)/293)^5.26
  return(P)
}

#'@title Calculating psychrometric constant
#'@param atmospheric_pressure atmospheric pressure (kPa).
#'@return A vector for Psychrometric constant (kPa/degree Celsius)
#'@examples cal_psychrometriCconstant(100.1235)
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_psychrometriCconstant<-function(atmospheric_pressure){
  gamma<-0.665*10^(-3)*atmospheric_pressure
  return(gamma)
}

#'@title calculating the mean daily air temperature
#'@param Tmax the daily maximum.The temperature is given in degree Celsius or Fahrenhei.
#'@param Tmin the daily minimum.The temperature is given in degree Celsius,or Fahrenhei.
#'@return A vector for the mean daily air temperature
#'@details It is only employed in the FAO Penman-Monteith equation to calculate
#'    the slope of the saturation vapor pressure curves and the impact of mean
#'    air density as the effect of temperature variations on the value of the
#'    climatic parameter is small in these cases. For standardization, Tmean for
#'    24-hour periods is defined as the mean of mean of the daily maximum and
#'    minimum temperatures rather than as the average of hourly temperature
#'    measurements.
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_TemMean<-function(Tmax,Tmin){
  Tmean<-(Tmax+Tmin)/2
  return(Tmean)
}


#'@title Calculating saturation vapour pressure
#'@description saturation vapour pressure at the air temperature T.
#'@param Tem air temperature (degrees Celsius).
#'@return A vector for saturation vapour pressure at the air temperature T (kPa).
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_saturationVapourPressure<-function(Tem){
  saturation_vapour_pressure<-0.6108*exp((17.27*Tem)/(Tem+237.3))
  return(saturation_vapour_pressure)
}


#'@title Actual vapour pressure derived from dewpoint temperature
#'@description As the dewpint temperature is the temperature to which the air
#'    needs to be cooled to make the air saturated, the actual vapour pressure
#'    is the saturation vapour pressure at the dewpoint temperature.
#'@param Tdew dew point temperature(degrees Celsius).
#'@return A vector for actual vapour pressure
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_ActualVapourPressure_from_dewPoint<-function(Tdew){
  ActualVapourPressure<-cal_saturationVapourPressure(Tdew)
  return(ActualVapourPressure)
}

#'@title Actual vapour pressure (ea) derived from psychrometric data
#'@description The actual vapour pressure can be determined from the difference
#'    between the dry and wet bulb temperatures, the so-called wet bulb depression.
#'@param Twet,Tdry wet bulb depression, with Tdry the dry bulb and Twet the wet bulb temperature (degrees Celsius).
#'@param P is the atmospheric pressure (kPa).
#'@param type psychrometer type ("Asmann type","natural ventilated","non-ventilated").
#'@return A vector for Actual vapour pressure (ea)
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_ActualVapourPressure_from_psychrometricData<-function(Twet,Tdry,P,type){
  # if(type=="Asmann type"){
  #   a_psy<-0.000662
  # }else if(type=="natural ventilated"){
  #   a_psy<-0.000800
  # }else if(type=="non-ventilated"){
  #   a_psy<-0.001200
  # }
  a_psy<-ifelse(type=="Asmann type",a_psy<-0.000662,
                ifelse(type=="natural ventilated",a_psy<-0.000800,
                       ifelse(type=="non-ventilated",a_psy<-0.001200,
                              NA)))

  gamma_psy<-a_psy*P
  ActualVapourPressure<-cal_saturationVapourPressure(Twet)-gamma_psy*(Tdry-Twet)
}



#'@title Calculating mean saturation vapour pressure
#'@description Due to the non-linearity of the above equation, the mean saturation
#'    vapour pressure for a day, week, decade or month should be computed as the
#'    mean between the saturation vapour pressure at the mean daily maximum and
#'    minimum air temperatures for that period.
#'@details Using mean air temperature instead of daily minimum and maximum
#'    temperatures results in lower estimates for the mean saturation vapour
#'    pressure. The corresponding vapour pressure deficit (a parameter expressing
#'    the evaporating power of the atmosphere) will also be smaller and the result
#'    will be some underestimation of the reference crop evapotranspiration.
#'    Therefore,the mean sauration vapour pressure should be calculated as the
#'    mean between the saturation vapour pressure at both the daily maximum
#'    and minimum air temperature.
#'@param Tmax the daily maximum  air temperature(degrees Celsius).
#'@param Tmin the daily  minimum air temperature(degrees Celsius).
#'@return A vector for mean saturation vapour pressure (es)
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_meanSaturationVapourPressure<-function(Tmax,Tmin){
  mean_saturation_vapour_pressure<-(cal_saturationVapourPressure(Tmax)+cal_saturationVapourPressure(Tmin))/2
  return(mean_saturation_vapour_pressure)
}


#'@title Actual vapour pressure derived from RHmax and RHmin
#'@description The actual vapour pressure can also be calculated from the
#'    relative humidity. Depending on the availability of the humidity data,
#'    different equations should be used.
#'@param Tmax daily maximum temperature (kPa).
#'@param Tmin daily minimum temperature (KPa).
#'@param RHmax maximum relative humidity  \%.
#'@param RHmin minimum relative humidity  \%.
#'@details For periods of a week, ten days or a month, RHmax and RHmin are
#'    obtained by dividing the sum of the daily values by the number of days
#'    in that period.
#'@return A vector for actual vapour pressure
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_ActualVapourPressure_from_RHmaxAndRHmin<-function(Tmax,Tmin,RHmax,RHmin){
  ActualVapourPressure<-(cal_saturationVapourPressure(Tmin)*(RHmax/100)+cal_saturationVapourPressure(Tmax)*(RHmin/100))/2
  return(ActualVapourPressure)
}


#'@title Calculating actual vapour pressure derived from RHmax
#'@description When using equipment where errors in estimating RHmin can be large,
#'  or when RH data integrity are in doubt, then one should use only RHmax.
#'@param Tmin daily minimum temperature (degrees Celsius).
#'@param RHmax maximum relative humidity (\%).
#'@return A vector for actual vapour pressure
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_ActualVapourPressure_from_RHmax<-function(Tmin,RHmax){
  ActualVapourPressure<-cal_saturationVapourPressure(Tmin)*(RHmax/100)
  return(ActualVapourPressure)
}


#'@title Calculating actual vapour pressure derived from RHmean
#'@description In the absence of RH max and RHmin, it can be used to estimate
#'    actual vapour pressure.
#'@param RHmean mean relative humidity(\%).
#'@param Tmax daily maximum temperature (degrees Celsius).
#'@param Tmin daily minimum temperature (degrees Celsius).
#'@return A vector for actual vapour pressure
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_ActualVapourPressure_from_RHmean<-function(RHmean,Tmax,Tmin){
  ActualVapourPressure<-(RHmean/100)*((cal_saturationVapourPressure(Tmax)+cal_saturationVapourPressure(Tmin))/2)
  return(ActualVapourPressure)
}

#'@title Calculating actual vapour pressure for hourly time step
#'@param Thr is average hourly temperature (degrees Celsius).
#'@param RHhr is average hourly relative humidity [\%].
#'@return  A vector for average hourly actual vapour pressure [kPa].
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_ActualVapourPressure_for_hourly<-function(Thr,RHhr){
  ea<-cal_saturationVapourPressure(Thr)*(RHhr/100)
  return(ea)
}

#'@title Calculating slope of saturation vapour pressure curve
#'@details In the FAO Penman-Monteith equation, where it occurs in the numerator
#'    and denominator, the slope of the vapour pressure curve is calculated
#'    using mean air temperature.
#'@param Tem is air temperature (degrees Celsius).
#'@return A vector for slope of saturation vapour pressure curve at air temperature T
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_slopeOfSaturationVapourPressureCurve<-function(Tem){
  slopeOfSaturationVapourPressureCurve<-(4098*(0.6108*exp((17.27*Tem)/(Tem+237.3))))/(Tem+237.3)^2
  return(slopeOfSaturationVapourPressureCurve)
}

#'@title Calculating relative humidity
#'@description The relative humidity (RH) expresses the degree of saturation of
#'    the air as a ratio of the actual (ea) to the saturation (eo(T)) vapour
#'    pressure at the same temperature (T).
#'@param ea actual saturation vapour pressure. From cal_ActualVapourPressure_for_*
#'@param e0 saturation vapour pressure. From cal_saturationVapourPressure()
#'@details Relative humidity is the ratio between the amount of water the ambient
#'    air actually holds and the amount it could hold at the same temperature.
#'    It is dimensionless and is commonly given as a percentage. Although the
#'    actual vapour pressure might be relatively constant throughout the day,
#'    the relative humidity fluctuates between a maximum near sunrise and a
#'    minimum around early afternoon (Figure 12). The variation of the relative
#'    humidity is the result of the fact that the saturation vapour pressure is
#'    determined by the air temperature. As the temperature changes during the
#'     day, the relative humidity also changes substantially.
#'@export
#'@return A vector for relative humidity \%
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_relativeHumidity<-function(ea,e0){
  RH<-100*(ea/e0)
  return(RH)
}



#--Radiation

#'@title Calculating solar declination with FAO56 method
#'@param J is the number of the day in the year between 1 (1 January) and 365 or 366 (31 December)
#'@return A vector for solar declination(Radian)
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.


cal_solarDeclination_in_FAO<-function(J){
  solarDeclination<-0.409*sin((2*pi/365)*J-1.39)
  return(solarDeclination)
}

#'@title Calculating sunset hour angle
#'@param lat latitude (Radian), positive for the northern hemisphere and
#'    negative for the southern hemisphere.
#'@param solar_declination solar declination(Radina).
#'@return A vector for sunset Hour Angle(Radian)
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_sunsetHourAngle<-function(lat,solar_declination){
  sunsetHourAngle<-acos(-tan(lat)*tan(solar_declination))
  return(sunsetHourAngle)
}

#'@title Calculating inverse relative distance Earth-sun
#'@param J is the number of the day in the year between 1 (1 January) and 365
#'    or 366 (31 December).
#'@return A vector for inverse relative distance Earth-sun (Radian)
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_inverseRelativeDistance_Earth_sun<-function(J){
  inverseRelativeDistance_Earth_sun<-1+0.033*cos((2*pi/365)*J)
  return(inverseRelativeDistance_Earth_sun)
}


#'@title Calculating extraterrestrial radiation for daily periods
#'@description The extraterrestrial radiation, Ra, for each day of the year and
#'    for different latitudes can be estimated from the solar constant, the
#'    solar declination and the time of the year.
#'@param J is the number of the day in the year between 1 (1 January) and 365
#'    or 366 (31 December).
#'@param lat latitude (Radian), positive for the northern hemisphere and
#'    negative for the southern hemisphere.
#'@return A vector for extraterrestrial radiation for daily(MJ m-2 day-1)
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_extraterrestrialRadiation_for_daily<-function(J,lat){
  Gsc<-0.0820#MJ m-2 day-1
  dr<-cal_inverseRelativeDistance_Earth_sun(J)
  solar_declination<-cal_solarDeclination_in_FAO(J)
  ws<-cal_sunsetHourAngle(lat,solar_declination)
  Ra<-24*60*0.082*dr*(ws*sin(lat)*sin(solar_declination)+cos(lat)*cos(solar_declination)*sin(ws))/pi
  #24*60*0.082*Q2*(S2*SIN($A$5)*SIN(R2)+COS($A$5)*COS(R2)*SIN(S2))/3.14
  return(Ra)
}

#'@title Calculating extraterrestrial radiation for hourly or shorter periods
#'@param lat latitude (radian), positive for the northern hemisphere and negative
#'    for the southern hemisphere.
#'@param J is the number of the day in the year between 1 (1 January) and 365
#'    or 366 (31 December).
#'@param t standard clock time at the midpoint of the period (hour). For example
#'    for a period between 14.00 and 15.00 hours, t = 14.5.
#'@param lz longitude of the centre of the local time zone (degrees west of Greenwich).
#'    For example, Lz = 75, 90, 105 and 120° for the Eastern, Central, Rocky Mountain
#'    and Pacific time zones (United States) and Lz = 0° for Greenwich, 330°
#'    for Cairo (Egypt),and 255° for Bangkok (Thailand), radian.
#'@param t1 length of the calculation period (hour)
#'@param lm longitude of the measurement site (degrees west of Greenwich) radian.
#'@return A vector for extraterrestrial Radiation (MJ m-2 hour-1)
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_extraterrestrialRadiation_for_shorter<-function(lat,J,t,lz,lm,t1){
  Gsc<-0.0820
  dr<-cal_inverseRelativeDistance_Earth_sun(J)
  b<-2*pi*(J-81)/364
  Sc<-0.1645*sin(2*b)-0.1255*cos(b)-0.025*sin(b)#The seasonal correction for solar time
  w<-(pi/12)*((t+0.06667*(lz-lm)+Sc)-12)#The solar time angle at midpoint of the period
  w1<-w-(pi*t1)/24
  w2<-w+(pi*t1)/24
  solar_declination<-cal_solarDeclination_in_FAO(J)
  Ra<-(12*60/pi)*Gsc*dr*((w2-w1)*sin(lat)*sin(solar_declination)+cos(lat)*cos(solar_declination)*(sin(w2)-sin(w1)))
  return(Ra)
}


#'@title Calculating Daylight hours
#'@param sunsetHourAngle is the sunset hour angle in radians from cal_sunsetHourAngle().
#'@return A vector for day light Hours
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_daylightHours<-function(sunsetHourAngle){
  N<-(24/pi)*sunsetHourAngle
  return(N)
}

#'@title Calculating Solar radiation
#'@description If the solar radiation, Rs, is not measured, it can be calculated
#'    with the Angstrom formula, which relates solar radiation to extraterrestrial
#'    radiation and relative sunshine duration. This is a shortwave radiation.
#'@param as  regression constant, expressing the fraction of extraterrestrial
#'    radiation reaching the earth on overcast days (n = 0).Default is 0.25.
#'@param bs as+bs is fraction of extraterrestrial radiation reaching the earth
#'    on clear days (n = N). Default is 0.50.
#'@param n actual duration of sunshine [hour].
#'@param N maximum possible duration of sunshine or daylight hours [hour].from cal_daylightHours()
#'@param Ra  extraterrestrial radiation [MJ m-2 day-1]. From cal_extraterrestrialRadiation_for_daily()
#'@return A vector for solar or shortwave radiation [MJ m-2 day-1]
#'@note Rs is expressed in the above equation in MJ m-2 day-1. The corresponding
#'    equivalent evaporation in mm day-1 is obtained by multiplying Rs by 0.408
#'    (Equation 20). Depending on atmospheric conditions (humidity, dust) and
#'    solar declination (latitude and month), the Angstrom values as and bs will
#'    vary. Where no actual solar radiation data are available and no calibration
#'    has been carried out for improved as and bs parameters, the values
#'    as = 0.25 and bs = 0.50 are recommended.
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_solarRadiation<-function(as=0.25,bs=0.50,n,N,Ra){
  Rs<-(as+bs*(n/N))*Ra
}


#'@title Calculating Solar radiation from actual duration of sunshine
#'@param as  regression constant, expressing the fraction of extraterrestrial
#'    radiation reaching the earth on overcast days (n = 0).Default is 0.25.
#'@param bs as+bs is fraction of extraterrestrial radiation reaching the earth
#'    on clear days (n = N). Default is 0.50.
#'@param Na actual duration of sunshine [hour].
#'@param Latitude latitude (angert).
#'@param J is the number of the day in the year between 1 (1 January) and 365
#'    or 366 (31 December).
#'@export
#'@return A vector for solar radiation(MJ m-2 d-1)
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_Rs_from_Na<-function(as=0.25,bs=0.50,Na,Latitude,J){

  #计算Ra
  Lat=convert_angert_to_radian(Latitude)
  Ra=cal_extraterrestrialRadiation_for_daily(J,Lat)
  #计算Na
  solar_declination<-cal_solarDeclination_in_FAO(J)

  sunsetHourAngle<-cal_sunsetHourAngle(Lat,solar_declination)
  N<-cal_daylightHours(sunsetHourAngle)
  Rs<-cal_solarRadiation(as=0.25,bs=0.50,Na,N,Ra)
  return(Rs)
}



#'@title Calculating clear sky solar radiation with as and bs
#'@description The calculation of the clear-sky radiation, Rso, when n = N,
#'    is required for computing net longwave radiation.
#'@param as,bs as+bs fraction of extraterrestrial radiation reaching the earth
#'    on clear-sky days (n = N).
#'@param Ra extraterrestrial radiation [MJ m-2 day-1].
#'    From cal_extraterrestrialRadiation_for_daily()
#'@return A vector for clear-sky solar radiation [MJ m-2 day-1].
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_skySolarRadiation_withas_bs<-function(as,bs,Ra){
  Rso<-(as+bs)*Ra
  return(Rso)
}


#'@title Calculating clear sky solar radiation with elevation
#'@description The calculation of the clear-sky radiation, Rso, when n = N,
#'    is required for computing net longwave radiation.
#'@param z station elevation above sea level [m].
#'@param Ra extraterrestrial radiation [MJ m-2 day-1].
#'    From cal_extraterrestrialRadiation_for_daily()
#'@return A vector for clear-sky solar radiation [MJ m-2 day-1].
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_skySolarRadiation_withas_elevation<-function(z,Ra){
  Rso<-(0.75+2*10^(-5)*z)*Ra
  return(Rso)
}

#'@title Calculating net solar (shortware radiation) Rns
#'@description The net shortwave radiation resulting from the balance between
#'    incoming and reflected solar radiation.
#'@param alpha albedo or canopy reflection coefficient, which is 0.23 for the
#'    hypothetical grass reference crop [dimensionless].
#'@param Rs the incoming solar radiation [MJ m-2 day-1]. From cal_solarRadiation()
#'@return  A vector for net solar or shortwave radiation [MJ m-2 day-1].
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_netSolarRadiation<-function(alpha,Rs){
  Rns<-(1-alpha)*Rs
  return(Rns)
}

#'@title Calculating net longwave radiation Rnl
#'@param TKmax  maximum absolute temperature during the 24-hour period [K].
#'@param TKmin  minimum absolute temperature during the 24-hour period [K].
#'@param ea actual vapour pressure [kPa].
#'@param Rs measured or calculated solar radiation [MJ m-2 day-1].
#'    From cal_solarRadiation().
#'@param Rso calculated clear-sky radiation [MJ m-2 day-1].
#'    From cal_skySolarRadiation_withas_bs() or cal_skySolarRadiation_withas_elevation().
#'@note The rate of longwave energy emission is proportional to the absolute
#'    temperature of the surface raised to the fourth power. This relation is
#'    expressed quantitatively by the Stefan-Boltzmann law. The net energy flux
#'    leaving the earth's surface is, however, less than that emitted and given
#'    by the Stefan-Boltzmann law due to the absorption and downward radiation
#'    from the sky. Water vapour, clouds, carbon dioxide and dust are absorbers
#'    and emitters of longwave radiation. Their concentrations should be known
#'    when assessing the net outgoing flux. As humidity and cloudiness play an
#'    important role, the Stefan-Boltzmann law is corrected by these two factors
#'    when estimating the net outgoing flux of longwave radiation. It is thereby
#'    assumed that the concentrations of the other absorbers are constant. An
#'    average of the maximum air temperature to the fourth power and the minimum
#'    air temperature to the fourth power is commonly used in the Stefan-Boltzmann
#'    equation for 24-hour time steps. The term (0.34-0.14*sqrt(ea)) expresses
#'    the correction for air humidity, and will be smaller if the humidity increases.
#'    The effect of cloudiness is expressed by (1.35 Rs/Rso - 0.35). The term
#'    becomes smaller if the cloudiness increases and hence Rs decreases. The
#'    smaller the correction terms, the smaller the net outgoing flux of longwave
#'    radiation.  Note that the Rs/Rso term in Equation 39 must be limited so that
#'    Rs/Rso <= 1.0. Where measurements of incoming and outgoing short and longwave
#'    radiation during bright sunny and overcast hours are available, calibration
#'    of the coefficients in Equation 39 can be carried out.
#'@return A vector for net outgoing longwave radiation [MJ m-2 day-1]
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.


cal_netLongwaveRadiation<-function(TKmax,TKmin,ea,Rs,Rso){
  delta<-4.903*10^(-9)#Stefan-Boltzmann constant MJ K-4 m-2 day-1
  Rnl<-delta*((TKmax^4+TKmin^4)/2)*(0.34-0.14*sqrt(ea))*(1.35*(Rs/Rso)-0.35)
  return(Rnl)
}

#'@title Calculating net radiation Rn
#'@description The net radiation (Rn) is the difference between the incoming
#'    net shortwave radiation (Rns) and the outgoing net longwave radiation (Rnl).
#'@param Rns incoming net shortwave radiation. From cal_netSolarRadiation().
#'@param Rnl outgoing net longwave radiation. From cal_netLongwaveRadiation().
#'@return A vector for net radiation
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_netRadiation<-function(Rns,Rnl){
  Rn<-Rns-Rnl
  return(Rn)
}


#'@title Calculating soil heat flux (G) for general
#'@description Complex models are available to describe soil heat flux. Because
#'    soil heat flux is small compared to Rn, particularly when the surface is
#'    covered by vegetation and calculation time steps are 24 hours or longer,
#'    a simple calculation procedure is presented here for long time steps,
#'    based on the idea that the soil temperature follows air temperature.
#'@param cs soil heat capacity [MJ m-3 degrees Celsius-1].
#'@param T1 air temperature at time i [degrees Celsius].
#'@param T0 air temperature at time i-1 [degrees Celsius].
#'@param delta_t length of time interval [day].
#'@param delta_z effective soil depth [m].
#'@note Complex models are available to describe soil heat flux. Because soil
#'    heat flux is small compared to Rn, particularly when the surface is
#'    covered by vegetation and calculation time steps are 24 hours or longer,
#'    a simple calculation procedure is presented here for long time steps,
#'    based on the idea that the soil temperature follows air temperature.
#'@return A vector for soil heat flux [MJ m-2 day-1]
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.
#'


cal_soilHeatFlux_general<-function(cs,T1,T0,delta_t,delta_z){
  G<-cs*((T1+T0)/delta_t)*delta_z
  return(G)
}



#'@title Calculating soil heat flux(G) for day/ten-day periods
#'@description As the magnitude of the day or ten-day soil heat flux beneath
#'    the grass reference surface is relatively small, it may be ignored .
#'@return A value for 0
#'@export
#'@references FAO Irrigation and drainage paper 56 (P54)

cal_soilHeatFlux_day<-function(){
  G<-0
  return(G)
}

#'@title Calculating soil heat flux(G) for monthly periods
#'@description When assuming a constant soil heat capacity of 2.1 MJ m-3 °C-1
#'    and an appropriate soil depth, cal_soilHeatFlux_general can be used to
#'    derive G for monthly periods.
#'@param T1 air temperature at time i [degrees Celsius].
#'@param T0 air temperature at time i-1 [degrees Celsius].
#'@param Tmonth2 Is the mean air temperature of next month know?
#'@return A vector for soil heat flux [MJ m-2 day-1]
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_soilHeatFlux_monthly<-function(T1,T0,Tmonth2=TRUE){
  G<-ifelse(Tmonth2==TRUE,G<-0.07*(T1-T0),G<-0.14*(T1-T0))
  return(G)
}

#'@title Calculating soil heat flux(G) for hourly/shorter periods
#'@description For hourly (or shorter) calculations, G beneath a dense cover
#'    of grass does not correlate well with air temperature.
#'@param Rn net radiation.From cal_netRadiation().
#'@param periods "daylight" or "nighttime".
#'@note Where the soil is warming, the soil heat flux G is positive. The amount
#'    of energy required for this process is subtracted from Rn when estimating
#'    evapotranspiration.
#'@return A vector for soil heat flux [MJ m-2 day-1]
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_soilHeatFlux_hourly<-function(Rn,periods){
  Ghr<-ifelse(periods=="daylight",Ghr<-0.1*Rn,
         ifelse(periods=="nighttime",Ghr<-0.5*Rn,
                NA))
  return(Ghr)
}

