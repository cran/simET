#'@title calculating reference evapotranspiration from Penman-Monteith method
#'@description The FAO Penman-Monteith method is maintained as the sole standard
#'    method for the computation of ETo from meteorological data.
#'@param delta  slope vapour pressure curve (kPa &deg;C). From cal_slopeOfSaturationVapourPressureCurve()
#'@param Rn net Radiation at the crop surface [MJ m-2 day-1]. From cal_netRadiation()
#'@param G soil heat flux density [MJ m-2 day-1].
#'@param gamma psychrometric constant (kPa &deg;C).
#'@param Tem air temperature at 2 m height [&deg;C].
#'@param u2 wind speed at 2 m height [m s-1].
#'@param es  saturation vapour pressure [kPa].
#'@param ea actual vapour pressure [kPa].
#'@return A vector for reference evapotranspiration [mm day-1].
#'@export
#'@note Ten-day or monthly time step :
#'
#'  Notwithstanding the non-linearity in the Penman-Monteith equation and some weather
#'  parameter methods, mean ten-day or monthly weather data can be used to compute the mean
#'  ten-day or monthly values for the reference evapotranspiration. The value of the reference
#'  evapotranspiration calculated with mean monthly weather data is indeed very similar to the
#'  average of the daily ETo values calculated with daily average weather data for that month.
#'
#'  When the soil is warming (spring) or cooling (autumn), the soil heat flux (G) for monthly
#'  periods may become significant relative to the mean monthly Rn. In these cases G cannot be
#'  ignored and its value should be determined from the mean monthly air temperatures of the
#'  previous and next month.
#'
#'  Daily time step:
#'
#'  Calculation of ETo with the Penman-Monteith equation on 24-hour time scales will generally
#'  provide accurate results.
#'
#'  As the magnitude of daily soil heat flux (G) beneath the reference grass surface is relatively
#'  small, it may be ignored for 24-hour time steps.
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_ET0_from_PM<-function(delta,Rn,G,gamma,Tem,u2,es,ea){
  # ET0<-((0.408*slopVapourPressureCurve*(netRadiation-soilHeartFluxDensity)+(psychrometricConstant*900*windSpeed*(saturationVapourPressure-actualVapourPressure)))/(airTemperature+273))/(slopVapourPressureCurve+slopVapourPressureCurve*(1+0.34*windSpeed))
  ET0<-(0.408*delta*(Rn-G)+gamma*(900/(Tem+273))*u2*(es-ea))/(delta+gamma*(1+0.34*u2))
  #(0.408*J2*(Z2-AA2)+K2*900*E2*P2/(C2+273))/(J2+K2*(1+0.34*E2))
  #0.408是单位的转换因子
  return(ET0)
}

#'@title Calculating reference evapotranspiration from Penman-Monteith for
#'    daily
#'@description Based on lat, z, J, Tmax, Tmin, n, RHmax, RHmin, windSpeed parameters,
#'    reference evapotranspiration was calculated by Penman-Monteith.
#'@param Latitude latitude (radian), positive for the northern hemisphere and
#'    negative for the southern hemisphere.
#'@param Altitude station elevation above sea level [m].
#'@param J is the number of the day in the year between 1 (1 January) and 365 or
#'    366 (31 December).
#'@param Tmax daily maximum air temperature (degrees Celsius).
#'@param Tmin daily minimum air temperature (degrees Celsius).
#'@param Rs Solar radiation [MJ m-2 d-1].
#'@param RHmean daily mean relative humidity \%.
#'@param Wind wind speed at 2 m height [m s-1].
#'@export
#'@return A vector for reference evapotranspiration (mm/day)
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.
#'@examples
#'  library(simET)
#'  data("FIalfalfa")
#'  names(FIalfalfa)
#'  Result_data<- dplyr::mutate(FIalfalfa,
#'                       ET0=cal_ET0_from_PM_for_daily(Latitude=Latitude,
#'                                                    Altitude=Altitude,
#'                                                    J=Julian,
#'                                                    Tmax=Tmax,
#'                                                    Tmin=Tmin,
#'                                                    Rs=Rs,
#'                                                    RHmean=RHmean,
#'                                                    Wind=Wind))
#'  names(Result_data)



cal_ET0_from_PM_for_daily<-function(Latitude,Altitude,J,Tmax,Tmin,Rs,RHmean,Wind){
  Tmean=(Tmax+Tmin)/2
  P=cal_atmosphericPressure(Altitude)
  Delta=cal_slopeOfSaturationVapourPressureCurve(Tmean)
  gamma=cal_psychrometriCconstant(P)
  es=cal_meanSaturationVapourPressure(Tmax ,Tmin)
  ea=cal_ActualVapourPressure_from_RHmean(RHmean,Tmax,Tmin)
  Deficit=es-ea
  # dr=cal_inverseRelativeDistance_Earth_sun(J)
  # Solar_D=cal_solarDeclination_in_FAO(J)
  Lat=convert_angert_to_radian(Latitude)
  # ws=cal_sunsetHourAngle(Lat,Solar_D)
  Ra=cal_extraterrestrialRadiation_for_daily(J,Lat)
  # Nmax=cal_daylightHours(ws)
  # Rs=cal_solarRadiation(as=0.25,bs=0.5,n=Na,N=Nmax,Ra=Ra)
  Rso=cal_skySolarRadiation_withas_elevation(z=Altitude,Ra=Ra)
  Rns=cal_netSolarRadiation(alpha=0.23,Rs=Rs)
  TKmax=convert_degreesCelsius_to_Fahrenheit(Tmax)
  TKmin=convert_degreesCelsius_to_Fahrenheit(Tmin)
  Rnl=cal_netLongwaveRadiation(TKmax,TKmin,ea,Rs,Rso)
  Rn=cal_netRadiation(Rns,Rnl)
  G=0
  ET0=cal_ET0_from_PM(Delta,Rn,G,gamma,Tmean,Wind,es,ea)
}

#'@title Calculating reference evapotranspiration from Penman-Monteith method
#'    for hourly time step
#'@details In areas where substantial changes in wind speed, dewpoint or
#'    cloudiness occur during the day, calculation of the ETo equation using
#'    hourly time steps is generally better than using 24-hour calculation
#'    time steps. Such weather changes can cause 24-hour means to misrepresent
#'    evaporative power of the environment during parts of the day and may
#'    introduce error into the calculations. However, under most conditions,
#'    application of the FAO Penman-Monteith equation with 24-hour data produces
#'    accurate results.
#'@param slopVapourPressureCurve saturation slope vapour pressure curve at Thr [kPa &deg;C].
#'@param netRadiation net radiation at the grass surface [MJ m-2 hour-1].
#'@param soilHeatFlux soil heat flux density [MJ m-2 hour-1].
#'@param psychrometricConstant psychrometric constant [kPa &deg;C].
#'@param meanHourlyTem mean hourly air temperature [&deg;C].
#'@param windSpeed  average hourly wind speed [m s-1].
#'@param saturationVapourPressure saturation vapour pressure at air temperature Thr [kPa].
#'@param actualVapourPressure  average hourly actual vapour pressure [kPa].
#'@return A vector for reference evapotranspiration [mm hour-1].
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.
#'@note With the advent of electronic, automated weather stations, weather
#'    data are increasingly reported for hourly or shorter periods. Therefore,
#'    in situations where calculations are computerized, the FAO Penman-Monteith
#'    equation can be applied on an hourly basis with good results. When applying
#'    the FAO Penman-Monteith equation on an hourly or shorter time scale, the
#'    equation and some of the procedures for calculating meteorological data
#'    should be adjusted for the smaller time step.
#'
#'  For the calculation of radiation parameters, see P74-75

cal_ET0_from_PM_for_hourly<-function(slopVapourPressureCurve,netRadiation,soilHeatFlux,psychrometricConstant,meanHourlyTem,windSpeed,saturationVapourPressure,actualVapourPressure){
  ET0<-(0.408* slopVapourPressureCurve*(netRadiation-soilHeatFlux)+ psychrometricConstant*(37/(meanHourlyTem+273))*windSpeed*(saturationVapourPressure-actualVapourPressure))/(slopVapourPressureCurve+psychrometricConstant*(1+0.34*windSpeed))
  return(ET0)
}


