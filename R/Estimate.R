
#'@title Estimating missing humidity data
#'@description Where humidity data are lacking or are of questionable quality,
#'    an estimate of actual vapour pressure, ea, can be obtained by assuming
#'    that dewpoint temperature (Tdew) is near the daily minimum temperature (Tmin).
#'    This statement implicitly assumes that at sunrise, when the air temperature
#'    is close to Tmin, that the air is nearly saturated with water vapour and
#'    the relative humidity is nearly 100 \%.
#'@param Tmin the minimum tem daily.
#'@note The relationship Tdew near Tmin holds for locations where the cover
#'    crop of the station is well watered. However, particularly for arid regions,
#'    the air might not be saturated when its temperature is at its minimum.
#'    Hence, Tmin might be greater than Tdew and a further calibration may be
#'    required to estimate dewpoint temperatures. In these situations, "Tmin" in
#'    the above equation may be better approximated by subtracting 2-3 degrees
#'    Celsius from Tmin. Appropriate correction procedures are given in Annex 6.
#'    In humid and subhumid climates, Tmin and Tdew measured in early morning
#'    may be less than Tdew measured during the daytime because of condensation
#'    of dew during the night.  After sunrise, evaporation of the dew will once
#'    again humidify the air and will increase the value measured for Tdew
#'    during the daytime. This phenomenon is demonstrated in Figure 5.4 of Annex 5.
#'    However, it is standard practice in 24-hour calculations of ETo to use
#'    Tdew measured or calculated durin early morning. The estimate for ea from
#'    Tmin should be checked. When the prediction by Equation 48 is validated
#'    for a region, it can be used for daily estimates of ea.
#'@export
#'@return A vector for humidity
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

estimate_ea<-function(Tmin){
  ea<-0.611*exp((17.27*Tmin)/(Tmin+237.3))
  return(ea)
}


#'@title Estimating solar radiation data derived from air temperature differences
#'@description The difference between the maximum and minimum air temperature
#'    is related to the degree of cloud cover in a location. Clear-sky conditions
#'    result in high temperatures during the day (Tmax) because the atmosphere
#'    is transparent to the incoming solar radiation and in low temperatures
#'    during the night (Tmin) because less outgoing longwave radiation is absorbed
#'    by the atmosphere. On the other hand, in overcast conditions, Tmax is
#'    relatively smaller because a significant part of the incoming solar radiation
#'    never reaches the earth's surface and is absorbed and reflected by the clouds.
#'    Similarly, Tmin will be relatively higher as the cloud cover acts as a
#'    blanket and decreases the net outgoing longwave radiation. Therefore, the
#'    difference between the maximum and minimum air temperature (Tmax - Tmin)
#'    can be used as an indicator of the fraction of extraterrestrial radiation
#'    that reaches the earth's surface. This principle has been utilized by
#'    Hargreaves and Samani to develop estimates of ETo using only air temperature
#'    data.
#'
#'@param Ra  extraterrestrial radiation [MJ m-2 d-1].
#'@param Tmax  maximum air temperature.
#'@param Tmin  minimum air temperature.
#'@param locations The adjustment coefficient kRs is empirical and differs for
#'    interior' or 'coastal' regions.
#'
#'@note The temperature difference method is recommended for locations where
#'    it is not appropriate to import radiation data from a regional station,
#'    either because homogeneous climate conditions do not occur, or because
#'    data for the region are lacking. For island conditions, the methodology
#'    of Equation 50 is not appropriate due to moderating effects of the
#'    surrounding water body. Caution is required when daily computations of
#'    ETo are needed. The advice given for Equation 49 fully applies. It is
#'    recommended that daily estimates of ETo that are based on estimated Rs be
#'    summed or averaged over a several-day period, such as a week, decade or
#'    month to reduce prediction error.
#'@return A vector for solar radiation
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

estimate_Rs_from_airTemDiff<-function(Ra,Tmax,Tmin,locations){
  # if(locations=="interior"){
  #   k_Rs<-0.16
  # }else if (locations=="coastal"){
  #   k_Rs<-0.19
  # }
  k_Rs<-ifelse(locations=="interior",k_Rs<-0.16,
               ifelse(locations=="coastal",k_Rs<-0.19,
                      NA))
  Rs<-k_Rs*sqrt(Tmax-Tmin)*Ra
}

#'@title Estimating solar radiation for island locations
#'@description For island locations, where the land mass has a width
#'    perpendicular to the coastline of 20 km or less, the air masses influencing
#'    the atmospheric conditions are dominated by the adjacent water body in all
#'    directions. The temperature method is not appropriate for this situation.
#'    Where radiation data from another location on the island are not available,
#'    a first estimate of the monthly solar average can be obtained from the
#'    empirical relation.
#'@param Ra  extraterrestrial radiation [MJ m-2 day-1].
#'@param b empirical constant, equal to 4 MJ m-2 day-1.
#'@note This relationship is only applicable for low altitudes (from 0 to 100 m).
#'    The empirical constant represents the fact that in island locations some
#'    clouds are usually present, thus making the mean solar radiation 4 MJ m-2 day-1
#'    below the nearly clear sky envelope (0.7 Ra). Local adjustment of the
#'    empirical constant may improve the estimation. The method is only appropriate
#'    for monthly calculations. The constant relation between Rs and Ra does not
#'    yield accurate daily estimates.
#'@export
#'@return A vector for solar radiation
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

estimate_Rs_for_islandLocations<-function(Ra,b=4){
  Rs<-0.7*Ra-b
  return(Rs)
}

#'@title Estimating ET0 with Tmax and Tmin
#'@description When solar radiation data, relative humidity data and/or wind speed
#'    data are missing, they should be estimated using the procedures presented
#'    in this section.  As an alternative, ETo can be estimated using the
#'    Hargreaves ETo equation.
#'@param Tmean mean temperature.
#'@param Tmax max temperature.
#'@param Tmin min temperature.
#'@param Ra extraterrestrial radiation [mm day-1].
#'@return A vector for reference evapotranspiration (mm day-1).
#'@export
#'@note  Units for both ETo and Ra in Equation 52 are mm day-1.  Equation 52
#'    should be verified in each new region by comparing with estimates by
#'    the FAO Penman-Monteith equation (Equation 6) at weather stations where
#'    solar radiation, air temperature, humidity, and wind speed are measured.
#'    If necessary, Equation 52 can be calibrated on a monthly or annual basis
#'    by determining empirical coefficients where ETo = a + b ETo Eq.52, where
#'    the Eq. 52 subscript refers to ETo predicted using Equation 52.  The
#'    coefficients a and b can be determined by regression analyses or by visual
#'    fitting.  In general, estimating solar radiation, vapor pressure and wind
#'    speed as described in Equations 48 to 51 and Table 4 and then utilizing
#'    these estimates in Equation 6 (the FAO Penman-Monteith equation) will
#'    provide somewhat more accurate estimates as compared to estimating ETo
#'    directly using Equation 52. This is due to the ability of the estimation
#'    equations to incorporate general climatic characteristics such as high or
#'    low wind speed or high or low relative humidity into the ETo estimate made
#'    using Equation 6. Equation 52 has a tendency to underpredict under high
#'    wind conditions (u2 > 3 m/s) and to overpredict under conditions of high
#'    relative humidity.

estimate_ET0_with_TmaxAndTmin<-function(Tmean,Tmax,Tmin,Ra){
  ET0<-0.0023*(Tmean+17.8)*(Tmax-Tmin)^0.5*Ra
  return(ET0)
}

#'@title Estimate LAI for alfalfa
#'@param hc is the vegetation height in meter. (in meter)
#'@export
#'@return A vector for leaf are index of alfalfa
#'@references Zhao C , Feng Z , Chen G . Soil water balance simulation of alfalfa (Medicago sativa L.)
#'  in the semiarid Chinese Loess Plateau[J]. Agricultural Water Management, 2004, 69(2):0-114.

estimate_LAI_for_alfalfa<-function(hc){
  LAI<-5.5+1.5*log(hc)
  return(LAI)
}

#'@title linear interpolation for vector
#'@description Linear interpolation is performed by using the values on both
#'    sides of the missing values.
#'@export
#'@return A interpolated vector
#'@param DataVector data vector.Note that the starting value of vector needs
#'    to be no missing value.


linear_interpolation<-function(DataVector){
  numNa<-c()
  for (i in 2:length(DataVector)) {
    if(!is.na(DataVector[i-1])&is.na(DataVector[i])&!is.na(DataVector[i+1])){
      pre0<-DataVector[i-1]
      nex0<-DataVector[i+1]
      if(pre0>nex0){
        DataVector[i]<-pre0-(abs(nex0-pre0)/2)
      }else{
        DataVector[i]<-pre0+(abs(nex0-pre0)/2)
      }
    }else if(!is.na(DataVector[i-1])&is.na(DataVector[i])&is.na(DataVector[i+1])){
      pre1<-DataVector[i-1]
      numNa<-c(numNa,i)
    }else if(is.na(DataVector[i-1])&is.na(DataVector[i])&is.na(DataVector[i+1])){
      numNa<-c(numNa,i)
    }else if(is.na(DataVector[i-1])&is.na(DataVector[i])&!is.na(DataVector[i+1])){
      nex1<-DataVector[i+1]
      if(pre1<nex1){
        numNa<-c(numNa,i)
        Aver1<-abs(nex1-pre1)/(length(numNa)+1)
        ValueNa<-pre1+Aver1
        for(f in numNa){
          DataVector[f]<-ValueNa
          ValueNa<-ValueNa+Aver1
        }
        numNa<-NULL
        ValueNa<-NULL
      }else{
        numNa<-c(numNa,i)
        Aver1<-abs(nex1-pre1)/(length(numNa)+1)
        ValueNa<-pre1-Aver1
        for(f in numNa){
          DataVector[f]<-ValueNa
          ValueNa<-ValueNa-Aver1
        }
        nex1<-NA
        pre1<-NA
        Aver1<-NA
        numNa<-c()
        ValueNa<-NA
      }
    }
  }
  return(DataVector)
}



#'@title Calculating the goodness-of-fit indicators between measured and simulated values
#'@param Sim The simualtion value of model.
#'@param Obs The observed value.
#'@export
#'@return A vector for the goodness-of-fit indicators
estimate_goodnessOfFit<-function(Sim,Obs){
  DF<-data.frame(Sim,Obs)
  DF<-tidyr::drop_na(DF)
  RMSE<-with(DF,sqrt(sum((Obs-Sim)^2)/nrow(DF)))
  R2<-with(DF,(sum((Obs-mean(Obs))*(Sim-mean(Sim)))/(sum((Obs-mean(Obs))^2)^0.5*sum((Sim-mean(Sim))^2)^0.5))^2)
  b0<-with(DF,sum(Sim*Obs)/sum(Obs^2))
  RRMSE<-with(DF,sqrt((sum(Obs-Sim)^2)/nrow(DF))/mean(Obs))
  RSR<-with(DF,sqrt((sum(Obs-Sim)^2)/nrow(DF))/sqrt(sum((Obs-mean(Obs))^2)))
  MRE<-with(DF,(100/nrow(DF))*sum(abs((Obs-Sim)/Obs)))
  ME<-with(DF,sum(Obs-Sim)/nrow(DF))
  MAE<-with(DF,sum(abs(Obs-Sim))/nrow(DF))
  PBIAS<-with(DF,100*(sum(Obs-Sim)/sum(Obs)))
  EF<-with(DF,1-(sum(Obs-Sim)^2)/(sum(Obs-mean(Obs))^2))
  Result<-c(RMSE,R2,b0,RRMSE,RSR,MRE,ME,MAE,PBIAS,EF)
  names(Result)<-c("RMSE","R2","b0","RRMSE","RSR","MRE","ME","MAE","PBIAS","EF")
  return(Result)
}


#'@title Show the results of different models
#'@param model_list List. Including output results of different models.
#'@param names Vector. Name of models.
#'@importFrom ggplot2 aes_
#'@export
#'@return A list for ggplot2 plot
compare_model_plot<-function(model_list,names){
  names(model_list)<-names

  Selcet_Sim_Obs<-function(modelRe){
    data<-modelRe$Result[,c("Julian","Sim_SoilWater","SoilWater")]
    return(data)
  }

  mergeData<-plyr::ldply(.data = model_list,.fun = Selcet_Sim_Obs,.id="Model")

  names(mergeData)
  p1<-ggplot2::ggplot(data = mergeData)+
    ggplot2::geom_line(aes_(x=~Julian,y=~Sim_SoilWater),size=1)+
    ggplot2::geom_point(aes_(x=~Julian,y=~SoilWater),size=2)+
    ggplot2::scale_color_manual(values = c("red","blue","black"))+
    ggplot2::labs(x="Julian",y="Soil water(mm)",color="Model")+
    ggplot2::facet_grid(.~Model)+
    ggplot2::theme_bw()+
    ggplot2::theme(panel.grid=ggplot2::element_blank(),
          legend.position = "top")

  p2<-ggplot2::ggplot(data = mergeData,aes_(x=~SoilWater,y=~Sim_SoilWater))+
    ggplot2::geom_point(size=2)+
    ggplot2::geom_smooth(method = "lm",formula = "y~x",se = F,color="black")+
    ggpmisc::stat_poly_eq(aes_(label = ~paste(..eq.label.., ..rr.label.., sep = '~~~~')), formula = y ~ x, parse = T) + #添加回归方程和调整R方
    ggplot2::scale_color_manual(values = c("red","blue","black"))+
    ggplot2::facet_grid(.~Model)+
    ggplot2::labs(x="Measure value (mm)",y="Simnlation value (mm)")+
    ggplot2::theme_bw()+
    ggplot2::theme(panel.grid=ggplot2::element_blank())

  result_list<-list(p1,p2)
  names(result_list)<-c("Dynamic","Relationship")
  return(result_list)
}
