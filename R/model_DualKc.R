#=================================================
#Model_DualKc2
#=================================================

#'@title Adjust the recommended Kc values at the middle and late stages
#'@param Kcb_table Recommended value of KC in FAO 56 at the middle and late stages
#'@param u2 wind speed at 2 m
#'@param RHmin Minimum relative humidity
#'@param h Plant height
#'@return A value for adjust Kc at middle and late stages
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.
Kcb_adj_for_DualKc<-function(Kcb_table,u2,RHmin,h){
  if((RHmin!=45 | u2!=2)&Kcb_table>=0.45){
    Kcb<-Kcb_table+(0.04*(u2-2)-0.004*(RHmin-45))*(h/3)^0.3
  }else{
    Kcb<-Kcb_table
  }
  return(Kcb)
}

#'@title An upper limit on the evaporation and transpiration from any cropped surface
#'@description  It is imposed to reflect the natural constraints placed on
#'    available energy represented by the energy balance difference Rn - G - H
#'@param u2 The wind speed at 2 m
#'@param RHmin Minimum relative humidity
#'@param h Plant height
#'@param Kcb Basal crop coefficient
#'@export
#'@return A vector for the upper limit on the evaporation and transpiration from any cropped surface
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.
cal_Kc_max_for_DualKc<-function(u2,RHmin,h,Kcb){
  limit_A <- 1.2+(0.04*(u2-2)-0.004*(RHmin-45))*(h/3)^0.3
  limit_B <- Kcb+0.05
  Kc_max <- ifelse(limit_A >= limit_B,
                   limit_A,limit_B)
  # Kc_max<-max(1.2+(0.04*(u2-2)-0.004*(RHmin-45))*(h/3)^0.3,(Kcb+0.05))
  return(Kc_max)
}


#'@title calculating total evaporable water
#'@description maximum depth of water that can be evaporated from the soil
#'    when the topsoil has been initially completely wetted.
#'    Estimated TEW for Kr calculation
#'@param FC Soil water content at field capacity, m3 m-3
#'@param WP Soil water content at wilting point, m3 m-3
#'@param Ze Depth of the surface soil layer that is subject to drying by way of evaporation, 0.10-0.15 m.
#'@export
#'@return A value for total evaporable water
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.
cal_TEW_for_DualKc<-function(FC,WP,Ze){
  TEW<-1000*(FC-0.5*WP)*Ze
  return(TEW)
}


#'@title Dimensionless evaporation reduction coefficient
#'@description  It dependent on the soil water depletion
#'    (cumulative depth of evaporation) from the topsoil
#'    layer (Kr = 1 when De,i-1 is equal or lesser than REW)
#'@param TEW maximum cumulative depth of evaporation (depletion) from
#'    the soil surface layer when Kr = 0 (TEW = total evaporable water)
#'@param REW cumulative depth of evaporation (depletion) at the end of stage 1
#'    (REW = readily evaporable water), mm
#'@param De cumulative depth of evaporation (depletion) from the soil surface
#'    layer at the end of day i-1 (the previous day),mm
#'@export
#'@return A value for evaporation reduction coefficient
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.
cal_Kr_for_DualKc<-function(TEW,REW,De){
  if(De<REW){
    Kr<-1
  }else if(De>REW&De<TEW){
    Kr<-(TEW-De)/(TEW-REW)
  }else if(De>=TEW){
    De<-TEW
    # Kr<-(TEW-De)/(TEW-REW)
    Kr<-0
  }
  return(Kr)
}


#'@title Deep percolation loss from the topsoil layer
#'@param P Precipitation
#'@param I Irrigation
#'@param Dei_start Depletion in the topsoil layer
#'@param fw Fraction of soil surface wetted by irrigation, 0.01-1
#'@export
#'@return A value for deep percolation loss from the topsoil layer
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.
cal_DPe_for_DualKc<-function(P,I,Dei_start,fw){#fw,
  DP0<-P+(I/fw)-Dei_start
  if(DP0>=0){
    DP<-DP0
  }else{
    DP<-0
  }
  return(DP)
}


#'@title Deep percolation loss from the root layer
#'@param P Precipitation
#'@param Irrigation Irrigation
#'@param ETa Actual evapotranspiration
#'@param Dri_start Depletion in the root layer
#'@export
#'@return A value for deep percolation loss from the root layer
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.
cal_DPr_for_DualKc<-function(P,Irrigation,ETa,Dri_start){
  DP0<-P+Irrigation-ETa-Dri_start
  if(DP0>=0){
    DP<-DP0
  }else{
    DP<-0
  }
  return(DP)
}
# cal_DPr<-function(DPe,Ta,Dri_start){
#   DP0<-DPe-Ta-Dri_start
#   if(DP0>=0){
#     DP<-DP0
#   }else{
#     DP<-0
#   }
#   return(DP)
# }


#'@title calculating the depletion in the  topsoil layer  at the end of the day
#'@description In fact, it performs water balance in a day
#'@param Dei_start Depletion in the topsoil layer
#'@param P Precipitation
#'@param I Irrigation
#'@param E Evaporation on day i, mm
#'@param Dep Deep percolation loss from the topsoil layer on day i
#'    if soil water content exceeds field capacity, mm
#'@param TEW Maximum cumulative depth of evaporation (depletion) from the topsoil layer
#'@export
#'@return A value for the depletion in the topsoil layer at the end of the day
cal_Dei_for_DualKc<-function(Dei_start,P,I,E,Dep,TEW){
  #I and E have been divided by fw and few in the formula parameters in the cycle.
  #Note that transpiration is not considered
  Dei_end0<-Dei_start-P-I+E+Dep
  if(Dei_end0>TEW){
    Dei_end<-TEW
  }else{
    Dei_end<-Dei_end0
  }
  return(Dei_end)
}

################################################################################
#'@title Simulation of evapotranspiration using dual crop coefficient method
#'@param data A data box. Contains the daily data required by the model. You can refer to the function create_modelData()
#'@param param A list. Contains additional parameters.
#'    list(Kini,Kmid,Kend,fw,rootDepth,Dei_start,Dri_start,FCe,WPe,Ze,REW,TAW,p,FCr,CR_param)
#'@note The stages of data should include all four stages.
#'  If a crop has multiple growth cycles, each cycle should include all four stages.
#'@return A list for the model result including a data frame of daily model result ,a list of plots, A data frame of summary data
#'@importFrom dplyr n
#'@importFrom ggplot2 aes_
#'@export
#'@examples
#'  library(simET)
#'  data("FIalfalfa")
#'  names(FIalfalfa)
#'  #--Model parameter
#'  Dparam_FI<-list(Kini=0.3,#Kcb for initial stage
#'                Kmid=1.15,#Kcb for mid-season stage
#'                Kend=1.1,#Kcb for late season stage
#'                DI=FALSE,#Is it drip irrigation?
#'                fw=1,# The fraction of the surface wetted
#'                rootDepth=1.2,#Maximum root depth
#'                Dei_start=0,#Initial depletion of evaporation layer
#'                Dri_start=35,#Initial depletion of root layer
#'                FCe=0.22,#Field capacity of evaporation layer
#'                WPe=0.15,#Wilting point of evaporation layer
#'                Ze=0.15,#Depth of the surface soil layer
#'                REW=6,#Readily evaporable water
#'                TAW=297,#Total available soil water of the root zone
#'                p=0.55,#Evapotranspiration depletion factor
#'                FCrmm=430,#Field capacity of root layer
#'                CR_param=c(430,-0.32,310,-0.16,-1.4,6.8,1.11,-0.98)
#'                  #Capillary rise model parameters
#'                )
#'  #--Run model
#'  Model_re_FI<-Model_DualKc(data = FIalfalfa,param = Dparam_FI)
#'  #--The Result data
#'  Model_re_FI$Result
#'  Model_re_FI$Plot
#'  #--The goodness Of Fit
#'  estimate_goodnessOfFit(Sim = Model_re_FI$Result$Sim_SoilWater,
#'                          Obs = Model_re_FI$Result$SoilWater)


Model_DualKc<-function(data,param){
  # data=SI_2019_T3
  # param=param_SI2019_T3
  Kini<-param$Kini
  Kmid<-param$Kmid
  Kend<-param$Kend
  DI<-param$DI
  fw<-param$fw
  rootDepth<-param$rootDepth
  Dei_start<-param$Dei_start
  Dri_start<-param$Dri_start
  FCe<-param$FCe
  WPe<-param$WPe
  Ze<-param$Ze
  REW<-param$REW
  TAW<-param$TAW
  p<-param$p
  FCrmm<-param$FCrmm
  CR_param<-param$CR_param
  #----Calculating ET0
  data_ET0<-data%>%
    dplyr::mutate(ET0=cal_ET0_from_PM_for_daily(Latitude=data$Latitude,
                                         Altitude=data$Altitude,
                                         J=data$Julian,
                                         Tmax=data$Tmax,
                                         Tmin=data$Tmin,
                                         Rs=data$Rs,
                                         RHmean=data$RHmean,
                                         Wind=data$Wind))

  #----2.Calculating the daily Kc under the standard condition
  # Kcb_adj_for_DualKc(Kcb_table = 1.1,u2=2.2,RHmin = 30,h=0.4)
  # Kcb_adj_for_DualKc(Kcb_table = 0.25,u2=2.2,RHmin = 30,h=0.4)

  #Prepare data for adjusting Kcb
  Data2_pre<-data_ET0%>%
    dplyr::group_by(Cut,Stage)%>%
    dplyr::summarise(
      Days=n(),
      u2=mean(data_ET0$Wind),
      RHmin=mean(data_ET0$RHmin),
      h=mean(data_ET0$Height)
    )

  Data2<-within(Data2_pre,{
      Kcb_table<-NA
      Kcb_table[Stage=="Ini"]<-Kini
      Kcb_table[Stage=="Mid"]<-Kmid#
      Kcb_table[Stage=="End"]<-Kend#使用标准的
      Cut<-factor(Cut,ordered = T)
      Stage<-factor(Stage,levels = c("Ini","Dev","Mid","End"))
    })
  # str(Data2)
  #The Kcb adjustment preparation data is divided into four periods
  Da<-split(Data2,f=Data2$Stage)# by=c("Ini","Dev","Mid","End")

  #Calculate the Kcb of each cut and each stage respectively
  Kcini<-with(Da[["Ini"]],mapply(function(Kcb_table){return(Kcb_table)}, Kcb_table))# Kcb ini 不调整
  Kcmid<-with(Da[["Mid"]],mapply(Kcb_adj_for_DualKc,Kcb_table,u2,RHmin,h))
  Kcend<-with(Da[["End"]],mapply(Kcb_adj_for_DualKc,Kcb_table,u2,RHmin,h))#分别复制KC计算函数
  KcDeve<-rep(NA,time=max(data$Cut))#Determine how many cuts
  Kc<-c(Kcini,KcDeve,Kcmid,Kcend)#Merge the separately calculated KC values

  Data_Kcb<-Data2%>%
    dplyr::arrange(Stage,Cut)%>%
    as.data.frame()%>%
    dplyr::mutate(Kc=Kc)%>%
    dplyr::arrange(Cut,Stage)

  Data_D<-Data_Kcb%>%#Calculating the duration days of each cut and stage
    dplyr::select("Cut","Stage","Days")

  Data_D<-tidyr::spread(Data_D,key = Stage,value=3)#key = Stage,value=Days
  Data_D<-Data_D%>%dplyr::mutate(End1=Data_D$End-1,
           End2=1)%>%
    dplyr::select(-"End")

  Data_D<-tidyr::gather(Data_D,key = "Stage",value = "Days",2:6)#Ini:End2
  Data_D<-within(Data_D,{
      Stage<-factor(Stage,ordered = T,levels = unique(Stage))
      Cut<-factor(Cut,ordered = T,levels = unique(Cut))
    })%>%
    dplyr::arrange(Cut,Stage)

  Data_K<-Data_Kcb%>%##Calculate the Kcb corresponding to each stage of each cut
    dplyr::select(Cut,Stage,Kc)
  Data_K<-tidyr::spread(Data_K,key = Stage,value=3)#value=Kc
  Data_K<-Data_K%>%dplyr::mutate(End1=NA,End2=Data_K$End)%>%
    dplyr::select(-"End")

  Data_K<-tidyr::gather(Data_K,key = "Stage",value = "Kc",2:6)#`Ini`:`End2`

  Data_K<-within(Data_K,{
      Stage<-factor(Stage,ordered = T,levels = unique(Stage))
      Cut<-factor(Cut,ordered = T,levels = unique(Cut))
    })%>%
    dplyr::arrange(Cut,Stage)
  Kcb<-Data_K[,"Kc"]

  ##Merge the Kcb corresponding to each stage of each cut, duration days and cumulative days
  Data_DK<-cbind(Data_D,Kcb)
  Data_DK<-Data_DK%>%dplyr::mutate(Daycum=cumsum(Data_DK$Days))
  Data_DK<-Data_DK%>%dplyr::mutate(Day2=data$Julian[1]+Data_DK$Daycum-1)


  #Create a data frame to record daily Kcb
  Df<-data.frame(Day=data$Julian[1]:Data_DK$Day2[nrow(Data_DK)])
  dailyKc<-within(Df,{
      Item<-NA
      Item<-rep(1:nrow(Data_DK),time=Data_DK$Days)
      Kc<-NA
      Kc<-rep(Data_DK$Kc,time=Data_DK$Days)
      stage<-rep(Data_DK$stage,time=Data_DK$Days)
    })

  #Linear interpolation for Kcb
  dailyKc$Kcb<-linear_interpolation(dailyKc$Kc)

  #Generate result data with KCB
  Re_data<-data_ET0
  Re_data$Kcb<-dailyKc$Kcb

  #----3.Calculate daily transpiration and evaporation
  #Calculate the water balance of evaporation layer and root layer,respectively

  Re_data_fc_fw<-Re_data%>%
    dplyr::mutate(fc=(1-exp(-0.7*Re_data$LAI))^(1+Re_data$Height))
  Re_data_fc_fw<-Re_data_fc_fw%>%dplyr::mutate(fw_irr=ifelse(DI==TRUE,(1-0.67*Re_data_fc_fw$fc)*fw,fw))


  #Create a new data frame and add it to the beginning of the original data
  DF_row<-matrix(ncol =ncol(Re_data_fc_fw))
  colnames(DF_row)<-names(Re_data_fc_fw)
  Result_ETc2_row<-rbind(DF_row,Re_data_fc_fw)
  #Add blank 18 columns for calculation
  DF_col<-matrix(ncol = 18,nrow = nrow(Result_ETc2_row))
  colnames(DF_col)<-c("Dei_start","Dri_start","Kr","Kc_max","fw","Ke","few","Ks","E_few","DPe","CR","DPr","Kcb_adj","E","Ta","Dei_end","Dri_end","ETc_adj")
  Result_ETc2_col<-cbind(Result_ETc2_row,DF_col)

  ##########参数##########
  New_data<-Result_ETc2_col
  New_data[1,"Dei_end"]<-Dei_start #Initial soil water in evaporation layer
  New_data[1,"Dri_end"]<-Dri_start#Initial soil water in the root layer below the evaporation layer
  New_data[1,"fw"]<-1 #FW = 1 before irrigation



  #Note: it is the water balance of two soil layers
  names(CR_param)<-c("a1","b1","a2","b2","a3","b3","a4","b4")
  for (row in 2:nrow(New_data)) {
    # print(row)
    # row=9
    #Initial value of two-layer soil layer
    New_data[row,"Dei_start"]<-New_data[row-1,"Dei_end"]
    New_data[row,"Dri_start"]<-New_data[row-1,"Dri_end"]

    #Calculate the water balance of evaporation layer
    New_data[row,"Kr"]<-cal_Kr_for_DualKc(TEW=cal_TEW_for_DualKc(FC=FCe,WP=WPe,Ze=Ze),REW=REW,De=New_data[row,"Dei_start"])#Kr是针对表层土壤的；参数：FC,WP.Ze,REW
    New_data[row,"Kc_max"]<-cal_Kc_max_for_DualKc(u2=New_data[row,"Wind"],RHmin=New_data[row,"RHmin"],h=New_data[row,"Height"],Kcb=New_data[row,"Kcb"])
    New_data[row,"fw"]<-ifelse(New_data[row,"Irrigation"]>0&New_data[row,"Precipitation"]<4&New_data[row,"Precipitation"]>=0,New_data[row,"fw_irr"],
                               ifelse(New_data[row,"Irrigation"]>0&New_data[row,"Precipitation"]>4,1,
                                      ifelse(New_data[row,"Irrigation"]==0&New_data[row,"Precipitation"]>4,1,
                                             New_data[row-1,"fw"])
                               )
    )
    New_data[row,"few"]<-ifelse(1-New_data[row,"fc"]<=New_data[row,"fw"],1-New_data[row,"fc"],New_data[row,"fw"])
    New_data[row,"Ke"]<-min(New_data[row,"Kr"]*(New_data[row,"Kc_max"]-New_data[row,"Kcb"]),New_data[row,"few"]*New_data[row,"Kc_max"])
    New_data[row,"E_few"]<-(New_data[row,"Ke"]*New_data[row,"ET0"])/New_data[row,"few"]
      #这句代码好像没有用
    New_data[row,"DPe"]<-cal_DPe_for_DualKc(P=New_data[row,"Precipitation"],I=New_data[row,"Irrigation"],Dei_start=New_data[row,"Dei_start"],fw=New_data[row,"fw"])
    New_data[row,"E"]<-New_data[row,"Ke"]*New_data[row,"ET0"]
    New_data[row,"Dei_end"]<-cal_Dei_for_DualKc(Dei_start=New_data[row,"Dei_start"],P=New_data[row,"Precipitation"],I=(New_data[row,"Irrigation"]/New_data[row,"fw"]),E=(New_data[row,"E"]/New_data[row,"few"]),Dep=New_data[row,"DPe"],
                                     TEW=cal_TEW_for_DualKc(FC=FCe,WP=WPe,Ze=Ze))

    #Calculate the water balance of root layer
    New_data[row,"Ks"]<-cal_WaterStressCoef(Dr=New_data[row,"Dri_start"], TAW=TAW, p=p)#参数 TAW,p
    New_data[row,"Kcb_adj"]<-New_data[row,"Ks"]*New_data[row,"Kcb"]
    New_data[row,"Ta"]<-New_data[row,"Kcb_adj"]*New_data[row,"ET0"]
    New_data[row,"CR"]<-cal_capillaryRise(CR_param[[1]],CR_param[[2]],CR_param[[3]],CR_param[[4]],CR_param[[5]],CR_param[[6]],CR_param[[7]],CR_param[[8]],Dw=New_data[row,"GroundwaterDepth"]-rootDepth,Wa=New_data[row,"Dri_start"],LAI=New_data[row,"LAI"],ETm=(New_data[row,"Kcb"]+New_data[row,"Ke"])*New_data[row,"ET0"])
    New_data[row,"ETc_adj"]<-(New_data[row,"Ke"]+New_data[row,"Kcb_adj"])*New_data[row,"ET0"]
    New_data[row,"DPr"]<-cal_DPr_for_DualKc(P=New_data[row,"Precipitation"],Irrigation=New_data[row,"Irrigation"],ETa=New_data[row,"ETc_adj"],Dri_start=New_data[row,"Dri_start"])
    New_data[row,"Dri_end"]<-New_data[row,"Dri_start"]-New_data[row,"Precipitation"]-New_data[row,"Irrigation"]-New_data[row,"CR"]+New_data[row,"ETc_adj"]+New_data[row,"DPr"]

  }

  #Calculate the simulated soil water of the maximum root layer
  New_data<-dplyr::mutate(New_data,Sim_SoilWater=FCrmm-New_data$Dri_end)
  # #Calculation model fitting performance
  # Obs_sim_data<-New_data[2:nrow(New_data),c("SoilWater","Dri_end")]%>%
  #   drop_na()%>%
  #   mutate(Sim_value=FCrmm-Dri_end)
  #
  # RMSE<-with(Obs_sim_data,sqrt(sum((SoilWater-Sim_value)^2)/nrow(Obs_sim_data)))
  # R2<-with(Obs_sim_data,(sum((SoilWater-mean(SoilWater))*(Sim_value-mean(Sim_value)))/(sum((SoilWater-mean(SoilWater))^2)^0.5*sum((Sim_value-mean(Sim_value))^2)^0.5))^2)
  # b0<-with(Obs_sim_data,sum(Sim_value*SoilWater)/sum(SoilWater^2))
  #
  # stat<-c(RMSE,R2,b0)
  # names(stat)<-c("RMSE","R2","b0")
  # # #R2
  # # R2<-with(Obs_sim_data,
  # #          1-sum((Sim_value-SoilWater)^2)/sum((SoilWater-mean(SoilWater))^2)
  # #          )
  # # stat<-c(RMSE)
  # # names(stat)<-c("RMSE")

  #Plot result
  #mytheme for plot
  mytheme<-ggplot2::theme_bw()+
    ggplot2::theme(panel.grid=ggplot2::element_blank())

  p_ET0<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~ET0))+
    # geom_point()+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="ET0 mm")+
    mytheme
  p_Kcb<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~Kcb))+
    # geom_point()+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="Kcb")+
    mytheme
  p_Ks<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~Ks))+
    # geom_point()+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="Ks")+
    mytheme
  p_Kcb_adj<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~Kcb_adj))+
    # geom_point()+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="Kcb_adj")+
    mytheme
  p_Ke<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~Ke))+
    # geom_point()+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="Ke")+
    mytheme
  p_E<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~E))+
    # geom_point(data =New_data ,aes(x=Julian,y=EVAmm))+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="E (mm)")+
    mytheme
  p_Ta<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~Ta))+
    # geom_point()+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="Ta (mm)")+
    mytheme
  p_CR<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~CR))+
    # geom_point()+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="CR (mm)")+
    mytheme
  p_DP<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~DPr))+
    # geom_point()+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="DP (mm)")+
    mytheme

  # p_ETc_adj<-ggplot(data =New_data ,aes(x=Julian,y=ETc_adj))+
  #   # geom_point()+
  #   geom_line()+
  #   labs(x="Julian",y="ETc_adj (mm)")+
  #   theme_bw()
  # p_SWC15<-ggplot()+
  #   geom_line(data =New_data ,aes(x=Julian,y=FCe*Ze*1000-Dei_start))+
  #   geom_point(data =New_data ,aes(x=Julian,y=SWC15*Ze*1000))+
  #   # geom_hline(yintercept = FCe*100+FCrmm,color="blue")+
  #   # geom_hline(yintercept = 197,color="red")+
  #   # geom_hline(yintercept = TAW/2+197,color="black")+
  #   labs(x="Julian",y="SWC15 (mm)")+
  #   theme_bw()
  p_SWC<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~(FCrmm-Dri_start)))+
    ggplot2::geom_point(data =New_data ,aes_(x=~Julian,y=~SoilWater))+
    ggplot2::geom_line()+
    # geom_hline(yintercept = FCe*100+FCrmm,color="blue")+
    # geom_hline(yintercept = 197,color="red")+
    # geom_hline(yintercept = TAW/2+197,color="black")+
    ggplot2::labs(x="Julian",y="SWC (mm)")+
    mytheme

  # p_w2<-ggplot(data =New_data ,aes(x=SoilWater,y=FCrmm-Dri_start))+
  #   geom_point()+
  #   # geom_line()+
  #   geom_abline(intercept = 0)+
  #   xlim(c(0,FCrmm))+
  #   ylim(c(0,FCrmm))+
  #   labs(x="Obs(mm)",y="Sim (mm)")+
  #   theme_bw()
  #summary result
  New_data<-New_data[-1,]
  summ<-New_data%>%
   dplyr::summarise(DP=sum(New_data$DPr),
              CR=sum(New_data$CR),
              E=sum(New_data$E),
              Ta=sum(New_data$Ta),
              Tc=sum(Kcb*New_data$ET0),
              ETc_adj=sum(New_data$ETc_adj))



  p_result<-ggpubr::ggarrange(p_ET0,p_Kcb,p_Ks,p_Kcb_adj,p_Ke,p_E,p_Ta,p_CR,p_DP,p_SWC,ncol = 2,nrow = 5)
  Result<-list(New_data,p_result,summ)
  names(Result)<-c("Result","Plot","Summary")
  return(Result)
}
















