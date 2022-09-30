#'@title Calculating Kcini value
#'@param Pmean is the average depth od infiltrated water per wetting events(mm)
#'@param ET0 mean ET0 during initial period(mm/day)
#'@param tw is the mean interval between wetting events(days)
#'@param type soil type:"coarse soil textures" and "medium and fine soil textures"
#'@param fw the fraction of surfaces wetted by irrigation or rain (0-1)
#'@export
#'@return A value for Kcini value
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

cal_Kcini_for_SingleKc<-function(Pmean,ET0,tw,type,fw){
  Eso<-1.15*ET0
  if(Pmean<=10){
    TEWcor<-10
    REWcor<-min(max(2.5,6/(ET0^0.5)),7)
    t1<-REWcor/Eso
    Kcini<-fw*((TEWcor-(TEWcor-REWcor)*exp((-(tw-t1)*Eso*(1+REWcor/(TEWcor-REWcor)))/(TEWcor)))/(tw*ET0))
  }else if(Pmean>=40&type=="coarse soil textures"){
    TEWcor<-min(15,7*(ET0^0.5))
    REWcor<-min(9,TEWcor-0.01)
    t1<-REWcor/Eso
    Kcini<-fw*((TEWcor-(TEWcor-REWcor)*exp((-(tw-t1)*Eso*(1+REWcor/(TEWcor-REWcor)))/(TEWcor)))/(tw*ET0))
  }else if(Pmean>=40&type=="medium and fine soil textures"){
    TEWcor<-min(28,13*(ET0^0.5))
    REWcor<-min(9,TEWcor-0.01)
    t1<-REWcor/Eso
    Kcini<-fw*((TEWcor-(TEWcor-REWcor)*exp((-(tw-t1)*Eso*(1+REWcor/(TEWcor-REWcor)))/(TEWcor)))/(tw*ET0))
  }else if(Pmean>10&Pmean<40&type=="coarse soil textures"){
    TEWcor29<-10
    REWcor29<-min(max(2.5,6/(ET0^0.5)),7)
    t1<-REWcor29/Eso
    Kcini29<-(TEWcor29-(TEWcor29-REWcor29)*exp((-(tw-t1)*Eso*(1+REWcor29/(TEWcor29-REWcor29)))/(TEWcor29)))/(tw*ET0)

    TEWcor30<-min(15,7*(ET0^0.5))
    REWcor30<-min(9,TEWcor30-0.01)
    t1<-REWcor30/Eso
    Kcini30<-(TEWcor30-(TEWcor30-REWcor30)*exp((-(tw-t1)*Eso*(1+REWcor30/(TEWcor30-REWcor30)))/(TEWcor30)))/(tw*ET0)

    Kcini=Kcini29*fw+((Pmean-10)/(40-10))*(Kcini30*fw-Kcini29*fw)
  }else if(Pmean>10&Pmean<40&type=="medium and fine soil textures"){
    TEWcor29<-10
    REWcor29<-min(max(2.5,6/(ET0^0.5)),7)
    t1<-REWcor29/Eso
    Kcini29<-(TEWcor29-(TEWcor29-REWcor29)*exp((-(tw-t1)*Eso*(1+REWcor29/(TEWcor29-REWcor29)))/(TEWcor29)))/(tw*ET0)

    TEWcor30<-min(28,13*(ET0^0.5))
    REWcor30<-min(9,TEWcor-0.01)
    t1<-REWcor30/Eso
    Kcini30<-(TEWcor30-(TEWcor30-REWcor30)*exp((-(tw-t1)*Eso*(1+REWcor30/(TEWcor30-REWcor30)))/(TEWcor30)))/(tw*ET0)

    Kcini=Kcini29*fw+((Pmean-10)/(40-10))*(Kcini30*fw-Kcini29*fw)
  }
  return(Kcini)
}


#'@title Crop coefficient for the mid-season stage
#'@description Typical values for the crop coefficient at the end of the late season growth stage, Kc end, are
#'  listed in Table 12 for various agricultural crops.For specific adjustment in climates where RHmin differs from 45 \% or where u2 is larger or
#'  smaller than 2.0 m/s.
#'@param Ktable value for Kc mid taken from Table 12
#'@param RHmine mean value for daily minimum relative humidity during the mid-season
#'  growth stage , for 20 <= RHmine <= 80
#'@param u2e mean value for daily wind speed at 2 m height over grass during the mid-
#'  season growth stage (mls), for 1  <= u2e <= 6
#'@param he mean plant height during the mid-season stage [m] for 0.1 m < h < 10 m
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.
#'@export
#'@return A value for Kcmid value

cal_Kcmid_for_singleKc<-function(RHmine,u2e,Ktable,he){
  if (RHmine!=45&u2e!=2) {
    Kcmid<-Ktable+(0.04*(u2e-2)-0.004*(RHmine-45))*(he/3)^0.3
    return(Kcmid)
  }else{
    stop("ERROR:RHmin!=45&u2!=2")
    }
}

#'@title Crop coefficient for the end of the late season stage
#'@description Typical values for the crop coefficient at the end of the late season growth stage, Kc end, are
#'  listed in Table 12 for various agricultural crops.
#'@param Ktable value for Kc mid taken from Table 12
#'@param RHmine mean value for daily minimum relative humidity during the mid-season
#'  growth stage , for 20 <= RHmine<= 80
#'@param u2e mean value for daily wind speed at 2 m height over grass during the mid-
#'  season growth stage (m/s), for 1  <= u2e <= 6
#'@param he mean plant height during the mid-season stage (m) for 0.1 m < h < 10 m
#'@note only applied when the tabulated values for Kc end exceed 0.45. The
#'  equation reduces the Kc end with increasing RHmin. This reduction in Kc end is characteristic
#'  of crops that are harvested 'green' or before becoming completely dead and dry (Kc end >= 0.45).
#'
#'  No adjustment is made when Kc end (Table) < 0.45 (Kc end = Kc end (Tab)). When
#'   crops are allowed to senesce and dry in the field (as evidenced by Kc end < 0.45), u2 and
#'  RHmin have less effect on Kc end and no adjustment is necessary.
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.
#'@export
#'@return A value for Kcend value

cal_Kcend_for_singleKc<-function(RHmine,u2e,Ktable,he){
  while (RHmine!=45&u2e!=2) {
    if(Ktable>=0.45){
      Kcend<-Ktable+(0.04*(u2e-2)-0.004*(RHmine-45))*(he/3)^0.3
    }else{
      Kcend<-Ktable
    }
    return(Kcend)
  }
  stop("ERROR:RHmin!=45&u2!=2")
}

#'@title Calculating deep percolation
#'@param P Precipitation
#'@param I Irrigation
#'@param ETa Actual evapotranspiration
#'@param Dri_start The depletion of root layer
#'@return A value for deep percolation
#'@export
cal_DP_for_singleKc<-function(P,I,ETa,Dri_start){
  DP0<-P+I-ETa-Dri_start
  if(DP0>=0){
    DP<-DP0
  }else{
    DP<-0
  }
  return(DP)
}

#========================================================================
#Simulation for evapotranspiration using single crop coefficient method

#'@title Simulation for evapotranspiration using single crop coefficient method
#'@param data A data box. Contains the daily data required by the model. You can
#'  refer to the function create_modelData()
#'@param param A list. Contains additional parameters.
#'@note The stages of data should include all four stages.
#'  If a crop has multiple growth cycles, each cycle should include all four stages.
#'@importFrom dplyr n
#'@importFrom ggplot2 aes_
#'@importFrom rlang .data
#'@export
#'@return A list for the model result including a data frame of daily model result ,a list of plots, A data frame of summary data
#'@examples
#'   library(simET)
#'   #--Data preparation
#'   data("FIalfalfa")
#'   #--Parameter preparation
#'   param_SingleKc<-list(Kc_mid=1.2,#Kcb for mid-season stage
#'                        Kc_end=1.15,#Kcb for late season stage
#'                        rootDepth=1.2,#Maximum root depth
#'                        #The soil type used for calculating
#'                        #Kc for initial stage
#'                        soil_type="coarse soil textures",
#'                        Dr_start=40,#Initial depletion of root layer
#'                        TAW=290,#Total available soil water of the root zone
#'                        p=0.55,#Evapotranspiration depletion factor
#'                        Field_capacity=420,#Field capacity of root layer
#'                        fw=1,#The fraction of the surface wetted
#'                        #Capillary rise model parameters
#'                        CR_param=c(420,-0.32,303,-0.16,-1.4,6.8,1.11,-0.98)
#'                        )
#'  #--Run model
#'  Re_SingleKc<- Model_single_Kc(data = FIalfalfa, param = param_SingleKc)
#'  #--The Result data
#'  Re_SingleKc$Result
#'  Re_SingleKc$Plot
#'  #--The goodness Of Fit
#'  estimate_goodnessOfFit(Sim = Re_SingleKc$Result$Sim_SoilWater,
#'                          Obs = Re_SingleKc$Result$SoilWater)



Model_single_Kc<-function(data,param){
  Kc_mid<-param$Kc_mid
  Kc_end<-param$Kc_end
  rootDepth<-param$rootDepth
  soil_type<-param$soil_type
  Dr_start<-param$Dr_start
  TAW<-param$TAW
  p<-param$p
  Field_capacity<-param$Field_capacity
  CR_param<-param$CR_param
  fw<-param$fw
  #----1.ET0
  data_ET0<-dplyr::mutate(data,
                   ET0=cal_ET0_from_PM_for_daily(Latitude=.data$Latitude,
                                         Altitude=.data$Altitude,
                                         J=.data$Julian,
                                         Tmax=.data$Tmax,
                                         Tmin=.data$Tmin,
                                         Rs=.data$Rs,
                                         RHmean=.data$RHmean,
                                         Wind=.data$Wind))

  #----2.Calculating the corresponding parameters
        #for calculating and adjusting the standard KC value
  Kcparam_data<-within(data_ET0,{
      Cut<-factor(Cut,ordered = T,levels = unique(Cut))
      Stage<-factor(Stage,ordered = T,levels = unique(Stage))
    })%>%
    dplyr::mutate(Water=.data$Precipitation+.data$Irrigation)%>%
      #Irrigation and Precipitation as wetting events.
      #Effective rainfall can be used for Precipitation
    dplyr::group_by(Cut,Stage)%>%
    dplyr::summarise(Days=dplyr::n(),
              SumWater=sum(.data$Water),
              ET0=mean(.data$ET0),
              nw=sum(.data$SumWater>0),
                #Setting greater than zero is valid, but it is not.
                #It should be a parameter here (0 can be replaced by a parameter)
              RHmine=mean(.data$RHmin),
              u2e=mean(.data$Wind),
              he=mean(.data$Height))%>%
      #Calculate the days of each stage,
      #cumulative rainfall and irrigation,
      #average ET0, wetting times, average minimum humidity,
      #wind speed and plant height
    dplyr::mutate(Pmean=(.data$SumWater)/.data$nw,#Average wetting days
           tw=.data$Days/(.data$nw+0.5))
  Kcparam_data<-within(Kcparam_data,{
      type<-rep(soil_type,time=nrow(Kcparam_data))
      Ktable<-NA
      Ktable[Stage=="Mid"]<-Kc_mid
      Ktable[Stage=="End"]<-Kc_end
    })


  #----3.Calculating the KC value of each stage
  #It is divided into four lists corresponding to their respective calculation KC functions
  Kcparam_data_stage<-split(Kcparam_data,c("Ini","Dev","Mid","End"))#Divide into four data frames

  Kcini<-with(Kcparam_data_stage[["Ini"]],mapply(cal_Kcini_for_SingleKc,Pmean,ET0,tw,type,fw))
  Kcmid<-with(Kcparam_data_stage[["Mid"]],mapply(cal_Kcmid_for_singleKc,RHmine, u2e, Ktable, he))
  Kcend<-with(Kcparam_data_stage[["End"]],mapply(cal_Kcend_for_singleKc,RHmine, u2e, Ktable, he))#分别复制KC计算函数
  KcDeve<-rep(NA,time=max(data$Cut))
    #Generate the corresponding number of KcDve
    #according to the number of cuts
  Kc<-c(Kcini,KcDeve,Kcmid,Kcend)
    #Merge the separately calculated KC values

  #Add the calculated KC value to the data frame, Kc_data.
  Kc_data<-Kcparam_data%>%
    dplyr::arrange(Stage,Cut)
  Kc_data<-within(Kc_data,{Kc=Kc})%>%
    dplyr::arrange(Cut,Stage)


  #----4.Calculate daily Kc value
  #Calculate the days of the corresponding stages
  Data_D<-Kc_data%>%dplyr::select("Cut","Stage","Days")
  Data_D<-tidyr::spread(Data_D,key = Stage,value="Days")%>%
    dplyr::mutate(End1=.data$End-1,
           End2=1)%>%
    dplyr::select(-"End")
  Data_D<-tidyr::gather(Data_D,key = Stage,value = "Days","Ini":"End2")
  Data_D<-within(Data_D,{
      Stage<-factor(Stage,ordered = T,levels = unique(Stage))
      Cut<-factor(Cut,ordered = T,levels = unique(Cut))
    })%>%
    dplyr::arrange(Cut,Stage)%>%
    as.data.frame()#Calculate the days of each cut and each stage

  #Calculate the Kc value of the corresponding stages
  Data_K<-Kc_data%>%
    dplyr::select(Cut,Stage,Kc)

  Data_K<-tidyr::spread(Data_K,key = Stage,value=Kc)%>%
    dplyr::mutate(End1=NA,
           End2=.data$End)%>%
    dplyr::select(-"End")
  Data_K<-tidyr::gather(Data_K,key = Stage,value = Kc,"Ini":"End2")
  Data_K<-within(Data_K,{
      Stage<-factor(Stage,ordered = T,levels = unique(Stage))
      Cut<-factor(Cut,ordered = T,levels = unique(Cut))
    })%>%
    dplyr::arrange(Cut,Stage)%>%
    as.data.frame()#Calculate the Kc value of each cut and each stage

  #merge the days data and Kc data
  Kc1<-Data_K$Kc
  Data_DK<-cbind(Data_D,Kc1)%>%
    dplyr::mutate(Daycum=cumsum(.data$Days),
           Day2=data$Julian[1]+.data$Daycum-1)

  Df<-data.frame(Day=data$Julian[1]:Data_DK$Day2[nrow(Data_DK)])
  #Create a data frame to record daily Kc
  dailyKc<-within(Df,{
      Item<-NA
      Item<-rep(1:nrow(Data_DK),time=Data_DK$Days)
      Kc<-NA
      Kc<-rep(Data_DK$Kc,time=Data_DK$Days)
      stage<-rep(Data_DK$Stage,time=Data_DK$Days)
    })

  #Linear interpolation Kc to obtain daily Kc value
  DataVector<-dailyKc$Kc
  dailyKc$Kc2<-linear_interpolation(DataVector)

  #Merge the calculated daily KC value into the data frame
  #with daily ET0 data, and calculate the ETc under the standard condition
  Kc<-linear_interpolation(DataVector)
  Result_ETc<-cbind(data_ET0,Kc)%>%
    dplyr::mutate(ETc=.data$ET0*Kc)

  #----5.Calculate the water stress coefficient and adjust ETc,
  #Calculating capillary rise and daily water balance
  # names(Result_ETc)
  Result_ETc2<-Result_ETc%>%
    dplyr::select("Year":"Cut","Precipitation","LAI":"ETc")

  #Add a new row
  # head(Result_ETc2)
  DF_row<-matrix(ncol =ncol(Result_ETc2))
  colnames(DF_row)<-names(Result_ETc2)
  Result_ETc2_row<-rbind(DF_row,Result_ETc2)
  #Add 5 new columns for calculation
  DF_col<-matrix(ncol = 6,nrow = nrow(Result_ETc2_row))
  colnames(DF_col)<-c("Dri_start","Ks","ETc_adj","CR","DP","Dri_end")
  Result_ETc2_col<-cbind(Result_ETc2_row,DF_col)

  New_data<-Result_ETc2_col
  New_data[1,"Dri_end"]<-Dr_start
  # head(New_data)
  names(CR_param)<-c("a1","b1","a2","b2","a3","b3","a4","b4")
  for (row in 2:nrow(New_data)) {
    # row<-2
    New_data[row,"Dri_start"]<-New_data[row-1,"Dri_end"]
    New_data[row,"Ks"]<-cal_WaterStressCoef(Dr=New_data[row,"Dri_start"], TAW, p)
    New_data[row,"ETc_adj"]<-New_data[row,"ETc"]*New_data[row,"Ks"]
    New_data[row,"CR"]<-cal_capillaryRise(CR_param[["a1"]],CR_param[["b1"]],CR_param[["a2"]],CR_param[["b2"]],CR_param[["a3"]],CR_param[["b3"]],CR_param[["a4"]],CR_param[["b4"]],Dw=New_data[row,"GroundwaterDepth"]-rootDepth,Wa=New_data[row,"Dri_start"],LAI=New_data[row,"LAI"],ETm=New_data[row,"ETc"])
    New_data[row,"DP"]<-cal_DP_for_singleKc(P=New_data[row,"Precipitation"],I=New_data[row,"Irrigation"],ETa=New_data[row,"ETc_adj"],Dri_start=New_data[row,"Dri_start"])
    New_data[row,"Dri_end"]<-New_data[row,"Dri_start"]-New_data[row,"Precipitation"]-New_data[row,"Irrigation"]-New_data[row,"CR"]+New_data[row,"ETc_adj"]+New_data[row,"DP"]
  }

  New_data<-dplyr::mutate(New_data,Sim_SoilWater=Field_capacity-.data$Dri_end)
  # names(New_data)
  #----6.plot
  mytheme<-ggplot2::theme_bw()+
    ggplot2::theme(panel.grid=ggplot2::element_blank())
  p_ET0<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~ET0))+
    # geom_point()+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="ET0 (mm)")+
    mytheme
  p_Kc<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~Kc))+
    # geom_point()+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="Kc")+
    mytheme
  p_ETc<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~ETc))+
    # geom_point()+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="ETc (mm)")+
    mytheme
  p_Ks<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~Ks))+
    # geom_point()+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="Ks")+
    mytheme
  p_ETc_adj<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~ETc_adj))+
    # geom_point()+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="ETc_adj (mm)")+
    mytheme
  p_CR<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~CR))+
    # geom_point()+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="CR (mm)")+
    mytheme
  p_DP<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~DP))+
    # geom_point()+
    ggplot2::geom_line()+
    ggplot2::labs(x="Julian",y="DP (mm)")+
    mytheme
  # sum(New_data$CR,na.rm = T)
  # sum(New_data$ETc_adj,na.rm = T)
  p_W<-ggplot2::ggplot(data =New_data ,aes_(x=~Julian,y=~(Field_capacity-Dri_start)))+
    #geom_point()+
    ggplot2::geom_line()+
    ggplot2::geom_point(aes_(x=~Julian,y=~SoilWater))+
    ggplot2::geom_hline(yintercept = Field_capacity,color="blue")+
    ggplot2::geom_hline(yintercept = Field_capacity-TAW,color="red")+#萎蔫含水量
    ggplot2::geom_hline(yintercept = TAW*p+(Field_capacity-TAW),color="yellow",linetype=2)+#萎蔫含水量
    ggplot2::labs(x="Julian",y="Soil water (mm)")+
    mytheme

  summ<-New_data[-1,]%>%
    dplyr::summarise(DP=sum(.data$DP),
              CR=sum(.data$CR),
              ETc_adj=sum(.data$ETc_adj))

  p_result<-ggpubr::ggarrange(p_ET0,p_Kc,p_ETc,p_Ks,p_ETc_adj,p_CR,p_DP,p_W,ncol = 2,nrow = 4)
  Result<-list(New_data,p_result,summ)
  names(Result)<-c("Result","Plot","Summary")
  return(Result)
}
