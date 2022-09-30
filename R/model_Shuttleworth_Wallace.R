#===================================================
#Model_SW
#===================================================

#-----energy available to crop and soil
#'@title Calculating net radiation available to the system(soil and crop)
#'@param Ac net radiation available to the crop (W/m2 ground)
#'@param As net radiation available to the soil(W/m2 ground)
#'@return A vector for net radiation available to the system(soil and crop) (W/m2 ground)
#'@export

cal_netRadiationForSystem<-function(Ac,As){
  A<-Ac+As
  return(A)
}

#'@title The canopy penetration probability for net radiation
#'@param kRn the canopy extinction coefficient for net radiation (taken as 0.3)
#'@param L the leaf area index (m2 leaf m-2 ground)
#'@return The vector for canopy penetration probability for net radiation
#'@export

cal_canopyPenetrationProbabilityForNetRadiation<-function(kRn,L){
  pRn<-exp(-kRn*L)
  return(pRn)
}

#'@title Calculating net radiation available to the crop
#'@param pRn canopy penetration probability for net radiation
#'@param Rn the net radiation (W/m2 ground).see cal_hourlyNetRadiation()
#'@return A vector for net radiation available to the crop (W/m2 ground)
#'@export

cal_netRadiationForCrop<-function(pRn,Rn){
  Ac<-(1-pRn)*Rn
  return(Ac)
}


#'@title Calculating net radiation available to the soil
#'@param pRn canopy penetration probability for net radiation
#'@param Rn the net radiation (W/m2 ground).see cal_hourlyNetRadiation()
#'@param G the soil heat flux(W/m2 ground)
#'@return A vector for net radiation available to the soil (W/m2 ground)
#'@export

cal_netRadiationForSoil<-function(pRn,Rn,G){
  As<-pRn*Rn-G
  return(As)
}



#'@title Calculating Soil/ground heat flux
#'@param pRn the canopy penetration probility for net radiation
#'@param Rn the net radiation (W/m2 ground)
#'@return A vector for the soil heat flux (W/m2 ground)
#'@export

cal_soilHeatFlux<-function(pRn,Rn){
  G<-0.3*pRn*Rn
  return(G)
}

#---- Wind speed profile

#'@title Calculating zero plane displacement height
#'@param h the plant height (m)
#'@return A vector for zero plane displacement height (m)
#'@export

cal_zeroPlaneHeight<-function(h){
  d<-0.64*h
  return(d)
}

#'@title Calculating the crop roughness length
#'@param h the plant height (m)
#'@return A vector for the crop roughness length(m)
#'@export

cal_cropRoughnessLength<-function(h){
  z0<-0.13*h
  return(z0)
}

#'@title Calculating friction velocity
#'@param k the von Karman constant (0.4)
#'@param u_zr the wind speed(m/s) at the reference height zr (m).
#'@param zr the height of the weather station(m).
#'@param d zero plane displacement height (m)
#'@param z0 the crop roughness length(m)
#'@return A vector for friction velocity(m/s)
#'@export

cal_frictionVelocity<-function(k=0.4,zr=2,u_zr,d,z0){
  u_<-k*u_zr/log((zr-d)/z0)
  return(u_)
}

#'@title Calculating wind speed above and within the canopies.
#'@param z the height (m)
#'@param h the plant height (m)
#'@param u_ the friction velocity (m/s)
#'@param k the von Karman constant (0.4)
#'@param d zero plane displacement height (m)
#'@param z0 the crop roughness length(m)
#'@param nu the wind speed extinction coefficient (taken as 2)
#'@return A vector for the wind speed (m/s) at height z(m)
#'@export

cal_windSpeed_Canopy<-function(z,h,u_,k=0.4,d,z0,nu=2){
  u_z<-ifelse(z>=h,(u_/k)*log((z-d)/z0),
              (u_/k)*log((h-d)/z0)*exp(-nu*(1-z/h)))
  return(u_z)
}


#----Flux resistances

#'@title Calculating eddy diffusivity at the canopy top
#'@param k the von Karman constant (0.4)
#'@param u_ the friction velocity (m/s)
#'@param h the plant height (m)
#'@return A vector for eddy diffusivity at the canopy top(m2/s)
#'@export

cal_eddyDiffusivity_Canopytop<-function(k=0.4,u_,h){
  Kh<-k*u_*h
  return(Kh)
}


#'@title Calculating eddy diffusivity at height z
#'@param Kh eddy diffusivity at the canopy top(m2/s)
#'@param nK the eddy diffusivity extinction coefficient (taken as 2)
#'@param z height(m)
#'@param h the plant height(m)
#'@return A vector for eddy diffusivity at height z (m2/s)
#'@export

cal_eddyDiffusivity_heightZ<-function(Kh,nK,z,h){
  K_z<-Kh*exp(-nK*(1-z/h))
  return(K_z)
}

#'@title Calculating soil surface to mean canopy flow
#'@description the resistance between the soil surface and the mean canopy flow (s m-1)
#'@param h the plant height(m)
#'@param nK the eddy diffusivity extinction coefficient (taken as 2)
#'@param Kh eddy diffusivity at the canopy top(m2/s)
#'@param zs0 is the soil surface roughnesslength (m). Note: for flat, tilled land, zs0 can be taken as 0.004 m.
#'@param z0 the crop roughness length (m)
#'@param d zero plane displacement height (m)
#'@return A vector for aerodynamic resistancessoil surface to mean canopy flow (s m-1)
#'@export

cal_soilSurfaceToMeanCanopyFlow<-function(h,nK,Kh,zs0,z0,d){
  n1<-(h*exp(nK))/(nK*Kh)
  n2<-exp(-nK*zs0/h)
  n3<-exp(-nK*(z0+d)/h)
  r_s_a<-n1*(n2-n3)
  return(r_s_a)
}


#'@title Calculating mean canopy flow to reference level
#'@param k the von Karman constant (0.4)
#'@param u_ the friction velocity (m s-1)
#'@param zr is the reference height (m).the height of the weather station(m).
#'@param d zero plane displacement height (m)
#'@param h the plant height(m)
#'@param nK the eddy diffusivity extinction coefficient (taken as 2)
#'@param z0 the crop roughness length (m)
#'@return A vector for mean canopy flow to reference level
#'@export

cal_meanCanopyFlowToReferenceLevel<-function(k=0.4,u_,zr,d,h,nK,z0){
  n4<-log((zr-d)/(h-d))/(k*u_)
  n5<-1/(nK*k*u_)
  n6<-exp(nK*(1-(z0+d)/h))-1
  r_a_a<-n4+n5*n6
  return(r_a_a)
}

#'@title Calculating bulk boundary layer resistance
#'@param nu the wind speed extinction coefficient (taken as 2)
#'@param u_h the wind speed at the canopy top (i.e., at plant height h) (m s-1)
#'@param w is the mean leaf width (m)
#'@param L leaf area index
#'@return A vector for bulk boundary layer resistance (s/m)
#'@export

cal_bulkBoundaryLayerResistance<-function(nu,u_h,w,L){
  n7<-0.012*L*(1-exp(-nu/2))
  n8=sqrt(u_h/w)
  r_c_a<-nu/(n7*n8)
  return(r_c_a)
}

#'@title Calculating canopy resistance
#'@param a1,a2 are empirical coefficients, dependent on the crop type.
#'@param It the total hourly solar irradiance (W m-2 ground)
#'@param L the leaf area index (m2 leaf m-2 ground)
#'@param Lmax  the maximum total leaf area index (m2 leaf m-2 ground)
#'@return A vector for canopy resistance(s/m)
#'@export

cal_canopyResistance<-function(a1,a2,It,L,Lmax){
  # r_st<-((a1+0.5*It)/(a2*(0.5*It)))# the leaf stomatal resistance(s/m)
  totalPAR<-ifelse((It*0.5)<=0.1,0.1,It*0.5)
  rst<-(a1+totalPAR)/(a2*totalPAR)
  r_c_s<-ifelse(L<=(0.5*Lmax), rst/L, rst/(0.5*Lmax))
  return(r_c_s)
}


#'@title Calculating soil surface resistance
#'@param tau soil tortuosity (taken as 2)
#'@param l is the dry soil layer thickness (taken as the first soil layer thickness) (m)
#'@param PHI_p is soil porosity
#'@param Dm_v the vapor diffusion coefficient in air (24.7  10-6 m2 s-1)
#'@param lambda_p the soil pore-size distribution index from the Brooks-Corey equation.
#'@param THETA_v_l volumetric soil water content (m-3 m-3) of the first soil layer
#'@param THETA_v_sat_l saturated soil water content (m3 m-3) of the first soil layer
#'@return A vector for the soil surface resistance (s m-1)
#'@export

cal_soilSurfaceResistance<-function(tau,l,PHI_p,Dm_v,lambda_p,THETA_v_l,THETA_v_sat_l){
  r_s_sdry<-(tau*l)/(PHI_p*Dm_v)
  r_s_s<-r_s_sdry*exp(-(1/lambda_p)*(THETA_v_l/THETA_v_sat_l))
  return(r_s_s)
}

#'@title Calculating total latent heat flux
#'@param DELTA the slope of the saturated vapor pressure curve (mbar K-1)
#'@param gamma is the psychometric constant (0.658 mbar K-1)
#'@param r_a_a the aerodynamic resistance between the mean canopy flow and reference height (s m-1)
#'@param r_c_a the bulk boundary layer resistance (s m-1)
#'@param r_s_a is the aerodynamic resistance between the soil and mean canopy flow (s m-1)
#'@param r_c_s the canopy resistance(s m-1)
#'@param r_s_s soil surface resistance (s m-1)
#'@param A energy available to the system (total)(W m-2 ground)
#'@param As energy available to soil (W m-2 ground)
#'@param Ac energy available to crop (W m-2 ground)
#'@param rho_cp is the volumetric heat capacity for air (1221.09 J m-3 K-1)
#'@param D the vapor pressure deficit (mbar)
#'@return A vector for the total latent heat flux (W m-2 ground)
#'@export

cal_totalLatentHeatFlux<-function(DELTA,gamma,r_a_a,r_c_a,r_s_a,r_c_s,r_s_s,A,rho_cp=1221.09,D,As,Ac){
  Ra<-(DELTA+gamma)*r_a_a
  Rc<-(DELTA+gamma)*r_c_a+gamma*r_c_s
  Rs<-(DELTA+gamma)*r_s_a+gamma*r_s_s
  Cc<-(1+(Rc*Ra)/(Rs*(Rc+Ra)))^(-1)
  Cs<-(1+(Rs*Ra)/(Rc*(Rs+Ra)))^(-1)
  n9<-(DELTA*A+(rho_cp*D-DELTA*r_c_a*As)/(r_a_a+r_c_a))
  n10<-(DELTA+gamma*(1+r_c_s/(r_a_a+r_c_a)))
  PMc<-n9/n10
  n11<-(DELTA*A+(rho_cp*D-DELTA*r_s_a*Ac)/(r_a_a+r_s_a))
  n12<-(DELTA+gamma*(1+r_s_s/(r_a_a+r_s_a)))
  PMs<-n11/n12
  lambda_ET<-Cc*PMc+Cs*PMs
  return(lambda_ET)
}

#'@title Calculating air vapor pressure deficit at the mean canopy
#'@param D the vapor pressure deficit (mbar)
#'@param r_a_a the aerodynamic resistance between the mean canopy flow and reference height (s m-1)
#'@param rho_cp the volumetric heat capacity for air (1221.09 J m-3 K-1)
#'@param DELTA is the slope of the saturated vapor pressure curve (mbar K-1)
#'@param A is the total energy available to the system (W m-2 ground)
#'@param gamma psychometric constant (0.658 mbar K-1 )
#'@param lambda_ET the total latent heat flux (W m-2 ground)
#'@note Knowing D0 is essential because this value is used to calculate the latent and sensible heat fluxes for the soil and
#'    crop components.
#'@return A vector for the vapor pressure deficit at the mean canopy flow (mbar)
#'@export

cal_airVaporPressureDeficit_meanCanopyflow<-function(D,r_a_a,rho_cp=1221.09,DELTA,A,gamma=0.658,lambda_ET){
  D0<-D+(r_a_a/rho_cp)*(DELTA*A-(DELTA+gamma)*lambda_ET)
  return(D0)
}

#'@title Calculating latent heat fluxes for soil
#'@param DELTA the slope of the saturated vapor pressure curve (mbar K-1)
#'@param As energy available to the soil (W m-2 ground)
#'@param rho_cp the volumetric heat capacity for air (1221.09 J m-3 K-1)
#'@param D0 the vapor pressure deficit at the mean canopy flow
#'@param r_s_a the aerodynamic resistance between the soil and mean canopy flow (s m-1)
#'@param gamma the psychometric constant (0.658 mbar K-1)
#'@param r_s_s soil surface resistance,  (s m-1)
#'@return A vector for soil latent heat fluxes (W m-2 ground)
#'@export

cal_latentHeatFluxesForSoil<-function(DELTA,As,rho_cp=1221.09,D0,r_s_a,gamma=0.658,r_s_s){
  lambda_ETs<-(DELTA*As+rho_cp*D0/r_s_a)/(DELTA+gamma*(r_s_s+r_s_a)/r_s_a)
  return(lambda_ETs)
}

#'@title Calculating latent heat fluxes for crop
#'@param DELTA the slope of the saturated vapor pressure curve (mbar K-1)
#'@param Ac energy available to the crop (W m-2 ground)
#'@param rho_cp the volumetric heat capacity for air (1221.09 J m-3 K-1)
#'@param D0 the vapor pressure deficit at the mean canopy flow
#'@param r_c_a is the bulk boundary layer resistance (s m-1)
#'@param gamma is the psychometric constant (0.658 mbar K-1)
#'@param r_c_s the canopy resistance (s m-1)
#'@return A vector for latent heat fluxes for crop (W m-2 ground)
#'@export

cal_latentHeatFluxesForCrop<-function(DELTA,Ac,rho_cp,D0,r_c_a,gamma,r_c_s){
  lambda_ETc<-(DELTA*Ac+rho_cp*D0/r_c_a)/(DELTA+gamma*(r_c_s+r_c_a)/r_c_a)
  return(lambda_ETc)
}

#'@title Calaulating sensible heat fluxes for soil
#'@param gamma is the psychometric constant (0.658 mbar K-1)
#'@param As energy available to the soil (W m-2 ground)
#'@param  r_s_s soil surface resistance,  (s m-1)
#'@param r_s_a is the aerodynamic resistance between the soil and mean anopy flow (s m-1);
#'@param rho_cp the volumetric heat capacity for air (1221.09 J m-3 K-1)
#'@param D0 the vapor pressure deficit at the mean canopy flow
#'@param DELTA  the slope of the saturated vapor pressure curve (mbar K-1)
#'@return A vector for soil sensible heat fluxes(W m-2 ground)
#'@export

cal_sensibleHeatFluxesForSoil<-function(gamma=0.659,As,r_s_s,r_s_a,rho_cp,D0,DELTA){
  Hs<-(gamma*As*(r_s_s+r_s_a)-rho_cp*D0)/(DELTA*r_s_a+gamma*(r_s_s+r_s_a))
  return(Hs)
}


#'@title Calculating sensible heat fluxes for crop
#'@param gamma is the psychometric constant (0.658 mbar K-1)
#'@param Ac energy available to the crop (W m-2 ground)
#'@param r_c_s the canopy resistance (s m-1)
#'@param r_c_a is the bulk boundary layer resistance (s m-1)
#'@param rho_cp the volumetric heat capacity for air (1221.09 J m-3 K-1)
#'@param D0 the vapor pressure deficit at the mean canopy flow
#'@param DELTA the slope of the saturated vapor pressure curve (mbar K-1)
#'@return A vector for sensible heat fluxes for crop
#'@export

cal_sensibleHeatFluxesForCrop<-function(gamma,Ac,r_c_s,r_c_a,rho_cp,D0,DELTA){
  Hc<-(gamma*Ac*(r_c_s+r_c_a)-rho_cp*D0)/(DELTA*r_c_a+gamma*(r_c_s+r_c_a))
  return(Hc)
}

#'@title Calculating canopy temperature
#'@param Hc crop sensible heat fluxes (W m-2)
#'@param r_c_a the bulk boundary layer resistance (s m-1)
#'@param Hs soil sensible heat fluxes (W m-2)
#'@param r_a_a the aerodynamic resistance between the mean canopy flow and reference level (s m-1)
#'@param rho_cp the volumetric heat capacity for air (1221.09 J m-3 K-1)
#'@param Tr Tr is the air temperature at reference level (Celsius degree). weather station.
#'@return A vector for the canopy (foliage) temperature (Celsius degree)
#'@export

cal_canopyTem<-function(Hc,r_c_a,Hs,r_a_a,rho_cp,Tr){
  Tf<-((Hc*r_c_a+(Hs+Hc)*r_a_a)/rho_cp)+Tr
}

#------------------
#'@title Calculating percolation for excess water
#'@param THETA_i_t0 the water amount of the day before in soil layer i (mm)
#'@param Pe_i_t1 the percolation of previous soil layer(mm)
#'@param THETA_sat_i soil saturation water amount (mm)
#'@return A value for percolation for excess water (mm)
#'@export

cal_percolationForExcessWater<-function(THETA_i_t0,Pe_i_t1,THETA_sat_i){
  pe<-ifelse(THETA_i_t0+Pe_i_t1<=THETA_sat_i,
             0,
             THETA_i_t0+Pe_i_t1-THETA_sat_i)
  return(pe)
}

#'@title Calculating the volumetric water content after redistribution
#'@param THETA_v_sat soil saturation water content (m3 m-3)
#'@param alpha empirical coefficient. 13 for homogenous soil, 13-16 for heterogeneous soil
#'@param Ksat saturated hydraulic conductitity
#'@param deltaT time step difference (day)
#'@param L the thickness(m) of soil layer i
#'@param THETA11_1 the volumetric water content before redistribution (m3 m-3)
#'@return A value for the volumetric water content after redistribution(m3 m-3)
#'@export
cal_afterRedistribution<-function(THETA_v_sat,alpha,Ksat,deltaT,L,THETA11_1){
  THETA_v_i_t<-THETA_v_sat-(THETA_v_sat/alpha)*log(((alpha*Ksat*deltaT)/(L*THETA_v_sat))+exp((alpha/THETA_v_sat)*(THETA_v_sat-THETA11_1)))
  return(THETA_v_i_t)
}

#'@title Calculating reduction factor for evaporation
#'@param THETA the water amount of the day before in stop soil layer  (mm)
#'@param THETA_sat soil saturation water content in stop oil layer (mm)
#'@return A value for reduction factor for evaporation
#'@export
cal_reductionFactorForE<-function(THETA,THETA_sat){
  R_D_e<-1/(1+(3.6073*(THETA/THETA_sat))^(-9.3172))
  return(R_D_e)
}

#'@title calculating reduction factor for transpiration
#'@param THETA_v_wp soil water content at wilting point (m3 m-3)
#'@param p a coefficient. 0.5 for C3 and 0.3 for C4 plant
#'@param THETA_v_sat soil saturation water content (m3 m-3)
#'@param THETA_v the water amount of the day before in root layer(m3 m-3)
#'@return A value for reduction factor for transpiration
#'@export
cal_reductionFactorForT<-function(THETA_v_wp,p,THETA_v_sat,THETA_v){
  THETA_v_cr<-THETA_v_wp+p*(THETA_v_sat-THETA_v_wp)
  R_D_t<-(THETA_v-THETA_v_wp)/(THETA_v_cr-THETA_v_wp)
  return(R_D_t)
}

#===============================================================================

#'@title Simulation of evapotranspiration using Shuttleworth-Wallace model
#'@param data A data box. Contains the daily data required by the model.
#'  You can refer to the function create_modelData()
#'@param param A list. Contains additional parameters.
#'@note The stages of data should include all four stages.
#'  If a crop has multiple growth cycles, each cycle should include all four stages.
#'@importFrom dplyr n
#'@importFrom ggplot2 aes_
#'@export
#'@importFrom rlang .data
#'@return A list for the model result including a data frame of daily model result
#'        ,a list of plots, A data frame of summary data
#'@examples
#'    library(simET)
#'    #--Data preparation
#'    data("FIalfalfa")
#'    #--Parameter preparation
#'    param_SW<-list(
#'        plant=list(
#'                  #the canopy extinction coefficient for net radiation
#'                  kRn=0.3,
#'                  alpha_plant=0.3,#Canopy reflectance
#'                  w=0.01,#Leaf width
#'                  Lmax=10,#Maximum leaf area index
#'                  a1=10,# Leaf stomatal resistance coefficients
#'                  a2=0.005,# Leaf stomatal resistance coefficients
#'                  p=0.5,#the param of reduction factor for T
#'                  rootDepth=1.2 #Maximum root depth
#'                  ),
#'        Soil=list(
#'                zs0=0.04,#The soil surface roughnesslength (m)
#'                tau2=2,#Soil tortuosity
#'                PHI_p=2,	 #Soil porosity
#'                #The soil pore-size distribution index
#'                #from the Brooks-Corey equation.
#'                lambda_p=0.18,
#'                l1=0.02,#Depth of the surface soil layer (m)
#'                l2=1.2,#Depth of the root layer (m)
#'                #Saturation water content of evaporation layer
#'                THETA_v_sat_1=0.36,
#'                THETA_v_sat_2=0.40, #Saturation water content of root layer
#'                THETA_start_1=0.2,#Initial water content of evaporation layer
#'                THETA_start_2=0.36,#Initial water content of root layer
#'                THETA_wp1=0.15,#Wilting point of evaporation layer
#'                THETA_wp2=0.15,#Wilting point of root layer
#'                #Empirical coefficient of evaporation layer.
#'                #13 for homogenous soil
#'                alpha1=14,
#'                alpha2=14,#Empirical coefficient of root layer.
#'                #Saturated hydraulic conductitity of evaporation layer
#'                Ksat_1=13.52,
#'                Ksat_2=0.02,#Saturated hydraulic conductitity of root layer
#'                #Capillary rise model parameters
#'                CR_param=c(430,-0.32,313,-0.16,-1.4,6.8,0.5,-0.98)
#'                ),
#'      Mete=list(
#'               nu=2,#The wind speed extinction coefficient
#'               nK=2, #The eddy diffusivity extinction coefficient(taken as 2)
#'               zr=2,	 #The reference height (m)
#'               ##The vapor diffusion coefficient in air (24.7 10-6 m2 s-1)
#'               Dm_v=24.7*10^(-6),
#'               deltaT=1 #Time step difference (day)
#'               )
#'          )
#'    #--Run model
#'    Re_SW<-Model_SW(data = FIalfalfa, param = param_SW)
#'    #--The Result data
#'    Re_SW$Result
#'    Re_SW$Plot
#'    #--The goodness Of Fit
#'    estimate_goodnessOfFit(Sim = Re_SW$Result$Sim_SoilWater,
#'                          Obs = Re_SW$Result$SoilWater)

Model_SW<-function(data,param){
  data$LAI[data$LAI==0]<-0.01
  data$Height[data$Height==0]<-5
  #----1.Parameter setting
  #Parameter of plant
  kRn=param$plant$kRn
  alpha_plant=param$plant$alpha_plant
  w=param$plant$w
  Lmax=param$plant$Lmax
  a1=param$plant$a1
  a2=param$plant$a2
  p=param$plant$p
  rootDepth=param$plant$rootDepth

  #Parameter of soil
  zs0=param$Soil$zs0
  tau2=param$Soil$tau2
  PHI_p=param$Soil$PHI_p
  lambda_p=param$Soil$lambda_p
  l1=param$Soil$l1
  l2=param$Soil$l2
  THETA_v_sat_1=param$Soil$THETA_v_sat_1
  THETA_v_sat_2=param$Soil$THETA_v_sat_2
  THETA_start_1=param$Soil$THETA_start_1
  THETA_start_2=param$Soil$THETA_start_2
  THETA_wp1=param$Soil$THETA_wp1
  THETA_wp2=param$Soil$THETA_wp2
  alpha1=param$Soil$alpha1
  alpha2=param$Soil$alpha2
  Ksat_1=param$Soil$Ksat_1
  Ksat_2=param$Soil$Ksat_2
  CR_param=param$Soil$CR_param

  #Meteorological parameters
  nu=param$Mete$nu
  nK=param$Mete$nK
  zr=param$Mete$zr
  Dm_v=param$Mete$Dm_v
  deltaT=param$Mete$deltaT

  #----2.Calculate meteorological parameters
  Mete_growth<-data%>%
    dplyr::mutate(Tmean=(.data$Tmax+.data$Tmin)/2,#unit 摄氏度
           Lat=convert_angert_to_radian(.data$Latitude),#unit 弧度
           TKmax=convert_degreesCelsius_to_Fahrenheit(.data$Tmax),#华氏度
           TKmin=convert_degreesCelsius_to_Fahrenheit(.data$Tmin),#华氏度
           ea=cal_ActualVapourPressure_from_RHmean(RHmean=.data$RHmean, Tmax=.data$Tmax, Tmin=.data$Tmin)*10,#KPa转化为mba
           Ra=cal_extraterrestrialRadiation_for_daily(J=.data$Julian,lat=.data$Lat),
           Rso=cal_skySolarRadiation_withas_elevation(z=.data$Altitude, Ra=.data$Ra),
           Rns=cal_netSolarRadiation(alpha=alpha_plant, Rs=.data$Rs),
           Rnl=cal_netLongwaveRadiation(TKmax=.data$TKmax, TKmin=.data$TKmin, ea=.data$ea, Rs=.data$Rs, Rso=.data$Rso),
           It=convert_Rad_unit(rad = .data$Rs,type =  "MJ_m2_day_to_W_m2"),
           Rn=convert_Rad_unit(rad = cal_netRadiation(Rns=.data$Rns, Rnl=.data$Rnl),type =  "MJ_m2_day_to_W_m2"),
           es=cal_meanSaturationVapourPressure(.data$Tmax ,.data$Tmin)*10,#KPa转化为mba
           D=.data$es-.data$ea,
           DELTA=cal_slopeOfSaturationVapourPressureCurve(.data$Tmean)*10#mbar degres-1
    )


  #----3.Calculate evapotranspiration and water balance
  # names(Mete_growth)
  baseData<-Mete_growth%>%
    dplyr::mutate(pRn=cal_canopyPenetrationProbabilityForNetRadiation(kRn = kRn,L=.data$LAI),
           G=cal_soilHeatFlux(pRn=.data$pRn, Rn=.data$Rn),
           d=cal_zeroPlaneHeight(h=.data$Height/100),
           z0=cal_cropRoughnessLength(h=.data$Height/100),
           u_=cal_frictionVelocity(k = 0.4, zr = zr, u_zr=.data$Wind, d=.data$d, z0=.data$z0),
           u_h=cal_windSpeed_Canopy(z=.data$Height/100, h=.data$Height/100, u_=.data$u_, k = 0.4, d=.data$d, z0=.data$z0, nu = nu),
           Kh=cal_eddyDiffusivity_Canopytop(k = 0.4, u_=.data$u_, h=.data$Height/100),
           r_s_a=cal_soilSurfaceToMeanCanopyFlow(h=.data$Height/100,nK=nK,Kh=.data$Kh,zs0=zs0,z0=.data$z0,d=.data$d),
           r_a_a=cal_meanCanopyFlowToReferenceLevel(k = 0.4, u_=.data$u_, zr=zr, d=.data$d, h=.data$Height/100, nK=nK, z0=.data$z0),
           r_c_a=cal_bulkBoundaryLayerResistance(nu=nu, u_h=.data$u_h, w=w,L=.data$LAI),
           r_c_s=cal_canopyResistance(a1=a1, a2=a2, It=.data$It, L=.data$LAI, Lmax=Lmax))

  #Add a new row to place the initial water
  DF_row<-matrix(ncol =ncol(baseData))
  colnames(DF_row)<-names(baseData)
  baseData_row<-rbind(DF_row,baseData)
  #Add some new column for the daily water balance
  DF_col<-matrix(ncol = 26,nrow = nrow(baseData_row))
  colnames(DF_col)<-c("r_s_s","As","Ac","A","lambda_ET","ET","D0","lambda_ETs","lambda_ETc","ETs","ETc",
                      "Pe1","pe1","pd1","THETA11_1","THETA1_1","Ea","THETA_1",
                      "Pe2","pe2","pd2","THETA11_2","THETA1_2","CR","Ta","THETA_2")
  baseData_row_col<-cbind(baseData_row,DF_col)

  #Set initial moisture water
  New_data<-baseData_row_col
  New_data[1,"THETA_1"]<-THETA_start_1*l1*1000
  New_data[1,"THETA_2"]<-THETA_start_2*l2*1000
  names(CR_param)<-c("a1","b1","a2","b2","a3","b3","a4","b4")

  names(New_data)
  for (row in 2:nrow(New_data)) {
    # row<-67
    New_data[row,"r_s_s"]<-cal_soilSurfaceResistance(tau=tau2, l=l1, PHI_p=0.42, Dm_v=Dm_v, lambda_p=lambda_p,
                                                     THETA_v_l=New_data[row-1,"THETA_1"]/(l1*1000),THETA_v_sat_l=THETA_v_sat_1)
    New_data[row,"As"]<-cal_netRadiationForSoil(pRn=New_data[row,"pRn"], Rn=New_data[row,"Rn"], G=New_data[row,"G"])
    New_data[row,"Ac"]<-cal_netRadiationForCrop(pRn=New_data[row,"pRn"], Rn=New_data[row,"Rn"])
    New_data[row,"A"]<-cal_netRadiationForSystem(Ac=New_data[row,"Ac"], As=New_data[row,"As"])
    New_data[row,"lambda_ET"]<-cal_totalLatentHeatFlux(DELTA=New_data[row,"DELTA"],
                                                       gamma=0.658,
                                                       r_a_a=New_data[row,"r_a_a"],
                                                       r_c_a=New_data[row,"r_c_a"],
                                                       r_s_a=New_data[row,"r_s_a"],
                                                       r_c_s=New_data[row,"r_c_s"],
                                                       r_s_s=New_data[row,"r_s_s"],
                                                       A=New_data[row,"A"],
                                                       rho_cp = 1221.09,
                                                       D=New_data[row,"D"],
                                                       As=New_data[row,"As"],
                                                       Ac=New_data[row,"Ac"])
    New_data[row,"ET"]<-convert_Rad_unit(rad=New_data[row,"lambda_ET"],type = "W_m2_to_MJ_m2_day")/2.454
    New_data[row,"D0"]<-cal_airVaporPressureDeficit_meanCanopyflow(D=New_data[row,"D"],
                                                                   r_a_a=New_data[row,"r_a_a"],
                                                                   rho_cp = 1221.09,
                                                                   DELTA=New_data[row,"DELTA"],
                                                                   A=New_data[row,"A"],
                                                                   gamma = 0.658,
                                                                   lambda_ET=New_data[row,"lambda_ET"])
    New_data[row,"lambda_ETs"]<-cal_latentHeatFluxesForSoil(DELTA=New_data[row,"DELTA"],
                                                            As=New_data[row,"As"],
                                                            rho_cp = 1221.09,
                                                            D0=New_data[row,"D0"],
                                                            r_s_a=New_data[row,"r_s_a"],
                                                            gamma = 0.658,
                                                            r_s_s=New_data[row,"r_s_s"])
    New_data[row,"lambda_ETc"]<-cal_latentHeatFluxesForCrop(DELTA=New_data[row,"DELTA"],
                                                            Ac=New_data[row,"Ac"],
                                                            rho_cp=1221.09,
                                                            D0=New_data[row,"D0"],
                                                            r_c_a=New_data[row,"r_c_a"],
                                                            gamma=0.658,
                                                            r_c_s=New_data[row,"r_c_s"])
    New_data[row,"ETc"]<-convert_Rad_unit(New_data[row,"lambda_ETc"],type = "W_m2_to_MJ_m2_day")/2.454
    New_data[row,"ETs"]<-convert_Rad_unit(New_data[row,"lambda_ETs"],type = "W_m2_to_MJ_m2_day")/2.454
    #topsoil layer
    New_data[row,"Pe1"]<-New_data[row,"Precipitation"]+New_data[row,"Irrigation"]
    New_data[row,"pe1"]<-cal_percolationForExcessWater(THETA_i_t0=New_data[row-1,"THETA_1"],
                                                       Pe_i_t1=New_data[row,"Pe1"],
                                                       THETA_sat_i=THETA_v_sat_1*l1*1000)
    New_data[row,"THETA11_1"]<-New_data[row-1,"THETA_1"]+New_data[row,"Pe1"]-New_data[row,"pe1"]
    New_data[row,"THETA1_1"]<-cal_afterRedistribution(THETA_v_sat=THETA_v_sat_1,
                                                      alpha=alpha1,
                                                      Ksat=Ksat_1,
                                                      deltaT=deltaT,
                                                      L=l1,
                                                      THETA11_1=New_data[row,"THETA11_1"]/(l1*1000))*l1*1000
    New_data[row,"pd1"]<-New_data[row,"THETA11_1"]-New_data[row,"THETA1_1"]
    New_data[row,"Ea"]<-New_data[row,"ETs"]*cal_reductionFactorForE(THETA=New_data[row-1,"THETA_1"],THETA_sat=THETA_v_sat_1*l1*1000)
    New_data[row,"THETA_1"]<-ifelse(New_data[row,"THETA1_1"]-New_data[row,"Ea"]*0.26>THETA_wp1*l1*1000,
                                    New_data[row,"THETA1_1"]-New_data[row,"Ea"]*0.26,
                                    THETA_wp1*l1*1000)
    #root layer
    New_data[row,"Pe2"]<-New_data[row,"pe1"]+New_data[row,"pd1"]
    New_data[row,"pe2"]<-cal_percolationForExcessWater(THETA_i_t0=New_data[row-1,"THETA_2"],
                                                       Pe_i_t1=New_data[row,"Pe2"],
                                                       THETA_sat_i=THETA_v_sat_2*l2*1000)
    New_data[row,"THETA11_2"]<-New_data[row-1,"THETA_2"]+New_data[row,"Pe2"]-New_data[row,"pe2"]
    New_data[row,"THETA1_2"]<-cal_afterRedistribution(THETA_v_sat=THETA_v_sat_2,
                                                      alpha=alpha2,
                                                      Ksat=Ksat_2,
                                                      deltaT=deltaT,
                                                      L=l2,
                                                      THETA11_1=New_data[row,"THETA11_2"]/(l2*1000))*l2*1000
    New_data[row,"pd2"]<-New_data[row,"THETA11_2"]-New_data[row,"THETA1_2"]
    New_data[row,"Ta"]<-New_data[row,"ETc"]*cal_reductionFactorForT(THETA_v_wp=THETA_wp2,
                                                                    p=0.5,
                                                                    THETA_v_sat=THETA_v_sat_2,
                                                                    THETA_v=New_data[row-1,"THETA_2"]/(l2*1000))
    New_data[row,"CR"]<-cal_capillaryRise(CR_param[["a1"]],CR_param[["b1"]],CR_param[["a2"]],CR_param[["b2"]],CR_param[["a3"]],CR_param[["b3"]],CR_param[["a4"]],CR_param[["b4"]],Dw=New_data[row,"GroundwaterDepth"]-rootDepth,Wa=New_data[row-1,"THETA_2"],LAI=New_data[row,"LAI"],ETm=New_data[row,"ET"])
    New_data[row,"THETA_2"]<-ifelse(New_data[row,"THETA_1"]==THETA_wp1*l1*1000,
                                    New_data[row,"THETA1_2"]-New_data[row,"Ea"]-New_data[row,"Ta"]+New_data[row,"CR"],
                                    New_data[row,"THETA1_2"]-New_data[row,"Ea"]*0.74-New_data[row,"Ta"]+New_data[row,"CR"])
    # print(row)
  }

  New_data<-dplyr::mutate(New_data,Sim_SoilWater=.data$THETA_2)
  #注意水分的单位
  # names(New_data)

  #计算模型拟合结果
  # Obs_sim_data<-New_data[2:nrow(New_data),c("SoilWater","Dri_end")]%>%
  #   drop_na()%>%
  #   mutate(Sim_value=FCrmm-Dri_end)
  # #RMSE
  # RMSE<-with(Obs_sim_data,sqrt(sum((SoilWater-Sim_value)^2)/nrow(Obs_sim_data)))
  # R2<-with(Obs_sim_data,(sum((SoilWater-mean(SoilWater))*(Sim_value-mean(Sim_value)))/(sum((SoilWater-mean(SoilWater))^2)^0.5*sum((Sim_value-mean(Sim_value))^2)^0.5))^2)
  # b0<-with(Obs_sim_data,sum(Sim_value*SoilWater)/sum(SoilWater^2))
  #
  # stat<-c(RMSE,R2,b0)
  # names(stat)<-c("RMSE","R2","b0")

  #Plot result
  mytheme<-ggplot2::theme_bw()+
    ggplot2::theme(panel.grid=ggplot2::element_blank())
  PET<-ggplot2::ggplot(data = New_data)+
    ggplot2::geom_line(aes_(x=~Julian,y=~ET),color="black")+
    ggplot2::labs(y="PET (mm)")+
    mytheme
  PE<-ggplot2::ggplot(data = New_data)+
    ggplot2::geom_line(aes_(x=~Julian,y=~ETs),color="black")+
    ggplot2::labs(y="PE (mm)")+
    mytheme
  PT<-ggplot2::ggplot(data = New_data)+
    ggplot2::geom_line(aes_(x=~Julian,y=~ETc),color="black")+
    ggplot2::labs(y="PT (mm)")+
    mytheme
  Ea<-ggplot2::ggplot(data = New_data)+
    ggplot2::geom_line(aes_(x=~Julian,y=~Ea),color="black")+
    # geom_point(aes(x=Julian,y=EVAmm),color="black")+
    ggplot2::labs(y="Ea (mm)")+
    mytheme
  Ta<-ggplot2::ggplot(data = New_data)+
    ggplot2::geom_line(aes_(x=~Julian,y=~Ta),color="black")+
    ggplot2::labs(y="Ta (mm)")+
    mytheme
  ETa<-ggplot2::ggplot(data = New_data)+
    ggplot2::geom_line(aes_(x=~Julian,y=~(Ea+Ta)),color="black")+
    ggplot2::labs(y="ETa (mm)")+
    mytheme
  CR<-ggplot2::ggplot(data = New_data)+
    ggplot2::geom_line(aes_(x=~Julian,y=~CR),color="black")+
    ggplot2::labs(y="CR (mm)")+
    mytheme
  DP<-ggplot2::ggplot(data = New_data)+
    ggplot2::geom_line(aes_(x=~Julian,y=~(pe2+pd2)),color="black")+
    ggplot2::labs(y="DP (mm)")+
    mytheme
  SoilWater<-ggplot2::ggplot(data = New_data)+
    ggplot2::geom_line(aes_(x=~Julian,y=~THETA_2),color="black")+
    ggplot2::geom_point(aes_(x=~Julian,y=~SoilWater),color="black")+
    ggplot2::labs(y="The soil water (mm)")+
    mytheme

  plot<-ggpubr::ggarrange(PET,PE,PT,Ea,Ta,ETa,DP,CR,SoilWater)

  Result_data<-New_data

  # PET<-sum(New_data$ET,na.rm = T)
  # PE<-sum(New_data$Ea,na.rm = T)
  # PT<-sum(New_data$Ta,na.rm = T)
  # ETa<-sum(New_data$Ea+New_data$Ta,na.rm = T)
  # CR<-sum(New_data$CR,na.rm = T)
  # DP<-sum(New_data$pe2+New_data$pd2,na.rm = T)
  # num<-c(PET,PE,PT,ETa,CR,DP)
  # names(num)<-c("PET","PE","PT","ETa","CR","DP")
  #summary result
  summ<-New_data[-1,]%>%
    dplyr::summarise(Ea=sum(Ea),
              Ta=sum(Ta),
              CR=sum(CR),
              DP=sum(.data$pe2+.data$pd2))%>%
    dplyr::mutate(ETa=Ea+Ta)


  Result<-list(Result_data,plot,summ)
  names(Result)<-c("Result","Plot","Summary")
  return(Result)
}
