
#'@title Calculating capillary rise
#'@param a1 Soil water storage to maximum root depth at field capacity(mm).
#'@param b1 A parameter.
#'@param a2 storage above the average between those at field capacity and the wilting point(mm).
#'@param b2 A parameter.
#'@param ETm potential crop evaporanspiration (mm/day),usually ETm=ETc(mm/d).
#'@param a3 A parameter.
#'@param b3 A parameter.6.7 for clay and silty clay loam soils, decreasing to 6.2 for loamy sands
#'@param a4 A parameter.4.6 for silty loam and silty clay loam soils,decreasing to 6.2 for loamy sands.
#'@param b4 A parameter.-0.65 for silty loam soils and decreasing to -2.5 for loamy sand soils.
#'@param Dw Groudwater depth below root zone(m).
#'@param Wa actual soil water storage in the root zone.
#'@param LAI Leaf area index.
#'@return The value for capillary Rise (mm/day).
#'@export
#'@references Liu Y, Pereira L S, Fernando R M. Fluxes through the bottom boundary of the root zone in silty soils:
#'  Parametric approaches to estimate groundwater contribution and percolation[J].
#'  Agricultural Water Management, 2006, 84(1):27-40.

cal_capillaryRise<-function(a1,b1=-0.17,a2,b2=-0.27,a3=-1.3,b3,a4,b4,Dw,Wa,LAI,ETm){
  #----1.Wc=Critical soil water storage(mm)
  Wc<-a1*Dw^b1

  #----2.Ws=Steady soil water storage(mm)
  if(Dw<=3){
    Ws<-a2*Dw^b2 #
  }else{
    Ws<-240 #mm
  }

  #----3.Dwc=critical groundwater depth
  if(ETm<=4){
    Dwc<-a3*ETm+b3
  }else{
    Dwc<-1.4
  }

  #----4.k=factor relating evapotranspiration with transpiration(no-dimensional)
  if(ETm<=4){
    k<-1-exp(-0.6*LAI)
  }else{
    k<-3.8/ETm
  }

  #----5.CRmax=potential capillary flux(mm/d)
  if(Dw<=Dwc){
    CRmax<-k*ETm
  }else{
    CRmax<-a4*Dw^b4
  }

  #----6.CR
  if(Wa<Ws){
    CR<-CRmax
  }else if(Wa>=Ws&Wa<=Wc){
    CR<-CRmax*((Wc-Wa)/(Wc-Ws))
  }else if(Wa>Wc){
    CR<-0
  }
  return(CR)
}


#'@title Calculating Deep percolation
#'@param a A water storage value comprised between WFc and Wa at saturation.
#'@param b b<-0.0173 for soils draining quickly.Otherwise b>-0.0173.
#'@param Wa actual soil water storage in the root zone (mm)
#'@param Wfc soil water storage to maximum root depth (Zr ) at field capacity (mm)
#'@param t time after an irrigation or rain that produced a storage above field capacity (days)
#'@return A vector for deep percolation(mm/day).
#'@export
#'@references Liu Y, Pereira L S, Fernando R M. Fluxes through the bottom boundary of the root zone in silty soils:
#'  Parametric approaches to estimate groundwater contribution and percolation[J].
#'  Agricultural Water Management, 2006, 84(1):27-40.

cal_DeepPercolation<-function(Wa,Wfc,a,b,t){
  while (Wa>Wfc) {
    DeepPercolation<-a*t^b
  }
  return(DeepPercolation)
}
