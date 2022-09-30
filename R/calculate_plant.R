#'@title Calculating water stress coefficient
#'@param Dr root zone depletion(mm).
#'@param TAW total available soil water in the root zone(mm).
#'@param p fraction of TAW that a crop can extract from the root zone without
#'    suffering water stress.
#'@export
#'@return A value for water stress coefficient which is a dimensionless transpiration reduction factor
#'    dependent on available soil water

cal_WaterStressCoef<-function(Dr,TAW,p){
  RAW<-p*TAW
  if(Dr>RAW){
    WaterStressCoef<-(TAW-Dr)/((1-p)*TAW)
  }else{
    WaterStressCoef<-1
  }
  return(WaterStressCoef)
}



