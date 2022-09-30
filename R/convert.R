#'@title Converting angert to radian
#'@description Converting the unit of angle in longitude and latitude into
#'    the unit of radian.
#'@param anger Longitude or dimension in Angle.
#'@return A vector for longitude or dimension in radian.
#'@export
#'@examples convert_angert_to_radian(98.8)

convert_angert_to_radian<-function(anger){
  radian<-(anger*pi/180)
  return(radian)
}

#'@title Convert date to day of year
#'@param Date is a date format data.
#'@return A vector for the day of year.
#'@export
#'@importFrom lubridate ymd
#'@importFrom lubridate year

convert_Date_to_dayofyear<-function(Date){
  Days<-as.numeric(Date-ymd(stringr::str_c(year(Date),1,1,sep="-"))+1)
  return(Days)
}


#'@title Convert degrees Celsius to Fahrenheit
#'@param degrees_Celsius temperature in degrees Celsius(°C).
#'@return A vector for atemperature in Fahrenheit(°F).
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

convert_degreesCelsius_to_Fahrenheit<-function(degrees_Celsius){
  Fahrenheit<-degrees_Celsius+273.16
  return(Fahrenheit)
}

#'@title Convert Fahrenheit to degrees Celsius
#'@param Fahrenheit temperature in Fahrenheit(°F).
#'@return A vector for temperature in degrees Celsius(°C)
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

convert_Fahrenheit_to_degreesCelsius<-function(Fahrenheit){
  degrees_Celsius<-Fahrenheit-273.16
  return(degrees_Celsius)
}

#'@title Convert radiation unit
#'@param rad Radiation data need to be converted from one unit to another unit.
#'@param type Used to specify how to convert.
#'@description Type has the following types:
#'  MJ_m2_day_to_J_cm2_day;
#'  MJ_m2_day_to_cal_cm2_day;
#'  MJ_m2_day_to_W_m2;
#'  cal_cm2_day_to_MJ_m2_day
#'  cal_cm2_day_to_J_cm2_day
#'  cal_cm2_day_to_W_m2
#'  cal_cm2_day_to_mm_day
#'  W_m2_to_MJ_m2_day
#'  W_m2_to_J_cm2_day
#'  W_m2_to_cal_cm2_day
#'  W_m2_to_mm_day
#'  mm_day_to_MJ_m2_day
#'  mm_day_to_J_cm2_day
#'  mm_day_to_cal_cm2_day
#'  mm_day_to_W_m2
#'@export
#'@return A vector for radiation converted unit
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

convert_Rad_unit<-function(rad,type){
  if(type=="MJ_m2_day_to_J_cm2_day"){
    newrad<-rad*100
  }else if (type=="MJ_m2_day_to_cal_cm2_day"){
    newrad<-rad*23.9
  }else if (type=="MJ_m2_day_to_W_m2"){
    newrad<-rad*11.6
  }else if (type=="MJ_m2_day_to_mm_day"){
    newrad<-rad*0.408
  }else if (type=="cal_cm2_day_to_MJ_m2_day"){
    newrad<-rad*4.1868*10^(-2)
  }else if (type=="cal_cm2_day_to_J_cm2_day"){
    newrad<-rad*4.1868
  }else if (type=="cal_cm2_day_to_W_m2"){
    newrad<-rad*0.485
  }else if (type=="cal_cm2_day_to_mm_day"){
    newrad<-rad*0.0171
  }else if (type=="W_m2_to_MJ_m2_day"){
    newrad<-rad*0.0864
  }else if (type=="W_m2_to_J_cm2_day"){
    newrad<-rad*8.64
  }else if (type=="W_m2_to_cal_cm2_day"){
    newrad<-rad*2.06
  }else if (type=="W_m2_to_mm_day"){
    newrad<-rad*0.035
  }else if (type=="mm_day_to_MJ_m2_day"){
    newrad<-rad*2.45
  }else if (type=="mm_day_to_J_cm2_day"){
    newrad<-rad*245
  }else if (type=="mm_day_to_cal_cm2_day"){
    newrad<-rad*58.5
  }else if (type=="mm_day_to_W_m2"){
    newrad<-rad*28.4
  }
  return(newrad)
}


#'@title Convert wind speed to the standard of 2m
#'@description For the calculation of evapotranspiration, wind. Speed measured
#'    at 2 m above the surface is required. To adjust wind speed data obtained
#'    from instruments placed at elevations other than the standard height of 2 m,
#'    a logarithmic wind speed profile may be used for measurements above a short grassed surface.
#'@param uz  measured wind speed at z m above ground surface [m s-1].
#'@param z height of measurement above ground surface [m].
#'@return  A vector for wind speed at 2 m above ground surface [m s-1].
#'@export
#'@references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M.
#'    FAO Irrigation and drainage paper No. 56. Rome: Food and Agriculture
#'    Organization of the United Nations, 1998.

convert_windSpeed_to_2m<-function(uz,z){
  u2<-uz*(4.87/log(67.8*z-5.42))
  return(u2)
}

#' #'@title Converting the every sheets in XLSX file to csv files
#' #'@param xlsx_file the xlsx file path
#' #'@export
#' #'@importFrom utils write.csv
#' #'@return No return value
#'
#' convert_xlsx_to_csv<-function(xlsx_file){
#'   # library(readxl)
#'   # library(stringr)
#'   #sheetsname
#'   Sheetnames<-readxl::excel_sheets(xlsx_file)
#'
#'   #xlsx2csv
#'
#'   for (i in Sheetnames) {
#'     data<-readxl::read_xlsx(xlsx_file,sheet=i)
#'     write.csv(data,file = stringr::str_c(i,".csv"),row.names = FALSE)
#'   }
#'
#' }


#----Create data file
#'@title Create a csv file or a dataframe in R to store the model data
#'@param TreNum Number. Need to generate how many treatment columun
#'@param rowNum Number. The number of row.
#'@details Latitude and Longitude use radians as units; Altitude use 'm' as units;
#'    Na use 'hour' as units;Tmax and Tmin use Celsius as units; Wind use m/s as
#'    units;RHmean and RHmin use percent sign as uits;Rs use MJ M-2 day-1 as units
#'    Height use cm as units;SoilWater and Irrigation use mm as units;
#'    GroundwaterDepth use cm as unit.
#'@export
#'@return A dataframe and a csv file (if to_CSV_file==TRUE)
#'@importFrom utils write.csv
#'@note The column of soilwater refers to the measured soil water (mm) in the
#'    maximum root layer, which is used to compare the difference between the
#'    simulated value and the measured value, which is an optional variable.
#'    The columns of Stage includes four stages: Ini, Development, Mid, End,
#'    which are mainly determined by LAI and growth status and can refer
#'    to Allen et al., (1998).

create_modelDF<-function(TreNum=1,rowNum=1){
  #Generate standard data frame
  Data<-matrix(ncol = 20+TreNum,nrow = rowNum)#How many treatment are there
  colnames(Data)<-c("Year",paste("Tre",1:TreNum,sep = ""),"Julian","Cut","Latitude","Longitude","Altitude",
                    "Precipitation","Na","Tmax","Tmin","Wind","RHmean","RHmin","Rs",
                    "Height","LAI","Stage","SoilWater","Irrigation","GroundwaterDepth")
  Data<-as.data.frame(Data)
  # #Export it as CSV file
  # if(to_CSV_file==TRUE){
  #   write.csv(Data,row.names = FALSE)
  # }
  return(Data)
}


