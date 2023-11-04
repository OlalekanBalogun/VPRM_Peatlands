##########################################################################VPRM_SIF###########################################################
# Set the working directory to the directory of the source file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(devtools)
#install_github('Timothy-W-Hilton/VPRMLandSfcModel')

#Load VPRMLandSfcModle as an R library/package
library('VPRMLandSfcModel')

#The vignette document 'VPRM_vignette', included with the package, provides several "getting started" usage examples
#RShowDoc('VPRM_vignette', package='VPRMLandSfcModel') #opens in browser as html doc

#The 'VPRM_vignette' files can also be found in the code folder under inst/doc as .html and .R files

#Read in csv files of tower obs, vegetation index and reflectances to assemble VPRM driver data
ARB_NARR_met.SIF <- read.csv("ARB-NARR-Final-ready.csv", stringsAsFactors = FALSE)
ARB_evi_modis.SIF <- read.csv("ARB-SIF-ready.csv", stringsAsFactors = FALSE)
ARB_refl_modis.SIF <- read.csv("ARB-Reflectance-ready.csv", stringsAsFactors = FALSE)
ARB_phen_modis.SIF <- read.csv("ARB-Phenology-MAX-ready.csv", stringsAsFactors = FALSE)

#Convert "Character" dates to "Chron" dates and times for each VPRM driver data
library(chron)

ARB_NARR_met.SIF.dates <- t(as.data.frame(strsplit(ARB_NARR_met.SIF[["Date"]],' '))) #t is (matrix transpose)
row.names(ARB_NARR_met.SIF.dates) = NULL
ARB_NARR_met.SIF ["Date"] <- chron(dates=ARB_NARR_met.SIF.dates[,1],times=ARB_NARR_met.SIF.dates[,2], format = c(dates = "m/d/y", times = "h:m:s"))
class(ARB_NARR_met.SIF[["Date"]])

ARB_evi_modis.SIF.dates <- t(as.data.frame(strsplit(ARB_evi_modis.SIF[["Date"]],' '))) #t is (matrix transpose)
row.names(ARB_evi_modis.SIF.dates) = NULL
ARB_evi_modis.SIF ["Date"] <- chron(dates=ARB_evi_modis.SIF.dates[,1],times=ARB_evi_modis.SIF.dates[,2], format = c(dates = "m/d/y", times = "h:m:s"))
class(ARB_evi_modis.SIF[["Date"]])

ARB_refl_modis.SIF.dates <- t(as.data.frame(strsplit(ARB_refl_modis.SIF[["Date"]],' '))) #t is (matrix transpose)
row.names(ARB_refl_modis.SIF.dates) = NULL
ARB_refl_modis.SIF ["Date"] <- chron(dates=ARB_refl_modis.SIF.dates[,1],times=ARB_refl_modis.SIF.dates[,2], format = c(dates = "m/d/y", times = "h:m:s"))
class(ARB_refl_modis.SIF[["Date"]])

####Melt and rearrange Phenology data and Convert "Character" dates to "Chron" dates and times for each VPRM driver data
ARB_phen.melted.SIF <- reshape2::melt(ARB_phen_modis.SIF, id.var='Date')
names(ARB_phen.melted.SIF) <- c("var", "phen", "date")
head(ARB_phen.melted.SIF)
tail(ARB_phen.melted.SIF)

ARB_phen.melted.SIF.dates <- t(as.data.frame(strsplit(ARB_phen.melted.SIF[["date"]],' '))) #t is (matrix transpose)
row.names(ARB_phen.melted.SIF.dates) = NULL
ARB_phen.melted.SIF ["date"] <- chron(dates=ARB_phen.melted.SIF.dates[,1],times=ARB_phen.melted.SIF.dates[,2], format = c(dates = "m/d/y", times = "h:m:s"))
class(ARB_phen.melted.SIF[["date"]])

###Subset the relevant columns for aggregation
phen.varselect.SIF<- c("date", "phen")
ARB_phen.subset.SIF <- ARB_phen.melted.SIF[phen.varselect.SIF]

ARB_phen_rearranged.SIF <- ARB_phen.subset.SIF[order(as.Date(ARB_phen.subset.SIF$date, format="%m/%d/%y")),]
head(ARB_phen_rearranged.SIF)
tail(ARB_phen_rearranged.SIF)

## determine the phenology phase for each tower observation date
phen_filled.SIF <- interp_phenology(ARB_phen_rearranged.SIF, ARB_NARR_met.SIF[['Date']])



#####Read csv file from original optimized run####
specify_path <- (dirname(rstudioapi::getActiveDocumentContext()$path))
filename_parameters <- paste(specify_path, "/ARB_Optimized_Site_Parameter_Values",".csv",sep="")
SIF_optimized_parameters <- read.csv(filename_parameters)


######Assemble VPRM driver data frame

## determine the phenology phase for each tower observation date
#phen_filled.SIF <- interp_phenology(ARB_phen, ARB_NARR_met.SIF[['date']])
## place the tower observations and MODIS data into a VPRM_driver_data object.
ARB_dd.SIF <- VPRM_driver_data(name_long="Attawapiskat ON",
                               name_short = "ARB",
                               lat=52.695,
                               lon=-83.9452,
                               PFT='WET',    ## Permanent Wetland
                               tower_date=ARB_NARR_met.SIF[['Date']],
                               NEE_obs=NA,
                               T=ARB_NARR_met.SIF[['TA']],
                               PAR=ARB_NARR_met.SIF[['PAR']],
                               date_nir = ARB_refl_modis.SIF[['Date']],
                               rho_nir=ARB_refl_modis.SIF[['NIR']],
                               date_swir = ARB_refl_modis.SIF[['Date']],
                               rho_swir = ARB_refl_modis.SIF[['SWIR']],
                               date_EVI = ARB_evi_modis.SIF[['Date']],
                               EVI=ARB_evi_modis.SIF[['SIF']],
                               phen=phen_filled.SIF)

## take a look at the result
print(head(as.data.frame(ARB_dd.SIF)))
print(tail(as.data.frame(ARB_dd.SIF)))


######Calculate VPRM NEE#####
#Calculating VPRM fluxes requires values for the parameters lambda, PAR0, alpha, and beta.
#The dataset VPRM_parameters, part of the VPRMLandSfcModel package, includes several parameters sets.

attach(SIF_optimized_parameters)

ARB_dd.SIF[['data']][['NEE_VPRM']] <- vprm_calc_NEE(
  ARB_dd.SIF, lambda=lambda, PAR_0=PAR_0, alpha=alpha, beta=beta)

######Calculate VPRM GPP#####
ARB_dd.SIF[['data']][['GPP_VPRM']] <- vprm_calc_GEE(
  ARB_dd.SIF, lambda=lambda, PAR_0=PAR_0)

######Calculate VPRM ER#####
ARB_dd.SIF[['data']][['ER_VPRM']] <- vprm_calc_R(
  ARB_dd.SIF, alpha=alpha, beta=beta)


###plot the VPRM NEE
library(ggplot2)
library(chron)



##################################################Diagnosis#############################################
VPRM_Fluxes.SIF <- data.frame(Date = ARB_dd.SIF[["data"]][["date"]], NEE_VPRM = ARB_dd.SIF[["data"]][["NEE_VPRM"]], GPP_VPRM = ARB_dd.SIF[["data"]][["GPP_VPRM"]], ER_VPRM = ARB_dd.SIF[["data"]][["ER_VPRM"]], Temperature = ARB_dd.SIF[["data"]][["T"]], PAR = ARB_dd.SIF[["data"]][["PAR"]])
#NEE.Diff.SIF.2 <- data.frame(Date = fig_dNEE.SIF[["data"]][["date"]], NEE_Residual = fig_dNEE.SIF[["data"]][["NEE_VPRM"]] - fig_dNEE.SIF[["data"]][["NEE_obs"]])

###Subset relevant columns for NARR meteorological variables needed for analysis
variables.NARR.met <- c("Shore", "WindSpeed", "Precipitation","Humidity", "Evaporation", "SoilMoisture")
NARR.met.subset <- ARB_NARR_met.SIF[variables.NARR.met]

###Combine NEE data and NARR met data##
VPRM_Fluxes.Met.All <- cbind(VPRM_Fluxes.SIF,NARR.met.subset)


#Calculate Vapour Pressure Deficit (VPD) from Temperature (degrees Celsius) and Relative Humidity (%)
#VPRM_Fluxes.Met.All <- VPRM_Fluxes.Met.Study
#Calculate Saturation Vapor Pressure (es): es = 0.6108 * exp(17.27 * T / (T + 237.3))
#VPRM_Fluxes.Met.All [["es"]] <- 0.6108 * exp(17.27 * VPRM_Fluxes.Met.All$Temperature / (VPRM_Fluxes.Met.All$Temperature + 237.3))
#Calculate Actual Vapor Pressure (ea): ea = RH / 100 * es 
#VPRM_Fluxes.Met.All [["ea"]] <- VPRM_Fluxes.Met.All$Humidity / 100 * VPRM_Fluxes.Met.All$es
#Calculate Vapor Pressure Deficit: VPD = es - ea
#VPRM_Fluxes.Met.All [["VPD_1"]] <- (VPRM_Fluxes.Met.All$es - VPRM_Fluxes.Met.All$ea)
#VPRM_Fluxes.Met.All$es <- NULL #Drop or remove column
#VPRM_Fluxes.Met.All$ea <- NULL #Drop or remove column
VPRM_Fluxes.Met.All [["VPD"]] <- ((0.6108 * exp(17.27 * VPRM_Fluxes.Met.All$Temperature / (VPRM_Fluxes.Met.All$Temperature + 237.3))) - (VPRM_Fluxes.Met.All$Humidity / 100 * (0.6108 * exp(17.27 * VPRM_Fluxes.Met.All$Temperature / (VPRM_Fluxes.Met.All$Temperature + 237.3)))))
#VPRM_Fluxes.Met.All [["VPD_2"]] <- ((610.78 * exp(VPRM_Fluxes.Met.All$Temperature / (VPRM_Fluxes.Met.All$Temperature + 238.3) * 17.2694)) /1000) * (1 - VPRM_Fluxes.Met.All$Humidity/100) #Divide by 1000 to convert SVP from pascals to kPa


####Exclude all other months (e.g. May and November) just in case they are in the dataset)####
VPRM_Fluxes.Met.Study <- VPRM_Fluxes.Met.All #Duplicate dataframe to avoid potential conflicts
VPRM_Fluxes.Met.Study ["Month"] <- as.data.frame(months.Date(VPRM_Fluxes.Met.Study$Date)) #Create and add months column
VPRM_Fluxes.Met.Study <- subset(VPRM_Fluxes.Met.Study, Month=="June" | Month=="July" | Month=="August" | Month=="September" | Month=="October")
VPRM_Fluxes.Met.Study$Month <- NULL #Drop or remove the months column

detach(SIF_optimized_parameters)


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

###############################################################PLOT GRAPHS##################################################################################

library(gdata)
require(reshape2)
require(ggplot2)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(scales)

##################################################################Daily Averaged Analysis################################################################
##Subset the relevant columns for aggregation and convert month names to month numbers

#Create and add new columns for day, month, year, and time; extracted from the original date column
VPRM_Fluxes.Daily <- VPRM_Fluxes.Met.Study #Duplicate dataframe to avoid potential conflicts

VPRM_Fluxes.Daily ["Day"] <- as.data.frame(day(VPRM_Fluxes.Daily$Date))
VPRM_Fluxes.Daily ["Month"] <- as.data.frame(months.Date(VPRM_Fluxes.Daily$Date))
VPRM_Fluxes.Daily ["Year"] <- as.data.frame(year(VPRM_Fluxes.Daily$Date))
#VPRM_Fluxes.Daily ["Hours"] <- as.data.frame(hours(VPRM_Fluxes.Daily$Date))
#VPRM_Fluxes.Daily ["Minutes"] <- as.data.frame(minutes(VPRM_Fluxes.Daily$Date))
#VPRM_Fluxes.Daily ["Seconds"] <- as.data.frame(seconds(VPRM_Fluxes.Daily$Date))
VPRM_Fluxes.Daily ["Time"] <- data.frame(substr(VPRM_Fluxes.Daily[["Date"]], start = 11, stop = 15))

variables.fluxes.daily <- c("NEE_VPRM", "GPP_VPRM", "ER_VPRM", "Day", "Month")
VPRM_Fluxes.Daily.subset <- VPRM_Fluxes.Daily[variables.fluxes.daily]
VPRM_Fluxes.Daily.subset$Month = match(VPRM_Fluxes.Daily.subset$Month,(month.name)) #Convert month names to month numbers before agg for proper ordering
VPRM_Fluxes.Daily.Mean <- aggregate(.~ Day + Month, data = VPRM_Fluxes.Daily.subset, mean, na.action = na.pass)
#VPRM_Fluxes.Daily.Mean <- aggregate(.~ Day + Month + Year, data = VPRM_Fluxes.Daily.subset, mean, na.action = na.pass)
#VPRM_Fluxes.Daily.Mean.Line <- VPRM_Fluxes.Daily.Mean #Duplicate dataframe to avoid conflict and use for other graphs

VPRM_Fluxes.Daily.Mean ["Char_Date"] <- paste( month.abb[VPRM_Fluxes.Daily.Mean$Month], VPRM_Fluxes.Daily.Mean$Day, sep="-" )
VPRM_Fluxes.Daily.Mean ["Date"] <- as.Date(VPRM_Fluxes.Daily.Mean$Char_Date, format="%b-%d")
#VPRM_Fluxes.Daily.Mean ["Date"] <- as.data.frame(paste(VPRM_Fluxes.Daily.Mean$Month, VPRM_Fluxes.Daily.Mean$Day, sep="-") %>% ymd() %>% as.Date()) #Combine dates column into single date
VPRM_Fluxes.Daily.Mean$Char_Date <- NULL #Drop or remove the character Date column 
VPRM_Fluxes.Daily.Mean$Day <- NULL #Drop or remove the Day column 
VPRM_Fluxes.Daily.Mean$Month <- NULL #Drop or remove the Month column
#VPRM_Fluxes.Daily.Mean$Year <- NULL #Drop or remove the Year column

##Melt to prepare data form ggplot multiple variables
VPRM_Fluxes.Daily.Mean.Melted <- reshape2::melt(VPRM_Fluxes.Daily.Mean, id.var='Date')
head(VPRM_Fluxes.Daily.Mean.Melted)
tail(VPRM_Fluxes.Daily.Mean.Melted)

#Write VPRM Daily Mean Fluxes to csv file
write.csv(VPRM_Fluxes.Daily.Mean, file = "VPRM_Fluxes_Daily_Mean.csv", row.names = F)

#Move the excel .csv files to the stats folder
#library(filesstrings)
#file.move("VPRM_Fluxes_Daily_Mean.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/ARB_FINAL/Stats")



##################################################################Annual Averaged Analysis################################################################
##Subset the relevant columns for aggregation
VPRM_Fluxes.Annual <- VPRM_Fluxes.Daily #Copy previously used dataframe from above to avoid potential conflicts
variables.fluxes.annual <- c("NEE_VPRM", "GPP_VPRM", "ER_VPRM", "Year")
VPRM_Fluxes.Annual.subset <- VPRM_Fluxes.Annual[variables.fluxes.annual]
VPRM_Fluxes.Annual.Mean <- aggregate(.~ Year, data = VPRM_Fluxes.Annual.subset, mean, na.action = na.pass)
#VPRM_Fluxes.Annual.Mean.Line <- VPRM_Fluxes.Annual.Mean #Duplicate dataframe to avoid conflict and use for other graphs
VPRM_Fluxes.Annual.Mean ["Year_Char"] <- as.character(VPRM_Fluxes.Annual.Mean$Year) #Year column is in class ordered or factor. Need to change to Numeric
VPRM_Fluxes.Annual.Mean$Year <- NULL #Drop or remove the original ordered/factor Year column
VPRM_Fluxes.Annual.Mean ["Year"] <- as.numeric(VPRM_Fluxes.Annual.Mean$Year_Char)
VPRM_Fluxes.Annual.Mean$Year_Char <- NULL #Drop or remove the Character Year column 


##Melt to prepare data form ggplot multiple variables
VPRM_Fluxes.Annual.Mean.Melted <- reshape2::melt(VPRM_Fluxes.Annual.Mean, id.var='Year')
head(VPRM_Fluxes.Annual.Mean.Melted)
tail(VPRM_Fluxes.Annual.Mean.Melted)



###########Create plots of the annual averaged flux for the growing season
#plot <- ggplot(data = NEE.residual.EVI.Annual, aes(x=Year, y=NEE_VPRM, group=1)) + geom_line() + geom_point(shape = 1, size = 3)
#plot <- plot + geom_line(data = NEE.residual.EVI.Annual, aes(y = NEE_obs), size = 1.2) + geom_point(aes(y = NEE_obs), shape = 15, size = 3)
plot <- ggplot(VPRM_Fluxes.Annual.Mean.Melted, aes(x=Year, y=value, group=variable)) + geom_line(aes(color=variable), size = 1) + geom_point(aes(color=variable), size = 3)
plot <- plot + xlab("Year") + ylab(expression(CO[2]~Fluxes~(mu*mol~m^{-2}~s^{-1}))) + ggtitle(paste("ARB Mean Annual Fluxes" ,"(June-October, 2000-2019)"))
plot <- plot + scale_x_continuous(breaks = seq(from=2000,to=2020,by=2))
#plot <- plot + scale_x_continuous(labels =c("June","July","August","September","October"))
plot <- plot + theme(plot.title = element_text(hjust = 0.5)) #center title
plot <- plot + theme(axis.text=element_text(size=11, face="bold"))
plot <- plot + theme(axis.title=element_text(size=11))
plot <- plot + theme(plot.title=element_text(size=11,face="bold"))
#plot <- plot + theme(legend.key.size = unit(1.0, "cm"),legend.key.width = unit(1.0,"cm")) 
plot <- plot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank())
plot <- plot + labs(shape = "") #Remove or specify Legend title
plot <- plot + theme(legend.title = element_blank()) #Remove Legend title
plot <- plot + theme(legend.text = element_text(size = 11, face = "bold"))
plot <- plot + theme(legend.position = "bottom")
plot <- plot + geom_hline(yintercept=0, color="black")
plot <- plot + scale_y_continuous(breaks=seq(-3,4.5,1.5)) #Sets the axis range and interval ("by")
plot <- plot + coord_cartesian(ylim = c(-3, 4.5), xlim = c(2000, 2020)) #Forces the limits to be set accordingly
#plot
#save graph to output as high resoultion, Windows compatible format .EMF
fileout <- paste("Line-Final-ARB Mean Annual Fluxes","_LineCycle",".emf",sep="")
x11(width=7, height=4) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
dev.cur()
plot
savePlot(fileout, "emf", device = dev.cur())
dev.off()
graphics.off()

#Write VPRM Annual Mean Fluxes to csv file
write.csv(VPRM_Fluxes.Annual.Mean, file = "VPRM_Fluxes_Annual_Mean.csv", row.names = F)

#Move the excel .csv files to the stats folder
library(filesstrings)
#file.move("VPRM_Fluxes_Annual_Mean.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/ARB_FINAL/Stats")



##########Convert Annual Averaged to g C m^-2#########################
##Duplicate data frame from annual averaged analysis above
VPRM_Fluxes.Annual.Mean.Melted.gC <- VPRM_Fluxes.Annual.Mean.Melted
VPRM_Fluxes.Annual.Mean.Melted.gC [["value_gC"]] <- VPRM_Fluxes.Annual.Mean.Melted.gC$value * (86400 * 153 * 12.011) / 1000000 #Convert from umols of CO2 m^-2 s^-1 to g C m^-2. Note that 153 is the number of days in the study period per year.
VPRM_Fluxes.Annual.Mean.Melted.gC$value <- NULL #Drop or remove the original value column 

##Melt to prepare data form ggplot multiple variables
head(VPRM_Fluxes.Annual.Mean.Melted.gC)
tail(VPRM_Fluxes.Annual.Mean.Melted.gC)



###########Create plots of the annual averaged flux in g C m-2 for the growing season
#plot <- ggplot(data = NEE.residual.EVI.Annual, aes(x=Year, y=NEE_VPRM, group=1)) + geom_line() + geom_point(shape = 1, size = 3)
#plot <- plot + geom_line(data = NEE.residual.EVI.Annual, aes(y = NEE_obs), size = 1.2) + geom_point(aes(y = NEE_obs), shape = 15, size = 3)
plot <- ggplot(VPRM_Fluxes.Annual.Mean.Melted.gC, aes(x=Year, y=value_gC, group=variable)) + geom_line(aes(color=variable), size = 1) + geom_point(aes(color=variable), size = 3)
plot <- plot + xlab("Year") + ylab(expression(Carbon~Fluxes~(g~C~m^{-2}))) + ggtitle(paste("ARB Annual Fluxes" ,"(June-October, 2000-2019)"))
plot <- plot + scale_x_continuous(breaks = seq(from=2000,to=2020,by=2))
#plot <- plot + scale_x_continuous(labels =c("June","July","August","September","October"))
plot <- plot + theme(plot.title = element_text(hjust = 0.5)) #center title
plot <- plot + theme(axis.text=element_text(size=11, face="bold"))
plot <- plot + theme(axis.title=element_text(size=11))
plot <- plot + theme(plot.title=element_text(size=11,face="bold"))
#plot <- plot + theme(legend.key.size = unit(1.0, "cm"),legend.key.width = unit(1.0,"cm")) 
plot <- plot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank())
plot <- plot + labs(shape = "") #Remove or specify Legend title
plot <- plot + theme(legend.title = element_blank()) #Remove Legend title
plot <- plot + theme(legend.text = element_text(size = 11, face = "bold"))
plot <- plot + theme(legend.position = "bottom")
plot <- plot + geom_hline(yintercept=0, color="black")
plot <- plot + scale_y_continuous(breaks=seq(-300,450,150)) #Sets the axis range and interval ("by")
plot <- plot + coord_cartesian(ylim = c(-300, 450), xlim = c(2000, 2020)) #Forces the limits to be set accordingly
#plot
#save graph to output as high resoultion, Windows compatible format .EMF
fileout <- paste("Line-Final-ARB Annual gC Fluxes","_LineCycle",".emf",sep="")
x11(width=7, height=4) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
dev.cur()
plot
savePlot(fileout, "emf", device = dev.cur())
dev.off()
graphics.off()

#Write VPRM Annual Mean Fluxes in g C to csv file
VPRM_Fluxes.Annual.Mean.gC <- dcast(VPRM_Fluxes.Annual.Mean.Melted.gC, Year ~ variable)

write.csv(VPRM_Fluxes.Annual.Mean.gC, file = "VPRM_Fluxes_Annual_gC.csv", row.names = F)

#Move the excel .csv files to the stats folder
#library(filesstrings)
#file.move("VPRM_Fluxes_Annual_gC.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/ARB_FINAL/Stats")






##########################################################Multiple Regression Analysis - GPP#############################################################

##########Daily################
variables.GPP_VPRM_MultipleReg.daily <- c("GPP_VPRM", "Temperature", "PAR", "Precipitation", "Evaporation", "SoilMoisture",  "VPD", "Day", "Month")
GPP_VPRM_MultipleReg.Daily.subset <- VPRM_Fluxes.Daily[variables.GPP_VPRM_MultipleReg.daily]
GPP_VPRM_MultipleReg.Daily.Mean <- aggregate(.~ Day + Month, data = GPP_VPRM_MultipleReg.Daily.subset, mean, na.action = na.pass)
GPP_VPRM_MultipleReg.Daily.Mean$Day <- NULL #Drop or remove the Day column
GPP_VPRM_MultipleReg.Daily.Mean$Month <- NULL #Drop or remove the Month column
GPP_VPRM_MultipleReg.Daily.Mean$Precipitation <- NULL #Drop or remove the mean Precipitation column. Sum needed
GPP_VPRM_MultipleReg.Daily.Mean$Evaporation <- NULL #Drop or remove the mean Evaporation column. Sum needed
temporary_SoilMoisture.Daily <- GPP_VPRM_MultipleReg.Daily.Mean$SoilMoisture #Create temporary dataframe for Soil Moisture for subsequent column rearrangement
temporary_VPD.Daily <- GPP_VPRM_MultipleReg.Daily.Mean$VPD #Create temporary dataframe for VPD for subsequent column rearrangement
GPP_VPRM_MultipleReg.Daily.Mean$SoilMoisture <- NULL #Drop or remove the Soil Moisture column. Will be re-added later.
GPP_VPRM_MultipleReg.Daily.Mean$VPD <- NULL #Drop or remove the VPD column. Will be re-added later.

variables.GPP_VPRM_MultipleReg.daily.2 <- c("Precipitation", "Evaporation", "Day", "Month", "Year")
GPP_VPRM_MultipleReg.Daily.subset.2 <- VPRM_Fluxes.Daily[variables.GPP_VPRM_MultipleReg.daily.2]
GPP_VPRM_MultipleReg.Daily.Sum <- aggregate(.~ Day + Month + Year, data = GPP_VPRM_MultipleReg.Daily.subset.2, sum, na.action = na.pass)
GPP_VPRM_MultipleReg.Daily.Sum$Year <- NULL #Drop or remove the Year column
GPP_VPRM_MultipleReg.Daily.Mean.2 <- aggregate(.~ Day + Month, data = GPP_VPRM_MultipleReg.Daily.Sum, mean, na.action = na.pass)
GPP_VPRM_MultipleReg.Daily.Mean ["Precipitation"] <- GPP_VPRM_MultipleReg.Daily.Mean.2$Precipitation
GPP_VPRM_MultipleReg.Daily.Mean ["Evaporation"] <- GPP_VPRM_MultipleReg.Daily.Mean.2$Evaporation
GPP_VPRM_MultipleReg.Daily.Mean ["SoilMoisture"] <- temporary_SoilMoisture.Daily #Re-add Soil Moisture column
GPP_VPRM_MultipleReg.Daily.Mean ["VPD"] <- temporary_VPD.Daily #Re-add VPD column

#Stepwise Regression
#Performs stepwise model selection by AIC.
library(MASS)
MultipleReg.model.Daily.GPP <- lm(GPP_VPRM ~., data = GPP_VPRM_MultipleReg.Daily.Mean)
MultipleReg.model.Daily.GPP <- stepAIC(MultipleReg.model.Daily.GPP, direction="both")
MultipleReg.model.Daily.GPP$anova # display results
summary(MultipleReg.model.Daily.GPP)
summary(MultipleReg.model.Daily.GPP)$coefficient


##########Annual###############
variables.GPP_VPRM_MultipleReg.annual <- c("GPP_VPRM", "Temperature", "PAR", "Precipitation", "Evaporation", "SoilMoisture",  "VPD", "Year")
GPP_VPRM_MultipleReg.Annual.subset <- VPRM_Fluxes.Annual[variables.GPP_VPRM_MultipleReg.annual]
GPP_VPRM_MultipleReg.Annual.Mean <- aggregate(.~ Year, data = GPP_VPRM_MultipleReg.Annual.subset, mean, na.action = na.pass)
GPP_VPRM_MultipleReg.Annual.Mean$Year <- NULL #Drop or remove the Day column
GPP_VPRM_MultipleReg.Annual.Mean$Precipitation <- NULL #Drop or remove the mean Precipitation column. Sum needed
GPP_VPRM_MultipleReg.Annual.Mean$Evaporation <- NULL #Drop or remove the mean Evaporation column. Sum needed
temporary_SoilMoisture.Annual <- GPP_VPRM_MultipleReg.Annual.Mean$SoilMoisture #Create temporary dataframe for Soil Moisture for subsequent column rearrangement
temporary_VPD.Annual <- GPP_VPRM_MultipleReg.Annual.Mean$VPD #Create temporary dataframe for VPD for subsequent column rearrangement
GPP_VPRM_MultipleReg.Annual.Mean$SoilMoisture <- NULL #Drop or remove the Soil Moisture column. Will be re-added later.
GPP_VPRM_MultipleReg.Annual.Mean$VPD <- NULL #Drop or remove the VPD column. Will be re-added later.

variables.GPP_VPRM_MultipleReg.annual.2 <- c("Precipitation", "Evaporation", "Year")
GPP_VPRM_MultipleReg.Annual.subset.2 <- VPRM_Fluxes.Annual[variables.GPP_VPRM_MultipleReg.annual.2]
GPP_VPRM_MultipleReg.Annual.Sum <- aggregate(.~ Year, data = GPP_VPRM_MultipleReg.Annual.subset.2, sum, na.action = na.pass)
GPP_VPRM_MultipleReg.Annual.Mean ["Precipitation"] <- GPP_VPRM_MultipleReg.Annual.Sum$Precipitation
GPP_VPRM_MultipleReg.Annual.Mean ["Evaporation"] <- GPP_VPRM_MultipleReg.Annual.Sum$Evaporation
GPP_VPRM_MultipleReg.Annual.Mean ["SoilMoisture"] <- temporary_SoilMoisture.Annual #Re-add Soil Moisture column
GPP_VPRM_MultipleReg.Annual.Mean ["VPD"] <- temporary_VPD.Annual #Re-add VPD column

#Stepwise Regression
#Performs stepwise model selection by AIC.
MultipleReg.model.Annual.GPP <- lm(GPP_VPRM ~., data = GPP_VPRM_MultipleReg.Annual.Mean)
MultipleReg.model.Annual.GPP <- stepAIC(MultipleReg.model.Annual.GPP, direction="both")
MultipleReg.model.Annual.GPP$anova # display results
summary(MultipleReg.model.Annual.GPP)
summary(MultipleReg.model.Annual.GPP)$coefficient



##########################################################Multiple Regression Analysis - ER#############################################################

##########Daily################
variables.ER_VPRM_MultipleReg.daily <- c("ER_VPRM", "Temperature", "PAR", "Precipitation", "Evaporation", "SoilMoisture",  "VPD", "Day", "Month")
ER_VPRM_MultipleReg.Daily.subset <- VPRM_Fluxes.Daily[variables.ER_VPRM_MultipleReg.daily]
ER_VPRM_MultipleReg.Daily.Mean <- aggregate(.~ Day + Month, data = ER_VPRM_MultipleReg.Daily.subset, mean, na.action = na.pass)
ER_VPRM_MultipleReg.Daily.Mean$Day <- NULL #Drop or remove the Day column
ER_VPRM_MultipleReg.Daily.Mean$Month <- NULL #Drop or remove the Month column
ER_VPRM_MultipleReg.Daily.Mean$Precipitation <- NULL #Drop or remove the mean Precipitation column. Sum needed
ER_VPRM_MultipleReg.Daily.Mean$Evaporation <- NULL #Drop or remove the mean Evaporation column. Sum needed
temporary_SoilMoisture.Daily <- ER_VPRM_MultipleReg.Daily.Mean$SoilMoisture #Create temporary dataframe for Soil Moisture for subsequent column rearrangement
temporary_VPD.Daily <- ER_VPRM_MultipleReg.Daily.Mean$VPD #Create temporary dataframe for VPD for subsequent column rearrangement
ER_VPRM_MultipleReg.Daily.Mean$SoilMoisture <- NULL #Drop or remove the Soil Moisture column. Will be re-added later.
ER_VPRM_MultipleReg.Daily.Mean$VPD <- NULL #Drop or remove the VPD column. Will be re-added later.

variables.ER_VPRM_MultipleReg.daily.2 <- c("Precipitation", "Evaporation", "Day", "Month", "Year")
ER_VPRM_MultipleReg.Daily.subset.2 <- VPRM_Fluxes.Daily[variables.ER_VPRM_MultipleReg.daily.2]
ER_VPRM_MultipleReg.Daily.Sum <- aggregate(.~ Day + Month + Year, data = ER_VPRM_MultipleReg.Daily.subset.2, sum, na.action = na.pass)
ER_VPRM_MultipleReg.Daily.Sum$Year <- NULL #Drop or remove the Year column
ER_VPRM_MultipleReg.Daily.Mean.2 <- aggregate(.~ Day + Month, data = ER_VPRM_MultipleReg.Daily.Sum, mean, na.action = na.pass)
ER_VPRM_MultipleReg.Daily.Mean ["Precipitation"] <- ER_VPRM_MultipleReg.Daily.Mean.2$Precipitation
ER_VPRM_MultipleReg.Daily.Mean ["Evaporation"] <- ER_VPRM_MultipleReg.Daily.Mean.2$Evaporation
ER_VPRM_MultipleReg.Daily.Mean ["SoilMoisture"] <- temporary_SoilMoisture.Daily #Re-add Soil Moisture column
ER_VPRM_MultipleReg.Daily.Mean ["VPD"] <- temporary_VPD.Daily #Re-add VPD column

#Stepwise Regression
#Performs stepwise model selection by AIC.
ER_VPRM_MultipleReg.Daily.Mean$Temperature <- NULL #Drop Temperature column ##FOR ER ONLY
MultipleReg.model.Daily.ER <- lm(ER_VPRM ~., data = ER_VPRM_MultipleReg.Daily.Mean)
MultipleReg.model.Daily.ER <- stepAIC(MultipleReg.model.Daily.ER, direction="both")
MultipleReg.model.Daily.ER$anova # display results
summary(MultipleReg.model.Daily.ER)
summary(MultipleReg.model.Daily.ER)$coefficient


##########Annual###############
variables.ER_VPRM_MultipleReg.annual <- c("ER_VPRM", "Temperature", "PAR", "Precipitation", "Evaporation", "SoilMoisture",  "VPD", "Year")
ER_VPRM_MultipleReg.Annual.subset <- VPRM_Fluxes.Annual[variables.ER_VPRM_MultipleReg.annual]
ER_VPRM_MultipleReg.Annual.Mean <- aggregate(.~ Year, data = ER_VPRM_MultipleReg.Annual.subset, mean, na.action = na.pass)
ER_VPRM_MultipleReg.Annual.Mean$Year <- NULL #Drop or remove the Day column
ER_VPRM_MultipleReg.Annual.Mean$Precipitation <- NULL #Drop or remove the mean Precipitation column. Sum needed
ER_VPRM_MultipleReg.Annual.Mean$Evaporation <- NULL #Drop or remove the mean Evaporation column. Sum needed
temporary_SoilMoisture.Annual <- ER_VPRM_MultipleReg.Annual.Mean$SoilMoisture #Create temporary dataframe for Soil Moisture for subsequent column rearrangement
temporary_VPD.Annual <- ER_VPRM_MultipleReg.Annual.Mean$VPD #Create temporary dataframe for VPD for subsequent column rearrangement
ER_VPRM_MultipleReg.Annual.Mean$SoilMoisture <- NULL #Drop or remove the Soil Moisture column. Will be re-added later.
ER_VPRM_MultipleReg.Annual.Mean$VPD <- NULL #Drop or remove the VPD column. Will be re-added later.

variables.ER_VPRM_MultipleReg.annual.2 <- c("Precipitation", "Evaporation", "Year")
ER_VPRM_MultipleReg.Annual.subset.2 <- VPRM_Fluxes.Annual[variables.ER_VPRM_MultipleReg.annual.2]
ER_VPRM_MultipleReg.Annual.Sum <- aggregate(.~ Year, data = ER_VPRM_MultipleReg.Annual.subset.2, sum, na.action = na.pass)
ER_VPRM_MultipleReg.Annual.Mean ["Precipitation"] <- ER_VPRM_MultipleReg.Annual.Sum$Precipitation
ER_VPRM_MultipleReg.Annual.Mean ["Evaporation"] <- ER_VPRM_MultipleReg.Annual.Sum$Evaporation
ER_VPRM_MultipleReg.Annual.Mean ["SoilMoisture"] <- temporary_SoilMoisture.Annual #Re-add Soil Moisture column
ER_VPRM_MultipleReg.Annual.Mean ["VPD"] <- temporary_VPD.Annual #Re-add VPD column

#Stepwise Regression
#Performs stepwise model selection by AIC.
ER_VPRM_MultipleReg.Annual.Mean$Temperature <- NULL #Drop Temperature column ##FOR ER ONLY
MultipleReg.model.Annual.ER <- lm(ER_VPRM ~., data = ER_VPRM_MultipleReg.Annual.Mean)
MultipleReg.model.Annual.ER <- stepAIC(MultipleReg.model.Annual.ER, direction="both")
MultipleReg.model.Annual.ER$anova # display results
summary(MultipleReg.model.Annual.ER)
summary(MultipleReg.model.Annual.ER)$coefficient



##########################################################Multiple Regression Analysis - NEE#############################################################

##########Daily################
variables.NEE_VPRM_MultipleReg.daily <- c("NEE_VPRM", "Temperature", "PAR", "Precipitation", "Evaporation", "SoilMoisture",  "VPD", "Day", "Month")
NEE_VPRM_MultipleReg.Daily.subset <- VPRM_Fluxes.Daily[variables.NEE_VPRM_MultipleReg.daily]
NEE_VPRM_MultipleReg.Daily.Mean <- aggregate(.~ Day + Month, data = NEE_VPRM_MultipleReg.Daily.subset, mean, na.action = na.pass)
NEE_VPRM_MultipleReg.Daily.Mean$Day <- NULL #Drop or remove the Day column
NEE_VPRM_MultipleReg.Daily.Mean$Month <- NULL #Drop or remove the Month column
NEE_VPRM_MultipleReg.Daily.Mean$Precipitation <- NULL #Drop or remove the mean Precipitation column. Sum needed
NEE_VPRM_MultipleReg.Daily.Mean$Evaporation <- NULL #Drop or remove the mean Evaporation column. Sum needed
temporary_SoilMoisture.Daily <- NEE_VPRM_MultipleReg.Daily.Mean$SoilMoisture #Create temporary dataframe for Soil Moisture for subsequent column rearrangement
temporary_VPD.Daily <- NEE_VPRM_MultipleReg.Daily.Mean$VPD #Create temporary dataframe for VPD for subsequent column rearrangement
NEE_VPRM_MultipleReg.Daily.Mean$SoilMoisture <- NULL #Drop or remove the Soil Moisture column. Will be re-added later.
NEE_VPRM_MultipleReg.Daily.Mean$VPD <- NULL #Drop or remove the VPD column. Will be re-added later.

variables.NEE_VPRM_MultipleReg.daily.2 <- c("Precipitation", "Evaporation", "Day", "Month", "Year")
NEE_VPRM_MultipleReg.Daily.subset.2 <- VPRM_Fluxes.Daily[variables.NEE_VPRM_MultipleReg.daily.2]
NEE_VPRM_MultipleReg.Daily.Sum <- aggregate(.~ Day + Month + Year, data = NEE_VPRM_MultipleReg.Daily.subset.2, sum, na.action = na.pass)
NEE_VPRM_MultipleReg.Daily.Sum$Year <- NULL #Drop or remove the Year column
NEE_VPRM_MultipleReg.Daily.Mean.2 <- aggregate(.~ Day + Month, data = NEE_VPRM_MultipleReg.Daily.Sum, mean, na.action = na.pass)
NEE_VPRM_MultipleReg.Daily.Mean ["Precipitation"] <- NEE_VPRM_MultipleReg.Daily.Mean.2$Precipitation
NEE_VPRM_MultipleReg.Daily.Mean ["Evaporation"] <- NEE_VPRM_MultipleReg.Daily.Mean.2$Evaporation
NEE_VPRM_MultipleReg.Daily.Mean ["SoilMoisture"] <- temporary_SoilMoisture.Daily #Re-add Soil Moisture column
NEE_VPRM_MultipleReg.Daily.Mean ["VPD"] <- temporary_VPD.Daily #Re-add VPD column

#Stepwise Regression
#Performs stepwise model selection by AIC.
MultipleReg.model.Daily.NEE <- lm(NEE_VPRM ~., data = NEE_VPRM_MultipleReg.Daily.Mean)
MultipleReg.model.Daily.NEE <- stepAIC(MultipleReg.model.Daily.NEE, direction="both")
MultipleReg.model.Daily.NEE$anova # display results
summary(MultipleReg.model.Daily.NEE)
summary(MultipleReg.model.Daily.NEE)$coefficient


##########Annual###############
variables.NEE_VPRM_MultipleReg.annual <- c("NEE_VPRM", "Temperature", "PAR", "Precipitation", "Evaporation", "SoilMoisture",  "VPD", "Year")
NEE_VPRM_MultipleReg.Annual.subset <- VPRM_Fluxes.Annual[variables.NEE_VPRM_MultipleReg.annual]
NEE_VPRM_MultipleReg.Annual.Mean <- aggregate(.~ Year, data = NEE_VPRM_MultipleReg.Annual.subset, mean, na.action = na.pass)
NEE_VPRM_MultipleReg.Annual.Mean$Year <- NULL #Drop or remove the Day column
NEE_VPRM_MultipleReg.Annual.Mean$Precipitation <- NULL #Drop or remove the mean Precipitation column. Sum needed
NEE_VPRM_MultipleReg.Annual.Mean$Evaporation <- NULL #Drop or remove the mean Evaporation column. Sum needed
temporary_SoilMoisture.Annual <- NEE_VPRM_MultipleReg.Annual.Mean$SoilMoisture #Create temporary dataframe for Soil Moisture for subsequent column rearrangement
temporary_VPD.Annual <- NEE_VPRM_MultipleReg.Annual.Mean$VPD #Create temporary dataframe for VPD for subsequent column rearrangement
NEE_VPRM_MultipleReg.Annual.Mean$SoilMoisture <- NULL #Drop or remove the Soil Moisture column. Will be re-added later.
NEE_VPRM_MultipleReg.Annual.Mean$VPD <- NULL #Drop or remove the VPD column. Will be re-added later.

variables.NEE_VPRM_MultipleReg.annual.2 <- c("Precipitation", "Evaporation", "Year")
NEE_VPRM_MultipleReg.Annual.subset.2 <- VPRM_Fluxes.Annual[variables.NEE_VPRM_MultipleReg.annual.2]
NEE_VPRM_MultipleReg.Annual.Sum <- aggregate(.~ Year, data = NEE_VPRM_MultipleReg.Annual.subset.2, sum, na.action = na.pass)
NEE_VPRM_MultipleReg.Annual.Mean ["Precipitation"] <- NEE_VPRM_MultipleReg.Annual.Sum$Precipitation
NEE_VPRM_MultipleReg.Annual.Mean ["Evaporation"] <- NEE_VPRM_MultipleReg.Annual.Sum$Evaporation
NEE_VPRM_MultipleReg.Annual.Mean ["SoilMoisture"] <- temporary_SoilMoisture.Annual #Re-add Soil Moisture column
NEE_VPRM_MultipleReg.Annual.Mean ["VPD"] <- temporary_VPD.Annual #Re-add VPD column

#Stepwise Regression
#Performs stepwise model selection by AIC.
MultipleReg.model.Annual.NEE <- lm(NEE_VPRM ~., data = NEE_VPRM_MultipleReg.Annual.Mean)
MultipleReg.model.Annual.NEE <- stepAIC(MultipleReg.model.Annual.NEE, direction="both")
MultipleReg.model.Annual.NEE$anova # display results
summary(MultipleReg.model.Annual.NEE)
summary(MultipleReg.model.Annual.NEE)$coefficient



######Arrange Outputs - Daily######
library(broom)
MultipleReg.model.Daily.GPP.transp <- t(cbind(tidy(MultipleReg.model.Daily.GPP), glance(MultipleReg.model.Daily.GPP)))
MultipleReg.model.Daily.GPP.transp[1] <- "GPP" #Rename

MultipleReg.model.Daily.ER.transp <- t(cbind(tidy(MultipleReg.model.Daily.ER), glance(MultipleReg.model.Daily.ER)))
MultipleReg.model.Daily.ER.transp[1] <- "ER" #Rename

MultipleReg.model.Daily.NEE.transp <- t(cbind(tidy(MultipleReg.model.Daily.NEE), glance(MultipleReg.model.Daily.NEE)))
MultipleReg.model.Daily.NEE.transp[1] <- "NEE" #Rename

#Now combine all fluxes
MultipleReg.model.Daily.All.Fluxes <- cbind(MultipleReg.model.Daily.GPP.transp, MultipleReg.model.Daily.ER.transp, MultipleReg.model.Daily.NEE.transp)


######Arrange Outputs - Annual######
MultipleReg.model.Annual.GPP.transp <- t(cbind(tidy(MultipleReg.model.Annual.GPP), glance(MultipleReg.model.Annual.GPP)))
MultipleReg.model.Annual.GPP.transp[1] <- "GPP" #Rename

MultipleReg.model.Annual.ER.transp <- t(cbind(tidy(MultipleReg.model.Annual.ER), glance(MultipleReg.model.Annual.ER)))
MultipleReg.model.Annual.ER.transp[1] <- "ER" #Rename

MultipleReg.model.Annual.NEE.transp <- t(cbind(tidy(MultipleReg.model.Annual.NEE), glance(MultipleReg.model.Annual.NEE)))
MultipleReg.model.Annual.NEE.transp[1] <- "NEE" #Rename

#Now combine all fluxes
MultipleReg.model.Annual.All.Fluxes <- cbind(MultipleReg.model.Annual.GPP.transp, MultipleReg.model.Annual.ER.transp, MultipleReg.model.Annual.NEE.transp)


#######Finally write to csv and move to stats folder######
library(filesstrings)
write.csv(MultipleReg.model.Daily.All.Fluxes, file = "Stepwise_Daily_All_Fluxes.csv")
write.csv(MultipleReg.model.Annual.All.Fluxes, file = "Stepwise_Annual_All_Fluxes.csv")
#file.move("Stepwise_Daily_All_Fluxes.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/ARB_FINAL/Stats")
#file.move("Stepwise_Annual_All_Fluxes.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/ARB_FINAL/Stats")



############################################################################################################################################

################################## Calculate Relative Importance for Each Predictor ########################################################
library(relaimpo)
my.plot.relimplmbooteval <- source("~/Ecosystem_Modelling/VPRM_POLAR/Chapter_4/Edited_plot_relimplmbooteval.R") #Load edited version of plot.relimplmbooteval from folder
plot.relimplmbooteval <- my.plot.relimplmbooteval$value #Replace original with edited version
#help(calc.relimp)
#plot.relimplmbooteval <- edit(plot.relimplmbooteval)
#plot.relimplm <- edit(plot.relimplm)

##GPP
site.name <- "(ARB)"
calc.relimp(MultipleReg.model.Annual.GPP, type=c("lmg"), rela=TRUE)
#Bootstrap Measures of Relative Importance (1000 samples)
boot.Annual.GPP <- boot.relimp(MultipleReg.model.Annual.GPP, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
#boot.Annual.GPP <- boot.relimp(MultipleReg.model.Annual.GPP, b = 1000, type = c("lmg","last", "first", "pratt"), rank = TRUE, diff = TRUE, rela = TRUE)
booteval.relimp(boot.Annual.GPP) # print result
boot.Extract.GPP <- booteval.relimp(boot.Annual.GPP) # Extract relevant ouput
boot.Output.GPP <- boot.Extract.GPP@mark
colnames(boot.Output.GPP)[1] <- "GPP" #Rename

##ER
calc.relimp(MultipleReg.model.Annual.ER, type=c("lmg"), rela=TRUE)
#Bootstrap Measures of Relative Importance (1000 samples)
boot.Annual.ER <- boot.relimp(MultipleReg.model.Annual.ER, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
#boot.Annual.ER <- boot.relimp(MultipleReg.model.Annual.ER, b = 1000, type = c("lmg","last", "first", "pratt"), rank = TRUE, diff = TRUE, rela = TRUE)
booteval.relimp(boot.Annual.ER) # print result
boot.Extract.ER <- booteval.relimp(boot.Annual.ER) # Extract relevant ouput
boot.Output.ER <- boot.Extract.ER@mark
colnames(boot.Output.ER)[1] <- "ER" #Rename

##NEE
calc.relimp(MultipleReg.model.Annual.NEE, type=c("lmg"), rela=TRUE)
#Bootstrap Measures of Relative Importance (1000 samples)
boot.Annual.NEE <- boot.relimp(MultipleReg.model.Annual.NEE, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
#boot.Annual.NEE <- boot.relimp(MultipleReg.model.Annual.NEE, b = 1000, type = c("lmg","last", "first", "pratt"), rank = TRUE, diff = TRUE, rela = TRUE)
booteval.relimp(boot.Annual.NEE) # print result
boot.Extract.NEE <- booteval.relimp(boot.Annual.NEE) # Extract relevant ouput
boot.Output.NEE <- boot.Extract.NEE@mark
colnames(boot.Output.NEE)[1] <- "NEE" #Rename

######Create Output - Annual Only
boot.Output.GPP.transp <- t(boot.Output.GPP)
sep.column.ER <- matrix(c(ER = "ER", "ER", "ER", "ER"))
boot.Output.ER.transp <- t(boot.Output.ER)
sep.column.NEE <- matrix(c(NEE = "NEE", "NEE", "NEE", "NEE"))
boot.Output.NEE.transp <- t(boot.Output.NEE)

boot.FINAL.Annual.Fluxes <- cbind(boot.Output.GPP.transp, sep.column.ER, boot.Output.ER.transp, sep.column.NEE, boot.Output.NEE.transp)

#######Finally write to csv and move to stats folder######
library(filesstrings)
write.csv(boot.FINAL.Annual.Fluxes, file = "boot_FINAL_Annual_Fluxes.csv")
#file.move("boot_FINAL_Annual_Fluxes.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/ARB_FINAL/Stats")


#save graph to output as high resolution, Windows compatible format .EMF
fileout <- paste("ARB_boot_Annual_GPP","_Barplot",".emf",sep="")
x11(width=7, height=3.33) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
dev.cur()
plot(booteval.relimp(boot.Annual.GPP,sort=TRUE))
savePlot(fileout, "emf", device = dev.cur())
dev.off()
graphics.off()


#save graph to output as high resolution, Windows compatible format .EMF
fileout <- paste("ARB_boot_Annual_ER","_Barplot",".emf",sep="")
x11(width=7, height=3.33) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
dev.cur()
plot(booteval.relimp(boot.Annual.ER,sort=TRUE))
savePlot(fileout, "emf", device = dev.cur())
dev.off()
graphics.off()


#save graph to output as high resolution, Windows compatible format .EMF
fileout <- paste("ARB_boot_Annual_NEE","_Barplot",".emf",sep="")
x11(width=7, height=3.33) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
dev.cur()
plot(booteval.relimp(boot.Annual.NEE,sort=TRUE))
savePlot(fileout, "emf", device = dev.cur())
dev.off()
graphics.off()

############################################################################################################################################



##################################################################Anomalies - Annual Fluxes################################################################
##Subset the relevant columns for aggregation
VPRM_Fluxes.Annual <- VPRM_Fluxes.Daily #Copy previously used dataframe from above to avoid potential conflicts
variables.fluxes.annual <- c("NEE_VPRM", "GPP_VPRM", "ER_VPRM", "Year")
VPRM_Fluxes.Annual.subset <- VPRM_Fluxes.Annual[variables.fluxes.annual]
VPRM_Fluxes.Annual.Mean <- aggregate(.~ Year, data = VPRM_Fluxes.Annual.subset, mean, na.action = na.pass)

VPRM_Fluxes.Annual.Mean.gC <- as.data.frame(VPRM_Fluxes.Annual.Mean[["Year"]])
names(VPRM_Fluxes.Annual.Mean.gC) <- c("Year")
VPRM_Fluxes.Annual.Mean.gC [["NEE_VPRM"]] <- VPRM_Fluxes.Annual.Mean$NEE_VPRM * (86400 * 153 * 12.011) / 1000000 #Convert from umols of CO2 m^-2 s^-1 to g C m^-2. Note that 153 is the number of days in the study period per year.
VPRM_Fluxes.Annual.Mean.gC [["GPP_VPRM"]] <- VPRM_Fluxes.Annual.Mean$GPP_VPRM * (86400 * 153 * 12.011) / 1000000 #Convert from umols of CO2 m^-2 s^-1 to g C m^-2. Note that 153 is the number of days in the study period per year.
VPRM_Fluxes.Annual.Mean.gC [["ER_VPRM"]] <- VPRM_Fluxes.Annual.Mean$ER_VPRM * (86400 * 153 * 12.011) / 1000000 #Convert from umols of CO2 m^-2 s^-1 to g C m^-2. Note that 153 is the number of days in the study period per year.

#Calculate mean (climatology) for the study period (20 years)
#VPRM_Fluxes.climatology <- subset(VPRM_Fluxes.Annual.Mean.gC, Year %in% 2000:2019)
VPRM_Fluxes.climatology.mean.gC.NEE <- mean(VPRM_Fluxes.Annual.Mean.gC$NEE_VPRM)
VPRM_Fluxes.climatology.mean.gC.GPP <- mean(VPRM_Fluxes.Annual.Mean.gC$GPP_VPRM)
VPRM_Fluxes.climatology.mean.gC.ER <- mean(VPRM_Fluxes.Annual.Mean.gC$ER_VPRM)

VPRM_Fluxes.climatology.mean.gC.NEE
VPRM_Fluxes.climatology.mean.gC.GPP
VPRM_Fluxes.climatology.mean.gC.ER

##Create data for Anomaly timeseries. Start by creating a duplicate copy of mean annual dataframe
VPRM_Fluxes.Annual.Mean.gC.new <- VPRM_Fluxes.Annual.Mean.gC

#Add new columns with the mean climatology
VPRM_Fluxes.Annual.Mean.gC.new["Climatology_NEE"] <- VPRM_Fluxes.climatology.mean.gC.NEE
VPRM_Fluxes.Annual.Mean.gC.new["Climatology_GPP"] <- VPRM_Fluxes.climatology.mean.gC.GPP
VPRM_Fluxes.Annual.Mean.gC.new["Climatology_ER"] <- VPRM_Fluxes.climatology.mean.gC.ER

#Calculate anomalies and add new columns
VPRM_Fluxes.Annual.anomaly.gC.NEE <- VPRM_Fluxes.Annual.Mean.gC.new$NEE_VPRM - VPRM_Fluxes.Annual.Mean.gC.new$Climatology_NEE
VPRM_Fluxes.Annual.anomaly.gC.GPP <- VPRM_Fluxes.Annual.Mean.gC.new$GPP_VPRM - VPRM_Fluxes.Annual.Mean.gC.new$Climatology_GPP
VPRM_Fluxes.Annual.anomaly.gC.ER <- VPRM_Fluxes.Annual.Mean.gC.new$ER_VPRM - VPRM_Fluxes.Annual.Mean.gC.new$Climatology_ER

VPRM_Fluxes.Annual.Mean.gC.new["Anomaly_NEE"] <- VPRM_Fluxes.Annual.anomaly.gC.NEE
VPRM_Fluxes.Annual.Mean.gC.new["Anomaly_GPP"] <- VPRM_Fluxes.Annual.anomaly.gC.GPP
VPRM_Fluxes.Annual.Mean.gC.new["Anomaly_ER"] <- VPRM_Fluxes.Annual.anomaly.gC.ER

VPRM_Fluxes.Annual.Mean.gC.new ["Year_Char"] <- as.character(VPRM_Fluxes.Annual.Mean.gC.new$Year) #Year column is in class ordered or factor. Need to change to Numeric
VPRM_Fluxes.Annual.Mean.gC.new$Year <- NULL #Drop or remove the original ordered/factor Year column
VPRM_Fluxes.Annual.Mean.gC.new ["Year"] <- as.numeric(VPRM_Fluxes.Annual.Mean.gC.new$Year_Char)
VPRM_Fluxes.Annual.Mean.gC.new$Year_Char <- NULL #Drop or remove the Character Year column 

##Subset the relevant columns to plot anomalies
variables.fluxes.annual.anomalies.gC <- c("Year", "Anomaly_NEE", "Anomaly_GPP", "Anomaly_ER")
VPRM_Fluxes.Annual.anomaly.gC.Final <- VPRM_Fluxes.Annual.Mean.gC.new[variables.fluxes.annual.anomalies.gC]

#############Statistical tests for Trends################
# Install and load the packages
#install.packages(c("car", "Kendall", "zyp", "lmtest"))
library(car)
library(Kendall)
library(zyp)
library(lmtest)

#Test for Autocorrelation - durbinWatsonTest
y1.NEE <- VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_NEE
x1 <- VPRM_Fluxes.Annual.anomaly.gC.Final$Year
dwtest(y1.NEE ~ x1)

y1.GPP <- VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_GPP
x1 <- VPRM_Fluxes.Annual.anomaly.gC.Final$Year
dwtest(y1.GPP ~ x1)

y1.ER <- VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_ER
x1 <- VPRM_Fluxes.Annual.anomaly.gC.Final$Year
dwtest(y1.ER ~ x1)

#Run the Mann-Kendall test
y.NEE <- VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_NEE
Trend.Annual.gC.NEE <- MannKendall(y.NEE)
Trend.Annual.gC.NEE

y.GPP <- VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_GPP
Trend.Annual.gC.GPP <- MannKendall(y.GPP)
Trend.Annual.gC.GPP

y.ER <- VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_ER
Trend.Annual.gC.ER <- MannKendall(y.ER)
Trend.Annual.gC.ER

#TheilSen Slope Estimation
y.NEE <- VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_NEE
x <- VPRM_Fluxes.Annual.anomaly.gC.Final$Year
Slope.Annual.gC.NEE <- zyp.sen(y.NEE~x)
Slope.Annual.gC.NEE

y.GPP <- VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_GPP
x <- VPRM_Fluxes.Annual.anomaly.gC.Final$Year
Slope.Annual.gC.GPP <- zyp.sen(y.GPP~x)
Slope.Annual.gC.GPP

y.ER <- VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_ER
x <- VPRM_Fluxes.Annual.anomaly.gC.Final$Year
Slope.Annual.gC.ER <- zyp.sen(y.ER~x)
Slope.Annual.gC.ER

#Subset and create objects for trend outputs to add to plots
Trend.Annual.gC.NEE
Slope.Annual.gC.NEE
Trend.Annual.gC.NEE.tau <- format(round(Trend.Annual.gC.NEE$tau[1], 2), nsmall = 2) #Round to two decimals
Trend.Annual.gC.NEE.pvalue <- format(round(Trend.Annual.gC.NEE$sl[1], 2), nsmall = 2) #Round to two decimals
Slope.Annual.gC.NEE.20 <- format(round(Slope.Annual.gC.NEE$coefficients[["x"]]*20, 2), nsmall = 2) #Round to two decimals

Trend.Annual.gC.GPP
Slope.Annual.gC.GPP
Trend.Annual.gC.GPP.tau <- format(round(Trend.Annual.gC.GPP$tau[1], 2), nsmall = 2) #Round to two decimals
Trend.Annual.gC.GPP.pvalue <- format(round(Trend.Annual.gC.GPP$sl[1], 2), nsmall = 2) #Round to two decimals
Slope.Annual.gC.GPP.20 <- format(round(Slope.Annual.gC.GPP$coefficients[["x"]]*20, 2), nsmall = 2) #Round to two decimals

Trend.Annual.gC.ER
Slope.Annual.gC.ER
Trend.Annual.gC.ER.tau <- format(round(Trend.Annual.gC.ER$tau[1], 2), nsmall = 2) #Round to two decimals
Trend.Annual.gC.ER.pvalue <- format(round(Trend.Annual.gC.ER$sl[1], 2), nsmall = 2) #Round to two decimals
Slope.Annual.gC.ER.20 <- format(round(Slope.Annual.gC.ER$coefficients[["x"]]*20, 2), nsmall = 2) #Round to two decimals


##Melt to prepare data form ggplot multiple variables
VPRM_Fluxes.Annual.anomaly.gC.Final.Melted <- reshape2::melt(VPRM_Fluxes.Annual.anomaly.gC.Final, id.var='Year')
head(VPRM_Fluxes.Annual.anomaly.gC.Final.Melted)
tail(VPRM_Fluxes.Annual.anomaly.gC.Final.Melted)

###########Create plots of the annual flux anomalies for the growing season
library(ggplot2)
library(ggpubr)
library(ggpmisc)
plot <- ggplot(VPRM_Fluxes.Annual.anomaly.gC.Final.Melted, aes(x=Year, y=value, group=variable)) + geom_bar(stat = "identity", aes(fill = variable), position = "dodge", width = 0.7)
plot <- plot + annotate("text", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.05, size = 3, fontface = "bold", label = (paste("NEE: Tau = ", Trend.Annual.gC.NEE.tau, ", p = ", Trend.Annual.gC.NEE.pvalue, ", TSA = ", Slope.Annual.gC.NEE.20, sep="")), color = "#f8766d")
plot <- plot + annotate("text", x = -Inf, y = Inf, vjust = 3.0, hjust = -0.05, size = 3, fontface = "bold", label = (paste("GPP: Tau = ", Trend.Annual.gC.GPP.tau, ", p = ", Trend.Annual.gC.GPP.pvalue, ", TSA = ", Slope.Annual.gC.GPP.20, sep="")), color = "#00ba38")
plot <- plot + annotate("text", x = -Inf, y = Inf, vjust = 4.5, hjust = -0.05, size = 3, fontface = "bold", label = (paste("ER: Tau = ", Trend.Annual.gC.ER.tau, ", p = ", Trend.Annual.gC.ER.pvalue, ", TSA = ", Slope.Annual.gC.ER.20, sep="")), color = "#619cff")
plot <- plot + xlab("Year") + ylab(expression(Carbon~Fluxes~(g~C~m^{-2}))) + ggtitle(paste("ARB Annual Flux Anomalies" ,"(June-October, 2000-2019)"))
plot <- plot + scale_x_continuous(breaks = seq(from=2000,to=2020,by=2))
#plot <- plot + scale_x_continuous(labels =c("June","July","August","September","October"))
plot <- plot + theme(plot.title = element_text(hjust = 0.5)) #center title
plot <- plot + theme(axis.text=element_text(size=11, face="bold"))
plot <- plot + theme(axis.title=element_text(size=11))
plot <- plot + theme(plot.title=element_text(size=11,face="bold"))
#plot <- plot + theme(legend.key.size = unit(1.0, "cm"),legend.key.width = unit(1.0,"cm")) 
plot <- plot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank())
plot <- plot + labs(shape = "") #Remove or specify Legend title
plot <- plot + theme(legend.title = element_blank()) #Remove Legend title
plot <- plot + theme(legend.text = element_text(size = 11, face = "bold"))
plot <- plot + theme(legend.position = "bottom")
plot <- plot + geom_hline(yintercept=0, color="black")
plot <- plot + scale_y_continuous(breaks=seq(-50,75,25)) #Sets the axis range and interval ("by")
plot <- plot + coord_cartesian(ylim = c(-50, 75), xlim = c(2000, 2020)) #Forces the limits to be set accordingly
plot.Anomalies.Annual <- plot

#save graph to output as high resoultion, Windows compatible format .EMF
fileout <- paste("Line-Final-ARB Annual Flux gC Anomalies","_Bar",".emf",sep="")
#x11(width=7, height=4) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
x11(width=7, height=4) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
dev.cur()
plot.Anomalies.Annual 
savePlot(fileout, "emf", device = dev.cur())
dev.off()
graphics.off()

library(filesstrings)
#Write VPRM Annual Anomalies to csv file and Move the excel .csv files to the stats folder
write.csv(VPRM_Fluxes.Annual.anomaly.gC.Final, file = "VPRM_Fluxes_Annual_anomaly_gC_Final.csv", row.names = F)
#file.move("VPRM_Fluxes_Annual_anomaly_gC_Final.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/ARB_FINAL/Stats")




##################################################################Anomalies - Annual Temperature################################################################
##Subset the relevant columns for aggregation
variables.Temperature.annual <- c("Temperature", "Year")
VPRM_Temperature.Annual.subset <- VPRM_Fluxes.Annual[variables.Temperature.annual]
VPRM_Temperature.Annual.Mean <- aggregate(.~ Year, data = VPRM_Temperature.Annual.subset, mean, na.action = na.pass)

#Calculate mean (climatology) for the study period (20 years)
#VPRM_Temperature.climatology <- subset(VPRM_Temperature.Annual.Mean, Year %in% 2000:2019)
VPRM_Temperature.climatology.mean <- mean(VPRM_Temperature.Annual.Mean$Temperature)

VPRM_Temperature.climatology.mean

##Create data for Anomaly timeseries. Start by creating a duplicate copy of mean annual dataframe
VPRM_Temperature.Annual.Mean.new <- VPRM_Temperature.Annual.Mean

#Add new columns with the mean climatology
VPRM_Temperature.Annual.Mean.new["Climatology_Temperature"] <- VPRM_Temperature.climatology.mean

#Calculate anomalies and add new columns
VPRM_Temperature.Annual.anomaly <- VPRM_Temperature.Annual.Mean.new$Temperature - VPRM_Temperature.Annual.Mean.new$Climatology_Temperature

VPRM_Temperature.Annual.Mean.new["Anomaly_Temperature"] <- VPRM_Temperature.Annual.anomaly

VPRM_Temperature.Annual.Mean.new ["Year_Char"] <- as.character(VPRM_Temperature.Annual.Mean.new$Year) #Year column is in class ordered or factor. Need to change to Numeric
VPRM_Temperature.Annual.Mean.new$Year <- NULL #Drop or remove the original ordered/factor Year column
VPRM_Temperature.Annual.Mean.new ["Year"] <- as.numeric(VPRM_Temperature.Annual.Mean.new$Year_Char)
VPRM_Temperature.Annual.Mean.new$Year_Char <- NULL #Drop or remove the Character Year column 

##Subset the relevant columns to plot anomalies
variables.Temperature.annual.anomalies <- c("Year", "Anomaly_Temperature")
VPRM_Temperature.Annual.anomaly.Final <- VPRM_Temperature.Annual.Mean.new[variables.Temperature.annual.anomalies]

#############Statistical tests for Trends################
# Install and load the packages
#install.packages(c("car", "Kendall", "zyp", "lmtest"))
library(car)
library(Kendall)
library(zyp)
library(lmtest)

#Test for Autocorrelation - durbinWatsonTest
y1.Temperature <- VPRM_Temperature.Annual.anomaly.Final$Anomaly_Temperature
x1 <- VPRM_Temperature.Annual.anomaly.Final$Year
dwtest(y1.Temperature ~ x1)

#Run the Mann-Kendall test
y.Temperature <- VPRM_Temperature.Annual.anomaly.Final$Anomaly_Temperature
Trend.Annual.Temperature <- MannKendall(y.Temperature)
Trend.Annual.Temperature

#TheilSen Slope Estimation
y.Temperature <- VPRM_Temperature.Annual.anomaly.Final$Anomaly_Temperature
x <- VPRM_Temperature.Annual.anomaly.Final$Year
Slope.Annual.Temperature <- zyp.sen(y.Temperature~x)
Slope.Annual.Temperature

#Subset and create objects for trend outputs to add to plots
Trend.Annual.Temperature
Slope.Annual.Temperature
Trend.Annual.Temperature.tau <- format(round(Trend.Annual.Temperature$tau[1], 2), nsmall = 2) #Round to two decimals
Trend.Annual.Temperature.pvalue <- format(round(Trend.Annual.Temperature$sl[1], 2), nsmall = 2) #Round to two decimals
Slope.Annual.Temperature.20 <- format(round(Slope.Annual.Temperature$coefficients[["x"]]*20, 2), nsmall = 2) #Round to two decimals

###########Create plots of the annual flux anomalies for the growing season
library(ggplot2)
library(ggpubr)
library(ggpmisc)
VPRM_Temperature.Annual.anomaly.Final ["pos"] <- (VPRM_Temperature.Annual.anomaly.Final$Anomaly_Temperature >=0) #'pos' is a new column specifying +ve or -ve values as TRUE or FALSE
plot <- ggplot(VPRM_Temperature.Annual.anomaly.Final, aes(x=Year, y=Anomaly_Temperature, group=1)) + geom_bar(stat = "identity", aes(fill = pos), position = "dodge", width = 0.7) + scale_fill_manual(values = c("#A9A9A9","black"), guide = FALSE) + geom_smooth(method = "lm", colour = "red", se = FALSE)
plot <- plot + annotate("text", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.05, size = 3, fontface = "bold", label = (paste("Tau = ", Trend.Annual.Temperature.tau, ", p = ", Trend.Annual.Temperature.pvalue, ", TSA = ", Slope.Annual.Temperature.20, sep="")), color = "black")
plot <- plot + xlab("Year") + ylab(expression("Temperature (\u00B0C)")) + ggtitle(paste("ARB Annual Temperature Anomalies" ,"(June-October, 2000-2019)"))
plot <- plot + scale_x_continuous(breaks = seq(from=2000,to=2020,by=2))
#plot <- plot + scale_x_continuous(labels =c("June","July","August","September","October"))
plot <- plot + theme(plot.title = element_text(hjust = 0.5)) #center title
plot <- plot + theme(axis.text=element_text(size=11, face="bold"))
plot <- plot + theme(axis.title=element_text(size=11))
plot <- plot + theme(plot.title=element_text(size=11,face="bold"))
#plot <- plot + theme(legend.key.size = unit(1.0, "cm"),legend.key.width = unit(1.0,"cm")) 
plot <- plot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank())
plot <- plot + labs(shape = "") #Remove or specify Legend title
plot <- plot + theme(legend.title = element_blank()) #Remove Legend title
plot <- plot + theme(legend.text = element_text(size = 11, face = "bold"))
plot <- plot + theme(legend.position = "bottom")
plot <- plot + geom_hline(yintercept=0, color="black")
plot <- plot + scale_y_continuous(breaks=seq(-2,2,1)) #Sets the axis range and interval ("by")
plot <- plot + coord_cartesian(ylim = c(-2, 2), xlim = c(2000, 2020)) #Forces the limits to be set accordingly
plot.Temperature.Anomalies.Annual <- plot

#save graph to output as high resoultion, Windows compatible format .EMF
fileout <- paste("Line-Final-ARB Annual Temperature Anomalies","_Bar",".emf",sep="")
x11(width=7, height=4) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
dev.cur()
plot.Temperature.Anomalies.Annual 
savePlot(fileout, "emf", device = dev.cur())
dev.off()
graphics.off()

#Write VPRM Annual Anomalies to csv file and Move the excel .csv files to the stats folder
write.csv(VPRM_Temperature.Annual.anomaly.Final, file = "VPRM_Temperature_Annual_anomaly_Final.csv", row.names = F)
#file.move("VPRM_Temperature_Annual_anomaly_Final.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/ARB_FINAL/Stats")



##################################################################Anomalies - Annual PAR################################################################
##Subset the relevant columns for aggregation
variables.PAR.annual <- c("PAR", "Year")
VPRM_PAR.Annual.subset <- VPRM_Fluxes.Annual[variables.PAR.annual]
VPRM_PAR.Annual.Mean <- aggregate(.~ Year, data = VPRM_PAR.Annual.subset, mean, na.action = na.pass)

#Calculate mean (climatology) for the study period (20 years)
#VPRM_PAR.climatology <- subset(VPRM_PAR.Annual.Mean, Year %in% 2000:2019)
VPRM_PAR.climatology.mean <- mean(VPRM_PAR.Annual.Mean$PAR)

VPRM_PAR.climatology.mean

##Create data for Anomaly timeseries. Start by creating a duplicate copy of mean annual dataframe
VPRM_PAR.Annual.Mean.new <- VPRM_PAR.Annual.Mean

#Add new columns with the mean climatology
VPRM_PAR.Annual.Mean.new["Climatology_PAR"] <- VPRM_PAR.climatology.mean

#Calculate anomalies and add new columns
VPRM_PAR.Annual.anomaly <- VPRM_PAR.Annual.Mean.new$PAR - VPRM_PAR.Annual.Mean.new$Climatology_PAR

VPRM_PAR.Annual.Mean.new["Anomaly_PAR"] <- VPRM_PAR.Annual.anomaly

VPRM_PAR.Annual.Mean.new ["Year_Char"] <- as.character(VPRM_PAR.Annual.Mean.new$Year) #Year column is in class ordered or factor. Need to change to Numeric
VPRM_PAR.Annual.Mean.new$Year <- NULL #Drop or remove the original ordered/factor Year column
VPRM_PAR.Annual.Mean.new ["Year"] <- as.numeric(VPRM_PAR.Annual.Mean.new$Year_Char)
VPRM_PAR.Annual.Mean.new$Year_Char <- NULL #Drop or remove the Character Year column 

##Subset the relevant columns to plot anomalies
variables.PAR.annual.anomalies <- c("Year", "Anomaly_PAR")
VPRM_PAR.Annual.anomaly.Final <- VPRM_PAR.Annual.Mean.new[variables.PAR.annual.anomalies]

#############Statistical tests for Trends################
# Install and load the packages
#install.packages(c("car", "Kendall", "zyp", "lmtest"))
library(car)
library(Kendall)
library(zyp)
library(lmtest)

#Test for Autocorrelation - durbinWatsonTest
y1.PAR <- VPRM_PAR.Annual.anomaly.Final$Anomaly_PAR
x1 <- VPRM_PAR.Annual.anomaly.Final$Year
dwtest(y1.PAR ~ x1)

#Run the Mann-Kendall test
y.PAR <- VPRM_PAR.Annual.anomaly.Final$Anomaly_PAR
Trend.Annual.PAR <- MannKendall(y.PAR)
Trend.Annual.PAR

#TheilSen Slope Estimation
y.PAR <- VPRM_PAR.Annual.anomaly.Final$Anomaly_PAR
x <- VPRM_PAR.Annual.anomaly.Final$Year
Slope.Annual.PAR <- zyp.sen(y.PAR~x)
Slope.Annual.PAR

#Subset and create objects for trend outputs to add to plots
Trend.Annual.PAR
Slope.Annual.PAR
Trend.Annual.PAR.tau <- format(round(Trend.Annual.PAR$tau[1], 2), nsmall = 2) #Round to two decimals
Trend.Annual.PAR.pvalue <- format(round(Trend.Annual.PAR$sl[1], 2), nsmall = 2) #Round to two decimals
Slope.Annual.PAR.20 <- format(round(Slope.Annual.PAR$coefficients[["x"]]*20, 2), nsmall = 2) #Round to two decimals

###########Create plots of the annual flux anomalies for the growing season
library(ggplot2)
library(ggpubr)
library(ggpmisc)
VPRM_PAR.Annual.anomaly.Final ["pos"] <- (VPRM_PAR.Annual.anomaly.Final$Anomaly_PAR >=0) #'pos' is a new column specifying +ve or -ve values as TRUE or FALSE
plot <- ggplot(VPRM_PAR.Annual.anomaly.Final, aes(x=Year, y=Anomaly_PAR, group=1)) + geom_bar(stat = "identity", aes(fill = pos), position = "dodge", width = 0.7) + scale_fill_manual(values = c("#A9A9A9","black"), guide = FALSE) + geom_smooth(method = "lm", colour = "red", se = FALSE)
plot <- plot + annotate("text", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.05, size = 3, fontface = "bold", label = (paste("Tau = ", Trend.Annual.PAR.tau, ", p = ", Trend.Annual.PAR.pvalue, ", TSA = ", Slope.Annual.PAR.20, sep="")), color = "black")
plot <- plot + xlab("Year") + ylab(expression(PAR~(mu*mol~m^{-2}~s^{-1}))) + ggtitle(paste("ARB Annual PAR Anomalies" ,"(June-October, 2000-2019)"))
plot <- plot + scale_x_continuous(breaks = seq(from=2000,to=2020,by=2))
#plot <- plot + scale_x_continuous(labels =c("June","July","August","September","October"))
plot <- plot + theme(plot.title = element_text(hjust = 0.5)) #center title
plot <- plot + theme(axis.text=element_text(size=11, face="bold"))
plot <- plot + theme(axis.title=element_text(size=11))
plot <- plot + theme(plot.title=element_text(size=11,face="bold"))
#plot <- plot + theme(legend.key.size = unit(1.0, "cm"),legend.key.width = unit(1.0,"cm")) 
plot <- plot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank())
plot <- plot + labs(shape = "") #Remove or specify Legend title
plot <- plot + theme(legend.title = element_blank()) #Remove Legend title
plot <- plot + theme(legend.text = element_text(size = 11, face = "bold"))
plot <- plot + theme(legend.position = "bottom")
plot <- plot + geom_hline(yintercept=0, color="black")
plot <- plot + scale_y_continuous(breaks=seq(-50,50,25)) #Sets the axis range and interval ("by")
plot <- plot + coord_cartesian(ylim = c(-50, 50), xlim = c(2000, 2020)) #Forces the limits to be set accordingly
plot.PAR.Anomalies.Annual <- plot

#save graph to output as high resoultion, Windows compatible format .EMF
fileout <- paste("Line-Final-ARB Annual PAR Anomalies","_Bar",".emf",sep="")
x11(width=7, height=4) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
dev.cur()
plot.PAR.Anomalies.Annual
savePlot(fileout, "emf", device = dev.cur())
dev.off()
graphics.off()

library(filesstrings)
#Write VPRM Annual Anomalies to csv file and Move the excel .csv files to the stats folder
write.csv(VPRM_PAR.Annual.anomaly.Final, file = "VPRM_PAR_Annual_anomaly_Final.csv", row.names = F)
#file.move("VPRM_PAR_Annual_anomaly_Final.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/ARB_FINAL/Stats")



##################################################################Anomalies - Annual Precip_Evap################################################################
##Subset the relevant columns for aggregation
variables.Precip_Evap.annual <- c("Precipitation", "Evaporation", "Year")
VPRM_Precip_Evap.Annual.subset <- VPRM_Fluxes.Annual[variables.Precip_Evap.annual]
VPRM_Precip_Evap.Annual.Mean <- aggregate(.~ Year, data = VPRM_Precip_Evap.Annual.subset, sum, na.action = na.pass) #Note that this is actually sum and not mean

#Calculate mean (climatology) for the study period (20 years)
#VPRM_Precipitation.climatology <- subset(VPRM_Precipitation.Annual.Mean, Year %in% 2000:2019)
VPRM_Precipitation.climatology.mean <- mean(VPRM_Precip_Evap.Annual.Mean$Precipitation)
VPRM_Evaporation.climatology.mean <- mean(VPRM_Precip_Evap.Annual.Mean$Evaporation)

VPRM_Precipitation.climatology.mean
VPRM_Evaporation.climatology.mean

##Create data for Anomaly timeseries. Start by creating a duplicate copy of mean annual dataframe
VPRM_Precip_Evap.Annual.Mean.new <- VPRM_Precip_Evap.Annual.Mean

#Add new columns with the mean climatology
VPRM_Precip_Evap.Annual.Mean.new["Climatology_Precipitation"] <- VPRM_Precipitation.climatology.mean
VPRM_Precip_Evap.Annual.Mean.new["Climatology_Evaporation"] <- VPRM_Evaporation.climatology.mean

#Calculate anomalies and add new columns
VPRM_Precipitation.Annual.anomaly <- VPRM_Precip_Evap.Annual.Mean.new$Precipitation - VPRM_Precip_Evap.Annual.Mean.new$Climatology_Precipitation
VPRM_Evaporation.Annual.anomaly <- VPRM_Precip_Evap.Annual.Mean.new$Evaporation - VPRM_Precip_Evap.Annual.Mean.new$Climatology_Evaporation

VPRM_Precip_Evap.Annual.Mean.new["Anomaly_Precipitation"] <- VPRM_Precipitation.Annual.anomaly
VPRM_Precip_Evap.Annual.Mean.new["Anomaly_Evaporation"] <- VPRM_Evaporation.Annual.anomaly

VPRM_Precip_Evap.Annual.Mean.new ["Year_Char"] <- as.character(VPRM_Precip_Evap.Annual.Mean.new$Year) #Year column is in class ordered or factor. Need to change to Numeric
VPRM_Precip_Evap.Annual.Mean.new$Year <- NULL #Drop or remove the original ordered/factor Year column
VPRM_Precip_Evap.Annual.Mean.new ["Year"] <- as.numeric(VPRM_Precip_Evap.Annual.Mean.new$Year_Char)
VPRM_Precip_Evap.Annual.Mean.new$Year_Char <- NULL #Drop or remove the Character Year column 

##Subset the relevant columns to plot anomalies
variables.Precip_Evap.annual.anomalies <- c("Year", "Anomaly_Precipitation", "Anomaly_Evaporation")
VPRM_Precip_Evap.Annual.anomaly.Final <- VPRM_Precip_Evap.Annual.Mean.new[variables.Precip_Evap.annual.anomalies]

#############Statistical tests for Trends################
# Install and load the packages
#install.packages(c("car", "Kendall", "zyp", "lmtest"))
library(car)
library(Kendall)
library(zyp)
library(lmtest)

#Test for Autocorrelation - durbinWatsonTest
y1.Precipitation <- VPRM_Precip_Evap.Annual.anomaly.Final$Anomaly_Precipitation
x1 <- VPRM_Precip_Evap.Annual.anomaly.Final$Year
dwtest(y1.Precipitation ~ x1)

y1.Evaporation <- VPRM_Precip_Evap.Annual.anomaly.Final$Anomaly_Evaporation
x1 <- VPRM_Precip_Evap.Annual.anomaly.Final$Year
dwtest(y1.Evaporation ~ x1)


#Run the Mann-Kendall test
y.Precipitation <- VPRM_Precip_Evap.Annual.anomaly.Final$Anomaly_Precipitation
Trend.Annual.Precipitation <- MannKendall(y.Precipitation)
Trend.Annual.Precipitation

y.Evaporation <- VPRM_Precip_Evap.Annual.anomaly.Final$Anomaly_Evaporation
Trend.Annual.Evaporation <- MannKendall(y.Evaporation)
Trend.Annual.Evaporation

#TheilSen Slope Estimation
y.Precipitation <- VPRM_Precip_Evap.Annual.anomaly.Final$Anomaly_Precipitation
x <- VPRM_Precip_Evap.Annual.anomaly.Final$Year
Slope.Annual.Precipitation <- zyp.sen(y.Precipitation~x)
Slope.Annual.Precipitation

y.Evaporation <- VPRM_Precip_Evap.Annual.anomaly.Final$Anomaly_Evaporation
x <- VPRM_Precip_Evap.Annual.anomaly.Final$Year
Slope.Annual.Evaporation <- zyp.sen(y.Evaporation~x)
Slope.Annual.Evaporation

#Subset and create objects for trend outputs to add to plots
Trend.Annual.Precipitation
Slope.Annual.Precipitation
Trend.Annual.Precipitation.tau <- format(round(Trend.Annual.Precipitation$tau[1], 2), nsmall = 2) #Round to two decimals
Trend.Annual.Precipitation.pvalue <- format(round(Trend.Annual.Precipitation$sl[1], 2), nsmall = 2) #Round to two decimals
Slope.Annual.Precipitation.20 <- format(round(Slope.Annual.Precipitation$coefficients[["x"]]*20, 2), nsmall = 2) #Round to two decimals

Trend.Annual.Evaporation
Slope.Annual.Evaporation
Trend.Annual.Evaporation.tau <- format(round(Trend.Annual.Evaporation$tau[1], 2), nsmall = 2) #Round to two decimals
Trend.Annual.Evaporation.pvalue <- format(round(Trend.Annual.Evaporation$sl[1], 2), nsmall = 2) #Round to two decimals
Slope.Annual.Evaporation.20 <- format(round(Slope.Annual.Evaporation$coefficients[["x"]]*20, 2), nsmall = 2) #Round to two decimals


##Melt to prepare data form ggplot multiple variables
VPRM_Precip_Evap.Annual.anomaly.Final.Melted <- reshape2::melt(VPRM_Precip_Evap.Annual.anomaly.Final, id.var='Year')
head(VPRM_Precip_Evap.Annual.anomaly.Final.Melted)
tail(VPRM_Precip_Evap.Annual.anomaly.Final.Melted)

###########Create plots of the annual flux anomalies for the growing season
library(ggplot2)
library(ggpubr)
library(ggpmisc)
VPRM_Precip_Evap.Annual.anomaly.Final.Melted ["pos"] <- (VPRM_Precip_Evap.Annual.anomaly.Final.Melted$value >=0) #'pos' is a new column specifying +ve or -ve values as TRUE or FALSE
plot <- ggplot(VPRM_Precip_Evap.Annual.anomaly.Final.Melted, aes(x=Year, y=value, group=variable)) + geom_bar(stat = "identity", aes(fill = pos), colour = "black", linetype = ifelse(VPRM_Precip_Evap.Annual.anomaly.Final.Melted$variable == "Anomaly_Precipitation", "solid", "dotted"), position = "dodge", width = 0.7) + scale_fill_manual(values = c("#A9A9A9","black"), guide = FALSE) + geom_smooth(method = "lm", colour = "red", aes(linetype = variable), se = FALSE, show.legend = TRUE)
plot <- plot + annotate("text", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.05, size = 3, fontface = "bold", label = (paste("Prec: Tau = ", Trend.Annual.Precipitation.tau, ", p = ", Trend.Annual.Precipitation.pvalue, ", TSA = ", Slope.Annual.Precipitation.20, sep="")), color = "black")
plot <- plot + annotate("text", x = -Inf, y = Inf, vjust = 3.0, hjust = -0.05, size = 3, fontface = "bold", label = (paste("Evap: Tau = ", Trend.Annual.Evaporation.tau, ", p = ", Trend.Annual.Evaporation.pvalue, ", TSA = ", Slope.Annual.Evaporation.20, sep="")), color = "black")
plot <- plot + scale_shape_manual(values=c(16, 1)) #Specify shape type for geom_point
plot <- plot + xlab("Year") + ylab(expression("P & E (mm)")) + ggtitle(paste("ARB Annual Total Precip & Evap Anomalies" ,"(June-October, 2000-2019)"))
plot <- plot + scale_x_continuous(breaks = seq(from=2000,to=2020,by=2))
plot <- plot + theme(plot.title = element_text(hjust = 0.5)) #center title
plot <- plot + theme(axis.text=element_text(size=11, face="bold"))
plot <- plot + theme(axis.title=element_text(size=11))
plot <- plot + theme(plot.title=element_text(size=11,face="bold"))
plot <- plot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank())
plot <- plot + theme(legend.title = element_blank()) #Remove Legend title
plot <- plot + theme(legend.text = element_text(size = 11, face = "bold"))
#plot <- plot + theme(legend.position = c(0.83, 0.85), legend.text = element_text(size = 8, face = "bold"), legend.key = element_rect(colour = "transparent", fill = "white")) + theme(legend.background=element_blank()) #Remove grey key bg and white legend bg
plot <- plot + theme(legend.position = "bottom")
plot <- plot + geom_hline(yintercept=0, color="black")
plot <- plot + scale_y_continuous(breaks=seq(-150,150,75)) #Sets the axis range and interval ("by")
plot <- plot + coord_cartesian(ylim = c(-150, 160), xlim = c(2000, 2020)) #Forces the limits to be set accordingly
plot.Precip_Evap.Anomalies.Annual <- plot

#save graph to output as high resoultion, Windows compatible format .EMF
fileout <- paste("Line-Final-ARB Annual Total Precip_Evap Anomalies","_Bar",".emf",sep="")
x11(width=7, height=4) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
dev.cur()
plot.Precip_Evap.Anomalies.Annual 
savePlot(fileout, "emf", device = dev.cur())
dev.off()
graphics.off()

#Write VPRM Annual Anomalies to csv file and Move the excel .csv files to the stats folder
write.csv(VPRM_Precip_Evap.Annual.anomaly.Final, file = "VPRM_Precip_Evap_Annual_anomaly_Final.csv", row.names = F)
#file.move("VPRM_Precip_Evap_Annual_anomaly_Final.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/ARB_FINAL/Stats")



##################################################################Anomalies - Annual SoilMoisture################################################################
##Subset the relevant columns for aggregation
variables.SoilMoisture.annual <- c("SoilMoisture", "Year")
VPRM_SoilMoisture.Annual.subset <- VPRM_Fluxes.Annual[variables.SoilMoisture.annual]
VPRM_SoilMoisture.Annual.Mean <- aggregate(.~ Year, data = VPRM_SoilMoisture.Annual.subset, mean, na.action = na.pass)

#Calculate mean (climatology) for the study period (20 years)
#VPRM_SoilMoisture.climatology <- subset(VPRM_SoilMoisture.Annual.Mean, Year %in% 2000:2019)
VPRM_SoilMoisture.climatology.mean <- mean(VPRM_SoilMoisture.Annual.Mean$SoilMoisture)

VPRM_SoilMoisture.climatology.mean

##Create data for Anomaly timeseries. Start by creating a duplicate copy of mean annual dataframe
VPRM_SoilMoisture.Annual.Mean.new <- VPRM_SoilMoisture.Annual.Mean

#Add new columns with the mean climatology
VPRM_SoilMoisture.Annual.Mean.new["Climatology_SoilMoisture"] <- VPRM_SoilMoisture.climatology.mean

#Calculate anomalies and add new columns
VPRM_SoilMoisture.Annual.anomaly <- VPRM_SoilMoisture.Annual.Mean.new$SoilMoisture - VPRM_SoilMoisture.Annual.Mean.new$Climatology_SoilMoisture

VPRM_SoilMoisture.Annual.Mean.new["Anomaly_SoilMoisture"] <- VPRM_SoilMoisture.Annual.anomaly

VPRM_SoilMoisture.Annual.Mean.new ["Year_Char"] <- as.character(VPRM_SoilMoisture.Annual.Mean.new$Year) #Year column is in class ordered or factor. Need to change to Numeric
VPRM_SoilMoisture.Annual.Mean.new$Year <- NULL #Drop or remove the original ordered/factor Year column
VPRM_SoilMoisture.Annual.Mean.new ["Year"] <- as.numeric(VPRM_SoilMoisture.Annual.Mean.new$Year_Char)
VPRM_SoilMoisture.Annual.Mean.new$Year_Char <- NULL #Drop or remove the Character Year column 

##Subset the relevant columns to plot anomalies
variables.SoilMoisture.annual.anomalies <- c("Year", "Anomaly_SoilMoisture")
VPRM_SoilMoisture.Annual.anomaly.Final <- VPRM_SoilMoisture.Annual.Mean.new[variables.SoilMoisture.annual.anomalies]

#############Statistical tests for Trends################
# Install and load the packages
#install.packages(c("car", "Kendall", "zyp", "lmtest"))
library(car)
library(Kendall)
library(zyp)
library(lmtest)

#Test for Autocorrelation - durbinWatsonTest
y1.SoilMoisture <- VPRM_SoilMoisture.Annual.anomaly.Final$Anomaly_SoilMoisture
x1 <- VPRM_SoilMoisture.Annual.anomaly.Final$Year
dwtest(y1.SoilMoisture ~ x1)

#Run the Mann-Kendall test
y.SoilMoisture <- VPRM_SoilMoisture.Annual.anomaly.Final$Anomaly_SoilMoisture
Trend.Annual.SoilMoisture <- MannKendall(y.SoilMoisture)
Trend.Annual.SoilMoisture

#TheilSen Slope Estimation
y.SoilMoisture <- VPRM_SoilMoisture.Annual.anomaly.Final$Anomaly_SoilMoisture
x <- VPRM_SoilMoisture.Annual.anomaly.Final$Year
Slope.Annual.SoilMoisture <- zyp.sen(y.SoilMoisture~x)
Slope.Annual.SoilMoisture

#Subset and create objects for trend outputs to add to plots
Trend.Annual.SoilMoisture
Slope.Annual.SoilMoisture
Trend.Annual.SoilMoisture.tau <- format(round(Trend.Annual.SoilMoisture$tau[1], 2), nsmall = 2) #Round to two decimals
Trend.Annual.SoilMoisture.pvalue <- format(round(Trend.Annual.SoilMoisture$sl[1], 2), nsmall = 2) #Round to two decimals
Slope.Annual.SoilMoisture.20 <- format(round(Slope.Annual.SoilMoisture$coefficients[["x"]]*20, 2), nsmall = 2) #Round to two decimals

###########Create plots of the annual flux anomalies for the growing season
library(ggplot2)
library(ggpubr)
library(ggpmisc)
VPRM_SoilMoisture.Annual.anomaly.Final ["pos"] <- (VPRM_SoilMoisture.Annual.anomaly.Final$Anomaly_SoilMoisture >=0) #'pos' is a new column specifying +ve or -ve values as TRUE or FALSE
plot <- ggplot(VPRM_SoilMoisture.Annual.anomaly.Final, aes(x=Year, y=Anomaly_SoilMoisture, group=1)) + geom_bar(stat = "identity", aes(fill = pos), position = "dodge", width = 0.7) + scale_fill_manual(values = c("#A9A9A9","black"), guide = FALSE) + geom_smooth(method = "lm", colour = "red", se = FALSE)
plot <- plot + annotate("text", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.05, size = 3, fontface = "bold", label = (paste("Tau = ", Trend.Annual.SoilMoisture.tau, ", p = ", Trend.Annual.SoilMoisture.pvalue, ", TSA = ", Slope.Annual.SoilMoisture.20, sep="")), color = "black")
plot <- plot + xlab("Year") + ylab(expression("Soil Moisture (mm)")) + ggtitle(paste("ARB Annual Soil Moisture Anomalies" ,"(June-October, 2000-2019)"))
plot <- plot + scale_x_continuous(breaks = seq(from=2000,to=2020,by=2))
#plot <- plot + scale_x_continuous(labels =c("June","July","August","September","October"))
plot <- plot + theme(plot.title = element_text(hjust = 0.5)) #center title
plot <- plot + theme(axis.text=element_text(size=11, face="bold"))
plot <- plot + theme(axis.title=element_text(size=11))
plot <- plot + theme(plot.title=element_text(size=11,face="bold"))
#plot <- plot + theme(legend.key.size = unit(1.0, "cm"),legend.key.width = unit(1.0,"cm")) 
plot <- plot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank())
plot <- plot + labs(shape = "") #Remove or specify Legend title
plot <- plot + theme(legend.title = element_blank()) #Remove Legend title
plot <- plot + theme(legend.text = element_text(size = 11, face = "bold"))
plot <- plot + theme(legend.position = "bottom")
plot <- plot + geom_hline(yintercept=0, color="black")
plot <- plot + scale_y_continuous(breaks=seq(-100,100,50)) #Sets the axis range and interval ("by")
plot <- plot + coord_cartesian(ylim = c(-100, 100), xlim = c(2000, 2020)) #Forces the limits to be set accordingly
plot.SoilMoisture.Anomalies.Annual <- plot

#save graph to output as high resoultion, Windows compatible format .EMF
fileout <- paste("Line-Final-ARB Annual Soil Moisture Anomalies","_Bar",".emf",sep="")
x11(width=7, height=4) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
dev.cur()
plot.SoilMoisture.Anomalies.Annual 
savePlot(fileout, "emf", device = dev.cur())
dev.off()
graphics.off()


#Write VPRM Annual Anomalies to csv file and Move the excel .csv files to the stats folder
write.csv(VPRM_SoilMoisture.Annual.anomaly.Final, file = "VPRM_SoilMoisture_Annual_anomaly_Final.csv", row.names = F)
#file.move("VPRM_SoilMoisture_Annual_anomaly_Final.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/ARB_FINAL/Stats")



##################################################################Anomalies - Annual VPD################################################################
##Subset the relevant columns for aggregation
variables.VPD.annual <- c("VPD", "Year")
VPRM_VPD.Annual.subset <- VPRM_Fluxes.Annual[variables.VPD.annual]
VPRM_VPD.Annual.Mean <- aggregate(.~ Year, data = VPRM_VPD.Annual.subset, mean, na.action = na.pass)

#Calculate mean (climatology) for the study period (20 years)
#VPRM_VPD.climatology <- subset(VPRM_VPD.Annual.Mean, Year %in% 2000:2019)
VPRM_VPD.climatology.mean <- mean(VPRM_VPD.Annual.Mean$VPD)

VPRM_VPD.climatology.mean

##Create data for Anomaly timeseries. Start by creating a duplicate copy of mean annual dataframe
VPRM_VPD.Annual.Mean.new <- VPRM_VPD.Annual.Mean

#Add new columns with the mean climatology
VPRM_VPD.Annual.Mean.new["Climatology_VPD"] <- VPRM_VPD.climatology.mean

#Calculate anomalies and add new columns
VPRM_VPD.Annual.anomaly <- VPRM_VPD.Annual.Mean.new$VPD - VPRM_VPD.Annual.Mean.new$Climatology_VPD

VPRM_VPD.Annual.Mean.new["Anomaly_VPD"] <- VPRM_VPD.Annual.anomaly

VPRM_VPD.Annual.Mean.new ["Year_Char"] <- as.character(VPRM_VPD.Annual.Mean.new$Year) #Year column is in class ordered or factor. Need to change to Numeric
VPRM_VPD.Annual.Mean.new$Year <- NULL #Drop or remove the original ordered/factor Year column
VPRM_VPD.Annual.Mean.new ["Year"] <- as.numeric(VPRM_VPD.Annual.Mean.new$Year_Char)
VPRM_VPD.Annual.Mean.new$Year_Char <- NULL #Drop or remove the Character Year column 

##Subset the relevant columns to plot anomalies
variables.VPD.annual.anomalies <- c("Year", "Anomaly_VPD")
VPRM_VPD.Annual.anomaly.Final <- VPRM_VPD.Annual.Mean.new[variables.VPD.annual.anomalies]

#############Statistical tests for Trends################
# Install and load the packages
#install.packages(c("car", "Kendall", "zyp", "lmtest"))
library(car)
library(Kendall)
library(zyp)
library(lmtest)

#Test for Autocorrelation - durbinWatsonTest
y1.VPD <- VPRM_VPD.Annual.anomaly.Final$Anomaly_VPD
x1 <- VPRM_VPD.Annual.anomaly.Final$Year
dwtest(y1.VPD ~ x1)

#Run the Mann-Kendall test
y.VPD <- VPRM_VPD.Annual.anomaly.Final$Anomaly_VPD
Trend.Annual.VPD <- MannKendall(y.VPD)
Trend.Annual.VPD

#TheilSen Slope Estimation
y.VPD <- VPRM_VPD.Annual.anomaly.Final$Anomaly_VPD
x <- VPRM_VPD.Annual.anomaly.Final$Year
Slope.Annual.VPD <- zyp.sen(y.VPD~x)
Slope.Annual.VPD

#Subset and create objects for trend outputs to add to plots
Trend.Annual.VPD
Slope.Annual.VPD
Trend.Annual.VPD.tau <- format(round(Trend.Annual.VPD$tau[1], 2), nsmall = 2) #Round to two decimals
Trend.Annual.VPD.pvalue <- format(round(Trend.Annual.VPD$sl[1], 2), nsmall = 2) #Round to two decimals
Slope.Annual.VPD.20 <- format(round(Slope.Annual.VPD$coefficients[["x"]]*20, 2), nsmall = 2) #Round to two decimals

###########Create plots of the annual flux anomalies for the growing season
library(ggplot2)
library(ggpubr)
library(ggpmisc)
VPRM_VPD.Annual.anomaly.Final ["pos"] <- (VPRM_VPD.Annual.anomaly.Final$Anomaly_VPD >=0) #'pos' is a new column specifying +ve or -ve values as TRUE or FALSE
plot <- ggplot(VPRM_VPD.Annual.anomaly.Final, aes(x=Year, y=Anomaly_VPD, group=1)) + geom_bar(stat = "identity", aes(fill = pos), position = "dodge", width = 0.7) + scale_fill_manual(values = c("#A9A9A9","black"), guide = FALSE) + geom_smooth(method = "lm", colour = "red", se = FALSE)
plot <- plot + annotate("text", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.05, size = 3, fontface = "bold", label = (paste("Tau = ", Trend.Annual.VPD.tau, ", p = ", Trend.Annual.VPD.pvalue, ", TSA = ", Slope.Annual.VPD.20, sep="")), color = "black")
plot <- plot + xlab("Year") + ylab(expression("VPD (kPa)")) + ggtitle(paste("ARB Annual VPD Anomalies" ,"(June-October, 2000-2019)"))
plot <- plot + scale_x_continuous(breaks = seq(from=2000,to=2020,by=2))
#plot <- plot + scale_x_continuous(labels =c("June","July","August","September","October"))
plot <- plot + theme(plot.title = element_text(hjust = 0.5)) #center title
plot <- plot + theme(axis.text=element_text(size=11, face="bold"))
plot <- plot + theme(axis.title=element_text(size=11))
plot <- plot + theme(plot.title=element_text(size=11,face="bold"))
#plot <- plot + theme(legend.key.size = unit(1.0, "cm"),legend.key.width = unit(1.0,"cm")) 
plot <- plot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank())
plot <- plot + labs(shape = "") #Remove or specify Legend title
plot <- plot + theme(legend.title = element_blank()) #Remove Legend title
plot <- plot + theme(legend.text = element_text(size = 11, face = "bold"))
plot <- plot + theme(legend.position = "bottom")
plot <- plot + geom_hline(yintercept=0, color="black")
plot <- plot + scale_y_continuous(breaks=seq(-0.10,0.10,0.05)) #Sets the axis range and interval ("by")
plot <- plot + coord_cartesian(ylim = c(-0.10, 0.10), xlim = c(2000, 2020)) #Forces the limits to be set accordingly
plot.VPD.Anomalies.Annual <- plot

#save graph to output as high resoultion, Windows compatible format .EMF
fileout <- paste("Line-Final-ARB Annual VPD Anomalies","_Bar",".emf",sep="")
x11(width=7, height=4) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
dev.cur()
plot.VPD.Anomalies.Annual
savePlot(fileout, "emf", device = dev.cur())
dev.off()
graphics.off()

#Write VPRM Annual Anomalies to csv file and Move the excel .csv files to the stats folder
write.csv(VPRM_VPD.Annual.anomaly.Final, file = "VPRM_VPD_Annual_anomaly_Final.csv", row.names = F)
#file.move("VPRM_VPD_Annual_anomaly_Final.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/ARB_FINAL/Stats")

cat("\n"); 
VPRMstate = "Final Model Run and Data Analysis Completed for ARB!"
cat("\n"); 
cat(VPRMstate);cat("\n"); cat("\n")

