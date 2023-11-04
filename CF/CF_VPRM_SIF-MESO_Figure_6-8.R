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
CF_NARR_met.SIF <- read.csv("CF-NARR-Final-ready.csv", stringsAsFactors = FALSE)
CF_evi_modis.SIF <- read.csv("CF-SIF-ready.csv", stringsAsFactors = FALSE)
CF_refl_modis.SIF <- read.csv("CF-Reflectance-ready.csv", stringsAsFactors = FALSE)
CF_phen_modis.SIF <- read.csv("CF-Phenology-MAX-ready.csv", stringsAsFactors = FALSE)

#Convert "Character" dates to "Chron" dates and times for each VPRM driver data
library(chron)

CF_NARR_met.SIF.dates <- t(as.data.frame(strsplit(CF_NARR_met.SIF[["Date"]],' '))) #t is (matrix transpose)
row.names(CF_NARR_met.SIF.dates) = NULL
CF_NARR_met.SIF ["Date"] <- chron(dates=CF_NARR_met.SIF.dates[,1],times=CF_NARR_met.SIF.dates[,2], format = c(dates = "m/d/y", times = "h:m:s"))
class(CF_NARR_met.SIF[["Date"]])

CF_evi_modis.SIF.dates <- t(as.data.frame(strsplit(CF_evi_modis.SIF[["Date"]],' '))) #t is (matrix transpose)
row.names(CF_evi_modis.SIF.dates) = NULL
CF_evi_modis.SIF ["Date"] <- chron(dates=CF_evi_modis.SIF.dates[,1],times=CF_evi_modis.SIF.dates[,2], format = c(dates = "m/d/y", times = "h:m:s"))
class(CF_evi_modis.SIF[["Date"]])

CF_refl_modis.SIF.dates <- t(as.data.frame(strsplit(CF_refl_modis.SIF[["Date"]],' '))) #t is (matrix transpose)
row.names(CF_refl_modis.SIF.dates) = NULL
CF_refl_modis.SIF ["Date"] <- chron(dates=CF_refl_modis.SIF.dates[,1],times=CF_refl_modis.SIF.dates[,2], format = c(dates = "m/d/y", times = "h:m:s"))
class(CF_refl_modis.SIF[["Date"]])

####Melt and rearrange Phenology data and Convert "Character" dates to "Chron" dates and times for each VPRM driver data
CF_phen.melted.SIF <- reshape2::melt(CF_phen_modis.SIF, id.var='Date')
names(CF_phen.melted.SIF) <- c("var", "phen", "date")
head(CF_phen.melted.SIF)
tail(CF_phen.melted.SIF)

CF_phen.melted.SIF.dates <- t(as.data.frame(strsplit(CF_phen.melted.SIF[["date"]],' '))) #t is (matrix transpose)
row.names(CF_phen.melted.SIF.dates) = NULL
CF_phen.melted.SIF ["date"] <- chron(dates=CF_phen.melted.SIF.dates[,1],times=CF_phen.melted.SIF.dates[,2], format = c(dates = "m/d/y", times = "h:m:s"))
class(CF_phen.melted.SIF[["date"]])

###Subset the relevant columns for aggregation
phen.varselect.SIF<- c("date", "phen")
CF_phen.subset.SIF <- CF_phen.melted.SIF[phen.varselect.SIF]

CF_phen_rearranged.SIF <- CF_phen.subset.SIF[order(as.Date(CF_phen.subset.SIF$date, format="%m/%d/%y")),]
head(CF_phen_rearranged.SIF)
tail(CF_phen_rearranged.SIF)

## determine the phenology phase for each tower observation date
phen_filled.SIF <- interp_phenology(CF_phen_rearranged.SIF, CF_NARR_met.SIF[['Date']])



#####Read csv file from original optimized run####
#####Read csv file from original optimized run####
specify_path <- (dirname(rstudioapi::getActiveDocumentContext()$path))
filename_parameters <- paste(specify_path, "/CF_Optimized_Site_Parameter_Values",".csv",sep="")
SIF_optimized_parameters <- read.csv(filename_parameters)


######Assemble VPRM driver data frame

## determine the phenology phase for each tower observation date
#phen_filled.SIF <- interp_phenology(CF_phen, CF_NARR_met.SIF[['date']])
## place the tower observations and MODIS data into a VPRM_driver_data object.
CF_dd.SIF <- VPRM_driver_data(name_long="Churchill MB",
                               name_short = "CF",
                               lat=58.6658,
                               lon=-93.83,
                               PFT='WET',    ## Permanent Wetland
                               tower_date=CF_NARR_met.SIF[['Date']],
                               NEE_obs=NA,
                               T=CF_NARR_met.SIF[['TA']],
                               PAR=CF_NARR_met.SIF[['PAR']],
                               date_nir = CF_refl_modis.SIF[['Date']],
                               rho_nir=CF_refl_modis.SIF[['NIR']],
                               date_swir = CF_refl_modis.SIF[['Date']],
                               rho_swir = CF_refl_modis.SIF[['SWIR']],
                               date_EVI = CF_evi_modis.SIF[['Date']],
                               EVI=CF_evi_modis.SIF[['SIF']],
                               phen=phen_filled.SIF)

## take a look at the result
print(head(as.data.frame(CF_dd.SIF)))
print(tail(as.data.frame(CF_dd.SIF)))


######Calculate VPRM NEE#####
#Calculating VPRM fluxes requires values for the parameters lambda, PAR0, alpha, and beta.
#The dataset VPRM_parameters, part of the VPRMLandSfcModel package, includes several parameters sets.

attach(SIF_optimized_parameters)

CF_dd.SIF[['data']][['NEE_VPRM']] <- vprm_calc_NEE(
  CF_dd.SIF, lambda=lambda, PAR_0=PAR_0, alpha=alpha, beta=beta)

######Calculate VPRM GPP#####
CF_dd.SIF[['data']][['GPP_VPRM']] <- vprm_calc_GEE(
  CF_dd.SIF, lambda=lambda, PAR_0=PAR_0)

######Calculate VPRM ER#####
CF_dd.SIF[['data']][['ER_VPRM']] <- vprm_calc_R(
  CF_dd.SIF, alpha=alpha, beta=beta)



##################################################Diagnosis#############################################
VPRM_Fluxes.SIF <- data.frame(Date = CF_dd.SIF[["data"]][["date"]], NEE_VPRM = CF_dd.SIF[["data"]][["NEE_VPRM"]], GPP_VPRM = CF_dd.SIF[["data"]][["GPP_VPRM"]], ER_VPRM = CF_dd.SIF[["data"]][["ER_VPRM"]], Temperature = CF_dd.SIF[["data"]][["T"]], PAR = CF_dd.SIF[["data"]][["PAR"]])
#NEE.Diff.SIF.2 <- data.frame(Date = fig_dNEE.SIF[["data"]][["date"]], NEE_Residual = fig_dNEE.SIF[["data"]][["NEE_VPRM"]] - fig_dNEE.SIF[["data"]][["NEE_obs"]])

###Subset relevant columns for NARR meteorological variables needed for analysis
variables.NARR.met <- c("Shore", "WindSpeed", "Precipitation","Humidity", "Evaporation", "SoilMoisture")
NARR.met.subset <- CF_NARR_met.SIF[variables.NARR.met]

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


#################################BOXPLOTS###########################
library(gdata)
require(reshape2)
library(ggplot2)



##Bar Graph (METHOD 2)
library(plyr)
library(dplyr)
library(Rmisc)
bar.variables <- c("Date", "GPP_VPRM", "ER_VPRM", "NEE_VPRM", "Shore")
bar.data.frame <- VPRM_Fluxes.Met.Study[bar.variables]
mydata.m.bar <- melt(bar.data.frame, id.var = c("Date", "Shore"))
mydata.m.bar_summary <- summarySE(mydata.m.bar, measurevar="value", groupvars=c("Shore","variable"))
plot <- ggplot(data = mydata.m.bar_summary, aes(x=variable, y=value, fill=Shore)) + geom_bar(stat = "identity", position = "dodge", width = 0.7)
plot <- plot + geom_errorbar(aes(ymin=value-se, ymax=value+se), size=.3, width=.1, position=position_dodge(0.7))
plot <- plot + xlab("Variables") + ylab(expression(CO[2]~Fluxes~(mu*mol~m^{-2}~s^{-1}))) + ggtitle(paste("CF Mean Fluxes " ,"(June-October, 2000-2019)"))
plot <- plot + scale_fill_manual(values = c("#FF3333", "#0099FF"))
plot <- plot + theme(plot.title = element_text(hjust = 0.5)) #center title
plot <- plot + theme(axis.text=element_text(size=11, face="bold"))
plot <- plot + theme(axis.title=element_text(size=11))
plot <- plot + theme(plot.title=element_text(size=11,face="bold"))
plot <- plot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank())
plot <- plot + labs(shape = "") #Remove or specify Legend title
plot <- plot + theme(legend.title = element_blank()) #Remove Legend title
plot <- plot + theme(legend.text = element_text(size = 11, face = "bold"))
plot <- plot + theme(legend.position = "right")
plot <- plot + geom_hline(yintercept=0, color="black")
plot <- plot + scale_y_continuous(breaks=seq(0,2,0.5), expand = c(0, 0)) #Sets the axis range and interval ("by")
plot <- plot + coord_cartesian(ylim = c(0, 2)) #Forces the limits to be set accordingly
#save graph to output as high resoultion, Windows compatible format .EMF
fileout <- paste("Meso-CF Mean Fluxes","_Bar",".emf",sep="")
x11(width=7, height=4) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
dev.cur()
plot
savePlot(fileout, "emf", device = dev.cur())
dev.off()
graphics.off()



#Write Wind Frequency mean and Anomalies to csv file and Move the excel .csv files to the stats folder
library(filesstrings)
write.csv(mydata.m.bar_summary, file = "Meso_CO2Fluxes_Summary_Stats.csv", row.names = F)
#file.move("Meso_CO2Fluxes_Summary_Stats.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/CF_MESO/Stats")



#####Calculate Student Tests and create output dataframe - All
library("dplyr")
library("ggpubr")
OFF.bar.data.frame <- subset(bar.data.frame, Shore=="OFF")
ON.bar.data.frame <- subset(bar.data.frame, Shore=="ON")



########################################################################################################################################################
##################################################################Diurnal Cycle Analysis################################################################

#Create and add new columns for day, month, year, and time; extracted from the original date column
VPRM_Fluxes.Diurnal <- VPRM_Fluxes.Met.Study #Duplicate dataframe to avoid potential conflicts
VPRM_Fluxes.Diurnal ["Day"] <- as.data.frame(days(VPRM_Fluxes.Diurnal$Date))
VPRM_Fluxes.Diurnal ["Month"] <- as.data.frame(months.Date(VPRM_Fluxes.Diurnal$Date))
VPRM_Fluxes.Diurnal ["Year"] <- as.data.frame(years(VPRM_Fluxes.Diurnal$Date))
#VPRM_Fluxes.Diurnal ["Hours"] <- as.data.frame(hours(VPRM_Fluxes.Diurnal$Date))
#VPRM_Fluxes.Diurnal ["Minutes"] <- as.data.frame(minutes(VPRM_Fluxes.Diurnal$Date))
#VPRM_Fluxes.Diurnal ["Seconds"] <- as.data.frame(seconds(VPRM_Fluxes.Diurnal$Date))
VPRM_Fluxes.Diurnal ["Time"] <- data.frame(substr(VPRM_Fluxes.Diurnal[["Date"]], start = 11, stop = 15))





# NEW LINES STARTS - MESO: WIND FREQUENCY --------------------------------------------------- ###########################################################

################BEGIN ANNUAL ANALYSES###
library(dplyr)
VPRM.freq.var <- c("Shore", "Year")
VPRM.freq.shore <- VPRM_Fluxes.Diurnal[VPRM.freq.var]
VPRM.freq.ON <- subset(VPRM.freq.shore, Shore=="ON")
VPRM.freq.OFF <- subset(VPRM.freq.shore, Shore=="OFF")

####Onshore
freq.year <- 2000:2019

for (year in freq.year) {
  
  i <- year
  
  cat("\n"); 
  #cat("\n"); cat("\n")
  startstate = paste("Simulation Beginning for Year", i, "!", sep="")
  cat(startstate);cat("\n"); cat("\n")
  
  
  run.out.current<-subset(VPRM.freq.ON, Year==i)
  count.freq <- count (run.out.current)
  assign(paste("out.",i,sep=""),count.freq)
  
  
}

ON.freq.df <- rbind(out.2000,out.2001,out.2002,out.2003,out.2004,out.2005,out.2006,out.2007,out.2008,out.2009,out.2010,out.2011,out.2012,out.2013,out.2014,out.2015,out.2016,out.2017,out.2018,out.2019)
colnames(ON.freq.df) <- c("ON")

####Offshore
freq.year <- 2000:2019

for (year in freq.year) {
  
  i <- year
  
  cat("\n"); 
  #cat("\n"); cat("\n")
  startstate = paste("Simulation Beginning for Year", i, "!", sep="")
  cat(startstate);cat("\n"); cat("\n")
  
  
  run.out.current<-subset(VPRM.freq.OFF, Year==i)
  count.freq <- count (run.out.current)
  assign(paste("out.",i,sep=""),count.freq)
  
  
}

OFF.freq.df <- rbind(out.2000,out.2001,out.2002,out.2003,out.2004,out.2005,out.2006,out.2007,out.2008,out.2009,out.2010,out.2011,out.2012,out.2013,out.2014,out.2015,out.2016,out.2017,out.2018,out.2019)
colnames(OFF.freq.df) <- c("OFF")

Year.df <- data.frame(c(2000:2019))
colnames(Year.df) <- c("Year")

####Combine Offshore and Onshore
final.freq.df <-  data.frame(Year.df,OFF.freq.df,ON.freq.df)

##Create new data for Percentage dataframe
##Start by creating a duplicate copy of final.freq dataframe
freq.percent <- final.freq.df

#Calculate and add new column with the total frequency
freq.percent["Total"] <- final.freq.df$OFF + final.freq.df$ON

#Calculate and add new column with the percentages
freq.percent["OFF.per"] <- (freq.percent$OFF/freq.percent$Total)*100
freq.percent["ON.per"] <- (freq.percent$ON/freq.percent$Total)*100
freq.percent

#Select OFF.per and ON.per columns only and Remove all other columns
select.per.variables <- c("Year", "OFF.per", "ON.per")
freq.percent.final <- freq.percent[select.per.variables]
freq.percent.final
#Rename colnames for meaningful Legend
names(freq.percent.final) <- c("Year","OFF","ON")
freq.percent.final





#####################################################################Anomalies###

##Create Separate dataframe for ON and OFF
extract.cols.ON <- c("Year", "ON.per")
ON.final.df <- freq.percent[extract.cols.ON]
#Rename colnames
names(ON.final.df) <- c("Year","ON")

extract.cols.OFF <- c("Year", "OFF.per")
OFF.final.df <- freq.percent[extract.cols.OFF]
#Rename colnames
names(OFF.final.df) <- c("Year","OFF")

ON.freq.climatology.mean <- mean(ON.final.df$ON)
OFF.freq.climatology.mean <- mean(OFF.final.df$OFF)

##print to console
ON.freq.climatology.mean
OFF.freq.climatology.mean

##Create data for Anomaly timeseries
##Start by creating a duplicate copy of mean annual dataframe
ON.mean.annual.new <- ON.final.df
OFF.mean.annual.new <- OFF.final.df

#Add new column with the mean climatology
ON.mean.annual.new["Climatology"] <- ON.freq.climatology.mean
OFF.mean.annual.new["Climatology"] <- OFF.freq.climatology.mean

#Calculate anomaly and add new column
ON.mean.annual.new["Anomaly"] <- ON.mean.annual.new$ON - ON.mean.annual.new$Climatology
OFF.mean.annual.new["Anomaly"] <- OFF.mean.annual.new$OFF - OFF.mean.annual.new$Climatology



####################################################STATS ANALYSES###
#R code for Statistical tests

# Install and load the packages
#install.packages("forecast")
library("forecast")

#install.packages(c("car", "Kendall", "zyp", "lmtest"))
library(car)
library(Kendall)
library(zyp)
library(lmtest)




################BEGIN SEPARATE MONTHLY ANALYSES###
library(dplyr)
VPRM.freq.var.month <- c("Shore", "Month", "Year")
VPRM.freq.shore.month <- VPRM_Fluxes.Diurnal[VPRM.freq.var.month]
VPRM.freq.ON.month <- subset(VPRM.freq.shore.month, Shore=="ON")
VPRM.freq.ON.month$Month = match(VPRM.freq.ON.month$Month,(month.name))
VPRM.freq.OFF.month <- subset(VPRM.freq.shore.month, Shore=="OFF")
VPRM.freq.OFF.month$Month = match(VPRM.freq.OFF.month$Month,(month.name))

VPRM.freq.ON.June <- subset(VPRM.freq.ON.month, Month=="6")
VPRM.freq.ON.July <- subset(VPRM.freq.ON.month, Month=="7")
VPRM.freq.ON.August <- subset(VPRM.freq.ON.month, Month=="8")
VPRM.freq.ON.September <- subset(VPRM.freq.ON.month, Month=="9")
VPRM.freq.ON.October <- subset(VPRM.freq.ON.month, Month=="10")

VPRM.freq.OFF.June <- subset(VPRM.freq.OFF.month, Month=="6")
VPRM.freq.OFF.July <- subset(VPRM.freq.OFF.month, Month=="7")
VPRM.freq.OFF.August <- subset(VPRM.freq.OFF.month, Month=="8")
VPRM.freq.OFF.September <- subset(VPRM.freq.OFF.month, Month=="9")
VPRM.freq.OFF.October <- subset(VPRM.freq.OFF.month, Month=="10")



###############################JUNE############################################################################
#####Onshore for June
VPRM.freq.year <- 2000:2019

for (year in VPRM.freq.year) {
  
  i <- year
  
  cat("\n"); 
  #cat("\n"); cat("\n")
  startstate = paste("Simulation Beginning for Year", i, "!", sep="")
  cat(startstate);cat("\n"); cat("\n")
  
  
  output.June.ON<-subset(VPRM.freq.ON.June, Year==i)
  count.freq <- count (output.June.ON)
  assign(paste("out.",i,sep=""),count.freq)
  
  
  
}

ON.VPRM.freq.df.June <- rbind(out.2000,out.2001,out.2002,out.2003,out.2004,out.2005,out.2006,out.2007,out.2008,out.2009,out.2010,out.2011,out.2012,out.2013,out.2014,out.2015,out.2016,out.2017,out.2018,out.2019)
colnames(ON.VPRM.freq.df.June) <- c("ON")

#####Offshore for June
VPRM.freq.year <- 2000:2019

for (year in VPRM.freq.year) {
  
  i <- year
  
  cat("\n"); 
  #cat("\n"); cat("\n")
  startstate = paste("Simulation Beginning for Year", i, "!", sep="")
  cat(startstate);cat("\n"); cat("\n")
  
  
  output.June.OFF<-subset(VPRM.freq.OFF.June, Year==i)
  count.freq <- count (output.June.OFF)
  assign(paste("out.",i,sep=""),count.freq)
  
  
}

OFF.VPRM.freq.df.June <- rbind(out.2000,out.2001,out.2002,out.2003,out.2004,out.2005,out.2006,out.2007,out.2008,out.2009,out.2010,out.2011,out.2012,out.2013,out.2014,out.2015,out.2016,out.2017,out.2018,out.2019)
colnames(OFF.VPRM.freq.df.June) <- c("OFF")

Year.df <- data.frame(c(2000:2019))
colnames(Year.df) <- c("Year")

#####Combine Offshore and Onshore for June
final.VPRM.freq.df.June <-  data.frame(Year.df,OFF.VPRM.freq.df.June,ON.VPRM.freq.df.June)

##Create new data for Percentage dataframe
##Start by creating a duplicate copy of final.freq dataframe
VPRM.freq.percent.June <- final.VPRM.freq.df.June

#Calculate and add new column with the total frequency
Total.col <- final.VPRM.freq.df.June$OFF + final.VPRM.freq.df.June$ON
VPRM.freq.percent.June["Total"] <- Total.col

#Calculate and add new column with the percentages
OFF.per.col <- (VPRM.freq.percent.June$OFF/VPRM.freq.percent.June$Total)*100
ON.per.col <- (VPRM.freq.percent.June$ON/VPRM.freq.percent.June$Total)*100
VPRM.freq.percent.June["OFF.per"] <- OFF.per.col
VPRM.freq.percent.June["ON.per"] <- ON.per.col
VPRM.freq.percent.June

#Select OFF.per and ON.per columns only and Remove all other columns
select.variables <- c("Year", "OFF.per", "ON.per")
VPRM.freq.percent.June.final <- VPRM.freq.percent.June[select.variables]
VPRM.freq.percent.June.final
#Rename colnames for meaningful Legend
names(VPRM.freq.percent.June.final) <- c("Year","OFF","ON")
VPRM.freq.percent.June.final



###############################JULY############################################################################
#####Onshore for July
VPRM.freq.year <- 2000:2019

for (year in VPRM.freq.year) {
  
  i <- year
  
  cat("\n"); 
  #cat("\n"); cat("\n")
  startstate = paste("Simulation Beginning for Year", i, "!", sep="")
  cat(startstate);cat("\n"); cat("\n")
  
  
  output.July.ON<-subset(VPRM.freq.ON.July, Year==i)
  count.freq <- count (output.July.ON)
  assign(paste("out.",i,sep=""),count.freq)
  
  
  
}

ON.VPRM.freq.df.July <- rbind(out.2000,out.2001,out.2002,out.2003,out.2004,out.2005,out.2006,out.2007,out.2008,out.2009,out.2010,out.2011,out.2012,out.2013,out.2014,out.2015,out.2016,out.2017,out.2018,out.2019)
colnames(ON.VPRM.freq.df.July) <- c("ON")

#####Offshore for July
VPRM.freq.year <- 2000:2019

for (year in VPRM.freq.year) {
  
  i <- year
  
  cat("\n"); 
  #cat("\n"); cat("\n")
  startstate = paste("Simulation Beginning for Year", i, "!", sep="")
  cat(startstate);cat("\n"); cat("\n")
  
  
  output.July.OFF<-subset(VPRM.freq.OFF.July, Year==i)
  count.freq <- count (output.July.OFF)
  assign(paste("out.",i,sep=""),count.freq)
  
  
}

OFF.VPRM.freq.df.July <- rbind(out.2000,out.2001,out.2002,out.2003,out.2004,out.2005,out.2006,out.2007,out.2008,out.2009,out.2010,out.2011,out.2012,out.2013,out.2014,out.2015,out.2016,out.2017,out.2018,out.2019)
colnames(OFF.VPRM.freq.df.July) <- c("OFF")

Year.df <- data.frame(c(2000:2019))
colnames(Year.df) <- c("Year")

#####Combine Offshore and Onshore for July
final.VPRM.freq.df.July <-  data.frame(Year.df,OFF.VPRM.freq.df.July,ON.VPRM.freq.df.July)

##Create new data for Percentage dataframe
##Start by creating a duplicate copy of final.freq dataframe
VPRM.freq.percent.July <- final.VPRM.freq.df.July

#Calculate and add new column with the total frequency
Total.col <- final.VPRM.freq.df.July$OFF + final.VPRM.freq.df.July$ON
VPRM.freq.percent.July["Total"] <- Total.col

#Calculate and add new column with the percentages
OFF.per.col <- (VPRM.freq.percent.July$OFF/VPRM.freq.percent.July$Total)*100
ON.per.col <- (VPRM.freq.percent.July$ON/VPRM.freq.percent.July$Total)*100
VPRM.freq.percent.July["OFF.per"] <- OFF.per.col
VPRM.freq.percent.July["ON.per"] <- ON.per.col
VPRM.freq.percent.July

#Select OFF.per and ON.per columns only and Remove all other columns
select.variables <- c("Year", "OFF.per", "ON.per")
VPRM.freq.percent.July.final <- VPRM.freq.percent.July[select.variables]
VPRM.freq.percent.July.final
#Rename colnames for meaningful Legend
names(VPRM.freq.percent.July.final) <- c("Year","OFF","ON")
VPRM.freq.percent.July.final



###############################AUGUST############################################################################
#####Onshore for August
VPRM.freq.year <- 2000:2019

for (year in VPRM.freq.year) {
  
  i <- year
  
  cat("\n"); 
  #cat("\n"); cat("\n")
  startstate = paste("Simulation Beginning for Year", i, "!", sep="")
  cat(startstate);cat("\n"); cat("\n")
  
  
  output.August.ON<-subset(VPRM.freq.ON.August, Year==i)
  count.freq <- count (output.August.ON)
  assign(paste("out.",i,sep=""),count.freq)
  
  
  
}

ON.VPRM.freq.df.August <- rbind(out.2000,out.2001,out.2002,out.2003,out.2004,out.2005,out.2006,out.2007,out.2008,out.2009,out.2010,out.2011,out.2012,out.2013,out.2014,out.2015,out.2016,out.2017,out.2018,out.2019)
colnames(ON.VPRM.freq.df.August) <- c("ON")

#####Offshore for August
VPRM.freq.year <- 2000:2019

for (year in VPRM.freq.year) {
  
  i <- year
  
  cat("\n"); 
  #cat("\n"); cat("\n")
  startstate = paste("Simulation Beginning for Year", i, "!", sep="")
  cat(startstate);cat("\n"); cat("\n")
  
  
  output.August.OFF<-subset(VPRM.freq.OFF.August, Year==i)
  count.freq <- count (output.August.OFF)
  assign(paste("out.",i,sep=""),count.freq)
  
  
}

OFF.VPRM.freq.df.August <- rbind(out.2000,out.2001,out.2002,out.2003,out.2004,out.2005,out.2006,out.2007,out.2008,out.2009,out.2010,out.2011,out.2012,out.2013,out.2014,out.2015,out.2016,out.2017,out.2018,out.2019)
colnames(OFF.VPRM.freq.df.August) <- c("OFF")

Year.df <- data.frame(c(2000:2019))
colnames(Year.df) <- c("Year")

#####Combine Offshore and Onshore for August
final.VPRM.freq.df.August <-  data.frame(Year.df,OFF.VPRM.freq.df.August,ON.VPRM.freq.df.August)

##Create new data for Percentage dataframe
##Start by creating a duplicate copy of final.freq dataframe
VPRM.freq.percent.August <- final.VPRM.freq.df.August

#Calculate and add new column with the total frequency
Total.col <- final.VPRM.freq.df.August$OFF + final.VPRM.freq.df.August$ON
VPRM.freq.percent.August["Total"] <- Total.col

#Calculate and add new column with the percentages
OFF.per.col <- (VPRM.freq.percent.August$OFF/VPRM.freq.percent.August$Total)*100
ON.per.col <- (VPRM.freq.percent.August$ON/VPRM.freq.percent.August$Total)*100
VPRM.freq.percent.August["OFF.per"] <- OFF.per.col
VPRM.freq.percent.August["ON.per"] <- ON.per.col
VPRM.freq.percent.August

#Select OFF.per and ON.per columns only and Remove all other columns
select.variables <- c("Year", "OFF.per", "ON.per")
VPRM.freq.percent.August.final <- VPRM.freq.percent.August[select.variables]
VPRM.freq.percent.August.final
#Rename colnames for meaningful Legend
names(VPRM.freq.percent.August.final) <- c("Year","OFF","ON")
VPRM.freq.percent.August.final



###############################SEPTEMBER############################################################################
#####Onshore for September
VPRM.freq.year <- 2000:2019

for (year in VPRM.freq.year) {
  
  i <- year
  
  cat("\n"); 
  #cat("\n"); cat("\n")
  startstate = paste("Simulation Beginning for Year", i, "!", sep="")
  cat(startstate);cat("\n"); cat("\n")
  
  
  output.September.ON<-subset(VPRM.freq.ON.September, Year==i)
  count.freq <- count (output.September.ON)
  assign(paste("out.",i,sep=""),count.freq)
  
  
  
}

ON.VPRM.freq.df.September <- rbind(out.2000,out.2001,out.2002,out.2003,out.2004,out.2005,out.2006,out.2007,out.2008,out.2009,out.2010,out.2011,out.2012,out.2013,out.2014,out.2015,out.2016,out.2017,out.2018,out.2019)
colnames(ON.VPRM.freq.df.September) <- c("ON")

#####Offshore for September
VPRM.freq.year <- 2000:2019

for (year in VPRM.freq.year) {
  
  i <- year
  
  cat("\n"); 
  #cat("\n"); cat("\n")
  startstate = paste("Simulation Beginning for Year", i, "!", sep="")
  cat(startstate);cat("\n"); cat("\n")
  
  
  output.September.OFF<-subset(VPRM.freq.OFF.September, Year==i)
  count.freq <- count (output.September.OFF)
  assign(paste("out.",i,sep=""),count.freq)
  
  
}

OFF.VPRM.freq.df.September <- rbind(out.2000,out.2001,out.2002,out.2003,out.2004,out.2005,out.2006,out.2007,out.2008,out.2009,out.2010,out.2011,out.2012,out.2013,out.2014,out.2015,out.2016,out.2017,out.2018,out.2019)
colnames(OFF.VPRM.freq.df.September) <- c("OFF")

Year.df <- data.frame(c(2000:2019))
colnames(Year.df) <- c("Year")

#####Combine Offshore and Onshore for September
final.VPRM.freq.df.September <-  data.frame(Year.df,OFF.VPRM.freq.df.September,ON.VPRM.freq.df.September)

##Create new data for Percentage dataframe
##Start by creating a duplicate copy of final.freq dataframe
VPRM.freq.percent.September <- final.VPRM.freq.df.September

#Calculate and add new column with the total frequency
Total.col <- final.VPRM.freq.df.September$OFF + final.VPRM.freq.df.September$ON
VPRM.freq.percent.September["Total"] <- Total.col

#Calculate and add new column with the percentages
OFF.per.col <- (VPRM.freq.percent.September$OFF/VPRM.freq.percent.September$Total)*100
ON.per.col <- (VPRM.freq.percent.September$ON/VPRM.freq.percent.September$Total)*100
VPRM.freq.percent.September["OFF.per"] <- OFF.per.col
VPRM.freq.percent.September["ON.per"] <- ON.per.col
VPRM.freq.percent.September

#Select OFF.per and ON.per columns only and Remove all other columns
select.variables <- c("Year", "OFF.per", "ON.per")
VPRM.freq.percent.September.final <- VPRM.freq.percent.September[select.variables]
VPRM.freq.percent.September.final
#Rename colnames for meaningful Legend
names(VPRM.freq.percent.September.final) <- c("Year","OFF","ON")
VPRM.freq.percent.September.final



###############################OCTOBER############################################################################
#####Onshore for October
VPRM.freq.year <- 2000:2019

for (year in VPRM.freq.year) {
  
  i <- year
  
  cat("\n"); 
  #cat("\n"); cat("\n")
  startstate = paste("Simulation Beginning for Year", i, "!", sep="")
  cat(startstate);cat("\n"); cat("\n")
  
  
  output.October.ON<-subset(VPRM.freq.ON.October, Year==i)
  count.freq <- count (output.October.ON)
  assign(paste("out.",i,sep=""),count.freq)
  
  
  
}

ON.VPRM.freq.df.October <- rbind(out.2000,out.2001,out.2002,out.2003,out.2004,out.2005,out.2006,out.2007,out.2008,out.2009,out.2010,out.2011,out.2012,out.2013,out.2014,out.2015,out.2016,out.2017,out.2018,out.2019)
colnames(ON.VPRM.freq.df.October) <- c("ON")

#####Offshore for October
VPRM.freq.year <- 2000:2019

for (year in VPRM.freq.year) {
  
  i <- year
  
  cat("\n"); 
  #cat("\n"); cat("\n")
  startstate = paste("Simulation Beginning for Year", i, "!", sep="")
  cat(startstate);cat("\n"); cat("\n")
  
  
  output.October.OFF<-subset(VPRM.freq.OFF.October, Year==i)
  count.freq <- count (output.October.OFF)
  assign(paste("out.",i,sep=""),count.freq)
  
  
}

OFF.VPRM.freq.df.October <- rbind(out.2000,out.2001,out.2002,out.2003,out.2004,out.2005,out.2006,out.2007,out.2008,out.2009,out.2010,out.2011,out.2012,out.2013,out.2014,out.2015,out.2016,out.2017,out.2018,out.2019)
colnames(OFF.VPRM.freq.df.October) <- c("OFF")

Year.df <- data.frame(c(2000:2019))
colnames(Year.df) <- c("Year")

#####Combine Offshore and Onshore for October
final.VPRM.freq.df.October <-  data.frame(Year.df,OFF.VPRM.freq.df.October,ON.VPRM.freq.df.October)

##Create new data for Percentage dataframe
##Start by creating a duplicate copy of final.freq dataframe
VPRM.freq.percent.October <- final.VPRM.freq.df.October

#Calculate and add new column with the total frequency
Total.col <- final.VPRM.freq.df.October$OFF + final.VPRM.freq.df.October$ON
VPRM.freq.percent.October["Total"] <- Total.col

#Calculate and add new column with the percentages
OFF.per.col <- (VPRM.freq.percent.October$OFF/VPRM.freq.percent.October$Total)*100
ON.per.col <- (VPRM.freq.percent.October$ON/VPRM.freq.percent.October$Total)*100
VPRM.freq.percent.October["OFF.per"] <- OFF.per.col
VPRM.freq.percent.October["ON.per"] <- ON.per.col
VPRM.freq.percent.October

#Select OFF.per and ON.per columns only and Remove all other columns
select.variables <- c("Year", "OFF.per", "ON.per")
VPRM.freq.percent.October.final <- VPRM.freq.percent.October[select.variables]
VPRM.freq.percent.October.final
#Rename colnames for meaningful Legend
names(VPRM.freq.percent.October.final) <- c("Year","OFF","ON")
VPRM.freq.percent.October.final





##Create Separate dataframe for ON and OFF for each month
extract.cols.ON <- c("Year", "ON")
ON.VPRM.freq.percent.June.final <- VPRM.freq.percent.June.final[extract.cols.ON]
extract.cols.OFF <- c("Year", "OFF")
OFF.VPRM.freq.percent.June.final <- VPRM.freq.percent.June.final[extract.cols.OFF]

extract.cols.ON <- c("Year", "ON")
ON.VPRM.freq.percent.July.final <- VPRM.freq.percent.July.final[extract.cols.ON]
extract.cols.OFF <- c("Year", "OFF")
OFF.VPRM.freq.percent.July.final <- VPRM.freq.percent.July.final[extract.cols.OFF]

extract.cols.ON <- c("Year", "ON")
ON.VPRM.freq.percent.August.final <- VPRM.freq.percent.August.final[extract.cols.ON]
extract.cols.OFF <- c("Year", "OFF")
OFF.VPRM.freq.percent.August.final <- VPRM.freq.percent.August.final[extract.cols.OFF]

extract.cols.ON <- c("Year", "ON")
ON.VPRM.freq.percent.September.final <- VPRM.freq.percent.September.final[extract.cols.ON]
extract.cols.OFF <- c("Year", "OFF")
OFF.VPRM.freq.percent.September.final <- VPRM.freq.percent.September.final[extract.cols.OFF]

extract.cols.ON <- c("Year", "ON")
ON.VPRM.freq.percent.October.final <- VPRM.freq.percent.October.final[extract.cols.ON]
extract.cols.OFF <- c("Year", "OFF")
OFF.VPRM.freq.percent.October.final <- VPRM.freq.percent.October.final[extract.cols.OFF]





#########################################################ANOMALIES CALCULATIONS#############################################################################

######################################June###
#Calculate mean for the entire study period (one single value)
ON.VPRM.freq.climatology.mean.June <- mean(ON.VPRM.freq.percent.June.final$ON)
OFF.VPRM.freq.climatology.mean.June <- mean(OFF.VPRM.freq.percent.June.final$OFF)

##print to console
ON.VPRM.freq.climatology.mean.June
OFF.VPRM.freq.climatology.mean.June

##Create data for Anomaly timeseries
##Start by creating a duplicate copy of mean annual dataframe
ON.mean.annual.new.June <- ON.VPRM.freq.percent.June.final
OFF.mean.annual.new.June <- OFF.VPRM.freq.percent.June.final

#Add new column with the mean climatology
ON.mean.annual.new.June["Climatology"] <- ON.VPRM.freq.climatology.mean.June
OFF.mean.annual.new.June["Climatology"] <- OFF.VPRM.freq.climatology.mean.June

#Calculate anomaly and add new column
ON.anomaly.result.June <- ON.mean.annual.new.June$ON - ON.mean.annual.new.June$Climatology
ON.mean.annual.new.June["Anomaly"] <- ON.anomaly.result.June

OFF.anomaly.result.June <- OFF.mean.annual.new.June$OFF - OFF.mean.annual.new.June$Climatology
OFF.mean.annual.new.June["Anomaly"] <- OFF.anomaly.result.June



######################################July###
#Calculate mean for the entire study period (one single value)
ON.VPRM.freq.climatology.mean.July <- mean(ON.VPRM.freq.percent.July.final$ON)
OFF.VPRM.freq.climatology.mean.July <- mean(OFF.VPRM.freq.percent.July.final$OFF)

##print to console
ON.VPRM.freq.climatology.mean.July
OFF.VPRM.freq.climatology.mean.July

#Calculate Anomaly for the study period
ON.anomaly.period.July <- ON.VPRM.freq.climatology.mean.July - ON.VPRM.freq.climatology.mean.July
OFF.anomaly.period.July <- OFF.VPRM.freq.climatology.mean.July - OFF.VPRM.freq.climatology.mean.July

##Create data for Anomaly timeseries
##Start by creating a duplicate copy of mean annual dataframe
ON.mean.annual.new.July <- ON.VPRM.freq.percent.July.final
OFF.mean.annual.new.July <- OFF.VPRM.freq.percent.July.final

#Add new column with the mean climatology
ON.mean.annual.new.July["Climatology"] <- ON.VPRM.freq.climatology.mean.July
OFF.mean.annual.new.July["Climatology"] <- OFF.VPRM.freq.climatology.mean.July

#Calculate anomaly and add new column
ON.anomaly.result.July <- ON.mean.annual.new.July$ON - ON.mean.annual.new.July$Climatology
ON.mean.annual.new.July["Anomaly"] <- ON.anomaly.result.July

OFF.anomaly.result.July <- OFF.mean.annual.new.July$OFF - OFF.mean.annual.new.July$Climatology
OFF.mean.annual.new.July["Anomaly"] <- OFF.anomaly.result.July



######################################August###
#Calculate mean for the entire study period (one single value)
ON.VPRM.freq.climatology.mean.August <- mean(ON.VPRM.freq.percent.August.final$ON)
OFF.VPRM.freq.climatology.mean.August <- mean(OFF.VPRM.freq.percent.August.final$OFF)

##print to console
ON.VPRM.freq.climatology.mean.August
OFF.VPRM.freq.climatology.mean.August

#Calculate Anomaly for the study period
ON.anomaly.period.August <- ON.VPRM.freq.climatology.mean.August - ON.VPRM.freq.climatology.mean.August
OFF.anomaly.period.August <- OFF.VPRM.freq.climatology.mean.August - OFF.VPRM.freq.climatology.mean.August

##Create data for Anomaly timeseries
##Start by creating a duplicate copy of mean annual dataframe
ON.mean.annual.new.August <- ON.VPRM.freq.percent.August.final
OFF.mean.annual.new.August <- OFF.VPRM.freq.percent.August.final

#Add new column with the mean climatology
ON.mean.annual.new.August["Climatology"] <- ON.VPRM.freq.climatology.mean.August
OFF.mean.annual.new.August["Climatology"] <- OFF.VPRM.freq.climatology.mean.August

#Calculate anomaly and add new column
ON.anomaly.result.August <- ON.mean.annual.new.August$ON - ON.mean.annual.new.August$Climatology
ON.mean.annual.new.August["Anomaly"] <- ON.anomaly.result.August

OFF.anomaly.result.August <- OFF.mean.annual.new.August$OFF - OFF.mean.annual.new.August$Climatology
OFF.mean.annual.new.August["Anomaly"] <- OFF.anomaly.result.August



######################################September###
#Calculate mean for the entire study period (one single value)
ON.VPRM.freq.climatology.mean.September <- mean(ON.VPRM.freq.percent.September.final$ON)
OFF.VPRM.freq.climatology.mean.September <- mean(OFF.VPRM.freq.percent.September.final$OFF)

##print to console
ON.VPRM.freq.climatology.mean.September
OFF.VPRM.freq.climatology.mean.September

#Calculate Anomaly for the study period
ON.anomaly.period.September <- ON.VPRM.freq.climatology.mean.September - ON.VPRM.freq.climatology.mean.September
OFF.anomaly.period.September <- OFF.VPRM.freq.climatology.mean.September - OFF.VPRM.freq.climatology.mean.September

##Create data for Anomaly timeseries
##Start by creating a duplicate copy of mean annual dataframe
ON.mean.annual.new.September <- ON.VPRM.freq.percent.September.final
OFF.mean.annual.new.September <- OFF.VPRM.freq.percent.September.final

#Add new column with the mean climatology
ON.mean.annual.new.September["Climatology"] <- ON.VPRM.freq.climatology.mean.September
OFF.mean.annual.new.September["Climatology"] <- OFF.VPRM.freq.climatology.mean.September

#Calculate anomaly and add new column
ON.anomaly.result.September <- ON.mean.annual.new.September$ON - ON.mean.annual.new.September$Climatology
ON.mean.annual.new.September["Anomaly"] <- ON.anomaly.result.September

OFF.anomaly.result.September <- OFF.mean.annual.new.September$OFF - OFF.mean.annual.new.September$Climatology
OFF.mean.annual.new.September["Anomaly"] <- OFF.anomaly.result.September



######################################October###
#Calculate mean for the entire study period (one single value)
ON.VPRM.freq.climatology.mean.October <- mean(ON.VPRM.freq.percent.October.final$ON)
OFF.VPRM.freq.climatology.mean.October <- mean(OFF.VPRM.freq.percent.October.final$OFF)

##print to console
ON.VPRM.freq.climatology.mean.October
OFF.VPRM.freq.climatology.mean.October

#Calculate Anomaly for the study period
ON.anomaly.period.October <- ON.VPRM.freq.climatology.mean.October - ON.VPRM.freq.climatology.mean.October
OFF.anomaly.period.October <- OFF.VPRM.freq.climatology.mean.October - OFF.VPRM.freq.climatology.mean.October

##Create data for Anomaly timeseries
##Start by creating a duplicate copy of mean annual dataframe
ON.mean.annual.new.October <- ON.VPRM.freq.percent.October.final
OFF.mean.annual.new.October <- OFF.VPRM.freq.percent.October.final

#Add new column with the mean climatology
ON.mean.annual.new.October["Climatology"] <- ON.VPRM.freq.climatology.mean.October
OFF.mean.annual.new.October["Climatology"] <- OFF.VPRM.freq.climatology.mean.October

#Calculate anomaly and add new column
ON.anomaly.result.October <- ON.mean.annual.new.October$ON - ON.mean.annual.new.October$Climatology
ON.mean.annual.new.October["Anomaly"] <- ON.anomaly.result.October

OFF.anomaly.result.October <- OFF.mean.annual.new.October$OFF - OFF.mean.annual.new.October$Climatology
OFF.mean.annual.new.October["Anomaly"] <- OFF.anomaly.result.October






#####################Plot Percentage Bargraph - Mean Monthly and Season####

bar.frame.frequency.OFF <- c(June = OFF.VPRM.freq.climatology.mean.June, July = OFF.VPRM.freq.climatology.mean.July, August = OFF.VPRM.freq.climatology.mean.August,
                             September = OFF.VPRM.freq.climatology.mean.September, October = OFF.VPRM.freq.climatology.mean.October, Season = OFF.freq.climatology.mean)

bar.frame.frequency.ON <- c(June = ON.VPRM.freq.climatology.mean.June, July = ON.VPRM.freq.climatology.mean.July, August = ON.VPRM.freq.climatology.mean.August,
                            September = ON.VPRM.freq.climatology.mean.September, October = ON.VPRM.freq.climatology.mean.October, Season = ON.freq.climatology.mean)

bar.frame.frequency.COMBINED <- cbind.data.frame(OFF = bar.frame.frequency.OFF, ON = bar.frame.frequency.ON)
bar.frame.frequency.COMBINED ["Month"] <- c("June", "July", "August", "September", "October", "Season")




# NEW LINES ENDS - MESO: WIND FREQUENCY --------------------------------------------------- ##############################################################




# NEW LINES STARTS - MESO: FLUX CYCLES --------------------------------------------------- ##############################################################

####################Diurnal Analysis###
###Subset the relevant columns for aggregation
variables.fluxes.diurnal <- c("GPP_VPRM", "ER_VPRM", "NEE_VPRM", "Time", "Shore")
VPRM_Fluxes.Diurnal.subset <- VPRM_Fluxes.Diurnal[variables.fluxes.diurnal]
VPRM_Fluxes.Diurnal.Mean <- aggregate(.~ Time + Shore, data = VPRM_Fluxes.Diurnal.subset, mean, na.action = na.pass)

##Melt to prepare data form ggplot multiple variables
VPRM_Fluxes.Diurnal.Mean.Melted <- reshape2::melt(VPRM_Fluxes.Diurnal.Mean, id.var=c('Time', 'Shore'))
head(VPRM_Fluxes.Diurnal.Mean.Melted)
tail(VPRM_Fluxes.Diurnal.Mean.Melted)


####################Daily Averaged Analysis###
##Subset the relevant columns for aggregation and convert month names to month numbers
VPRM_Fluxes.Daily <- VPRM_Fluxes.Diurnal #Copy previously used dataframe from above to avoid potential conflicts
variables.fluxes.daily <- c("GPP_VPRM", "ER_VPRM", "NEE_VPRM", "Day", "Month", "Shore")
VPRM_Fluxes.Daily.subset <- VPRM_Fluxes.Daily[variables.fluxes.daily]
VPRM_Fluxes.Daily.subset$Month = match(VPRM_Fluxes.Daily.subset$Month,(month.name)) #Convert month names to month numbers before agg for proper ordering
VPRM_Fluxes.Daily.Mean <- aggregate(.~ Day + Month + Shore, data = VPRM_Fluxes.Daily.subset, mean, na.action = na.pass)

library(tidyverse)
library(lubridate)
VPRM_Fluxes.Daily.Mean ["Char_Date"] <- paste( month.abb[VPRM_Fluxes.Daily.Mean$Month], VPRM_Fluxes.Daily.Mean$Day, sep="-" )
VPRM_Fluxes.Daily.Mean ["Date"] <- as.Date(VPRM_Fluxes.Daily.Mean$Char_Date, format="%b-%d")
#VPRM_Fluxes.Daily.Mean ["Date"] <- as.data.frame(paste(VPRM_Fluxes.Daily.Mean$Month, VPRM_Fluxes.Daily.Mean$Day, sep="-") %>% ymd() %>% as.Date()) #Combine dates column into single date
VPRM_Fluxes.Daily.Mean$Char_Date <- NULL #Drop or remove the character Date column 
VPRM_Fluxes.Daily.Mean$Day <- NULL #Drop or remove the Day column 
VPRM_Fluxes.Daily.Mean$Month <- NULL #Drop or remove the Month column
#VPRM_Fluxes.Daily.Mean$Year <- NULL #Drop or remove the Year column

##Melt to prepare data form ggplot multiple variables
VPRM_Fluxes.Daily.Mean.Melted <- reshape2::melt(VPRM_Fluxes.Daily.Mean, id.var=c('Date', 'Shore'))
head(VPRM_Fluxes.Daily.Mean.Melted)
tail(VPRM_Fluxes.Daily.Mean.Melted)


####################Monthly Averaged Analysis###
##Subset the relevant columns for aggregation and convert month names to month numbers
VPRM_Fluxes.Monthly <- VPRM_Fluxes.Diurnal #Copy previously used dataframe from above to avoid potential conflicts
variables.fluxes.monthly <- c("GPP_VPRM", "ER_VPRM", "NEE_VPRM", "Month", "Shore")
VPRM_Fluxes.Monthly.subset <- VPRM_Fluxes.Monthly[variables.fluxes.monthly]
VPRM_Fluxes.Monthly.subset$Month = match(VPRM_Fluxes.Monthly.subset$Month,(month.name)) #Convert month names to month numbers before agg for proper ordering
VPRM_Fluxes.Monthly.Mean <- aggregate(.~ Month + Shore, data = VPRM_Fluxes.Monthly.subset, mean, na.action = na.pass)

##Melt to prepare data form ggplot multiple variables
VPRM_Fluxes.Monthly.Mean.Melted <- reshape2::melt(VPRM_Fluxes.Monthly.Mean, id.var=c('Month', 'Shore'))
head(VPRM_Fluxes.Monthly.Mean.Melted)
tail(VPRM_Fluxes.Monthly.Mean.Melted)


####################Annual Averaged Analysis###
##Subset the relevant columns for aggregation
VPRM_Fluxes.Annual <- VPRM_Fluxes.Diurnal #Copy previously used dataframe from above to avoid potential conflicts
variables.fluxes.annual <- c("GPP_VPRM", "ER_VPRM", "NEE_VPRM", "Year", "Shore")
VPRM_Fluxes.Annual.subset <- VPRM_Fluxes.Annual[variables.fluxes.annual]
VPRM_Fluxes.Annual.Mean <- aggregate(.~ Year + Shore, data = VPRM_Fluxes.Annual.subset, mean, na.action = na.pass)
#VPRM_Fluxes.Annual.Mean.Line <- VPRM_Fluxes.Annual.Mean #Duplicate dataframe to avoid conflict and use for other graphs
VPRM_Fluxes.Annual.Mean ["Year_Char"] <- as.character(VPRM_Fluxes.Annual.Mean$Year) #Year column is in class ordered or factor. Need to change to Numeric
VPRM_Fluxes.Annual.Mean$Year <- NULL #Drop or remove the original ordered/factor Year column
VPRM_Fluxes.Annual.Mean ["Year"] <- as.numeric(VPRM_Fluxes.Annual.Mean$Year_Char)
VPRM_Fluxes.Annual.Mean$Year_Char <- NULL #Drop or remove the Character Year column 

##Melt to prepare data form ggplot multiple variables
VPRM_Fluxes.Annual.Mean.Melted <- reshape2::melt(VPRM_Fluxes.Annual.Mean, id.var=c('Year', 'Shore'))
head(VPRM_Fluxes.Annual.Mean.Melted)
tail(VPRM_Fluxes.Annual.Mean.Melted)

#Write VPRM Annual Mean Fluxes to csv file
write.csv(VPRM_Fluxes.Annual.Mean, file = "Meso_Fluxes_Annual_Mean.csv", row.names = F)

#Move the excel .csv files to the stats folder
library(filesstrings)
#file.move("Meso_Fluxes_Annual_Mean.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/CF_MESO/Stats")



##########Convert Annual Averaged to g C m^-2#########################
##Duplicate data frame from annual averaged analysis above
VPRM_Fluxes.Annual.Mean.Melted.gC <- VPRM_Fluxes.Annual.Mean

#VPRM_Fluxes.Annual.Mean.Melted.gC [["value_gC"]] <- VPRM_Fluxes.Annual.Mean.Melted.gC$value * (86400 * 153 * 12.011) / 1000000 #Convert from umols of CO2 m^-2 s^-1 to g C m^-2. Note that 153 is the number of days in the study period per year.

VPRM_Fluxes.Annual.Mean.Melted.gC[1:20,2:4] <- VPRM_Fluxes.Annual.Mean.Melted.gC[1:20,2:4] * (86400 * 153 * (bar.frame.frequency.COMBINED[["Season","OFF"]]/100) * 12.011) / 1000000 #Conversion that includes weighting by Wind freq percentage
VPRM_Fluxes.Annual.Mean.Melted.gC[21:40,2:4] <- VPRM_Fluxes.Annual.Mean.Melted.gC[21:40,2:4] * (86400 * 153 * (bar.frame.frequency.COMBINED[["Season","ON"]]/100) * 12.011) / 1000000 #Conversion that includes weighting by Wind freq percentage

##Melt to prepare data form ggplot multiple variables
VPRM_Fluxes.Annual.Mean.Melted.gC <- reshape2::melt(VPRM_Fluxes.Annual.Mean.Melted.gC, id.var=c('Year', 'Shore'))
names(VPRM_Fluxes.Annual.Mean.Melted.gC)[names(VPRM_Fluxes.Annual.Mean.Melted.gC) == 'value'] <- 'value_gC' #Rename last column
head(VPRM_Fluxes.Annual.Mean.Melted.gC)
tail(VPRM_Fluxes.Annual.Mean.Melted.gC)


#Write VPRM Annual Mean Fluxes in g C to csv file
VPRM_Fluxes.Annual.Mean.gC <- dcast(VPRM_Fluxes.Annual.Mean.Melted.gC, Year + Shore ~ variable)

write.csv(VPRM_Fluxes.Annual.Mean.gC, file = "Meso_Fluxes_Annual_gC.csv", row.names = F)

#Move the excel .csv files to the stats folder
library(filesstrings)
#file.move("Meso_Fluxes_Annual_gC.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/CF_MESO/Stats")





##############################################################################STUDENT T-TESTS ANALYSIS###############################################################################################################################
####################Create output dataframe - Diurnal
library("dplyr")
library("ggpubr")
OFF.VPRM_Fluxes.Diurnal.Mean <- subset(VPRM_Fluxes.Diurnal.Mean, Shore=="OFF")
ON.VPRM_Fluxes.Diurnal.Mean <- subset(VPRM_Fluxes.Diurnal.Mean, Shore=="ON")




####################Create output dataframe - Daily
OFF.VPRM_Fluxes.Daily.Mean <- subset(VPRM_Fluxes.Daily.Mean, Shore=="OFF")
ON.VPRM_Fluxes.Daily.Mean <- subset(VPRM_Fluxes.Daily.Mean, Shore=="ON")


####################Create output dataframe - Annual
OFF.VPRM_Fluxes.Annual.Mean <- subset(VPRM_Fluxes.Annual.Mean, Shore=="OFF")
ON.VPRM_Fluxes.Annual.Mean <- subset(VPRM_Fluxes.Annual.Mean, Shore=="ON")


####################Create output dataframe - Annual_gc
OFF.VPRM_Fluxes.Annual.Mean.gC <- subset(VPRM_Fluxes.Annual.Mean.gC, Shore=="OFF")
ON.VPRM_Fluxes.Annual.Mean.gC <- subset(VPRM_Fluxes.Annual.Mean.gC, Shore=="ON")



# NEW LINES ENDS - MESO: FLUX CYCLES --------------------------------------------------- ##############################################################






# NEW LINES STARTS - MESO: ANOMALIES ---------------------------------------- #######################################################################





##################################################################Anomalies - Annual Fluxes################################################################
##Subset the relevant columns for aggregation
VPRM_Fluxes.Annual <- VPRM_Fluxes.Diurnal #Copy previously used dataframe from above to avoid potential conflicts
variables.fluxes.annual <- c("GPP_VPRM", "ER_VPRM", "NEE_VPRM", "Year", "Shore")
VPRM_Fluxes.Annual.subset <- VPRM_Fluxes.Annual[variables.fluxes.annual]
VPRM_Fluxes.Annual.Mean <- aggregate(.~ Year + Shore, data = VPRM_Fluxes.Annual.subset, mean, na.action = na.pass)

VPRM_Fluxes.Annual.Mean.gC <- as.data.frame(VPRM_Fluxes.Annual.Mean[c("Year", "Shore")])
VPRM_Fluxes.Annual.Mean.gC [["NEE_VPRM"]] <- VPRM_Fluxes.Annual.Mean$NEE_VPRM
VPRM_Fluxes.Annual.Mean.gC [["GPP_VPRM"]] <- VPRM_Fluxes.Annual.Mean$GPP_VPRM
VPRM_Fluxes.Annual.Mean.gC [["ER_VPRM"]] <- VPRM_Fluxes.Annual.Mean$ER_VPRM

#VPRM_Fluxes.Annual.Mean.gC <- as.data.frame(VPRM_Fluxes.Annual.Mean[c("Year", "Shore")])
#VPRM_Fluxes.Annual.Mean.gC [["NEE_VPRM"]] <- VPRM_Fluxes.Annual.Mean$NEE_VPRM * (86400 * 153 * 12.011) / 1000000 #Convert from umols of CO2 m^-2 s^-1 to g C m^-2. Note that 153 is the number of days in the study period per year.
#VPRM_Fluxes.Annual.Mean.gC [["GPP_VPRM"]] <- VPRM_Fluxes.Annual.Mean$GPP_VPRM * (86400 * 153 * 12.011) / 1000000 #Convert from umols of CO2 m^-2 s^-1 to g C m^-2. Note that 153 is the number of days in the study period per year.
#VPRM_Fluxes.Annual.Mean.gC [["ER_VPRM"]] <- VPRM_Fluxes.Annual.Mean$ER_VPRM * (86400 * 153 * 12.011) / 1000000 #Convert from umols of CO2 m^-2 s^-1 to g C m^-2. Note that 153 is the number of days in the study period per year.

VPRM_Fluxes.Annual.Mean.gC[1:20,3:5] <- VPRM_Fluxes.Annual.Mean.gC[1:20,3:5] * (86400 * 153 * (bar.frame.frequency.COMBINED[["Season","OFF"]]/100) * 12.011) / 1000000 #Conversion that includes weighting by Wind freq percentage
VPRM_Fluxes.Annual.Mean.gC[21:40,3:5] <- VPRM_Fluxes.Annual.Mean.gC[21:40,3:5] * (86400 * 153 * (bar.frame.frequency.COMBINED[["Season","ON"]]/100) * 12.011) / 1000000 #Conversion that includes weighting by Wind freq percentage

#Calculate mean (climatology) for the study period (20 years) - OFF
OFF.VPRM_Fluxes.Annual.Mean.gC <- subset(VPRM_Fluxes.Annual.Mean.gC, Shore=="OFF")
#OFF.VPRM_Fluxes.climatology <- subset(OFF.VPRM_Fluxes.Annual.Mean.gC, Year %in% 2000:2019)
OFF.VPRM_Fluxes.climatology.mean.gC.NEE <- mean(OFF.VPRM_Fluxes.Annual.Mean.gC$NEE_VPRM)
OFF.VPRM_Fluxes.climatology.mean.gC.GPP <- mean(OFF.VPRM_Fluxes.Annual.Mean.gC$GPP_VPRM)
OFF.VPRM_Fluxes.climatology.mean.gC.ER <- mean(OFF.VPRM_Fluxes.Annual.Mean.gC$ER_VPRM)

OFF.VPRM_Fluxes.climatology.mean.gC.NEE
OFF.VPRM_Fluxes.climatology.mean.gC.GPP
OFF.VPRM_Fluxes.climatology.mean.gC.ER

##Create data for Anomaly timeseries. Start by creating a duplicate copy of mean annual dataframe
OFF.VPRM_Fluxes.Annual.Mean.gC.new <- OFF.VPRM_Fluxes.Annual.Mean.gC

#Add new columns with the mean climatology
OFF.VPRM_Fluxes.Annual.Mean.gC.new["Climatology_NEE"] <- OFF.VPRM_Fluxes.climatology.mean.gC.NEE
OFF.VPRM_Fluxes.Annual.Mean.gC.new["Climatology_GPP"] <- OFF.VPRM_Fluxes.climatology.mean.gC.GPP
OFF.VPRM_Fluxes.Annual.Mean.gC.new["Climatology_ER"] <- OFF.VPRM_Fluxes.climatology.mean.gC.ER

#Calculate anomalies and add new columns
OFF.VPRM_Fluxes.Annual.anomaly.gC.NEE <- OFF.VPRM_Fluxes.Annual.Mean.gC.new$NEE_VPRM - OFF.VPRM_Fluxes.Annual.Mean.gC.new$Climatology_NEE
OFF.VPRM_Fluxes.Annual.anomaly.gC.GPP <- OFF.VPRM_Fluxes.Annual.Mean.gC.new$GPP_VPRM - OFF.VPRM_Fluxes.Annual.Mean.gC.new$Climatology_GPP
OFF.VPRM_Fluxes.Annual.anomaly.gC.ER <- OFF.VPRM_Fluxes.Annual.Mean.gC.new$ER_VPRM - OFF.VPRM_Fluxes.Annual.Mean.gC.new$Climatology_ER

OFF.VPRM_Fluxes.Annual.Mean.gC.new["Anomaly_NEE"] <- OFF.VPRM_Fluxes.Annual.anomaly.gC.NEE
OFF.VPRM_Fluxes.Annual.Mean.gC.new["Anomaly_GPP"] <- OFF.VPRM_Fluxes.Annual.anomaly.gC.GPP
OFF.VPRM_Fluxes.Annual.Mean.gC.new["Anomaly_ER"] <- OFF.VPRM_Fluxes.Annual.anomaly.gC.ER

OFF.VPRM_Fluxes.Annual.Mean.gC.new ["Year_Char"] <- as.character(OFF.VPRM_Fluxes.Annual.Mean.gC.new$Year) #Year column is in class ordered or factor. Need to change to Numeric
OFF.VPRM_Fluxes.Annual.Mean.gC.new$Year <- NULL #Drop or remove the original ordered/factor Year column
OFF.VPRM_Fluxes.Annual.Mean.gC.new ["Year"] <- as.numeric(OFF.VPRM_Fluxes.Annual.Mean.gC.new$Year_Char)
OFF.VPRM_Fluxes.Annual.Mean.gC.new$Year_Char <- NULL #Drop or remove the Character Year column 

##Subset the relevant columns to plot anomalies
OFF.variables.fluxes.annual.anomalies.gC <- c("Year", "Shore", "Anomaly_GPP", "Anomaly_ER", "Anomaly_NEE")
#OFF.variables.fluxes.annual.anomalies.gC <- c("Year", "Shore", "Anomaly_NEE", "Anomaly_GPP", "Anomaly_ER")
OFF.VPRM_Fluxes.Annual.anomaly.gC.Final <- OFF.VPRM_Fluxes.Annual.Mean.gC.new[OFF.variables.fluxes.annual.anomalies.gC]

#############Statistical tests for OFF.Trends################
# Install and load the packages
#install.packages(c("car", "Kendall", "zyp", "lmtest"))
library(car)
library(Kendall)
library(zyp)
library(lmtest)

#Test for Autocorrelation - durbinWatsonTest
y1.NEE <- OFF.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_NEE
x1 <- OFF.VPRM_Fluxes.Annual.anomaly.gC.Final$Year
dwtest(y1.NEE ~ x1)

y1.GPP <- OFF.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_GPP
x1 <- OFF.VPRM_Fluxes.Annual.anomaly.gC.Final$Year
dwtest(y1.GPP ~ x1)

y1.ER <- OFF.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_ER
x1 <- OFF.VPRM_Fluxes.Annual.anomaly.gC.Final$Year
dwtest(y1.ER ~ x1)

#Run the Mann-Kendall test
y.NEE <- OFF.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_NEE
OFF.Trend.Annual.gC.NEE <- MannKendall(y.NEE)
OFF.Trend.Annual.gC.NEE

y.GPP <- OFF.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_GPP
OFF.Trend.Annual.gC.GPP <- MannKendall(y.GPP)
OFF.Trend.Annual.gC.GPP

y.ER <- OFF.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_ER
OFF.Trend.Annual.gC.ER <- MannKendall(y.ER)
OFF.Trend.Annual.gC.ER

#TheilSen Slope Estimation
y.NEE <- OFF.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_NEE
x <- OFF.VPRM_Fluxes.Annual.anomaly.gC.Final$Year
OFF.Slope.Annual.gC.NEE <- zyp.sen(y.NEE~x)
OFF.Slope.Annual.gC.NEE

y.GPP <- OFF.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_GPP
x <- OFF.VPRM_Fluxes.Annual.anomaly.gC.Final$Year
OFF.Slope.Annual.gC.GPP <- zyp.sen(y.GPP~x)
OFF.Slope.Annual.gC.GPP

y.ER <- OFF.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_ER
x <- OFF.VPRM_Fluxes.Annual.anomaly.gC.Final$Year
OFF.Slope.Annual.gC.ER <- zyp.sen(y.ER~x)
OFF.Slope.Annual.gC.ER

#Subset and create objects for trend outputs to add to plots
OFF.Trend.Annual.gC.NEE
OFF.Slope.Annual.gC.NEE
OFF.Trend.Annual.gC.NEE.tau <- format(round(OFF.Trend.Annual.gC.NEE$tau[1], 2), nsmall = 2) #Round to two decimals
OFF.Trend.Annual.gC.NEE.pvalue <- format(round(OFF.Trend.Annual.gC.NEE$sl[1], 2), nsmall = 2) #Round to two decimals
OFF.Slope.Annual.gC.NEE.20 <- format(round(OFF.Slope.Annual.gC.NEE$coefficients[["x"]]*20, 2), nsmall = 2) #Round to two decimals

OFF.Trend.Annual.gC.GPP
OFF.Slope.Annual.gC.GPP
OFF.Trend.Annual.gC.GPP.tau <- format(round(OFF.Trend.Annual.gC.GPP$tau[1], 2), nsmall = 2) #Round to two decimals
OFF.Trend.Annual.gC.GPP.pvalue <- format(round(OFF.Trend.Annual.gC.GPP$sl[1], 2), nsmall = 2) #Round to two decimals
OFF.Slope.Annual.gC.GPP.20 <- format(round(OFF.Slope.Annual.gC.GPP$coefficients[["x"]]*20, 2), nsmall = 2) #Round to two decimals

OFF.Trend.Annual.gC.ER
OFF.Slope.Annual.gC.ER
OFF.Trend.Annual.gC.ER.tau <- format(round(OFF.Trend.Annual.gC.ER$tau[1], 2), nsmall = 2) #Round to two decimals
OFF.Trend.Annual.gC.ER.pvalue <- format(round(OFF.Trend.Annual.gC.ER$sl[1], 2), nsmall = 2) #Round to two decimals
OFF.Slope.Annual.gC.ER.20 <- format(round(OFF.Slope.Annual.gC.ER$coefficients[["x"]]*20, 2), nsmall = 2) #Round to two decimals


##Melt to prepare data form ggplot multiple variables
OFF.VPRM_Fluxes.Annual.anomaly.gC.Final.Melted <- reshape2::melt(OFF.VPRM_Fluxes.Annual.anomaly.gC.Final, id.var=c('Year', 'Shore'))
head(OFF.VPRM_Fluxes.Annual.anomaly.gC.Final.Melted)
tail(OFF.VPRM_Fluxes.Annual.anomaly.gC.Final.Melted)






#Calculate mean (climatology) for the study period (20 years) - ON
ON.VPRM_Fluxes.Annual.Mean.gC <- subset(VPRM_Fluxes.Annual.Mean.gC, Shore=="ON")
#ON.VPRM_Fluxes.climatology <- subset(ON.VPRM_Fluxes.Annual.Mean.gC, Year %in% 2000:2019)
ON.VPRM_Fluxes.climatology.mean.gC.NEE <- mean(ON.VPRM_Fluxes.Annual.Mean.gC$NEE_VPRM)
ON.VPRM_Fluxes.climatology.mean.gC.GPP <- mean(ON.VPRM_Fluxes.Annual.Mean.gC$GPP_VPRM)
ON.VPRM_Fluxes.climatology.mean.gC.ER <- mean(ON.VPRM_Fluxes.Annual.Mean.gC$ER_VPRM)

ON.VPRM_Fluxes.climatology.mean.gC.NEE
ON.VPRM_Fluxes.climatology.mean.gC.GPP
ON.VPRM_Fluxes.climatology.mean.gC.ER

##Create data for Anomaly timeseries. Start by creating a duplicate copy of mean annual dataframe
ON.VPRM_Fluxes.Annual.Mean.gC.new <- ON.VPRM_Fluxes.Annual.Mean.gC

#Add new columns with the mean climatology
ON.VPRM_Fluxes.Annual.Mean.gC.new["Climatology_NEE"] <- ON.VPRM_Fluxes.climatology.mean.gC.NEE
ON.VPRM_Fluxes.Annual.Mean.gC.new["Climatology_GPP"] <- ON.VPRM_Fluxes.climatology.mean.gC.GPP
ON.VPRM_Fluxes.Annual.Mean.gC.new["Climatology_ER"] <- ON.VPRM_Fluxes.climatology.mean.gC.ER

#Calculate anomalies and add new columns
ON.VPRM_Fluxes.Annual.anomaly.gC.NEE <- ON.VPRM_Fluxes.Annual.Mean.gC.new$NEE_VPRM - ON.VPRM_Fluxes.Annual.Mean.gC.new$Climatology_NEE
ON.VPRM_Fluxes.Annual.anomaly.gC.GPP <- ON.VPRM_Fluxes.Annual.Mean.gC.new$GPP_VPRM - ON.VPRM_Fluxes.Annual.Mean.gC.new$Climatology_GPP
ON.VPRM_Fluxes.Annual.anomaly.gC.ER <- ON.VPRM_Fluxes.Annual.Mean.gC.new$ER_VPRM - ON.VPRM_Fluxes.Annual.Mean.gC.new$Climatology_ER

ON.VPRM_Fluxes.Annual.Mean.gC.new["Anomaly_NEE"] <- ON.VPRM_Fluxes.Annual.anomaly.gC.NEE
ON.VPRM_Fluxes.Annual.Mean.gC.new["Anomaly_GPP"] <- ON.VPRM_Fluxes.Annual.anomaly.gC.GPP
ON.VPRM_Fluxes.Annual.Mean.gC.new["Anomaly_ER"] <- ON.VPRM_Fluxes.Annual.anomaly.gC.ER

ON.VPRM_Fluxes.Annual.Mean.gC.new ["Year_Char"] <- as.character(ON.VPRM_Fluxes.Annual.Mean.gC.new$Year) #Year column is in class ordered or factor. Need to change to Numeric
ON.VPRM_Fluxes.Annual.Mean.gC.new$Year <- NULL #Drop or remove the original ordered/factor Year column
ON.VPRM_Fluxes.Annual.Mean.gC.new ["Year"] <- as.numeric(ON.VPRM_Fluxes.Annual.Mean.gC.new$Year_Char)
ON.VPRM_Fluxes.Annual.Mean.gC.new$Year_Char <- NULL #Drop or remove the Character Year column 

##Subset the relevant columns to plot anomalies
ON.variables.fluxes.annual.anomalies.gC <- c("Year", "Shore", "Anomaly_GPP", "Anomaly_ER", "Anomaly_NEE")
#ON.variables.fluxes.annual.anomalies.gC <- c("Year", "Shore", "Anomaly_NEE", "Anomaly_GPP", "Anomaly_ER")
ON.VPRM_Fluxes.Annual.anomaly.gC.Final <- ON.VPRM_Fluxes.Annual.Mean.gC.new[ON.variables.fluxes.annual.anomalies.gC]

#############Statistical tests for ON.Trends################
# Install and load the packages
#install.packages(c("car", "Kendall", "zyp", "lmtest"))
library(car)
library(Kendall)
library(zyp)
library(lmtest)

#Test for Autocorrelation - durbinWatsonTest
y1.NEE <- ON.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_NEE
x1 <- ON.VPRM_Fluxes.Annual.anomaly.gC.Final$Year
dwtest(y1.NEE ~ x1)

y1.GPP <- ON.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_GPP
x1 <- ON.VPRM_Fluxes.Annual.anomaly.gC.Final$Year
dwtest(y1.GPP ~ x1)

y1.ER <- ON.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_ER
x1 <- ON.VPRM_Fluxes.Annual.anomaly.gC.Final$Year
dwtest(y1.ER ~ x1)

#Run the Mann-Kendall test
y.NEE <- ON.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_NEE
ON.Trend.Annual.gC.NEE <- MannKendall(y.NEE)
ON.Trend.Annual.gC.NEE

y.GPP <- ON.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_GPP
ON.Trend.Annual.gC.GPP <- MannKendall(y.GPP)
ON.Trend.Annual.gC.GPP

y.ER <- ON.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_ER
ON.Trend.Annual.gC.ER <- MannKendall(y.ER)
ON.Trend.Annual.gC.ER

#TheilSen Slope Estimation
y.NEE <- ON.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_NEE
x <- ON.VPRM_Fluxes.Annual.anomaly.gC.Final$Year
ON.Slope.Annual.gC.NEE <- zyp.sen(y.NEE~x)
ON.Slope.Annual.gC.NEE

y.GPP <- ON.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_GPP
x <- ON.VPRM_Fluxes.Annual.anomaly.gC.Final$Year
ON.Slope.Annual.gC.GPP <- zyp.sen(y.GPP~x)
ON.Slope.Annual.gC.GPP

y.ER <- ON.VPRM_Fluxes.Annual.anomaly.gC.Final$Anomaly_ER
x <- ON.VPRM_Fluxes.Annual.anomaly.gC.Final$Year
ON.Slope.Annual.gC.ER <- zyp.sen(y.ER~x)
ON.Slope.Annual.gC.ER

#Subset and create objects for trend outputs to add to plots
ON.Trend.Annual.gC.NEE
ON.Slope.Annual.gC.NEE
ON.Trend.Annual.gC.NEE.tau <- format(round(ON.Trend.Annual.gC.NEE$tau[1], 2), nsmall = 2) #Round to two decimals
ON.Trend.Annual.gC.NEE.pvalue <- format(round(ON.Trend.Annual.gC.NEE$sl[1], 2), nsmall = 2) #Round to two decimals
ON.Slope.Annual.gC.NEE.20 <- format(round(ON.Slope.Annual.gC.NEE$coefficients[["x"]]*20, 2), nsmall = 2) #Round to two decimals

ON.Trend.Annual.gC.GPP
ON.Slope.Annual.gC.GPP
ON.Trend.Annual.gC.GPP.tau <- format(round(ON.Trend.Annual.gC.GPP$tau[1], 2), nsmall = 2) #Round to two decimals
ON.Trend.Annual.gC.GPP.pvalue <- format(round(ON.Trend.Annual.gC.GPP$sl[1], 2), nsmall = 2) #Round to two decimals
ON.Slope.Annual.gC.GPP.20 <- format(round(ON.Slope.Annual.gC.GPP$coefficients[["x"]]*20, 2), nsmall = 2) #Round to two decimals

ON.Trend.Annual.gC.ER
ON.Slope.Annual.gC.ER
ON.Trend.Annual.gC.ER.tau <- format(round(ON.Trend.Annual.gC.ER$tau[1], 2), nsmall = 2) #Round to two decimals
ON.Trend.Annual.gC.ER.pvalue <- format(round(ON.Trend.Annual.gC.ER$sl[1], 2), nsmall = 2) #Round to two decimals
ON.Slope.Annual.gC.ER.20 <- format(round(ON.Slope.Annual.gC.ER$coefficients[["x"]]*20, 2), nsmall = 2) #Round to two decimals

##Melt to prepare data form ggplot multiple variables
ON.VPRM_Fluxes.Annual.anomaly.gC.Final.Melted <- reshape2::melt(ON.VPRM_Fluxes.Annual.anomaly.gC.Final, id.var=c('Year', 'Shore'))
head(ON.VPRM_Fluxes.Annual.anomaly.gC.Final.Melted)
tail(ON.VPRM_Fluxes.Annual.anomaly.gC.Final.Melted)

MESO_Fluxes.Annual.anomaly.gC.COMBINED <- rbind(OFF.VPRM_Fluxes.Annual.anomaly.gC.Final, ON.VPRM_Fluxes.Annual.anomaly.gC.Final)

##Melt to prepare data form ggplot multiple variables
MESO_Fluxes.Annual.anomaly.gC.COMBINED.Melted <- reshape2::melt(MESO_Fluxes.Annual.anomaly.gC.COMBINED, id.var=c('Year', 'Shore'))
head(MESO_Fluxes.Annual.anomaly.gC.COMBINED.Melted)
tail(MESO_Fluxes.Annual.anomaly.gC.COMBINED.Melted)



###########Create plots of the annual averaged flux for the growing season
library(ggplot2)
library(ggpmisc)
detach("package:ggpubr", unload=TRUE) #Detach 'ggpubr' before loading 'egg' for the first time, to avoid error in ggarrange
#install.packages('egg', dependencies = TRUE)
library(egg) #For function tag_facet to annotate individual facets #Make sure to load 'ggpubr' only after 'egg'. Otherwise, error in ggarrange
library(ggpubr) #IMPORTANT: Make sure to load 'ggpubr' only after 'egg'. Otherwise, error in ggarrange

###########Graph Type 2 - Bar Graph (For Annual Only)
#library(ggplot2)
#library(ggpmisc)
#detach("package:ggpubr", unload=TRUE) #Detach 'ggpubr' before loading 'egg' for the first time, to avoid error in ggarrange
#install.packages('egg', dependencies = TRUE)
#library(egg) #For function tag_facet to annotate individual facets #Make sure to load 'ggpubr' only after 'egg'. Otherwise, error in ggarrange
#library(ggpubr) #IMPORTANT: Make sure to load 'ggpubr' only after 'egg'. Otherwise, error in ggarrange
plot <- ggplot(MESO_Fluxes.Annual.anomaly.gC.COMBINED.Melted, aes(x=Year, y=value, group=Shore)) + geom_bar(stat = "identity", aes(fill = Shore), position = "dodge", width = 0.7) + geom_smooth(method = "lm", aes(colour = Shore), se = FALSE, show.legend = FALSE)
plot <- plot + scale_fill_manual(values = c("#FF3333", "#0099FF"))
plot <- plot + scale_color_manual(values = c("#FF3333", "#0099FF"))
plot <- plot + scale_linetype_manual(values = c("dashed", "dotted", "solid"))
plot <- plot + facet_grid(variable ~ .) + theme(strip.text.y = element_text(size=10, face="bold"))

label.OFF.GPP = (paste("GPP: Tau = ", OFF.Trend.Annual.gC.GPP.tau, ", p = ", OFF.Trend.Annual.gC.GPP.pvalue, ", TSA = ", OFF.Slope.Annual.gC.GPP.20, sep=""))
label.ON.GPP = (paste("GPP: Tau = ", ON.Trend.Annual.gC.GPP.tau, ", p = ", ON.Trend.Annual.gC.GPP.pvalue, ", TSA = ", ON.Slope.Annual.gC.GPP.20, sep=""))
label.OFF.ER = (paste("ER: Tau = ", OFF.Trend.Annual.gC.ER.tau, ", p = ", OFF.Trend.Annual.gC.ER.pvalue, ", TSA = ", OFF.Slope.Annual.gC.ER.20, sep=""))
label.ON.ER = (paste("ER: Tau = ", ON.Trend.Annual.gC.ER.tau, ", p = ", ON.Trend.Annual.gC.ER.pvalue, ", TSA = ", ON.Slope.Annual.gC.ER.20, sep=""))
label.OFF.NEE = (paste("NEE: Tau = ", OFF.Trend.Annual.gC.NEE.tau, ", p = ", OFF.Trend.Annual.gC.NEE.pvalue, ", TSA = ", OFF.Slope.Annual.gC.NEE.20, sep=""))
label.ON.NEE = (paste("NEE: Tau = ", ON.Trend.Annual.gC.NEE.tau, ", p = ", ON.Trend.Annual.gC.NEE.pvalue, ", TSA = ", ON.Slope.Annual.gC.NEE.20, sep=""))
my_tag.OFF <- c(label.OFF.GPP, label.OFF.ER, label.OFF.NEE)
my_tag.ON <- c(label.ON.GPP, label.ON.ER, label.ON.NEE) 

plot <- tag_facet(plot, x = -Inf, y = Inf, vjust = 1.7, hjust = -0.05, open = "", close = "", size = 3, tag_pool = c(my_tag.OFF), color = "#FF3333") #For OFF
plot <- tag_facet(plot, x = -Inf, y = Inf, vjust = 3.2, hjust = -0.05, open = "", close = "", size = 3, tag_pool = c(my_tag.ON), color = "#0099FF") #For ON
plot <- plot + theme(strip.text = element_text(size=10, face="bold")) + theme(strip.background =element_rect(fill="lightgrey"))
#plot <- plot + facet_grid(. ~ variable) + theme(strip.text.x = element_text(size=10, face="bold"))
#plot <- plot + scale_color_manual(values = c("#00BA38", "#619CFF")) #New line for ONLY-SIF plot to set appropriate colors##################
plot <- plot + xlab("Year") + ylab(expression(Carbon~Fluxes~(g~C~m^{-2}))) + ggtitle(paste("CF Annual Flux Anomalies for OFF & ON" ,"(June-October, 2000-2019)"))
plot <- plot + scale_x_continuous(breaks = seq(from=2000,to=2020,by=2))
plot <- plot + theme(plot.title = element_text(hjust = 0.5)) #center title
plot <- plot + theme(axis.text=element_text(size=10, face="bold"))
plot <- plot + theme(axis.title=element_text(size=10))
plot <- plot + theme(plot.title=element_text(size=10,face="bold"))
#plot <- plot + theme(legend.key.size = unit(1.0, "cm"),legend.key.width = unit(1.0,"cm")) 
plot <- plot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank())
#plot <- plot + labs(shape = "") #Remove or specify Legend title
plot <- plot + theme(legend.title = element_blank()) #Remove Legend title
plot <- plot + theme(legend.text = element_text(size = 10, face = "bold"))
plot <- plot + theme(legend.position = "bottom")
plot <- plot + geom_hline(yintercept=0, color="black")
plot <- plot + scale_y_continuous(breaks=seq(-40,80,40)) #Sets the axis range and interval ("by")
plot <- plot + coord_cartesian(ylim = c(-40, 80), xlim = c(2000, 2020)) #Forces the limits to be set accordingly
plot.Anomalies.Annual_Bar.Meso <- plot


#save graph to output as high resoultion, Windows compatible format .EMF
fileout <- paste("Line-Meso-CF Annual Flux gC Anomalies","_LineCycle",".emf",sep="")
#x11(width=7, height=4) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
x11(width=7, height=5) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
dev.cur()
plot.Anomalies.Annual_Bar.Meso
savePlot(fileout, "emf", device = dev.cur())
dev.off()
graphics.off()

library(filesstrings)
#Write All (June-Annual) VPRM Anomalies to csv file and Move the excel .csv files to the stats folder
write.csv(MESO_Fluxes.Annual.anomaly.gC.COMBINED, file = "Meso-ALL_Annual_Fluxes_Annual_anomaly_gC_Combined.csv", row.names = F)
#file.move("Meso-ALL_June_Annual_Fluxes_Annual_anomaly_gC_Combined.csv", "~/Ecosystem_Modelling/VPRM_POLAR/MODELRUN_FINAL/CF_MESO/Stats")

# NEW LINES ENDS - MESO: ANOMALIES ---------------------------------------- #######################################################################



cat("\n"); 
VPRMstate = "Final Model Run and Data Analysis Completed for CF!"
cat("\n"); 
cat(VPRMstate);cat("\n"); cat("\n")

