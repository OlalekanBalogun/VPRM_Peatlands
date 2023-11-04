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
CF_tower_obs.SIF <- read.csv("CF-Tower-NARR-Obs-ready.csv", stringsAsFactors = FALSE)
CF_evi_modis.SIF <- read.csv("CF-SIF-ready.csv", stringsAsFactors = FALSE)
CF_refl_modis.SIF <- read.csv("CF-Reflectance-ready.csv", stringsAsFactors = FALSE)
CF_phen_modis.SIF <- read.csv("CF-Phenology-MAX-ready.csv", stringsAsFactors = FALSE)

#Convert "Character" dates to "Chron" dates and times for each VPRM driver data
library(chron)

CF_tower_obs.SIF.dates <- t(as.data.frame(strsplit(CF_tower_obs.SIF[["Date"]],' '))) #t is (matrix transpose)
row.names(CF_tower_obs.SIF.dates) = NULL
CF_tower_obs.SIF ["Date"] <- chron(dates=CF_tower_obs.SIF.dates[,1],times=CF_tower_obs.SIF.dates[,2], format = c(dates = "m/d/y", times = "h:m:s"))
class(CF_tower_obs.SIF[["Date"]])

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
phen_filled.SIF <- interp_phenology(CF_phen_rearranged.SIF, CF_tower_obs.SIF[['Date']])



#####Read csv file from original optimized run####
specify_path <- (dirname(rstudioapi::getActiveDocumentContext()$path))
filename_parameters <- paste(specify_path, "/CF_Optimized_Site_Parameter_Values",".csv",sep="")
SIF_optimized_parameters <- read.csv(filename_parameters)


######Assemble VPRM driver data frame

## determine the phenology phase for each tower observation date
#phen_filled.SIF <- interp_phenology(CF_phen, CF_tower_obs.SIF[['date']])
## place the tower observations and MODIS data into a VPRM_driver_data object.
CF_dd.SIF <- VPRM_driver_data(name_long="Churchill MB",
                               name_short = "CF",
                               lat=58.6658,
                               lon=-93.83,
                               PFT='WET',    ## Permanent Wetland
                               tower_date=CF_tower_obs.SIF[['Date']],
                               NEE_obs=CF_tower_obs.SIF[['FC']],
                               T=CF_tower_obs.SIF[['TA']],
                               PAR=CF_tower_obs.SIF[['PAR']],
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


######Calculate VPRM NEE
#Calculating VPRM fluxes requires values for the parameters lambda, PAR0, alpha, and beta.
#The dataset VPRM_parameters, part of the VPRMLandSfcModel package, includes several parameters sets.

## plot the VPRM NEE
attach(SIF_optimized_parameters)

CF_dd.SIF[['data']][['NEE_VPRM']] <- vprm_calc_NEE(
  CF_dd.SIF, lambda=lambda, PAR_0=PAR_0, alpha=alpha, beta=beta)

##################################################Diagnosis#############################################
NEE.Diff.SIF <- data.frame(Date = CF_dd.SIF[["data"]][["date"]], NEE_VPRM = CF_dd.SIF[["data"]][["NEE_VPRM"]], NEE_obs = CF_dd.SIF[["data"]][["NEE_obs"]], NEE_Residual = CF_dd.SIF[["data"]][["NEE_VPRM"]] - CF_dd.SIF[["data"]][["NEE_obs"]])
#NEE.Diff.SIF.2 <- data.frame(Date = fig_dNEE.SIF[["data"]][["date"]], NEE_Residual = fig_dNEE.SIF[["data"]][["NEE_VPRM"]] - fig_dNEE.SIF[["data"]][["NEE_obs"]])

##Subset only rows with both EC obs and VPRM NEE data values (Residuals were calculated only at time stamps where both quantities were available-Hilton 2014)
NEE.residual.SIF <- NEE.Diff.SIF[complete.cases(NEE.Diff.SIF["NEE_Residual"]),]
#vr.test <- subset(NEE.Diff.SIF, !is.na(NEE_Residual)) Alternate Method


####Exclude all other months (e.g. May and November) just in case they are in the dataset)####
NEE.residual.SIF ["Month"] <- as.data.frame(months.Date(NEE.residual.SIF$Date)) #Create and add months column
NEE.residual.SIF <- subset(NEE.residual.SIF, Month=="June" | Month=="July" | Month=="August" | Month=="September" | Month=="October")
NEE.residual.SIF$Month <- NULL #Drop or remove the months column


detach(SIF_optimized_parameters)


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

###############################################################PLOT GRAPHS##################################################################################

##################################################################Diurnal Cycle Analysis################################################################

#Create and add new columns for day, month, year, and time; extracted from the original date column
NEE.residual.SIF ["Day"] <- as.data.frame(days(NEE.residual.SIF$Date))
NEE.residual.SIF ["Month"] <- as.data.frame(months.Date(NEE.residual.SIF$Date))
NEE.residual.SIF ["Year"] <- as.data.frame(years(NEE.residual.SIF$Date))
#NEE.residual.SIF ["Hours"] <- as.data.frame(hours(NEE.residual.SIF$Date))
#NEE.residual.SIF ["Minutes"] <- as.data.frame(minutes(NEE.residual.SIF$Date))
#NEE.residual.SIF ["Seconds"] <- as.data.frame(seconds(NEE.residual.SIF$Date))
NEE.residual.SIF ["Time"] <- data.frame(substr(NEE.residual.SIF[["Date"]], start = 11, stop = 15))

###Subset the relevant columns for aggregation
variables.daily.SIF <- c("NEE_VPRM", "NEE_obs", "Time")
NEE.residual.SIF.subset <- NEE.residual.SIF[variables.daily.SIF]
NEE.residual.SIF.Daily <- aggregate(.~ Time, data = NEE.residual.SIF.subset, mean, na.action = na.pass)
NEE.residual.Line.Daily <- NEE.residual.SIF.Daily #Duplicate dataframe to avoid conflict and use for other graphs

##Melt to prepare data form ggplot multiple variables
NEE.residual.Line.Daily.Melted <- reshape2::melt(NEE.residual.Line.Daily, id.var='Time')
head(NEE.residual.Line.Daily.Melted)
tail(NEE.residual.Line.Daily.Melted)



###########Create plots of the mean diurnal cycle or variation for the growing season
library(ggplot2)
#library(ggpubr) #First install the package ggpubr #Required for the function ggarrange
#plot <- ggplot(data = NEE.residual.EVI.Daily, aes(x=Time, y=NEE_VPRM, group=1)) + geom_line() + geom_point(shape = 1, size = 3)
#plot <- plot + geom_line(data = NEE.residual.EVI.Daily, aes(y = NEE_obs), size = 1.2) + geom_point(aes(y = NEE_obs), shape = 15, size = 3)
plot <- ggplot(NEE.residual.Line.Daily.Melted, aes(x=Time, y=value, group=variable)) + geom_line(aes(color=variable), size = 1) + geom_point(aes(color=variable), size = 3)
plot <- plot + scale_color_manual(values = c("#00BA38", "#619CFF")) #New line for ONLY-SIF plot to set appropriate colors##################
plot <- plot + xlab("Time") + ylab(expression(NEE~(mu*mol~m^{-2}~s^{-1}))) + ggtitle(paste("CF Mean Diurnal NEE" ,"(June-October, 2007-2011*)"))
#plot <- plot + scale_x_discrete(breaks = c("00:00","01:00","02:00","03:00","04:00","05:00","06:00","07:00","08:00","09:00","10:00","11:00","12:00","13:00","14:00","15:00","16:00","17:00","18:00","19:00","20:00","21:00","22:00","23:00"))
#plot <- plot + scale_x_discrete(breaks = c("00:00","02:00","04:00","06:00","08:00","10:00","12:00","14:00","16:00","18:00","20:00","22:00"))
plot <- plot + scale_x_discrete(breaks = c("00:00","03:00","06:00","09:00","12:00","15:00","18:00","21:00"))
plot <- plot + theme(plot.title = element_text(hjust = 0.5)) #center title
plot <- plot + theme(axis.text=element_text(size=12,face="bold"))
plot <- plot + theme(axis.title=element_text(size=12))
plot <- plot + theme(plot.title=element_text(size=11,face="bold"))
#plot <- plot + theme(legend.key.size = unit(1.0, "cm"),legend.key.width = unit(1.0,"cm")) 
plot <- plot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank())
#plot <- plot + labs(shape = "") #Remove or specify Legend title
plot <- plot + theme(legend.title = element_blank()) #Remove Legend title
plot <- plot + theme(legend.text = element_text(size = 11, face = "bold"))
plot <- plot + theme(legend.position = "bottom")
plot <- plot + geom_hline(yintercept=0, color="black")
plot <- plot + scale_y_continuous(breaks=seq(-1.5,1.5,0.5)) #Sets the axis range and interval ("by")
plot <- plot + coord_cartesian(ylim = c(-1.5, 1.5)) #Forces the limits to be set accordingly
#plot
#save graph to output as high resoultion, Windows compatible format .EMF
fileout <- paste("Line-CF Mean Diurnal NEE","_LineCycle",".emf",sep="")
x11(width=6, height=4) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
dev.cur()
plot
savePlot(fileout, "emf", device = dev.cur())
dev.off()
graphics.off()








########################Seasonal Cycle Analysis##########################################################
##SIF VALID - Subset the relevant columns for aggregation and convert month names to month numbers
variables.monthly.SIF <- c("NEE_VPRM", "NEE_obs", "Month")
NEE.residual.SIF.subset <- NEE.residual.SIF[variables.monthly.SIF]
NEE.residual.SIF.subset$Month = match(NEE.residual.SIF.subset$Month,(month.name)) #Convert month names to month numbers before agg for proper ordering
NEE.residual.SIF.Monthly <- aggregate(.~ Month, data = NEE.residual.SIF.subset, mean, na.action = na.pass)
NEE.residual.Line.Monthly <- NEE.residual.SIF.Monthly #Duplicate dataframe to avoid conflict and use for other graphs

##Melt to prepare data form ggplot multiple variables
NEE.residual.Line.Monthly.Melted <- reshape2::melt(NEE.residual.Line.Monthly, id.var='Month')
head(NEE.residual.Line.Monthly.Melted)
tail(NEE.residual.Line.Monthly.Melted)



###########Create plots of the mean monthly growing season NEE variation
#plot <- ggplot(data = NEE.residual.EVI.Monthly, aes(x=Month, y=NEE_VPRM, group=1)) + geom_line() + geom_point(shape = 1, size = 3)
#plot <- plot + geom_line(data = NEE.residual.EVI.Monthly, aes(y = NEE_obs), size = 1.2) + geom_point(aes(y = NEE_obs), shape = 15, size = 3)
plot <- ggplot(NEE.residual.Line.Monthly.Melted, aes(x=Month, y=value, group=variable)) + geom_line(aes(color=variable), size = 1) + geom_point(aes(color=variable), size = 3)
plot <- plot + scale_color_manual(values = c("#00BA38", "#619CFF")) #New line for ONLY-SIF plot to set appropriate colors##################
plot <- plot + xlab("Month") + ylab(expression(NEE~(mu*mol~m^{-2}~s^{-1}))) + ggtitle(paste("CF Mean Monthly NEE" ,"(June-October, 2007-2011*)"))
plot <- plot + scale_x_continuous(labels =c("June","July","August","September","October"))
#plot <- plot + scale_x_continuous(labels=c("6" = "June", "7" = "July", "8" = "August", "9" = "September", "10" = "October"))
plot <- plot + theme(plot.title = element_text(hjust = 0.5)) #center title
plot <- plot + theme(axis.text=element_text(size=12,face="bold"))
plot <- plot + theme(axis.title=element_text(size=12))
plot <- plot + theme(plot.title=element_text(size=11,face="bold"))
#plot <- plot + theme(legend.key.size = unit(1.0, "cm"),legend.key.width = unit(1.0,"cm")) 
plot <- plot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank())
plot <- plot + labs(shape = "") #Remove or specify Legend title
plot <- plot + theme(legend.title = element_blank()) #Remove Legend title
plot <- plot + theme(legend.text = element_text(size = 11, face = "bold"))
plot <- plot + theme(legend.position = "bottom")
plot <- plot + geom_hline(yintercept=0, color="black")
plot <- plot + scale_y_continuous(breaks=seq(-1.5,1.0,0.5)) #Sets the axis range and interval ("by")
plot <- plot + coord_cartesian(ylim = c(-1.5, 1.2)) #Forces the limits to be set accordingly
#plot
#save graph to output as high resoultion, Windows compatible format .EMF
fileout <- paste("Line-NARR-CF Mean Monthly NEE","_LineCycle",".emf",sep="")
x11(width=6, height=4) #default dim is 7 #Open new device #2 (note: no #1 because #1 is null)
dev.cur()
plot
savePlot(fileout, "emf", device = dev.cur())
dev.off()
graphics.off()




cat("\n"); 
VPRMstate = "Model Run and Data Analysis Completed for CF!"
cat("\n"); 
cat(VPRMstate);cat("\n"); cat("\n")

