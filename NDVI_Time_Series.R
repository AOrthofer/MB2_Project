install.packages('devtools')
install_github('loicdtx/bfastSpatial')
library(devtools)
library(raster)
library(rgdal)
library(RStoolbox)

# define input and output directories
workingDirectory <- 'C:/Users/Anna/Desktop/Time_series_data/Landsat_data'
outputDirectory <- 'C:/Users/Anna/Desktop/Time_series_data/Results/'
shapefile <- 'C:/Users/Anna/Desktop/Time_series_data/Shapefile'


#set working directory
setwd(workingDirectory)

#list all files in directory
fileNames <- list.files(path = getwd(), pattern = '*MTL.txt')

for (fileName in fileNames) {
  #read meta Data
  MetaData <- readMeta(fileName) 
  #read raster data into stack
  rawData <- stackMeta(MetaData)
  
  #crop raster
  #load subset shapefile
  subset <- readOGR(shapefile, 'Study_area')
  # crop image
  dataCropped <- crop(rawData, subset)
  
  #apply cloud mask
  #set invalid data to NAs
  dataCropped[dataCropped == 0] <- NA
  #assign wavelength ranges to respective bands
  switch(MetaData$SATELLITE,
         LANDSAT5 = {
           blue = 'B1_dn'
           red = 'B3_dn'
           nir = 'B4_dn'
           tir = 'B6_dn'
         },
         LANDSAT8 = {
           blue = 'B2_dn'
           red = 'B4_dn'
           nir = 'B5_dn'
           tir = 'B10_dn'
         },
         {
           print(MetaData$SATELLITE)
           abort('unknown satellite')
         }
  )
  #calulate cloud mask by using the quality assessment layer
  qaLayer <- stackMeta(MetaData, category = 'qa')
  qaLayerCrop <- crop(qaLayer, subset)
  
  #shift bits to test for clouds and cloud shadows, cloud bits: 4-6; high confidence; shadows bits: 7-8, high confidence
  qaLayerCloudMask <- calc(qaLayerCrop, fun = function(x) bitwShiftR(x, 4)) ##shift values by 4 bits
  qaLayerCloudMask <- calc(qaLayerCloudMask, fun = function(x) bitwAnd(x, 7)) ##mask first 3 bits; if result is 7 then there are clouds
  qaLayerShadowMask <- calc(qaLayerCrop, fun = function(x) bitwShiftR(x, 7)) ##shift values by 7 bits
  qaLayerShadowMask <- calc(qaLayerShadowMask, fun = function(x) bitwAnd(x, 3)) ##mask first 2 bits, if result is 3 then there is cloud shadow
  
  #produce binary cloud mask
  qaLayerCloudMask[qaLayerCloudMask != 7] <- 0 # all pixels without clouds set to 0
  qaLayerCloudMask[qaLayerCloudMask == 7] <- 1 # all pixels with clouds set to 1
  qaLayerShadowMask[qaLayerShadowMask != 3] <- 0 # all pixels without shadows set to 0
  qaLayerShadowMask[qaLayerShadowMask == 3] <- 1 # all pixels with shadows set to 1
  #set cloud and shadow covered pixels to NA
  maskedImage <- mask(dataCropped, qaLayerCloudMask, maskvalue = 1, updatvalue = NA)
  maskedImage <- mask(maskedImage, qaLayerShadowMask, maskvalue = 1, updatevalue = NA)
  
  
  #calculate NDVI
  NDVI <- spectralIndices(maskedImage, red = red, nir = nir, indices = 'NDVI')
  
  #save raster
  resultFileName <- paste(outputDirectory, 'NDVI_', substr(MetaData$ACQUISITION_DATE, 1, 4),'_', MetaData$SATELLITE, '.tif', sep = '')
  writeRaster(NDVI, resultFileName, format ='GTiff', overwrite = TRUE)
} 


#Time Series Analysis!
#stack NDVI files
NDVIPath <- 'C:/Users/Anna/Desktop/Time_series_data/Results'
allNDVIFiles <- list.files(NDVIPath, full.names = TRUE, pattern='*.tif$')
LandsatNDVI <- stack(allNDVIFiles)

#calculate the mean of all NDVI layers
meanNDVI <- cellStats((LandsatNDVI), stat ='mean', na.rm = TRUE)

#create a dataframe with years and mean values which are set to NA
year <- c(2000:2017)
mean <- c((1:18)*NA)
dfNDVI <- data.frame(year, mean)

#Filling the empty dataframe with the mean NDVI values
for (name in names(meanNDVI)){
  print(name)
  year <- strtoi(substr(name, 6, 9))
  print(meanNDVI[name])
  dfNDVI[dfNDVI$year == year, "mean"] <- meanNDVI[name] 
  
}

#Visualize results by plotting them 
library(ggplot2)
ggplot(data = na.omit(dfNDVI), 
       aes(x = year, y = mean))+ 
  ylim(0.37,0.6)+
  labs(x = 'Year', y = 'mean NDVI') + 
  ggtitle('NDVI Development 2000 - 2017') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_line(colour = 'springgreen3', size = 1) + 
  geom_point(colour = 'springgreen3', size = 3) +
  scale_x_continuous(breaks = (dfNDVI$year)) +
  theme(panel.grid.minor =   element_blank(),
        panel.grid.major = element_line(colour = "grey",size=0.4))+
  theme(panel.background = element_rect(fill = 'grey95', colour = 'black'))



#Plot the overview map
library(maptools)
crs <- LandsatNDVI$NDVI_2003_LANDSAT5@crs
brazil <- readShapeSpatial('C:/Users/Anna/Desktop/Time_series_data/Shapefile/BRA_adm1.shp', proj4string=crs)
para = brazil[brazil$ID_1 == 14,]
plot0 <- plot(brazil, col = 'grey')
plot <- plot(para, col = 'lightgreen', add = TRUE)
scalebar(1000, xy = click(), type = 'bar', below = 'Kilometers')
title(main =' Brazil - State of Par?')
box(which = 'plot', lty = 'solid')

