library(sf) # Simple Features for R
library(tmap) # Thematic Maps
library(spdep) # Spatial Dependence: Weighting Schemes, Statistics
library(tidyverse) # Easily Install and Load the 'Tidyverse'
library(viridis) # Default Color Maps from 'matplotlib'
library(corrplot) # Visualization of a Correlation Matrix
library(Hmisc) # Harrell Miscellaneous
library(PerformanceAnalytics) # Econometric Tools for Performance and Risk Analysis
library(spgwr) # Geographically Weighted Regression
library(car) # Companion to Applied Regression
library(tidyverse)

shp04 <- st_read('/Users/gebeilei/OneDrive\ -\ University\ of\ Southampton/ucl/CASA0005/statistical-gis-boundaries-london/ESRI/MSOA_2004_London_High_Resolution.shp')
shp11 <- st_read('/Users/gebeilei/OneDrive\ -\ University\ of\ Southampton/ucl/CASA0005/statistical-gis-boundaries-london/ESRI/MSOA_2011_London_gen_MHW.shp')
ggplot(shp04) +
  geom_sf()
ggplot(shp11) +
  geom_sf()

df.2013to2014 <- read_csv('/Users/gebeilei/OneDrive\ -\ University\ of\ Southampton/ucl/CASA0005/2013-2014.csv') 
df.2011to2012 <- read_csv('/Users/gebeilei/OneDrive\ -\ University\ of\ Southampton/ucl/CASA0005/2011-2012.csv')

library(zoo)
df.2013to2014[,6:13] <- na.aggregate(df.2013to2014[,6:13])
df.2011to2012[,6:13] <- na.aggregate(df.2011to2012[,6:13])

df.2013to2014[,6:9] <- df.2013to2014[,6:9]*100
df.2011to2012[,6:9] <- df.2011to2012[,6:9]*100

shp.1314 <- left_join(shp04, df.2013to2014, by = c("MSOA_CODE" = "Code"))
shp.1112 <- left_join(shp04, df.2011to2012, by = c("MSOA_CODE" = "Code"))
shp.1314.2 <- left_join(shp11, df.2013to2014, by = c("MSOA11CD" = "Code"))
#shp.1314 <- merge(shp04, df.2013to2014, by.x='MSOA_CODE', by.y = 'Code')
#shp.1112 <- merge(shp, df.2011to2012, by.x='MSOA11CD', by.y = 'Code')

ggplot(shp.1112) +
  geom_sf(aes(fill=obese4_5))
ggplot(shp.1314) +
  geom_sf(aes(fill=obese4_5))

df.1112 <- df.2011to2012 %>% select(obese4_5,obese10_11,weight4_5,weight10_11,TWI,NWI,NIBHC,NIAHC)
chart.Correlation(df.1112, histogram=TRUE, pch=19, method = 'pearson') # see the correlation
df.1314 <- df.2013to2014 %>% select(obese4_5,obese10_11,weight4_5,weight10_11,TWI,NWI,NIBHC,NIAHC)
chart.Correlation(df.1314, histogram=TRUE, pch=19, method = 'pearson') # see the correlation
# selecting variance
model.1112 <- lm(obese4_5 ~ TWI, data = df.1112) # OLS model
summary(model.1112)
model.1112 <- lm(obese4_5 ~ NWI, data = df.1112) # OLS model
summary(model.1112)
model.1112 <- lm(obese4_5 ~ NIBHC, data = df.1112) # OLS model
summary(model.1112)
model.1112 <- lm(obese4_5 ~ NIAHC, data = df.1112) # OLS model
summary(model.1112)
model.1314 <- lm(obese4_5 ~ TWI, data = df.1314) # OLS model
summary(model.1314)
model.1314 <- lm(obese4_5 ~ NWI, data = df.1314) # OLS model
summary(model.1314)
model.1314 <- lm(obese4_5 ~ NIBHC, data = df.1314) # OLS model
summary(model.1314)
model.1314 <- lm(obese4_5 ~ NIAHC, data = df.1314) # OLS model
summary(model.1314)
model.1112S <- lm(obese10_11 ~ NIAHC, data = df.1112) # OLS model
summary(model.1112S)
model.1314S <- lm(obese10_11 ~ NIAHC, data = df.1314) # OLS model
summary(model.1314S)

shp.1112$ols_residuals_1112 <- residuals(model.1112)
shp.1314$ols_residuals_1314 <- residuals(model.1314)

#shp.point$ols_residuals <- residuals(model)
ggplot(data=shp.1112) +
  geom_sf(aes(fill=ols_residuals_1112))
ggplot(data=shp.1314) +
  geom_sf(aes(fill=ols_residuals_1314))

set.ZeroPolicyOption(TRUE)
neighbours <- poly2nb(shp04)
#shp.sp <- as(shp04,'Spatial')
#coords <- st_coordinates(st_centroid(st_geometry(shp04)))
plot(st_geometry(shp04))
plot(neighbours, coordinates(as(shp04, "Spatial")), add=TRUE, col="blue")

listw <- nb2listw(neighbours, style="C")
moran.test(shp.1112$ols_residuals_1112, listw = listw)
moran.test(shp.1314$ols_residuals_1314, listw = listw)

moran.test(shp.1112$obese4_5, listw = listw)
moran.test(shp.1314.2$obese4_5, listw = listw)

moran.test(shp.1112$obese10_11, listw = listw)
moran.test(shp.1314.2$obese10_11, listw = listw)

shp1112.sp <- as(shp.1112, 'Spatial')
GWRbandwidth.1112 <- gwr.sel(obese4_5 ~ NIAHC, data=shp1112.sp, adapt=T, method="cv") #calculate kernel bandwidth
GWRModel.1112 <- gwr(obese4_5 ~ TWI, data=shp1112.sp, adapt=GWRbandwidth.1112,
                     hatmatrix = TRUE,se.fit = TRUE) #run the gwr model
GWRModel.1112 #print the results of the model


#shp.1314.2 <- left_join(shp11, df.2013to2014, by = c("MSOA11CD" = "Code"))
shp1314.sp <- as(shp.1314.2, 'Spatial')
GWRbandwidth.1314 <- gwr.sel(obese4_5 ~ NIAHC, data=shp1314.sp, adapt=T, method="cv") #calculate kernel bandwidth
GWRModel.1314 <- gwr(obese4_5 ~ NIAHC, data=shp1314.sp, adapt=GWRbandwidth.1314,
                     hatmatrix = TRUE,se.fit = TRUE) #run the gwr model
GWRModel.1314 #print the results of the model




