maindir<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/TEFF/TEFF.distribution/analysis/0.load.shp.add.layer"

setwd(maindir)

library(raster)
library(maptools)
library(rgdal)

#load shape of sampling points
shpdir<-"../../data/GIS"
tef<-readOGR(shpdir, "teff.distribution")
#extract extension and attributes
minmax<-extent(tef)
att<-tef@data 

#get additional data
et0<-getData('GADM' , country="ETH", level=0)

#get bioclim at max resolution using tef extension to get the right tile
bioc <- getData("worldclim",var="bio",lon = minmax[1], lat=minmax[3],res=0.5)

#make a sample plot
plot(bioc[[1]], xlim=c(32, 49), ylim=c(1, 19))
plot(et0, add=T)
plot(tef, add=T)

#extract values
values <- extract(bioc,tef)
biovar <- cbind.data.frame(att[,1],coordinates(tef),values)

#save object for 
save(att, biovar, file="attributes.and.bioclim.tef.Rdata")
