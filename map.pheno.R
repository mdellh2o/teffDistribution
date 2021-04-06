options(stringsAsFactors = F)

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/TEFF/TEFF.distribution/analysis/2.map.phenotypes"
setwd(wd)

library(ggplot2)
library(gstat)
library(sp)
library(maptools)

#get in db
tab<-read.delim("../../data/geo.pheno.final.txt", header=T)
rownames(tab)<-tab[,2]
tab<-tab[,-2]

#extract coordinates
pheno<-tab[,c(10,9, 70:ncol(tab))]
colnames(pheno)[1:2]<-c("x", "y")

#remove incomplete observations
pheno<-na.omit(pheno)

#make it a spatial object
coordinates(pheno) = ~x+y
#check it by plotting
plot(pheno)

# first get the range in data
x.range <- as.integer(range(pheno@coords[,1]))
y.range <- as.integer(range(pheno@coords[,2]))

#here we will use a chosen amount of decimal degrees:
dd<-0.05
grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=dd), y=seq(from=y.range[1], to=y.range[2], by=dd))

## convert grid to SpatialPixel class
coordinates(grd) <- ~ x+y
gridded(grd) <- TRUE

#project points on grid
plot(grd)
points(pheno, pch=20, col='red', cex=1)

#do the interpolation!
idw<-idw(formula=SC ~ 1, locations=pheno, newdata=grd)
#get it out in a dataframe
idw.output=as.data.frame(idw)
names(idw.output)[1:3]<-c("long","lat","predicted")

#PLOT IT
plot<-ggplot(data=idw.output,aes(x=long,y=lat))
layer1<-c(geom_tile(data=idw.output,aes(fill=predicted)))
plot+layer1+scale_fill_gradient(low="white", high="brown")+coord_equal()

#try out Kringing
variog<-variogram(SC~1, locations=pheno, data=pheno)
plot(variog)
variog

range=1
sill=2
nugget=2.2

model.variog<-vgm(psill=sill, model="Lin", nugget=nugget, range=range)

fit.variog<-fit.variogram(variog, model.variog)
plot(variog, fit.variog)

#do the krig
krig<-krige(formula=SC ~ 1, locations=pheno, newdata=grd, model=model.variog)

#fix it 
krig.output=as.data.frame(krig)
names(krig.output)[1:3]<-c("long","lat","var1.pred")

#plot it
plot<-ggplot(data=krig.output,aes(x=long,y=lat))
layer1<-c(geom_tile(data=krig.output,aes(fill=var1.pred)))
plot+layer1+scale_fill_gradient(low="#FEEBE2", high="#7A0177")+coord_equal()
