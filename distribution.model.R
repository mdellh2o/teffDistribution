options(stringsAsFactors = F)
maindir<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/TEFF/TEFF.distribution/analysis/5.spatial.modelling"
setwd(maindir)

#USED PACKAGES
library(dismo)
library(maptools)
library(raster)
library(measurements)
library(patchwork)
library(car)


#get previous data
load("../3.paper.analysis/diversity.data.Rdata")

#get coordinate data
gps<-read.delim("../../data/coordR.txt")
gps$Lat = measurements::conv_unit(gps$Lat, from = 'deg_dec_min', to = 'dec_deg')
gps$Long = measurements::conv_unit(gps$Long, from = 'deg_dec_min', to = 'dec_deg')

coord<-gps[,2:3]
coord<-apply(coord, 2, as.numeric)
rownames(coord)<-gps[,1]

#reduce coordinates to those used for bioclim and heno analyses
coord<-coord[rownames(coord) %in% rownames(bc),]

#get ethiopian extension
eth <- getData("GADM", country="ETH", level=0)
eth1 <- getData("GADM", country="ETH", level=1)

#get environmental data
currentEnv=getData("worldclim", var="bio", res=2.5)
futureEnv=getData('CMIP5', var='bio', res=2.5, rcp=45, model='HE', year=70)
names(futureEnv)=names(currentEnv)

#crop climate data to study area
model.extent<-extent(eth)
modelEnv=crop(currentEnv,model.extent)
modelFutureEnv=crop(futureEnv, model.extent)

#save current data
save(modelEnv, modelFutureEnv, eth, eth1, coord, file="distribution.model.step.1.Rdata")


##########restart from here
maindir<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/TEFF/TEFF.distribution/analysis/5.spatial.modelling"
setwd(maindir)

library(raster)
library(rgdal)
library(maps)
library(mapdata)
library(dismo) 
library(rJava) 
library(maptools)
library(corrplot)
library(jsonlite)
library(ggplot2)
library(patchwork)
library(car)
library(BiodiversityR)
library(psych)
library(RStoolbox)
library(mclust)

load(file="distribution.model.step.1.Rdata")

#create spatial points
sp<-SpatialPoints(coord)

#do some plotting to check that everything is right
plot(modelEnv[["bio1"]]/10, main="Annual Mean Temperature")
points(sp, pch="+", cex=0.2)

plot(modelFutureEnv[["bio1"]]/10, main="Annual Mean Temperature")
points(sp, pch="+", cex=0.2)

#create a set of random points in the same number of presence points
bg <- randomPoints(modelEnv, nrow(coord))

##########
#check redundancy in bioclim variables using the current values
modelEnv<-stack(modelEnv)
vif <- ensemble.VIF(
  x = modelEnv,
  a = data.frame(sp),
  VIF.max = 10,
  keep = NULL,
  layer.drops = "bio7",
  factors = NULL,
  dummy.vars = NULL
)

tokeep<-names(vif$VIF.final)

###drop colinear variables from environmental datasets
redbc<-which(names(modelEnv) %in% tokeep)
modelEnv<-modelEnv[[redbc]]
modelEnv<-brick(modelEnv)

modelFutureEnv<-stack(modelFutureEnv)
modelFutureEnv<-modelFutureEnv[[redbc]]
modelFutureEnv<-brick(modelFutureEnv)

#set up the model
############under construction##############
k=10
fold <- kfold(coord, k=k) # add an index that makes five random groups of observations
teftest <- coord[fold == 1, ] # hold out one fifth as test data
teftrain <- coord[fold != 1, ] # the other four fifths are training data
e <- list()
for (i in 1:k) {
  teftest <- coord[fold == i, ] # hold out one fifth as test data
  teftrain <- coord[fold != i, ] # the other four fifths are training data
  tef.me <- bioclim(modelEnv, teftrain)
  e[[i]] <- evaluate(p=teftest, a=bg,  model=tef.me, x=modelEnv)
}

pdf("standard.roc.pdf")
  plot(e[[10]], 'ROC')
dev.off()

#consolidate validation results
auc <- sapply(e, function(x){x@auc})
mean(auc)
box<-ggplot() + geom_boxplot(aes(y=auc), fill="gray60") + theme_bw() + labs(y="AUC", title = "cross validation") +   theme(axis.title.x=element_blank(),
                                                                                          axis.text.x=element_blank(),
                                                                                          axis.ticks.x=element_blank())
box
ggsave("boxplot.AUC.pdf", height = 6, width=3)

#see spatial sorting bias issue
#https://rspatial.org/sdm/5_sdm_models.html

#create a subset of absence points
s <- sample(nrow(bg), 0.5 * nrow(bg))
back_train <- bg[-s, ]
back_test <- bg[s, ]

#check sampling bias
sb <- ssb(teftest, back_test, teftrain)
sb[,1] / sb[,2]
#very close to 0, suggesting BIG sorting bias

#fix this subsampling the sample
#select pairs of points from the two sets that have a similar distance to their nearest point in another set of points
i <- pwdSample(teftest, back_test, teftrain, n=100, tr=0.1, warn=T)
pres_test_pwd <- teftest[!is.na(i[,1]), ]
dim(pres_test_pwd)
back_test_pwd <- back_test[na.omit(as.vector(i)), ]
dim(back_test_pwd)

sb2 <- ssb(pres_test_pwd, back_test_pwd, teftrain)
sb2[,1]/ sb2[,2]
e
e<-evaluate(tef.me, p=pres_test_pwd, a=back_test_pwd, x=modelEnv)

pdf("fixed.roc.pdf")
  plot(e, 'ROC')
dev.off()

##################################
#run the real model using all data available

tef.me <- bioclim(modelEnv, coord)
plot(tef.me)
response(tef.me)

#evaluate the model
e<-evaluate(tef.me, p=coord, a=bg, x=modelEnv)

# Determine minimum threshold for "presence"
bc.threshold <- threshold(x = e, stat = "equal_sens_spec")
bc.threshold

#predict suitability to current conditions
tef.pred <- predict(tef.me, modelEnv)

#subset to values above threshold
tef.pred[tef.pred < bc.threshold] <- NA

pdf("current.suitability.pdf", height=6, width=6)
  plot(tef.pred, main="Current suitability", xlab="Lon", ylab="Lat", 
       bty="n", box=FALSE, yaxt="n",xaxt="n")
  lines(eth, col="gray50")
  lines(eth1, col="gray50", lty=2)
  points(sp, pch="+", cex=0.2)
  axis(2)
  axis(1)
dev.off()

#get future prediction
tef.2070 = predict(tef.me, modelFutureEnv)

#subset to values above threshold
tef.2070[tef.2070 < bc.threshold] <- NA

pdf("future.suitability.pdf", height=6, width=6)
  plot(tef.2070, main="Projected suitability", xlab="Lon", ylab="Lat", 
       bty="n", box=FALSE, yaxt="n",xaxt="n")
  lines(eth, col="gray50")
  lines(eth1, col="gray50", lty=2)
  points(sp, pch="+", cex=0.2)
  axis(2)
  axis(1)
dev.off()

#plot the predicted change in suitability
#with this computation, 0 means unchanged suitability. Positive values are an improvement, negative values are pejorative
tef.change=tef.2070-tef.pred
#set color
## Make a vector with n colors
ncol <- 100
cols <- RColorBrewer:::brewer.pal(11,"PuOr")  # OR c("purple","white","orange")  
rampcols <- colorRampPalette(colors = cols, space="Lab")(ncol)
rampcols[(ncol/2)] <- rgb(t(col2rgb("white")), maxColorValue=256) 

## Make a vector with n+1 breaks
rampbreaks <- seq(-1, 1, length.out = ncol)

#make arguments for the legend plotting
arg<-list(at=c(-1,-0.5,0,0.5,1), labels=c(-1,-0.5,0,0.5,1))


pdf("change.suitability.pdf", height=6, width=6)
  plot(tef.change, main="Projected change in suitability", bty="n", 
       box=FALSE, col=rampcols, breaks=rampbreaks, axis.args=arg, yaxt="n",xaxt="n")
  lines(eth,  col="gray50")
  lines(eth1, col="gray50", lty=2)
  points(sp, pch="+", cex=0.2)
  axis(2)
  axis(1)
dev.off()

#see how the current locations of tef will look like in the future
change = extract(tef.change, coord)
pdf("suitability variation.pdf")
  hist(change, main="", col="gray90", xlab="Suitability variation")
  abline(v=0, col="red", lty=2)
dev.off()

#get correlation between rate change and variation 
current = extract(modelEnv, coord)
forcor<-cbind(current, change)
corz<-corr.test(forcor, method ="spearman", use="pairwise.complete")

outr<-corz$r[1:(ncol(corz$r)-1),ncol(corz$r)]
outp<-corz$p[1:(ncol(corz$r)-1),ncol(corz$r)]

out<-cbind(outr,outp)
write.table(out, file="change.correlation.txt", sep="\t", row.names=F, quote=F)

#########################
#make publishable plots

#current conditions
mid<-0
current<-ggplot() +
        geom_raster(data = tef.pred, aes(x=x, y=y,  fill = layer )) + 
        scale_fill_gradient2(midpoint = mid, low = "white",
                             high = "#0e8562", space = "Lab", limits=c(0,0.8), na.value = "white") +
       #geom_point(data=sp, aes(x=long, y=lat), shape=3, size=0.5, col="gray50", alpha=0.5)+
        geom_polygon(data=eth, aes(x=long, y=lat), 
                     fill=NA,color="grey50") +
        geom_polygon(data=eth1, aes(x=long, y=lat, group=group), 
                     fill=NA,color="grey50", linetype = "dashed") +
        coord_quickmap() +
        theme_classic() + labs(fill = "Suitability", title="Current suitability", x = "Long", y="Lat")
current
ggsave("current.suitability.hires.pdf", width = 12, height = 10)

current<-current+labs(tag="A)")

#future conditions
mid<-0
future<-ggplot() +
    geom_raster(data = tef.2070, aes(x=x, y=y,  fill = layer )) + 
    scale_fill_gradient2(midpoint = mid, low = "white",
                         high = "#0e8562", space = "Lab" , limits=c(0,0.8), na.value = "white") +
    #geom_point(data=sp, aes(x=long, y=lat), shape=3, size=0.5, col="gray50", alpha=0.5)+
    geom_polygon(data=eth, aes(x=long, y=lat), 
                 fill=NA,color="grey50") +
    geom_polygon(data=eth1, aes(x=long, y=lat, group=group), 
                 fill=NA,color="grey50", linetype = "dashed") +
    coord_quickmap() +
    theme_classic() + labs(fill = "Suitability", title="Projected suitability", x = "Long", y="Lat")
future
ggsave("future.suitability.hires.pdf", width = 12, height = 10)

future<-future + labs(tag="B)")

#variation in conditions
mid<-0
tefchange<-ggplot() +
  geom_raster(data = tef.change, aes(x=x, y=y,  fill = layer )) + 
  scale_fill_gradient2(midpoint = mid, low = "red", mid = "white",
                       high = "blue", space = "Lab", limits=c(-0.7,0.7), na.value = "white") +
  #geom_point(data=sp, aes(x=long, y=lat), shape=3, size=0.5, col="gray50", alpha=0.5)+
  geom_polygon(data=eth, aes(x=long, y=lat), 
               fill=NA,color="grey50") +
  geom_polygon(data=eth1, aes(x=long, y=lat, group=group), 
               fill=NA,color="grey50", linetype = "dashed") +
  coord_quickmap() +
  theme_classic() + labs(fill = "Variation", title="Change in suitability", x = "Long", y="Lat")
tefchange
ggsave("change.suitability.hires.pdf", width = 12, height = 10)

tefchange<-tefchange + labs(tag="C)")
  
#histogram of variation
changehist<-ggplot() + geom_density(aes(x=change), color="gray20", fill="gray90") +
                geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
                theme_classic() + labs(tag="D)", x="Change in suitability", title="Variation in suitabilty in sampling points")
changehist  

#put everything together

fig4<-(current | future) / (tefchange | changehist)
fig4

ggsave("fig.4.pdf", width = 12, height = 10)     

############## check change in suitability in different phenotypic classes
#get phenotypes
load("../3.paper.analysis/diversity.data.Rdata")
names(qq)
#select only few qualitative phenotypes
phenolist<-c("Mdatec", "Gfpc", "Ptype", "Hdatec", "Plc")

#fix typos
qq[qq=="l"]<-1

#get coordinate data
gps<-read.delim("../../data/coordR.txt")
gps$Lat = measurements::conv_unit(gps$Lat, from = 'deg_dec_min', to = 'dec_deg')
gps$Long = measurements::conv_unit(gps$Long, from = 'deg_dec_min', to = 'dec_deg')
coord<-gps[,2:3]
coord<-apply(coord, 2, as.numeric)
rownames(coord)<-gps[,1]

#create subsets according to phentoypic classes
phenoclass<-list()
for (p in 1:length(phenolist)){
  tmppheno<-phenolist[p]
  tmp<-data.frame(cbind(rownames(qq), qq[,tmppheno]))
  #get rownames for each phenotypic class
  rntokeep<-list()
  classes<-unique(tmp[,2])
  #get change in suitability for phenotypic classes
  for(i in 1:length(classes)){
    hit<-tmp[which(tmp[,2]==classes[i]),1]
    coordtokeep<-coord[which(rownames(coord) %in% hit),]
    #here we get change values
    rntokeep[[i]]<-data.frame(class=rep(classes[i], nrow(coordtokeep)), change=extract(tef.change, coordtokeep))
  }
  names(rntokeep)<-classes
  phenoclass[[p]]<-rntokeep
}#for p
names(phenoclass)<-phenolist
lapply(phenoclass, names)

#now derive dataframes in tall format
changeclass<-list()
for (i in 1:length(phenoclass)){
  changeclass[[i]]<- do.call(rbind, phenoclass[[i]])
  rownames(changeclass[[i]])<-NULL
}#for i

names(changeclass)<-names(phenoclass)

#make plots and models
plts<-list()
mods<-list()
for (i in 1:length(changeclass)){
  plts[[i]]<- ggplot(data = changeclass[[i]], aes(x=class, y=change, fill=class))+geom_boxplot() + labs(title=names(changeclass)[i])
  mod1<-aov(change ~ class , data=changeclass[[i]])
  mods[[i]]<-summary(mod1)
}
names(mods)<-names(phenoclass)
mods

##################
#perform clustering on phenotypes

qt<-apply(qt,2,as.numeric)
qtsub<-qt[, c("Hdate", "GfP", "PL")]
BIC<-mclustBIC(qtsub)
BIC
plot(BIC)
mod1<-Mclust(qtsub, x = BIC)
summary(mod1, parameters = TRUE)

#writ outputs
plot(mod1)
dfout<-data.frame(ID=names(mod1$classification), group=mod1$classification)
write.table(dfout, file="clustering.output.txt", sep="\t", quote=F, row.names=F)

#get groups and coordinates in each group
grps<-mod1$classification
gpslist<-list()
for (i in 1:max(grps)){
  tmp<-i
  hit<-as.numeric(names(grps)[which(grps==i)])
  gpslist[[i]]<-coord[which(rownames(coord) %in% hit),]
}

#run the SDM for each group indipendently to get change in suitability
changelist<-list()
for (i in 1:length(gpslist)){
  tmp<-gpslist[[i]]
  #train the model in points
  tef.tmp <- bioclim(modelEnv, tmp)
  #model current conditions
  tef.tmp.current <- predict(tef.tmp, modelEnv)
  #model future conditions
  tef.tmp.2070 <- predict(tef.tmp, modelFutureEnv)
  #get change in suitability
  tmpchange<-tef.tmp.2070-tef.tmp.current
  changelist[[i]]<-data.frame(cluster=i, change=extract(tmpchange, tmp))
}#for i
changetall<-do.call(rbind, changelist)
head(changetall)
changetall[,1]<-as.factor(changetall[,1])

#plot and test differences
mod1<-lm(change ~ cluster , data=changetall)
summary(mod1)

p1<-ggplot(data=changetall, aes(x=cluster, y = change, group=cluster, col=cluster)) + geom_boxplot() + ylim(-0.6,0.6)+
  geom_hline(yintercept = mean(changetall$change), linetype = 2, col="red") +
  geom_hline(yintercept = 0, linetype = 2, col="black") +
  stat_compare_means(method = "anova", label.x= 2, label.y = 0.6)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", hide.ns = TRUE, label.y = 0.4) + 
  theme_bw()
p1 + theme(legend.position = "none")

ggsave("model.by.class.pdf", width = 12, height = 10)     



############## merge these info with diveristy
#get previous data
load("../3.paper.analysis/diversity.data.Rdata")

#sort appropriately all datasets
coord<-coord[order(match(rownames(coord),rownames(bc))),]
stopifnot(all(rownames(qt)==rownames(bc)))
stopifnot(all(rownames(qt) == rownames(coord)))

#get some correlations with bioclim PCA
cor.test(bc.pca$x[,1], change)
cor.test(bc.pca$x[,2], change)
cor.test(bc.pca$x[,3], change)

#get some correlations with bioclim PCA
cor.test(pcaqt$x[,1], change)
cor.test(pcaqt$x[,2], change)
cor.test(pcaqt$x[,3], change)


#try to make a meaningful plot
df<-data.frame(bc.pca$x[,1:3],change, bc, qq, qt)
df<-data.frame(apply(df, 2, as.numeric))
cr<-cor(df, method="spearman")
corrplot(cr, "pie")
sort(cr["change",])
head(df)

pdf("suitability_change_PCA.pdf", height=4, width=8)
  p1<-ggplot(aes(x=PC1, y=change, col=Altitude), data=df)+ geom_point()  +
      scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07")) # + geom_smooth()
  p1 + geom_hline(yintercept=0, linetype="dashed", color = "red") + theme_bw()
dev.off()


#Try 3d map
alt <- getData('alt', country='ETH')
adm <- getData('GADM', country='ETH', level=0)


p<-persp(alt, exp=0.2,phi=35, xlab="Longitude", ylab="Latitude", zlab="Elevation")

