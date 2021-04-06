options(stringsAsFactors = FALSE)

maindir<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/TEFF/TEFF.distribution/analysis/1.data.analysis"
setwd(maindir)
library(RColorBrewer)
library(corrplot)
library(plyr)
library(factoextra)
library(patchwork)

#####FOCUS ON THE ONES HAVING SAMPLING INFO##########

load("../0.load.shp.add.layer/attributes.and.bioclim.tef.Rdata")

#check pedigree data
#fix regions
att[,"Region"]<-as.character(att[,"Region"])
unique(att[,"Region"])
att[,"Region"]<-tolower(att[,"Region"])
att[grep("tig", att[,"Region"]), "Region"]<-"Tigray"
att[grep("snnp", att[,"Region"]), "Region"]<-"SNNP"
att[grep("am", att[,"Region"]), "Region"]<-"Amhara"
att[grep("benishangul", att[,"Region"]), "Region"]<-"Benishangul"
att[grep("soma", att[,"Region"]), "Region"]<-"Somali"
att[grep("orom", att[,"Region"]), "Region"]<-"Oromia"
att[grep("^na$", att[,"Region"]), "Region"]<-"NA"
unique(att[,"Region"])
att[,"Region"]<-factor(att[,"Region"])

#fix zones
att[,"Zone"]<-as.character(att[,"Zone"])
unique(att[,"Zone"])
att[grep("JIGJIGA", att[,"Zone"]), "Zone"]<-"Jigjiga"
att[grep("Benchi ", att[,"Zone"]), "Zone"]<-"Benchi Maji"
att[grep("Kembata Alab ", att[,"Zone"]), "Zone"]<-"Kembata Alab Tembaro"
unique(att[,"Zone"])
att[,"Zone"]<-factor(att[,"Zone"])

#order the dataframe
att<-att[order(att[,"Region"]),]

#add a column with colors according to regions
#set colors
colz<-c("#ffab7a","#261666","gray",
			"#8dda5d","#ff86eb","#9f5f00","#689eff")
coldf<-data.frame(unique(att[,"Region"]), colz)
names(coldf)<-c("region", "col")

#use the col idx to associate a color to each sample
att$colz<-mapvalues(att[,"Region"], from=coldf[,1], to=coldf[,2])

#fix classes
att[,10:27]<-apply(att[,10:27], 2, function(x) as.numeric(as.character(x)))

####################
#start looking at some distributions by region
#first, remove NAs
regatt<-att[-grep("NA", att[, "Region"]),]
regatt<-droplevels(regatt)

#make some boxplots
boxplot(Altitude_1 ~ Region, data=regatt, col=unique(as.character(regatt$colz)), las=2, cex.axis=0.8)
boxplot(Mdate ~ Region, data=regatt, col=unique(as.character(regatt$colz)), las=2, cex.axis=0.8)
#etc

###################
#forget about regions and get data altogether
alldat<-cbind.data.frame(att, biovar)
colz<-as.character(alldat[,28])
alldat<-alldat[,c(1,3,11:27,8,32:ncol(alldat))]
alldat[,3:ncol(alldat)]<-apply(alldat[,3:ncol(alldat)], 2, function(x) as.numeric(as.character(x)))
colnames(alldat)<-sub("_27$", "", colnames(alldat))

alldat<-alldat[order(alldat[,2]),]

#extract groups of data
qutr<-alldat[,c("Edate", "Hdate", "Mdate","GfilP", "PL")]
quatr<-alldat[,c("Edatec", "Hdatec", "Mdatec","Gfpc", "Plc","PCHdate", "PCMdate", "SBP", "Ptype", "seed.color")]
pheno<-cbind(qutr, quatr)
bioc<-alldat[,20:ncol(alldat)]
#fix temperature indexes reducing by 10
if(max(bioc[,2], na.rm=T)>100){
	bioc[2:12]<-bioc[2:12]*0.1		
}#if


#correlate among traits
cort<-cor(qutr, use = "complete.obs")
pdf("trait.correlations.pdf")
	corrplot(cort, cl.length=3, tl.col="black")
dev.off()

#correlate among bioclim
corb<-cor(bioc, use = "complete.obs")
pdf("bioclim.correlations.pdf")
	corrplot(corb, cl.length=3, tl.col="black")
dev.off()


#################BIOCLIM PCA########################
pcab<-prcomp(na.omit(bioc), scale=T)
eig<-(pcab$sdev)^2
vars<-eig*100/sum(eig)
vars<-round(vars,2)

pdf("bioclim.PCA.pdf")
	plot(pcab$x[,1], pcab$x[,2], col=colz, pch = 21, cex=1, 
		xlab=paste("PC1 ", vars[1], "%", sep=""), ylab=paste("PC2 ", vars[2], "%", sep=""))
	legend("topleft", legend=coldf[,1], col=coldf[,2], pch=21, cex=0.8, bty = "n")
dev.off()

#get bioclim correlations with pc axis
cors<-cor(pcab$x[,1:3], na.omit(bioc))
pdf("bioclim.PCA.correlations.pdf")
	corrplot(cors, cl.length=3, tl.col="black")
dev.off()

#plot pairs and PCA variances on axes
pdf("bioclim.PCApairs.pdf")
	pairs(pcab$x[,1:4], col=colz, pch=20)
dev.off()

pdf("bioclim.PCA.variance.barplot.pdf")
	par(mar=c(1,3,1,1))
	barplot(vars, ylab="Proportion of Variance Explained", xlab="PC")
dev.off()

#make biplot
pdf("bioclim.PCA.biplot.pdf")
	fviz_pca_biplot(pcab, repel = TRUE, label="var", title = "",
                col.var = "navyblue", # Variables color
                col.ind = "gray90"  # Individuals color
                )
dev.off()

###################################


#################Phenotype PCA########################
pcap<-prcomp(na.omit(pheno), scale=T)
eig<-(pcap$sdev)^2
vars<-eig*100/sum(eig)
vars<-round(vars,2)

pdf("phenotype.PCA.pdf")
	plot(pcap$x[,1], pcap$x[,2], col=colz, pch = 21, cex=1, 
		xlab=paste("PC1 ", vars[1], "%", sep=""), ylab=paste("PC2 ", vars[2], "%", sep=""))
	legend("topleft", legend=coldf[,1], col=coldf[,2], pch=21, cex=0.8, bty = "n")
dev.off()

#get phenotype correlations with pc axis
cors<-cor(pcap$x[,1:3], na.omit(pheno))
pdf("phenotype.PCA.correlations.pdf")
	corrplot(cors, cl.length=3, tl.col="black")
dev.off()

#plot pairs and PCA variances on axes
pdf("phenotype.PCApairs.pdf")
	pairs(pcap$x[,1:4], col=colz, pch=20)
dev.off()

pdf("phenotype.PCA.variance.barplot.pdf")
	par(mar=c(1,3,1,1))
	barplot(vars, ylab="Proportion of Variance Explained", xlab="PC")
dev.off()

#make biplot
pdf("phenotype.PCA.biplot.pdf")
	fviz_pca_biplot(pcap, repel = TRUE, label="var", title = "",
                col.var = "navyblue", # Variables color
                col.ind = "gray90"  # Individuals color
                )
dev.off()

###################################


###########################
#get bioclim distribution
#set qantities
qnt<-c("masl", "Annual Mean Temp (°C)", "Mean Diurnal Range (°C)", "Isothermality (°C)", "Temperature Seasonality (°C)",
		"Max Temperature of Warmer Month (°C)", "Min Temperature of Coldest Month (°C)", "Temperature Annual Range (°C)",
		"Mean Temperature of Wettest Quarter (°C)", "Mean Temperature of Coldest Quarter (°C)", "Annual Precipitation (mm)",
		"Precipitation of Wettest Month (mm)", "Precipitation of Driest Month (mm)", "Precipitation Seasonality (mm)",
		"Precipitation of Wettest Quarter (mm)", "Precipitation of Driest Quarter (mm)", "Precipitation of Warmest Quarter (mm)",
		"Precipitation of Coldest Quarter (mm)")

for(i in 1:ncol(bioc)){
	curname<-colnames(bioc[i])
	png(paste(curname, "density.png", sep="."))	
		d<-density(na.omit(bioc[,i]))
		plot(d, main=curname, xlab=qnt[i], col="black")	
		points(y=rep(0, length(bioc[,i])), x=bioc[,i], col=colz, pch="|")
		legend("topright", legend=coldf[,1], col=coldf[,2], pch=20, cex=0.8, bty = "n")
	dev.off()
}

#make a complete picture
# pdf("bioclim.density.plots.pdf")	
	# par(mfrow=c(4,3),	 # 2x2 layout
		# oma = c(5, 5, 0, 0), # two rows of text at the outer left and bottom margin
		# mar = c(5, 5, 0, 0), # space for one row of text at ticks and to separate plots
		# mgp = c(2, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
		# xpd = NA)         	
	# for(i in 1:10){
		# curname<-colnames(bioc[i])
		# d<-density(na.omit(bioc[,i]))
		# plot(d, main=curname, xlab=qnt[i], col="black")	
		# points(y=rep(0, length(bioc[,i])), x=bioc[,i], col=colz, pch="|")
		# legend("topright", legend=coldf[,1], col=coldf[,2], pch=20, cex=0.8, bty = "n")
	# }
# dev.off()



############################
#check for differential distribution of traits on the basis of bioclim

phebio<-cor(qutr, bioc, use="complete.obs")

pdf("heatmap.quantitative.phenotypes.bioc.pdf")
	heatmap(phebio)
dev.off()


#save all
save.image(file="data.analysis.step.1.Rdata")


#############################
###restart from here
options(stringsAsFactors = FALSE)
maindir<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/TEFF/TEFF.distribution/analysis/1.data.analysis"
setwd(maindir)
#library(MuMIn)

load("data.analysis.step.1.Rdata")

#start removing all NAs that may be harmful down the road
dim(alldat)
alldatr<-na.omit(alldat)


#set a global model for a target phenotype
#start by setting formulas to be used in a list
flist<-list()
for (i in 3:19){
	flist[[i-2]]<-as.formula(paste(colnames(alldatr)[i], "~", paste(colnames(alldatr)[20:ncol(alldatr)], collapse = "+"),sep = ""))
}
names(flist)<-colnames(alldatr)[3:19]

#set a list of models and run it!
modlist<-list()
for (i in 1:length(flist)){
	modlist[[i]] <- lm(flist[[i]], data = alldatr)
	#step.model <- stepAIC(modlist[[i]], direction = "both", trace = FALSE)
}
names(modlist)<-names(flist)

#extract only significant hits from each model
siglist<-list()
for (i in 1:length(modlist)){
	tmpmod<-summary(modlist[[i]])
	siglist[[i]]<-data.frame(tmpmod$coef[tmpmod$coef[,4] <= .05, 4])
	siglist[[i]]<-cbind(siglist[[i]], rep(tmpmod$r.squared, nrow(siglist[[i]])))
}
names(siglist)<-names(modlist)

pout<-do.call("rbind", siglist)
pout<-cbind(rownames(pout), rownames(pout), pout)
pout[,1]<-sub("\\..*$", "", pout[,1])
pout[,2]<-sub("^.*\\.", "", pout[,2])
colnames(pout)<-c("trait", "bioclim", "pval", "R^2")
rownames(pout)<-NULL


