#setting main directory
wd<-("C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/TEFF/TEFF.distribution/analysis/3.paper.analysis")
setwd(wd)

#keep the correct packages
library(ggplot2)
library(ggfortify) 
library(sp)
library(corrplot)
library(vegan)
library(CCA)
library(mclust)
library(GGally)
library(psych)
library(patchwork)
library(RColorBrewer)

##1. GETTING THE CORRECT DATASET TYPE
#load the table with all data (biocl and pheno)
tab<-read.delim("../../data/Bioclimatic.final.final.txt")
dim(tab)

#keep only wanted informations
tokeep<-c("Accession", "Region", "Sub_Region", "Lat","Long", "Altitude", "ph_ho_60cm",
          "ph_kcl_60cm","Edate","Hdate","Mdate","GfP","PL","NV",
          "NVc","Plc","Gfpc", "PCMdate", "SBP","Ptype", "SC", "PCHdate","Edatec", 
          "Hdatec","Mdatec","Bio_1" ,  "Bio_2" ,  "Bio_3",  
          "Bio_4"  , "Bio_5"  , "Bio_6",  "Bio_7"  ,
          "Bio_8"  , "Bio_9"  , "Bio_10" , "Bio_11" ,
          "Bio_12" , "Bio_13" , "Bio_14" , "Bio_15" ,
          "Bio_16",  "Bio_17" , "Bio_18",  "Bio_19" )
tb<-tab[, tokeep]
tab<-na.omit(tb)

dim(tb)

#get bioclim
biovar<-c("^Alt", "^ph_", "^Bio")
hit<-grep(paste0(biovar, collapse="|"), colnames(tab))
bc<-tab[,hit]
rownames(bc)<-tab[,"Accession"]
head(bc)
dim(bc)

#get phenos
qtvar<-c("^Edate$","^Hdate$","^Mdate$","^GfP$","^PL$")
hit<-grep(paste0(qtvar, collapse="|"), colnames(tab))
qt<-tab[,hit]
rownames(qt)<-tab[,"Accession"]
head(qt)
dim(qt)

#Check that everything is allrigth
stopifnot(all(rownames(bc)==rownames(qt)))

#get correlations done
forcor<-cbind(bc, qt)

#cr<-cor(csdonly[,which(colnames(csdonly) %in% forcor)], method ="spearman", use="pairwise.complete")
cr<-corr.test(forcor, method ="spearman", use="pairwise.complete")

pdf("correlations.pdf")
  corrplot(cr$r, p.mat = cr$p, insig ="blank", "shade", tl.col="black", tl.cex = 0.7)
dev.off()


#set colors
colz<-rev(brewer.pal(n = 6, name = "Set2"))

#perform PCA on bioclim
bc.pca <- prcomp(bc, center = TRUE,scale. = TRUE)
plt<- autoplot(bc.pca, data = tab, colour = "Region", title = "Climate", alpha = 0.8,  size = 1, loadings=TRUE,
               loadings.label = TRUE, loadings.colour ='gray',
               loadings.label.size = 3, loadings.label.color = "black") 
plt1<-plt + scale_color_manual(values=colz) +	theme_bw() + labs(tag="A)")
plt1

pdf("bioclim.pca.pdf", height=4, width=6)
  plt1
dev.off()

#correlate BC PCA with original BC
pcscores<-bc.pca$x[,1:3]
crbc<-cor(pcscores, bc, "complete.obs", method="spearman")
pdf("correlation.bioclim.PCA.pdf")
  corrplot(crbc, "color", tl.col="gray10", tl.cex = 0.7)
dev.off()

#make correlations b/w bio PC and quantitative traits
tocor<-cbind(pcscores, qt)

#get significances
pval<-list()
for(i in 1:3){
  pvalpc<-c()
  for(j in 1:5){
    ct<-cor.test(tocor[,i], tocor[,j+3])
    pvalpc[j]<-ct$p.value
  }
  names(pvalpc)<-names(tocor)[4:8]
  pval[[i]]<-pvalpc
}
names(pval)<-c("PC1", "PC2", "PC3")
lapply(pval, function(x) x<5e-2)

crpc<-ggpairs(tocor)
pdf("correlation.PC.bioclim.pdf")
  crpc
dev.off()

#pca on phenotypes now
pcaqt<-prcomp(qt, center = TRUE,scale. = TRUE)
plt<- autoplot(pcaqt, data = tab, colour = "Region", title = "Phenotypes", alpha = 0.6,  size = 1, loadings=TRUE,
               loadings.label = TRUE, loadings.colour ='gray',
               loadings.label.size = 3, loadings.label.color = "black") 
plt2<-plt + scale_color_manual(values=colz) +	theme_bw()  + labs(tag="B)")
plt2

pdf("phenotype.pca.pdf", height=4, width=6)
  plt2
dev.off()

#correlate trait PCA with original traits
pcscores<-pcaqt$x[,1:3]
crbc<-cor(pcscores, qt, "complete.obs", method="spearman")
pdf("correlation.trait.PCA.pdf")
  corrplot(crbc, "color", tl.col="gray10", tl.cex = 0.7)
dev.off()

#create figure two for the manuscript
combined <- plt1 / plt2  
combined + plot_layout(guides = 'collect')
  
ggsave("fig.2.pdf", width = 6, height = 9)
ggsave("fig.2.png", width = 6, height = 9)


#obtain tef ecological niche
#we assume that the tails of current teff distribution represent marginal environemnts
hist(bc[,"bio1"])
summary(bc[,"bio8"])
summary(bc[,"bio12"])
summary(qt[,"Mdate"])

quantile(bc[,"bio1"], .90)


#get qualitative info
qqvar<-c("NV", "Plc", "Gfpc", "PCMdate", "SBP", "Ptype", "SC", "PCHdate", "Edatec", "Hdatec", "Mdatec")
hit<-grep(paste0(qqvar, collapse="|"), colnames(tab))
qq<-tab[,hit]
rownames(qq)<-tab[,"Accession"]
head(qq)

#do some test
boxplot(bc[,"Altitude"] ~ as.factor(qq[,"SC"]))
anova<-aov(bc[,"Altitude"] ~ as.factor(qq[,"SC"]))
summary(anova)

boxplot(bc[,"Altitude"] ~ as.factor(qq[,"SC"]))
mod1<-lm(bc[,"Altitude"] ~ as.factor(qq[,"SC"]))
summary(mod1)
#plot(mod1)

plot(pcabcs[,"PC1"] ~ qt[,"Hdate"])
mod1<-lm(pcabcs[,"PC1"] ~ qt[,"Hdate"])
summary(mod1)


#save relevant data
save(bc, qt, qq, tab, bc.pca, pcaqt, file="diversity.data.Rdata")




#############NOW ON THE FULL DATASET###################
full<-read.delim("../../data/full.collection.phenotypes.passport.txt", header=T)
full<-na.omit(full)
#fix region names
full[full=="Oromiya"]<-"Oromia"
full[full=="Amara"]<-"Amhara"
full[full=="Tigay"]<-"Tigray"
full[full=="Tigeray"]<-"Tigray"
full[full=="Benishangul & Gumuz"]<-"Benish. Gumuz"
full[full=="SNNPR"]<-"SNNP"
unique(full[,"Region"])

#clean the dataset
rownames(full)<-full[,1]
fpheno<-full[,c(11,13,15,17,20)]


#pca on phenotypes now
pcaqtf<-prcomp(fpheno, center = TRUE,scale. = TRUE)
plt<- autoplot(pcaqtf, data = full, colour = "Region", alpha = 0.5, loadings=TRUE,
               loadings.label = TRUE, loadings.colour ='gray',
               loadings.label.size = 3, loadings.label.color = "black") 

pdf("phenotype.pca.full.dataset.pdf", height=4, width=6)
  plt	+	theme_bw()
dev.off()

#correlate trait PCA with original traits
pcscores<-pcaqt$x[,1:3]
crbc<-cor(pcscores, qt, "complete.obs", method="spearman")
pdf("correlation.trait.PCA.full.dataset.pdf")
  corrplot(crbc, "color", tl.col="gray10", tl.cex = 0.7)
dev.off()








#####################################

  pcabc<-prcomp(bc, scale=T)
  
  autoplot(pcabc, data=qt, loadings=TRUE,
           loadings.label = TRUE, loadings.colour ='blue',
           loadings.label.size = 3, col=as.numeric(tab[,"Region"]))
  
  
  autoplot(pcabc, data=qt, loadings=TRUE,
           loadings.label = TRUE, loadings.colour ='blue',
           loadings.label.size = 3, col=as.numeric(tab[,"SC"]))
  
  ggbiplot(pcabc,ellipse=TRUE,choices=c(3,4),
           labels=rownames(pcabc), groups=Region)
  
  summary(pcabc)
  pcabcs<-pcabc$x[,1:3]
  head(pcabcs)
  cr<-cor(pcabcs,bc)
  corrplot(cr)
  corrplot(cr,"pie")

#create a tab with biocl e PCbio
  biopc<-merge(bc,pcabcs)
  head(biopc)

#perform PCA on pheno
pcaqt<-prcomp(qt, scale=T)
autoplot(pcaqt, data=qt, loadings=TRUE,  
         loadings.label = TRUE, loadings.colour ='blue', 
         loadings.label.size = 3,col=as.numeric(tab[,"SC"]))  
summary(pcaqt)
pcaqts<-pcaqt$x[,1:3]
cr<-cor(pcscores,qt)
corrplot(cr)
corrplot(cr,"pie")

#check correlations bw traits
cr<-cor(qt,pcabcs)
corrplot(cr,"pie")

#do cca
cca<-cc(bc, qt) 
plt.cc(cca, type="v", var.label=T)



#####################
#CCA MAKES NOT SENSE; biclim variabity EXPLAINS ONLY 3% of the pneotypic variance

#do a CCA setting up a formula
form<-as.formula(paste("qt ~ ", paste(colnames(bc), collapse=" +")))
cc<-cca(form, data=bc) 

#do a CCA on bioclim PCs
#cc<-cca(qt ~ PC1+PC2+PC3, data=data.frame(bc.pca$x)) 
#plot(cc, scaling=3)

#check redundancy across variables
vif.cca(cc)

#drop redundant variables
drop<-drop1(cc, test="permutation")
#keep only significant variables
sign<-which(drop$`Pr(>F)`<= 0.05)
tokeep<-rownames(drop)[sign]

#re-run the model with significant variables
form1<-as.formula(paste("qt ~ ", paste(tokeep, collapse=" +")))
cc<-cca(form1, data=bc) 
vif.cca(cc)

#get TOTAL variance explained
cc$CCA$tot.chi/cc$tot.chi

#perform some plotting
barplot(cc$CA$eig/cc$tot.chi, names.arg = 1:cc$CA$rank, cex.names = 0.5, ylab="Proportion of variance explained", xlab="CCA axis")
plot(cc, scaling=3, col?)
#significance by variables
anova(cc, by="margin") 
#significance of the total model
anova(cc)

