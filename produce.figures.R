#setting main directory
wd<-("C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/TEFF/TEFF.distribution/analysis/3.paper.analysis")
setwd(wd)

library(ggpubr)
library(treemapify)
library(vcd)
library(ggmosaic)

df<-read.delim("../../publication/figures/Table.S2.pheno.txt")
#fix na in SC and NV
df[which(is.na(df[, "SC"])), "SC"]<-3
df[which(df[,"Ptype"]=="l"), "Ptype"]<-NA

#get figures
mod<-lm(df[,"Altitude"] ~ as.factor(df[,"SC"]))
summary(mod)
cols<-colorRampPalette(c("brown4", "bisque3"))

boxplot(df[,"Altitude"] ~ as.factor(df[,"SC"]), col=cols(5),
          ylab="masl", xlab="Seed color")

my_comparisons <- list( c("1", "5"), c("1", "3"), c("3", "5") )

p1<-ggboxplot(df, x="SC", y = "Altitude", color="SC", add = "jitter", legend = "none", ylim=c(1000, 3500)) +
  geom_hline(yintercept = mean(df$Altitude, na.rm=T), linetype = 2) +
  stat_compare_means(method = "anova", label.y = 3500)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", hide.ns = TRUE, label.y = 3200) +
  color_palette(cols(5)) # + theme_gray()

pdf("seed.color.pdf")
  print(p1)
dev.off()


dfred<-df[df[,"Region"] %in% c("Amara", "Oromiya", "SNNPR", "Tigray"),]
#work with PCAs

#get bioclim
biovar<-c("^Alt", "^ph_", "^bio.*2000$")
hit<-grep(paste0(biovar, collapse="|"), colnames(df))
bc<-df[,hit]
rownames(bc)<-df[,"Accession"]
head(bc)
dim(bc)

bc.pca <- prcomp(na.omit(bc), center = TRUE,scale. = TRUE)
bcvar<-bc.pca$x[,1:10]

dfbc<-df[which(df[,"Accession"] %in% rownames(bcvar)),]
stopifnot(all(dfbc[,"Accession"] == rownames(bcvar)))

df2<-cbind(dfbc, bcvar)
df3<-na.omit(df2)

p2<-ggboxplot(df3, x="Ptype", y = "PC2", color="Ptype", add = "jitter", legend = "none", ylim=c(-8,12)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_compare_means(method = "anova", label.y = 12)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", hide.ns = TRUE, label.y = 9) 

pdf("panicle.type.pdf")
  print(p2)
dev.off()


dfred<-df2[df2[,"Region"] %in% c("Amara", "Oromiya", "SNNPR", "Tigray"),]
dfred<-droplevels(dfred)

ggplot(dfred) +
  aes(x = Hdatec, fill = factor(PCHdate)) +
  facet_wrap(~ Region) +
  rotate_x_text(angle = 45) +
  geom_bar(position = "fill")

mosaicplot(~ Region + PCHdate, data = dfred, shade = TRUE, main="", las=2)
vcd::mosaic(~ Region + PCHdate, data = dfred, shade = TRUE, main="", las=2)

pdf("PCHdate.by.region.pdf")  
  vcd::assoc(~ PCHdate +  Region, data = dfred, shade = TRUE, main="",
                labeling= labeling_border(rot_labels = c(0,90,90,0)))
dev.off()
