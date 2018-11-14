###############################################
# Scripts from Paraphelidium transcriptome paper 
# (Torruella et al. Nat Comm Biol 2018)
# Xavier Grau-Bové (xavier.graubove@gmail.com)
# https://github.com/xgrau/paraphelidium2018
###############################################

# libraries
library(reshape2)
library(tidyr)
library(ape)
library(gplots)
library(viridis)

#### Input ####

input   = c("primary_transposed.txt") # assumes a wide-style matrix with headers, rows=genes/OGs & cols=sps 
outcode = "out"                       # code for output (final format: inputname_outcode_XXXX)

# Colors
hop     = colorRampPalette(interpolate="l",c("white","deepskyblue","dodgerblue3","dodgerblue4"))
hopviri = viridis(51,option = "A",begin=0.25,end=0.99,direction = -1)

# Load matrix
mi = read.table(input,header=T,sep="\t")

# Data formatting
# Set gene names as row names and remove gene names column
rownames(mi) = as.vector(mi[,1])
mi = mi[,-1]

mi[is.na(mi)] = 0                     # change missing to 0 (tweak if appropriate)
mi = as.matrix(mi)                    # format matrix
mi = mi[rowSums(mi)>0,colSums(mi)>0]  # retain genes that are present in >=1 sps; and sps that have >=1 genes
write.table(mi, file = "primary_transposed.txt", sep = "\t")
mi = t(mi)                            # transposa la matriu per a que row=sp, col=gen
mi[mi>0] = 1                          # Convert data to abs/pres (0/1)



#### Analysis ####

# PCA
pic = prcomp(mi)
pic$variancefraction = pic$sdev^2/sum(pic$sdev^2)

# Distances & clustering
# genes: euclidean distance + Ward clustering
mic.dist = dist(t(mi),method = "euclidean")   
mic.clus = hclust(mic.dist,method="ward.D2")
# genes: 1-pearson correlation + Ward clustering
mir.cor  = cor(t(mi),method = 'pe')
mir.dist = as.dist(1-mir.cor)
mir.clus = hclust(mir.dist,method="ward.D2")

# PCoA
pio = pcoa(mir.dist)
pio12 = as.data.frame(pio$vectors[,c(1,2)])
pio12$sp = rownames(pio12)


#### Plots ####

# open pdf: PCA & PCoA
pdf(paste(input,outcode,"reddim-colors.pdf",sep="."),height=6.5,width =6)

# PCA
# Plot PCA (components 1&2)
plot(pic$x[,c(1,2)],col="red",main="PCA")
text(pic$x[,c(1,2)],labels=factor(rownames(mi)),pos=3,col="black")
# Plot PCA (variance explained)
plot(pic$variancefraction,type="b",col="red",main="PCA variance")
text(pic$variancefraction,labels=signif(pic$variancefrac,digits=3),pos=3,col="black")


# Plot PCoA
# Color species per group
pio12$grup = ifelse(pio12$sp %in% c("Amp_sp","Ant_lo","Enc_cu","Enc_in","Ent_bi","Mit_da","Nem_pa","Nos_ce","Par_sa","Roz_al","Par_tr"),"Opisthosporidia",
                     ifelse(pio12$sp %in% c("Spi_pu","Gan_pr","Bat_de","All_ma","Cat_an","Bla_br","Mor_ve","Cop_ci","Ust_ma","Sch_po"),"Fungi",
                            ifelse(pio12$sp %in% c("Fon_al","Par_at"),"Nuclearid",
                                   ifelse(pio12$sp %in% c("Sal_ro","Mon_br","Cap_ow","Chr_pe","Cor_li","Ich_ho","Sph_ar"),"Holozoa",
                                          ifelse(pio12$sp %in% c("The_tr"),"Apusomonad",
                                                 ifelse(pio12$sp %in% c("Aca_ca","Dic_pu","Pol_pa","Ent_hi"),"Amoebozoa",
                                                        ifelse(pio12$sp %in% c("Tox_go","Pla_fa","Try_br","Lei_ma","Nae_gr","Phy_in"),"Other","Uncl")))))))
dpio12$grup=as.factor(pio12$grup)

# PCoA (components 1&2)
plot(x=pio12$Axis.1,y=pio12$Axis.2,col=pio12$grup,main="PCoA",
     xlab=paste("Axis 1",signif(pio$values$Relative_eig[1]*100,digits=3),"% variance"),
     ylab=paste("Axis 2",signif(pio$values$Relative_eig[2]*100,digits=3),"% variance"))
text(pio$vectors[,c(1,2)],labels=factor(rownames(mi)),pos=3,col="gray")
legend("bottomleft",legend = levels(pio12$grup),col=1:nlevels(pio12$grup),pch=1,cex=1)

# Plot PCoA (variance explained)
plot(pio$values$Relative_eig,col="red",main="PCoA relative eigenvalues",type="b",xlab="Axis rank",ylab="Fraction variance explained")
text(pio$values$Relative_eig,labels=signif(pio$values$Relative_eig,digits=3),pos=3,col="black")
lines(pio$values$Broken_stick[pio$values$Broken_stick>.01],
      col="gray",main="PCoA variance broken stick",type="b")

dev.off()



# open pdf: heatmaps
# heatmap: sps-to-sps 
pdf(paste(input,outcode,"heatmaps.pdf",sep="."),height=8,width=8)
mir.cor[mir.cor<0] = NA 
heatmap.2(mir.cor,
          col =hopviri,
          dendrogram = "both",
          symkey=F,
          trace="none",
          scale="none",
          useRaster=T,
          keysize = 1,
          main = "pairwise sample correlation",
          symm=F,
          Colv=as.dendrogram(mir.clus),
          Rowv=as.dendrogram(mir.clus),
          margins=c(10,10),
          na.color = "grey90",
          key.title = "corr (positive only")

# heatmap: gene counts per species 
heatmap.2(mi,
          col =hop(2),
          dendrogram = "both",
          symkey=F,
          trace="none",
          scale="none",
          useRaster=T,
          keysize = 1,
          main = paste(ncol(mi),"variables in",nrow(mi),"samples"),
          symm=F,
          Colv=as.dendrogram(mic.clus),
          Rowv=as.dendrogram(mir.clus),
          margins=c(10,10),
          na.color = "grey90",
          key.title = "gene counts")


dev.off()




