##################
##### WGCNA ######
##################
library(permute)
library(vegan)
library(MASS)
library(Hmisc)
library(ggbiplot)
library(nlme)
library(ggplot2)
library(ggrepel)
library(psych)
library(xlsx)
library(gridExtra)
library(phyloseq)
library(WGCNA)
library(survival)
require(DESeq2)

#############################################################################
##### Phyloseq, but reducing data further than other diversity analyses #####
#############################################################################

NUT <- read.table("GLCW_metadata_clean.txt", header = TRUE, row.names = 1, sep ="\t")
OTU <- read.table("OTU_table_singdoubrem.txt", header = TRUE, row.names = 1, sep="\t")
TAX <- read.table("OTU_taxonomy_singdoubrem_cols.txt", header = TRUE, row.names = 1, sep = "\t")

rownames(OTU) <- paste0("OTU", 1:nrow(OTU))
rownames(OTU)

TAX <- as.matrix(TAX, rownames.force = NA)
rownames(TAX) <- paste0("OTU", 1:nrow(TAX))
rownames(TAX)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)

physeq = phyloseq(OTU,TAX)

META = sample_data(NUT)
rownames(META) <- sample_names(physeq)

META = sample_data(META)

# get all the data in a phyloseq instance, or whatever
ALL = phyloseq(OTU,TAX,META)
ALL

### Prepping data to remove rare or erronious OTUs #####

# We need to de-noise the data by plotting the number of reads on a curve and look for the inflection point

at.least.n.in.m <- function(x, n, m){
  all(x[x>0]>=n)&length(x[x>0])>=m
}
counts<- rep(0,10)
for (i in 1:length(counts)){
  rows.to.keep<- apply(otu_table(ALL, taxa_are_rows = TRUE), 1,at.least.n.in.m, n=i, m=2)
  counts[i]<-sum(rows.to.keep)
}

plot(1:10, counts, xlab= 'Min sequences in 2 samples', ylab= 'Number of taxa remaining')

### Inflection point was 2 ###

# Filter taxa that arent seen more than twice in greater than 10% of the data.
CUT2<-filter_taxa(ALL, function(x) sum(x > 2) > (0.1*length(x)), TRUE)
CUT2
write.csv(tax_table(CUT2), "GLCW_seqnames.csv")
write.csv(otu_table(CUT2), "GLCW_seqOTU.csv")

### Deseq2 Normalization
### convert the counts to integer mode and group data to normalize together(e.g. region)
dds_WGCNA <- phyloseq_to_deseq2(CUT2, ~ region)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans_WGCNA = apply(counts(dds_WGCNA), 1, gm_mean)
dds_WGCNA = estimateSizeFactors(dds_WGCNA, geoMeans=geoMeans_WGCNA)
dds_WGCNA = estimateDispersions(dds_WGCNA)

### variance stabilizing transformation
#Make it a OTU again for phyloseq (no tax or sample data)
vst_WGCNA <- varianceStabilizingTransformation(dds_WGCNA, blind=FALSE)
vstMat_WGCNA <- assay(vst_WGCNA)
vstMat_WGCNA[vstMat_WGCNA<0]<-0
vst.otu.WGCNA <- otu_table(vstMat_WGCNA, taxa_are_rows=TRUE)
write.csv(otu_table(vst.otu.WGCNA), "vst.otu.WGCNA.csv")

######################
#### BEGIN WGCNA #####
######################

#### Maybe try this next time... ####
vst <- read.csv("vst.otu.WGCNA.csv")
dim(vst)
names(vst)
otu_WGCNA <- as.data.frame(t(vst))
names(otu_WGCNA) = vst$X
rownames(otu_WGCNA) = names(vst)
dim(otu_WGCNA)
names(otu_WGCNA)


NMDS_META <- read.table("GLCW_NMDS_metadata_clean.txt", header = T, row.names = 1)


# Because we are using RStudio, we have to disable threads
# As consequence of this, maybe it would be better to do this step in regular ol' R
disableWGCNAThreads()
options(stringsAsFactors = FALSE)
## Identify beta to ensure scale free topology
powers = c(seq(4,10,by=1), seq(12,20, by=2))

memory.size(max = TRUE)
memory.limit(size = 15000)
memory.limit()

pst <- pickSoftThreshold(otu_WGCNA, powerVector=powers, blockSize = 7562, verbose=2)

# Plot the results of soft thresholds if you wish:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(pst$fitIndices[,1], pst$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(pst$fitIndices[,1], pst$fitIndices[,5], labels=powers, cex=cex1,col="red")

#check the powerEstimate to choose a threshold
pst

# Create adjacency matrix by raising OTU matrix by beta and identify subnetworks (modules)
otu_WGCNA2 <- as.matrix(otu_WGCNA[2:75,1:7562])
mode(otu_WGCNA2)
class(otu_WGCNA2) <- "numeric"
mode(otu_WGCNA2)

# Check that the network ensures scale-free topology at that power
# R should be close to 1 (R > 0.8, I believe), should see a straight line.
##### scaleFreePlot #####
# here we define the adjacency matrix using soft thresholding with beta=4
ADJ1=abs(cor(otu_WGCNA2,use="p"))^4
# When you have relatively few genes (<5000) use the following code
#k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k=softConnectivity(datE=otu_WGCNA2,power=4)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
scaleFreeFitIndex(k)

#R^2 of 0.87, this suggests we meet the assumption of scale-free topol.

# power of 4 chosen based on powerEstimate from 'pst'
net = blockwiseModules(otu_WGCNA2, power=4, minModuleSize=30, maxBlockSize = 7562,
                       corType = "pearson", saveTOMs = TRUE, 
                       saveTOMFileBase = "blockwiseTOM", pamStage=FALSE, verbose=5)
# Plot the dendrogram
moduleLabels = net$colors
moduleColors = net$colors
MEs = net$MEs
geneTree = net$dendrograms[[1]]
pdf("plotDendro.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
# Save data

# Identify Eigenvalues for subnetworks by sample
nPop = ncol(otu_WGCNA2)
nSamples = nrow(otu_WGCNA2)
MEsO = moduleEigengenes(otu_WGCNA2, moduleColors)$eigengenes
MEs = orderMEs(MEsO)
save(MEs, moduleLabels, moduleColors, geneTree,file = "Module-networkConstruction-auto.RData")
# Save data
write.csv(file="Module_eigen_values.csv",MEs)
write.csv(file="Module_composition.csv",net$colors)

##Correlate Eigenvalues to metadata and create heatmap
# Optional rename of OM to NUTR
colnames(NMDS_META)[10] <- "NUTR"

moduleTraitCor = cor(MEs, NMDS_META[,c(4:12)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), " (",signif(moduleTraitPvalue, 1), ")"
                   , sep = "")
dim(textMatrix) = dim(moduleTraitCor)
pdf("Correlation.pdf",width=12,height=8)
par(mar = c(12, 12, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(NMDS_META[,c(4:12)]),yLabels = names(MEs)
               ,ySymbols = names(MEs),colorLabels = FALSE,colors = blueWhiteRed(50)
               ,textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,
               cex.lab = 0.5,zlim = c(-1,1),main = paste("Module-trait relationships"))
dev.off()

##example heatmap for presentations and such
#pdf("Correlation_example.pdf",width=12,height=8)
#par(mar = c(12, 12, 3, 3))
#labeledHeatmap(Matrix = moduleTraitCor[101:112,1:8],xLabels = names(NMDS_META[,3:10]),yLabels = names(MEs)
#               ,ySymbols = names(MEs),colorLabels = FALSE,colors = blueWhiteRed(50)
#               ,textMatrix = textMatrix[101:112,1:8],setStdMargins = FALSE,cex.text = 0.5,
#               cex.lab = 0.5,zlim = c(-1,1),main = paste("Module-trait relationships"))
#dev.off()


## Now make a plot for specific module <-> trait (metadata component) pairings
# This allows us to explore the structure of submodule OTU correlations with a given metadata component
# Here we will use "env variable" as trait and "color" as module
# First get the links between all modules and this trait

###OC-orange
parameter<-"NUTR"
weight <- as.data.frame(NMDS_META[,parameter])
names(weight) = "NUTR"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(otu_WGCNA2, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(otu_WGCNA2, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
# Then look at the specific module of interest
module<-"orange"
column = match(module, modNames);
moduleGenes = moduleColors==module
pdf(paste(module,"-vs-",parameter,".pdf"))
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]),xlab = paste("Module Membership in", module, "module"),ylab = paste("Population significance for ",parameter),main = paste("Module membership vs. population significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

###S-turquoise
parameter<-"S"
weight <- as.data.frame(NMDS_META[,parameter])
names(weight) = "S"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(otu_WGCNA2, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(otu_WGCNA2, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
# Then look at the specific module of interest
module<-"turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module
pdf(paste(module,"-vs-",parameter,".pdf"))
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]),xlab = paste("Module Membership in", module, "module"),ylab = paste("Population significance for ",parameter),main = paste("Module membership vs. population significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

###CNrat-pink
parameter<-"CNrat"
weight <- as.data.frame(NMDS_META[,parameter])
names(weight) = "CNrat"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(otu_WGCNA2, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(otu_WGCNA2, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
# Then look at the specific module of interest
module<-"pink"
column = match(module, modNames);
moduleGenes = moduleColors==module
pdf(paste(module,"-vs-",parameter,".pdf"))
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]),xlab = paste("Module Membership in", module, "module"),ylab = paste("Population significance for ",parameter),main = paste("Module membership vs. population significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

###NO3-orange
parameter<-"NO3"
weight <- as.data.frame(NMDS_META[,parameter])
names(weight) = "NO3"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(otu_WGCNA2, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(otu_WGCNA2, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
# Then look at the specific module of interest
module<-"orange"
column = match(module, modNames);
moduleGenes = moduleColors==module
pdf(paste(module,"-vs-",parameter,".pdf"))
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]),xlab = paste("Module Membership in", module, "module"),ylab = paste("Population significance for ",parameter),main = paste("Module membership vs. population significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

##Incorporate OTU taxonomy info to identify specific OTUs within a module, coalesce correlative data by OTU for later use in network analysis
annot2 = read.csv("GLCW_seqnames.csv", colClasses = "character", header = TRUE)
colnames(annot2)[1] <- "OTU"
dim(annot2)
names(annot2)
probes2 = names(otu_WGCNA)
probes2annot2 = match(probes2, annot2$OTU)
# The following is the number or probes without annotation:
sum(is.na(probes2annot2))
# Should return 0.
# Create the starting data frame
geneInfo0 = data.frame(otu_orig = probes2,
                       OTU = annot2$OTU[probes2annot2],
                       Domain = annot2$Domain[probes2annot2],
                       Phylum = annot2$Phylum[probes2annot2],
                       Class = annot2$Class[probes2annot2],
                       Order = annot2$Order[probes2annot2],
                       Family = annot2$Family[probes2annot2],
                       Genus = annot2$Genus[probes2annot2],
                       names = annot2$OTU[probes2annot2],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance according to metadata component
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.NUTR)); # <- change geneInfo$'X'
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = paste("OTUInfo_",module))

# Get the similarity between eigenvalues and weight
MET = orderMEs(cbind(MEs, weight))
pdf("Adjacency_trait_ME_NUTR.pdf") # <- This changes with env. variable
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene dendrogram with ",parameter), marDendro = c(0,4,2,0),plotHeatmaps = FALSE)
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene adjacency heatmap with ",parameter), marHeatmap = c(3,4,2,2),plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

################################################################################
#           Moving to machine learning PLS and VIP scores                      #
################################################################################

library("pls")
source("VIP.R")
# metadata component (e.g. S, designated "weight" here as above) is the same as before, we just replicate the row names for pls
#the below r2 threshold changes with user input, good correlation above "0.x"
#this is important, because VIP scores don't really mean anything without a good correlation
th_r2<-0.3
subnetwork<-otu_WGCNA[2:75,moduleGenes]
row.names(weight)<-row.names(subnetwork)
weight <- as.matrix(weight)
subnetwork <- as.matrix(subnetwork)
class(weight) <- "numeric"
class(subnetwork) <- "numeric"
pls_result<-plsr(weight ~ subnetwork, validation="LOO",method="oscorespls")

r2_vector<-R2(pls_result)
max<-0
max_comp<--1
for (j in 1:length(r2_vector$val)){
  if(r2_vector$val[j]>th_r2){         # We will only look at the PLS if the correlation is better than 0.3
    if(r2_vector$val[j]>max){
      max<-r2_vector$val[j]
      max_comp<-r2_vector$comp[j]
    }
  }
}
print(paste(" the max r2 is ",max," corresponding to comp ",max_comp,sep="",".pdf"))

if(max==0){
   		print ("No good correlation, we stop here")
} else{
  print("Good correlation, we check the VIP!")
}
  # Checking the VIP
  output<-paste("VIP_values_with_",parameter,sep="")
  vip_result<-VIP(pls_result)
  vip_components<-sort(vip_result[max_comp,],decreasing=TRUE)[1:69] # <- this value changes with subnetwork, check dim of 'subnetwork' file
  for (i in 1:69){ # <- this value changes as the line above
    cat(paste("Rank ",i," we have ",names(vip_components[i])," with a VIP of ",vip_components[i],"\n",sep=""),file=output,append=TRUE)
  }
  weight_2 <- as.data.frame(weight[!is.na(weight)])
  df<-data.frame(x=weight_2,y=pls_result$validation$pred[,,max_comp])
  colnames(df)<-c("x","y")
  pdf(paste("measured_vs_predicted_",module,"-vs-",parameter,".pdf"))
  ggplot(data=df) + geom_point(aes(x=x,y=y)) + geom_smooth(aes(x=x,y=y),method=lm) + xlab("Measured") + ylab("Predicted") + ggtitle(paste("Comparison of ",parameter," measured vs predicted for module ",module)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"))
  dev.off()
  # Establish the correlation between predicted and modeled
  # This is the data to report with the figure (R2, CI, signif, etc.)
  cor.test(df$x,df$y)

## Identify node centrality based on co-occurence data for each OTU in the module

TOM = TOMsimilarityFromExpr(otu_WGCNA2, power = 4, corType = "pearson");
# Select submodule of interest based on high correlation and signficance
module<-"orange"; # <- this changes with module color being currently explored
# Select module probes
probes = names(otu_WGCNA)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
write.csv(modTOM, paste("nodeconnections_",parameter,".csv"))

# Number of cells above 0.25 threshhold <-- this number is flexible and should change with your data
x<- as.data.frame(rowSums(modTOM > 0.25))
write.csv(x, paste("nodes_",parameter,".csv"))


# make scatter hive plots
# You will need to make the Nodeworksheet by combining OTUinfo table for submodule of interest, VIP socres, and Nodeconnections. See workflow for more details.

hive<- read.csv("NUTR_orange_nodeworksheet.csv", header=T) #(Figure 5-8)
hive$OTU<-factor(hive$OTU, levels = unique(hive$OTU))
p<- ggplot(hive, aes(x= hive$connectivity, y= hive$GS.NUTR)) +
  geom_point(aes(size = hive$VIP, colour = hive$Phylum, alpha = 0.5)) +
  scale_size_area(max_size= 10) +
  scale_color_brewer(type = qual, palette = "Set1", direction = -1) +
  theme_bw() +
  scale_alpha(guide=FALSE) +
  geom_text_repel(aes(label = hive$names), force = 3) +
  labs(x="Node centrality", y="Correlation to NUTR", color="Phylum",
       size = "VIP")
p

hive_S<- read.csv("S_turquoise_nodeworksheet.csv", header=T) #(Figure 5-8)
hive_S$OTU<-factor(hive_S$OTU, levels = unique(hive_S$OTU))
S<- ggplot(hive_S, aes(x= hive_S$connectivity, y= hive_S$GS.S)) +
  geom_point(aes(size = hive_S$VIP, colour = hive_S$Phylum, alpha = 0.5)) +
  scale_size_area(max_size= 10) +
  scale_color_brewer(palette = "Spectral", direction = -1) +
  theme_bw() +
  scale_alpha(guide=FALSE) +
  geom_text_repel(aes(label = hive_S$names), force = 3) +
  labs(x="Node centrality", y="Correlation to S", color="Phylum",
       size = "VIP")
S

hive_CN<- read.csv("CNrat_pink_nodeworksheet.csv", header=T) #(Figure 5-8)
hive_CN$OTU<-factor(hive_CN$OTU, levels = unique(hive_CN$OTU))
CN<- ggplot(hive_CN, aes(x= hive_CN$connectivity, y= hive_CN$GS.CNrat)) +
  geom_point(aes(size = hive_CN$VIP, colour = hive_CN$Phylum, alpha = 0.5)) +
  scale_size_area(max_size= 10) +
  scale_color_brewer(type = qual, palette = "Set1", direction = -1) +
  theme_bw() +
  scale_alpha(guide=FALSE) +
  geom_text_repel(aes(label = hive_CN$names), force = 3) +
  labs(x="Node centrality", y="Correlation to C:N", color="Phylum",
       size = "VIP")
CN

# PLS of NO3 doesn't correlate well (R2 < 0.3), thus, VIP scores are pretty meaningless,
#so we're just going to run correlations of individual taxa with NO3 instead

#import otu table from computer to return to non-deseq format
otu_tab_corr <- read.csv("vst.otu.WGCNA.csv", header = TRUE, row.names = 1)
corr_meta <- as.matrix(NMDS_META[,c(4:6,8,10)])

t_otu_tab_corr <- t(otu_tab_corr)

#Spearman's rank correlations
Spearcorr_GL <- rcorr(t_otu_tab_corr, corr_meta[,1:5], type = "spearman")
Spearcorr_GL
write.csv(Spearcorr_GL$r, "Spearcorr_GL_r.csv")
write.csv(Spearcorr_GL$P, "Spearcorr_GL_P.csv")

#edit spreadsheet in excel and import into R
# included in the below table are retained taxa with p < 0.001 correlation
# to Temp, pH, DO, DOC, or a combination of them. If p > 0.001, the value
# for that R is 0

# abundance graph
library(ggplot2)
library(scales)
require(reshape2)
require(plyr)

otu_bar_dat <- read.csv("NO3_corr_fig_dat_GLCW.csv", header = TRUE)

otu_bar_dat$Specificity <- factor(otu_bar_dat$Specificity, levels=unique(otu_bar_dat$Specificity))

ggplot(otu_bar_dat, aes(Specificity, direction)) +
  geom_col(size = .25, aes(fill = Phylum, col = corrnum)) + 
  facet_grid(env ~ .) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, size = 12),
        axis.title.x=element_blank()) +
  guides(col = FALSE) +
  ylab("Abundance") +
  scale_color_gradient(low = "white", high = "black") +
  geom_segment(aes(y = 0, yend = 0, x = 0, xend = 43))

###TESTS###
sampleTree = hclust(dist(otu_WGCNA2), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()



################################################
### Export to Cytoscape for Network Building ###
################################################

# Select modules
modules = "orange";
# Select module probes
inModule = is.finite(match(moduleColors, modules));
modProbes = probes2[inModule];
modGenes = annot2$gene_symbol[match(modProbes, annot2$OTU)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot2$OTU, annot2$OTU))
