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
library(multcomp)
require(DESeq2)
require(ggmap)
require(tidyverse)
library(RColorBrewer)

memory.size(max = TRUE)
memory.limit(size = 15000)
memory.limit()

setwd("~/Grad School/Central Michigan/Dissertation/Great Lakes Coastal Wetlands/Stats")

##########################
### Build a region map ###
##########################

#Here create a text file with all the coords needs a column of lat and long in decimal degrees
Sample_coords<- read.csv("GLCW_coordinates_region.csv", header=TRUE, row.names = 1)
map1<-get_googlemap(center=c(-85,44), zoom=7, size=c(420,640), maptype="satellite")
MapMiLakes<-ggmap(map1)+
  geom_point(data=Sample_coords, aes(x=long, y=lat,
                                     colour = rownames(Sample_coords)),
             size=5, show.legend = TRUE)+
  scale_color_brewer(type = "qual", palette = "Set1") +
  ylab("Latitude")+
  xlab("Longitude")+
  theme(legend.direction = 'vertical', 
        legend.position = 'right')+
  guides(colour=guide_legend(title="Region"))
MapMiLakes
ggsave("GLCW_Map.pdf", MapMiLakes)

GLstates <- map_data("state", region = c("arkansas", "new york", "pennsylvania", "delaware"))
OH <- map_data("state", region = "ohio")
OH

westMI <- subset(MI, subregion %in% "south")

ggplot() + geom_polygon(data = GLstates, aes(x=long, y = lat, group = group)) + 
  coord_fixed(1.2)

#use these colors throughout the rest of the figures to correspond to each region
g <- ggplot_build(MapMiLakes)
unique(g$data[[4]]["colour"])

###colour
#1 #F8766D
#2 #A3A500
#3 #00BF7D
#4 #00B0F6
#5 #E76BF3

####################
##### Site PCA #####
####################

GLCWsitePCA <- read.table("GLCW_site_dat.txt", header=TRUE, row.names = 1, sep = "\t")
names(GLCWsitePCA)

# transform pH
pH <- 10^-GLCWsitePCA[,6]
GLCWsitePCA <- cbind(GLCWsitePCA[,c(1:5,7:12)], pH)

# pearson correlation analysis
# check to see if any variables are autocorrelated
rcorr(as.matrix(GLCWsitePCA[,4:12]), type="pearson")

# List (corr > 0.7 and sig < 0.01 level)
# DO represents Redox_pot
# Temp represents Cond and TDS
GLCWsitePCA <- GLCWsitePCA[-c(6:7,9)]

# apply PCA
GLCWsite.pca <- prcomp(GLCWsitePCA[,4:9], center = TRUE, scale. = TRUE)

# print method
print(GLCWsite.pca)

# plot method
plot(GLCWsite.pca, type = "lines")

# summary method
summary(GLCWsite.pca)

#create PCA plot in ggplot

wetland <- as.factor(GLCWsitePCA$Wetland)
region <- as.factor(GLCWsitePCA$Region)
GLPCA <- ggbiplot(GLCWsite.pca, obs.scale = 1, var.scale = 1, 
                  groups = GLCWsitePCA$Region, ellipse = TRUE, 
                  circle = FALSE, varname.adjust = 1) +
  scale_color_discrete() +
  geom_point(aes(colour = wetland,
                 shape = region),size=3) +
  theme(legend.direction = 'vertical', 
        legend.position = 'right') +
  #geom_text_repel(aes(label = wetland)) +
  theme_classic()
GLPCA

#Finding point coordinates on biplot
names(GLchem.pca)
GLchem.pca$x


### EDIT: We've #'d out code pertaining to pH from hereforth, as we discovered
### pH could not be reliably assessed with the soil test implemented as soil was
### not tested in situ, influencing the accuracy of the test.

######################
##### Line plots #####
######################

# Create Line Chart

#Dean's WD
setwd("~/Grad School/Central Michigan/Dissertation/Great Lakes Coastal Wetlands/Stats")

#read in data
GLchem <- read.table("GLchemLineGraphs.txt", header = TRUE, sep="\t")

#create line plots for individual environmental variables

#GLphPlot <- ggplot(data=GLchem, aes(x = GLchem$pH, y = GLchem$depth, group = GLchem$wetland)) +
#  scale_y_reverse(lim = c(6,0)) +
#  geom_path(aes(color = GLchem$region), show.legend = FALSE) +
#  geom_point(aes(shape = GLchem$wetland, color = GLchem$region),show.legend = FALSE) +
#  scale_color_brewer(guide = FALSE, name ="Region",
#                       labels=c("BA","ESBT","LE","NSB","WSB"), type = "qual",
#                     palette = "Set1") +
#  scale_shape_manual(guide = FALSE, name="Wetland",
#                       labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB"), values = c(0,1,2,3,4,5,6,15,16,17)) +
#  scale_linetype_discrete(name ="Wetland",
#                          labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB")) +
#  labs(y = "Depth (cm)", x = "pH") +
#  theme_classic()
#GLphPlot

GLPPlot <- ggplot(data=GLchem, aes(x = GLchem$brayP..ppm., y = GLchem$depth, group = GLchem$wetland)) +
  scale_y_reverse(lim = c(6,0)) +
  geom_path(aes(color = GLchem$region), show.legend = FALSE) +
  geom_point(aes(shape = GLchem$wetland, color = GLchem$region),show.legend = FALSE) +
  scale_color_brewer(guide = FALSE, name ="Region",
                     labels=c("BA","ESBT","LE","NSB","WSB"), type = "qual",
                     palette = "Set1") +
  scale_shape_manual(guide = FALSE, name="Wetland",
                     labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB"), values = c(0,1,2,3,4,5,6,15,16,17)) +
  scale_linetype_discrete(name ="Wetland",
                          labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB")) +
  labs(y = "Depth (cm)", x = "P (ppm)") +
  theme_classic()
  #theme(axis.title.y=element_blank(),
  #    axis.text.y=element_blank(),
  #    axis.ticks.y=element_blank())
GLPPlot

GLSPlot <- ggplot(data=GLchem, aes(x = GLchem$sulfur..ppm., y = GLchem$depth, group = GLchem$wetland)) +
  scale_y_reverse(lim = c(6,0)) +
  geom_path(aes(color = GLchem$region), show.legend = FALSE) +
  geom_point(aes(shape = GLchem$wetland, color = GLchem$region),show.legend = FALSE) +
  scale_color_brewer(guide = FALSE, name ="Region",
                     labels=c("BA","ESBT","LE","NSB","WSB"), type = "qual",
                     palette = "Set1") +
  scale_shape_manual(guide = FALSE, name="Wetland",
                     labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB"), values = c(0,1,2,3,4,5,6,15,16,17)) +
  scale_linetype_discrete(name ="Wetland",
                          labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB")) +
  labs(y = "Depth (cm)", x = "S (ppm)") +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
GLSPlot

GLOMPlot <- ggplot(data=GLchem, aes(x = GLchem$organic.matter...., y = GLchem$depth, group = GLchem$wetland)) +
  scale_y_reverse(lim = c(6,0)) +
  geom_path(aes(color = GLchem$region), show.legend = FALSE) +
  geom_point(aes(shape = GLchem$wetland, color = GLchem$region),show.legend = FALSE) +
  scale_color_brewer(guide = FALSE, name ="Region",
                     labels=c("BA","ESBT","LE","NSB","WSB"), type = "qual",
                     palette = "Set1") +
  scale_shape_manual(guide = FALSE, name="Wetland",
                     labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB"), values = c(0,1,2,3,4,5,6,15,16,17)) +
  scale_linetype_discrete(name ="Wetland",
                          labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB")) +
  labs(y = "Depth (cm)", x = "OM (%)") +
  theme_classic()+
  theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
GLOMPlot

GLOCPlot <- ggplot(data=GLchem, aes(x = GLchem$organic.carbon...., y = GLchem$depth, group = GLchem$wetland)) +
  scale_y_reverse(lim = c(6,0)) +
  geom_path(aes(color = GLchem$region), show.legend = FALSE) +
  geom_point(aes(shape = GLchem$wetland, color = GLchem$region),show.legend = FALSE) +
  scale_color_brewer(guide = FALSE, name ="Region",
                     labels=c("BA","ESBT","LE","NSB","WSB"), type = "qual",
                     palette = "Set1") +
  scale_shape_manual(guide = FALSE, name="Wetland",
                     labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB"), values = c(0,1,2,3,4,5,6,15,16,17)) +
  scale_linetype_discrete(name ="Wetland",
                          labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB")) +
  labs(y = "Depth (cm)", x = "OC (%)") +
  theme_classic()
  #theme(axis.title.y=element_blank(),
  #      axis.text.y=element_blank(),
  #      axis.ticks.y=element_blank())
GLOCPlot

GLTNPlot <- ggplot(data=GLchem, aes(x = GLchem$tn...., y = GLchem$depth, group = GLchem$wetland)) +
  scale_y_reverse(lim = c(6,0)) +
  geom_path(aes(color = GLchem$region), show.legend = FALSE) +
  geom_point(aes(shape = GLchem$wetland, color = GLchem$region),show.legend = FALSE) +
  scale_color_brewer(guide = FALSE, name ="Region",
                     labels=c("BA","ESBT","LE","NSB","WSB"), type = "qual",
                     palette = "Set1") +
  scale_shape_manual(guide = FALSE, name="Wetland",
                     labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB"), values = c(0,1,2,3,4,5,6,15,16,17)) +
  scale_linetype_discrete(name ="Wetland",
                          labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB")) +
  labs(y = "Depth (cm)", x = "TN (%)") +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
GLTNPlot

GLCNPlot <- ggplot(data=GLchem, aes(x = GLchem$c.n, y = GLchem$depth, group = GLchem$wetland)) +
  scale_y_reverse(lim = c(6,0)) +
  geom_path(aes(color = GLchem$region), show.legend = FALSE) +
  geom_point(aes(shape = GLchem$wetland, color = GLchem$region),show.legend = FALSE) +
  scale_color_brewer(guide = FALSE, name ="Region",
                     labels=c("BA","ESBT","LE","NSB","WSB"), type = "qual",
                     palette = "Set1") +
  scale_shape_manual(guide = FALSE, name="Wetland",
                     labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB"), values = c(0,1,2,3,4,5,6,15,16,17)) +
  scale_linetype_discrete(name ="Wetland",
                          labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB")) +
  labs(y = "Depth (cm)", x = "C:N") +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
GLCNPlot

GLNO3Plot <- ggplot(data=GLchem, aes(x = GLchem$no3..ppm., y = GLchem$depth, group = GLchem$wetland)) +
  scale_y_reverse(lim = c(6,0)) +
  geom_path(aes(color = GLchem$region), show.legend = FALSE) +
  geom_point(aes(shape = GLchem$wetland, color = GLchem$region),show.legend = FALSE) +
  scale_color_brewer(guide = FALSE, name ="Region",
                     labels=c("BA","ESBT","LE","NSB","WSB"), type = "qual",
                     palette = "Set1") +
  scale_shape_manual(guide = FALSE, name="Wetland",
                     labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB"), values = c(0,1,2,3,4,5,6,15,16,17)) +
  scale_linetype_discrete(name ="Wetland",
                          labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB")) +
  labs(y = "Depth (cm)", x = "[NO3-] (ppm)") +
  theme_classic()
  #theme(axis.title.y=element_blank(),
  #      axis.text.y=element_blank(),
  #      axis.ticks.y=element_blank())
GLNO3Plot

GLNH4Plot <- ggplot(data=GLchem, aes(x = GLchem$nh4..ppm., y = GLchem$depth, group = GLchem$wetland)) +
  scale_y_reverse(lim = c(6,0)) +
  geom_path(aes(color = GLchem$region), show.legend = TRUE) +
  geom_point(aes(shape = GLchem$wetland, color = GLchem$region),show.legend = TRUE) +
  scale_color_brewer(name ="Region",
                     labels=c("BA","ESBT","LE","NSB","WSB"), type = "qual",
                     palette = "Set1") +
  scale_shape_manual(name="Wetland",
                     labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB"), values = c(0,1,2,3,4,5,6,15,16,17)) +
  scale_linetype_discrete(name ="Wetland",
                          labels=c("BA","ESBTA","ESBTB","ESBTC","LEC","LED","NSBA","NSBC","WSBA","WSBB")) +
  labs(y = "Depth (cm)", x = "[NH4+] (ppm)") +
  theme_classic() +
  theme(legend.position = "right", legend.box = "horizontal",panel.background = element_rect(fill = FALSE), 
        legend.key = element_rect(fill = FALSE)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  guides(shape = guide_legend(ncol = 2))
GLNH4Plot

# comboplot

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

GLconcLegend<-g_legend(GLNH4Plot)

GLChemPlot <- grid.arrange(arrangeGrob(
                               GLPPlot + theme(legend.position="none"),
                               GLSPlot + theme(legend.position="none"),
                               GLOMPlot + theme(legend.position="none"),
                               GLOCPlot + theme(legend.position="none"),
                               GLTNPlot + theme(legend.position="none"),
                               GLCNPlot + theme(legend.position="none"),
                               GLNO3Plot + theme(legend.position="none"),
                               GLNH4Plot + theme(legend.position="none"),
                               GLconcLegend,
                               nrow=3),ncol=1)


###############
##### PCA #####
###############

GLchemPCA <- read.table("GLchemPCAclean.txt",header=TRUE, sep="\t")
names(GLchemPCA)

# Arcsin(sqrt(x/100)) transformation for percentages
GLchemArcsin <- asin(sqrt(GLchemPCA[,7:9]/100))

# Log transform ratios
CNrat <- log(GLchemPCA[,10])

# change pH to [H+]
pH <- 10^-(GLchemPCA[,4])

# combine transformations into new table

GLchemtrans <- cbind(GLchemPCA[,c(1:3,5:6,11:12)], GLchemArcsin[,1:3], CNrat, pH)
#GLchemtrans.log <- (cbind(GLchemtrans[1:3],log(1 + GLchemtrans[,4:12])))
write.csv(GLchemtrans, "GLchemtrans.csv")

# pearson correlation analysis
# check to see if any variables are autocorrelated
rcorr(as.matrix(GLchemtrans[,4:12]), type="pearson")

# List (corr > 0.7 and sig < 0.001 level)
# OM represents NH4, OC, pH, & TN, as it was the most highly correlated to all vars.
# OM renamed as "[H+] + NUTR"
# S highly correlated with OM, OC, & NH4, but not NO3, so the cheese stands alone
# NO3 also did not correlate with S, so it was also left alone
# **WARNING** lots more variables sig <.01, but corr < 0.7 **WARNING**
# create new table with NH4, pH, OC, & TN removed (pH because it's useless now)

GLchemtranscorrs <- cbind(GLchemtrans[,c(1:6,8,11)])
colnames(GLchemtranscorrs)[7] <- "NUTR"


# apply PCA

GLchem.pca <- prcomp(GLchemtranscorrs[,4:8], center = TRUE, scale. = TRUE)


# print method
print(GLchem.pca)

# plot method
plot(GLchem.pca, type = "lines")

# summary method
summary(GLchem.pca)

#create PCA plot in ggplot

depth <- as.factor(GLchemtranscorrs$depth)
GLPCA <- ggbiplot(GLchem.pca, obs.scale = 1, var.scale = 1,
                  varname.adjust = 1, shape = GLchemtranscorrs$depth,
                  color = GLchemtranscorrs$region, fill = GLchemtranscorrs$region) +
  geom_point(aes(shape = factor(GLchemtranscorrs$depth),
                 color = factor(GLchemtranscorrs$region),
                 fill = factor(GLchemtranscorrs$region)),size=3) +
  scale_shape_manual(name = "Depth", values = c(21,22,24), labels = c("Top", "Middle", "Bottom")) +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Region") +
  scale_fill_brewer(type = "qual", palette = "Set1", guide = FALSE) +
  theme(legend.direction = 'vertical', 
               legend.position = 'right', aspect.ratio = 1) +
  stat_ellipse(aes(fill=factor(GLchemtranscorrs$region)),
               geom="polygon", level=0.95, alpha=0.2) +
  #coord_fixed(ratio = 1, xlim = c(-2,5), ylim = c(-2,2.5), expand = TRUE) +
  theme_classic()# +
GLPCA

#Finding point coordinates on biplot
names(GLchem.pca)
GLchem.pca$x


# MANOVA
PCA_manova <- manova(GLchem.pca$x ~ region*depth, data = GLchemtranscorrs)
summary(PCA_manova)

PCA_perMANOVA <- adonis(GLchem.pca$x ~ region*depth*wetland,
                        data = GLchemtranscorrs, method = "euclidean")
PCA_perMANOVA
pairwise.adonis(GLchem.pca$x, factors = GLchemtranscorrs$region)


###########################
###########################
##### Alpha Diversity #####
###########################
###########################

alpha <- read.table("alpha_div.txt", header = TRUE, row.names = 1, sep = "\t")
alpha

### rarefaction curves
alpha_otu <- read.table("GREATLAKES.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.pick.an.unique_list.shared")
alpha_otu2 <- alpha_otu[,-c(1,3)]
talpha_OTU <- t(alpha_otu2)
write.csv(talpha_OTU,"talpha_OTU.csv")
alpha_otu <- read.csv("talpha_OTU.csv", header = TRUE, row.names = 1)
alpha_otut <- t(alpha_otu)
alpharare <- rarecurve(alpha_otut, step = 20000, col = 2:75, label = FALSE, ylab = "OTUs")
abline(v = 48226)


### Chao1
#Region
#GL.chao.lme=lme(fixed=chao~Depth*Region, random=~1|factor(Wetland), data=alpha)
#GL_chao_lme_anova <- anova(GL.chao.lme)
#GL_chao_lme_anova

#region is significant, rerun w/out Depth since depth and interactions were not
#signficant. Also, we can run pairwise tests without glht getting pissed...
GL.chao.lme=lme(fixed=chao~Region, random=~1|factor(Wetland), data=alpha)
GL_chao_lme_anova <- anova(GL.chao.lme)
GL_chao_lme_anova

#Wetland
#GL.chao.lm.wetland=lm(chao~Wetland*Depth, data=alpha)
#GL_chao_wetland_anova <- anova(GL.chao.lm.wetland)
#GL_chao_wetland_anova
#wetland is significant, depth and interactions are not. Do same as above
GL.chao.lm.wetland=lm(chao~Wetland, data=alpha)
GL_chao_wetland_anova <- anova(GL.chao.lm.wetland)
GL_chao_wetland_anova


### npshannon region & depth

GL.npshannon.lme=lme(fixed=npshannon~Depth*Region, random=~1|factor(Wetland), data=alpha)
GL_npshannon_lme_anova <- anova(GL.npshannon.lme)
GL_npshannon_lme_anova
#no significance

#GL.npshannon.lm.wetland=lm(npshannon~Wetland*Depth, data=alpha)
#GL_npshannon_wetland_anova <- anova(GL.npshannon.lm.wetland)
#GL_npshannon_wetland_anova
#only wetland is significant. Do same as above.
GL.npshannon.lm.wetland=lm(npshannon~Wetland, data=alpha)
GL_npshannon_wetland_anova <- anova(GL.npshannon.lm.wetland)
GL_npshannon_wetland_anova


# Tukey Tests for significant results
Tukey_chao_lme <- summary(glht(GL.chao.lme, linfct = mcp(Region = "Tukey")),test = adjusted(type = "bonferroni"))
Tukey_chao_lme
a <- plot(cld(Tukey_chao_lme))

Tukey_chao_wetland <- summary(glht(GL.chao.lm.wetland, linfct = mcp(Wetland = "Tukey")),test = adjusted(type = "bonferroni"))
Tukey_chao_wetland
b <- plot(cld(Tukey_chao_wetland), las = 2)

Tukey_npshannon_lme <- summary(glht(GL.npshannon.lm.wetland, linfct = mcp(Wetland = "Tukey")),test = adjusted(type = "bonferroni"))
Tukey_npshannon_lme
confint(GL.npshannon.lm.wetland)
c <- plot(cld(Tukey_npshannon_lme), las = 2)

#layout(matrix(c(1,2,3,2), 2, 2, byrow = TRUE), heights = 1) couldn't get this damned thing to work
par(mfrow = c(1,3))
par(mar=c(5,5,8,1), xpd = NA)#oma=c(1,1,.0001,1))
plot(cld(Tukey_chao_lme, level = 0.01),xlab= '',ylab= '',las = 2, col = brewer.pal(5,"Set1"))
title(ylab="chao1", line=4)
mtext(paste0("(a)"), side = 3, adj = 0.05, 
      line = -1.3)
par(mar=c(5,0,8,1), xpd = NA)#oma=c(1,1,.0001,1))
plot(cld(Tukey_chao_wetland, level = 0.01),xlab= '',ylab='',yaxt = "n", las = 2, col = brewer.pal(5,"Set1")[c(1,2,2,2,3,3,4,4,5,5)])
mtext(paste0("(b)"), side = 3, adj = 0.05, 
        line = -1.3)
par(mar=c(5,4,8,1), xpd = NA)#oma=c(1,1,.0001,1))
plot(cld(Tukey_npshannon_lme, level = 0.01),xlab= '',las = 2., col = brewer.pal(5,"Set1")[c(1,2,2,2,3,3,4,4,5,5)])
mtext(paste0("(c)"), side = 3, adj = 0.05, 
      line = -1.3)

dev.off()

### Alpha - Environ correlations
corr_meta_GL <- read.csv("alpha_corr_meta.csv", row.names = 1)
corr_meta_GL_vals <- corr_meta_GL[3:8]
corr_meta_GL_vals <- as.matrix(corr_meta_GL_vals)
alphacorr_GL <- rcorr(as.matrix(alpha[,4:5]), corr_meta_GL_vals, type = "spearman")
alphacorr_GL

##########################
##########################
##### Beta Diversity #####
##########################
##########################

#gltransposed <- t(read.table("GREATLAKES.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.abund.shared"))
#write.table(gltransposed, "gl_t_singdoubrem.txt")

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

#################################################
##### Normalizing and transforming the data #####
#################################################

### Deseq2 Normalization
### convert the counts to integer mode and group data to normalize together(e.g. region)
dds_region <- phyloseq_to_deseq2(ALL, ~ region)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans_region = apply(counts(dds_region), 1, gm_mean)
dds_region = estimateSizeFactors(dds_region, geoMeans=geoMeans_region)
dds_region = estimateDispersions(dds_region)

### variance stabilizing transformation
#Make it a OTU again for phyloseq (no tax or sample data)
vst_region <- varianceStabilizingTransformation(dds_region, blind=FALSE)
vstMat_region <- assay(vst_region)
vstMat_region[vstMat_region<0]<-0
vst.otu.region <- otu_table(vstMat_region, taxa_are_rows=TRUE)
write.csv(otu_table(vst.otu.region), "vst.otu.region.csv")

###Create a distance matrix 
region_Dist <- phyloseq::distance(vst.otu.region, method= 'bray')

##NMDS of all samples
NMDS_region <- ordinate(vst.otu.region, method= 'NMDS', distance = region_Dist, formula = NULL)

## plot quick NMDS for comparisons
NMDS_region_plot <- plot_ordination(ALL, NMDS_region, color="region", label='depth')
NMDS_region_plot + geom_point(size=1) + theme_bw()


##################################################################
##### Publication-level NMDS figures and envfit correlations #####
##################################################################

##### Region NMDS
##### envfit depth
# NMDS_META <- read.table("GLCW_NMDS_metadata.txt", header = T, row.names = 1)
NMDS_META <- read.table("GLCW_NMDS_metadata_clean.txt", header = T, row.names = 1)
colnames(NMDS_META)[10] <- "NUTR"
ef_region <- envfit(NMDS_region, NMDS_META2[,c(3:6,8,10)], permu=999, na.rm = TRUE)
ef_region

vectors_reg<-as.data.frame(ef_region$vectors$arrows*sqrt(ef_region$vectors$r))
pvals_reg=ef_region$vectors$pvals
r_reg=ef_region$vectors$r
envfit_region<-cbind(vectors_reg,pvals_reg, r_reg)
envfit_region
#vectors are too long to fit on plot, transform with /2
envfit_region2 <- cbind(envfit_region[,1:2]/2, envfit_region[,3:4])
envfit_region2

NMDS_region_pub = data.frame(MDS1 = NMDS_region$points[,1], MDS2 = NMDS_region$points[,2])
NMDS_region_pub <- cbind(NMDS_region_pub,NMDS_META[,1:3])
NMDS_region_pub
glNMDSplotReg <- ggplot(NMDS_region_pub, aes(x=NMDS_region_pub$MDS1, y=NMDS_region_pub$MDS2)) +
  geom_point(aes(fill = factor(NMDS_region_pub$region),
                 shape = NMDS_region_pub$depth,
                 color = NMDS_region_pub$region),
             size=5) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(legend.direction="vertical",
        legend.position="right") +
  scale_shape_manual(name="Depth",
                     labels=c("Top", "Middle", "Bottom"),
                     values = c(21,22,24)) +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Region") +
  scale_fill_brewer(type = "qual", palette = "Set1", guide = FALSE) +
  geom_text(label="Stress = 0.148", x = -0.4, y = -0.4, size = 4) +
  coord_equal(ylim = 0, xlim = 0) +
  guides(fill = FALSE) +
  theme_classic()# +  
  #guides(fill = guide_legend(override.aes= list(colour = c("white","grey","black")))) +
  #theme(legend.key = element_rect(fill = "grey90"))
#glNMDSplotReg <- glNMDSplotReg + geom_text_repel(aes(label = NMDS_region_pub$depth))
glNMDSplotReg <- glNMDSplotReg + 
  geom_segment(data = envfit_region2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=envfit_region2, aes(x=NMDS1*1.05, y=NMDS2*1.05, label=rownames(envfit_region2)),size=4)
glNMDSplotReg

###################
##### ANOSIMS #####
###################

### REGION vst ANOSIMs
# testing region effect with ANOSIM

region_anosim <- anosim(region_Dist, NMDS_META$region, permutations = 999, 
                        distance = "bray", strata = NMDS_META$depth)
summary(region_anosim)
plot(region_anosim)

# testing wetland effect with ANOSIM
wetland_anosim <- anosim(region_Dist, NMDS_META$wetland, permutations = 999, distance = "bray")
summary(wetland_anosim)
plot(wetland_anosim)

# testing depth effect with ANOSIM 
depth_anosim <- anosim(region_Dist, NMDS_META$depth, permutations = 999, distance = "bray",
                       strata = NMDS_META$region)
summary(depth_anosim)
plot(depth_anosim)


##############
### ADONIS ###
##############

adonis_GL <- adonis(region_Dist ~ region*depth*wetland, data = NMDS_META)
adonis_GL


## Pairwise perMANOVA for above tests
BA_v_ESBT_META <- NMDS_META[1:35,]
BA_v_ESBT_anosim <- anosim(BA_v_ESBT_Dist, BA_v_ESBT_META$region, permutations = 999, 
                        distance = "bray", strata = BA_v_ESBT_META$depth)
summary(BA_v_ESBT_anosim)
plot(BA_v_ESBT_anosim)

vst_BA_v_ESBT <- otu_table(vstMat_region[,1:35], taxa_are_rows=TRUE)
BA_v_ESBT_Dist <- phyloseq::distance(vst_BA_v_ESBT, method= 'bray')

rdist_mat <- as.matrix(region_Dist)







##start copy here for function pairwise.adonis()

pairwise.adonis <- function(x,factors, sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
} 

#How can I do PerMANOVA pairwise contrasts in R?. Available from: https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R [accessed Jun 23, 2017].

pairwise.adonis(rdist_mat, factors = NMDS_META$region)

pairwise.adonis(rdist_mat, factors = NMDS_META$wetland)

pairwise.adonis(rdist_mat, factors = NMDS_META$depth)

## Pairwise wilcoxon rank sum tests for above tests
#dist <- as.matrix(region_Dist)
#dist

#pwrst_GL_wetland <- pairwise.wilcox.test(dist,NMDS_META$wetland, p.adj = "holm", conf.int = TRUE, conf.level = 0.05)
#pwrst_GL_wetland

#pwrst_GL_region <- pairwise.wilcox.test(dist,NMDS_META$region, p.adj = "holm", conf.int = TRUE, conf.level = 0.05)
#pwrst_GL_region

#pwrst_GL_depth <- pairwise.wilcox.test(dist,NMDS_META$depth, p.adj = "holm", conf.int = TRUE, conf.level = 0.05)
#pwrst_GL_depth
#######################
### Beta-dispersion ###
#######################

betadisp_region <- betadisper(region_Dist, NMDS_META$region)
betadisp_region
anova(betadisp_region)
permutest(betadisp_region)
plot(betadisp_region)
boxplot(betadisp_region)
TukeyHSD(betadisp_region)
plot(TukeyHSD(betadisp_region))
betadisp_region$eig/sum(betadisp_region$eig)
# No significant differences found post-hoc

betadisp_depth <- betadisper(region_Dist, NMDS_META$depth)
betadisp_depth
anova(betadisp_depth)
plot(betadisp_depth)
#no significance in beta-dispersion between depths

betadisp_wetland <- betadisper(region_Dist, NMDS_META$wetland)
betadisp_wetland
anova(betadisp_wetland)
permutest(betadisp_wetland)
plot(betadisp_wetland)
tuk_betadisp_wetl <- TukeyHSD(betadisp_wetland)
tuk_betadisp_wetl
betadisp_wetland$eig/sum(betadisp_wetland$eig)


boxplot(betadisp_wetland,col = brewer.pal(5,"Set1")[c(1,2,2,2,3,3,4,4,5,5)])#, las = 2)
lines(c(3,5), y= c(.455,.455))
lines(c(3,3),y = c(.455, .445))
lines(c(5,5),y = c(.455, .445))
lines(c(5,8),y = c(.475, .475))
lines(c(5,5),y = c(.475, .465))
lines(c(8,8),y = c(.475, .465))
lines(c(5,10),y = c(.495, .495))
lines(c(5,5),y = c(.495, .485))
lines(c(10,10),y = c(.495, .485))
text(c(4,6.5,7.6), y = c(.465,.485,.505),"*",cex = 2)



###########
### CCA ###
###########

#OTU table must have no rows with rowSums = 0
#otu_tab_cca <- read.csv("vst.otu.CMUBS_ALL.csv", header = TRUE, row.names = 1)
otu_tab_0less <- vst.otu.region[rowSums(vst.otu.region)!=0,]
otu_tab_0less_t <- t(otu_tab_0less)
otu_tab_0less_t <- as.matrix(otu_tab_0less_t)
CCA_meta <- as.data.frame(NMDS_META[,c(3,5,7,9)])
#both otu table and metadata must be 'numeric' and same rowlength
mode(CCA_meta)
mode(otu_tab_0less_t)
class(otu_tab_0less_t)
dim(otu_tab_0less)
dim(CCA_meta)

CMUBS_CCA <- cca(otu_tab_0less_t ~ depth + S + NUTR + pH, data = CCA_meta)
CMUBS_CCA
plot(CMUBS_CCA)

anova_CCA <- anova(CMUBS_CCA, alpha = 0.05, beta = 0.01, step=100, perm.max=9999)
anova_CCA
anova_terms <- anova(CMUBS_CCA, by = "terms", permu = 999, cutoff = 0.001)
anova_terms
anova_axes <- anova(CMUBS_CCA, by = "axis", permu = 999, cutoff = 0.001)
anova_axes
CCApermtest <- permutest(CMUBS_CCA, permutations = 999,
                         model = "reduced")
CCApermtest

#extract vectors
CCAvecs <- as.data.frame(cbind(CMUBS_CCA$CCA$biplot[,1], CMUBS_CCA$CCA$biplot[,2]))

CCApub <- data.frame(CC1 = CMUBS_CCA$CCA$u[,1], CC2 = CMUBS_CCA$CCA$u[,2])
CCApub <- cbind(CCApub,CCA_meta,CMUBS_META[,1:3])
CCAplot <- ggplot(CCApub, aes(x=CCApub$CC1, y=CCApub$CC2)) +
  geom_point(aes(fill = CCApub$Habitat,
                 shape = CCApub$Lake),
             size=5) +
  xlab("CCA1 (31.28%)") +
  ylab("CC2(28.56%)") +
  theme(legend.direction="vertical",
        legend.position="right") +
  theme_bw() +
  scale_shape_manual(values = c(21,22,23,24), aes(fill=Lake),
                     labels=c("BL", "FL", "LG", "LM"),
                     guide = FALSE) +
  scale_fill_manual(name="Habitat",
                    values = c(1,"white")) +
  coord_fixed()
CCAplot <- CCAplot + 
  geom_segment(data = CCAvecs, aes(x=0, xend=V1*2, y=0, yend=V2*2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=CCAvecs, aes(x=V1*2, y=V2*2, label=rownames(CCAvecs)),size=4)
CCAplot <- CCAplot + geom_text_repel(aes(label = CCApub$Time))

CCAplot


##############################
##### Pointless Heatmaps #####
##############################

# Create a new phyloseq file with the normalized OTUs, sample data, and taxonomy
NorReg = phyloseq(vst.otu.region, TAX, META)

# Heatmap for the 50 most abundant bacteria (vst region)
bact_region<- subset_taxa(NorReg, Domain == 'Bacteria(100)')
bact_region<- prune_taxa(names(sort(taxa_sums(bact_region), TRUE) [1:50]), bact_region)
plot_heatmap(bact_region)

# Heatmap for the 50 most abundant archaea (vst region)
arc_region<- subset_taxa(NorReg, Domain == 'Archaea(100)')
arc_region<- prune_taxa(names(sort(taxa_sums(arc_region), TRUE) [1:50]), arc_region)
plot_heatmap(arc_region)