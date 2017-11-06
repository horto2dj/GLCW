NMDS_META <- read.table("GLCW_NMDS_metadata_clean.txt", header = T, row.names = 1)

##########
### BA ###
##########

### Create new NMDS metadata so that depth contains integers
depthint <- as.numeric(NMDS_META$depth)
NMDS_META2 <- NMDS_META[,-c(3)]
NMDS_META2 <- cbind(NMDS_META2[,1:2],depthint,NMDS_META2[3:11])
colnames(NMDS_META2)[3] <- "depth"
colnames(NMDS_META2)[10] <- "NUTR"

###Create a distance matrix 
region_Dist_BA <- phyloseq::distance(vst.otu.region[,1:9], method= 'bray')

##NMDS of all samples
NMDS_region_BA <- ordinate(vst.otu.region[,1:9], method= 'NMDS', distance = region_Dist_BA, formula = NULL)

## plot quick NMDS for comparisons
NMDS_region_plot_BA <- plot_ordination(ALL, NMDS_region_BA, color="region", label='depth')
NMDS_region_plot_BA + geom_point(size=1) + theme_bw()

##### envfit depth
ef_BA <- envfit(NMDS_region_BA, NMDS_META2[1:9,c(3:6,8,10)], permu=999, na.rm = TRUE)
ef_BA

vectors_reg<-as.data.frame(ef_BA$vectors$arrows*sqrt(ef_BA$vectors$r))
pvals_reg=ef_BA$vectors$pvals
r_reg=ef_BA$vectors$r
envfit_BA<-cbind(vectors_reg,pvals_reg, r_reg)
envfit_BA
envfit_BA <- envfit_BA[c(1:3,6),]
#vectors are too long to fit on plot, transform with /2
envfit_BA2 <- cbind(envfit_BA[,1:2]/3, envfit_BA[,3:4]/3)
envfit_BA2

NMDS_region_pub_BA = data.frame(MDS1 = NMDS_region_BA$points[,1], MDS2 = NMDS_region_BA$points[,2])
NMDS_region_pub_BA <- cbind(NMDS_region_pub_BA,NMDS_META[1:9,1:3])
NMDS_region_pub_BA
glNMDSplotReg_BA <- ggplot(NMDS_region_pub_BA, aes(x=NMDS_region_pub_BA$MDS1, y=NMDS_region_pub_BA$MDS2)) +
  geom_point(aes(fill = factor(NMDS_region_pub_BA$region),
                 shape = NMDS_region_pub_BA$depth,
                 color = NMDS_region_pub_BA$wetland),
             size=5) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(legend.direction="vertical",
        legend.position="right") +
  scale_shape_manual(name="Depth",
                     labels=c("Top", "Middle", "Bottom"),
                     values = c(21,22,24)) +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Site") +
  scale_fill_brewer(type = "qual", palette = "Set1", guide = FALSE) +
  geom_text(label="Stress = 0.010", x = -0.2, y = -0.4, size = 4) +
  coord_equal(ylim = 0, xlim = 0) +
  guides(fill = FALSE, color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  geom_text(label="a. BA", x = -0.4, y = 0.45, size = 4) +
  theme_classic()# +  
#guides(fill = guide_legend(override.aes= list(colour = c("white","grey","black")))) +
#theme(legend.key = element_rect(fill = "grey90"))
#glNMDSplotReg <- glNMDSplotReg + geom_text_repel(aes(label = NMDS_region_pub$depth))
glNMDSplotReg_BA <- glNMDSplotReg_BA + 
  geom_segment(data = envfit_BA2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=envfit_BA2, aes(x=NMDS1*1.05, y=NMDS2*1.05, label=rownames(envfit_BA2)),size=4)
glNMDSplotReg_BA

############
### ESBT ###
############

###Create a distance matrix 
region_Dist_ESBT <- phyloseq::distance(vst.otu.region[,10:35], method= 'bray')

##NMDS of all samples
NMDS_region_ESBT <- ordinate(vst.otu.region[,10:35], method= 'NMDS', 
                             distance = region_Dist_ESBT, formula = NULL)
NMDS_region_ESBT
## plot quick NMDS for comparisons
NMDS_region_plot_ESBT <- plot_ordination(ALL, NMDS_region_ESBT, 
                                       color="wetland", label='depth')
NMDS_region_plot_ESBT + geom_point(size=1) + theme_bw()

##### envfit depth
ef_ESBT <- envfit(NMDS_region_ESBT, NMDS_META2[10:35,c(3:6,8,10)], permu=999, 
                  na.rm = TRUE)
ef_ESBT

vectors_reg<-as.data.frame(ef_ESBT$vectors$arrows*sqrt(ef_ESBT$vectors$r))
pvals_reg=ef_ESBT$vectors$pvals
r_reg=ef_ESBT$vectors$r
envfit_ESBT<-cbind(vectors_reg,pvals_reg, r_reg)
envfit_ESBT
#keep only those vectors which significantly correlated to structure
envfit_ESBT <- envfit_ESBT[c(1:3,5:6),]
#vectors are too long to fit on plot, transform with /2
envfit_ESBT2 <- cbind(envfit_ESBT[,1:2]/3, envfit_ESBT[,3:4]/3)
envfit_ESBT2

NMDS_region_pub_ESBT = data.frame(MDS1 = NMDS_region_ESBT$points[,1],
                                  MDS2 = NMDS_region_ESBT$points[,2])
NMDS_region_pub_ESBT <- cbind(NMDS_region_pub_ESBT,NMDS_META[10:35,1:3])
NMDS_region_pub_ESBT
glNMDSplotReg_ESBT <- ggplot(NMDS_region_pub_ESBT, aes(x=NMDS_region_pub_ESBT$MDS1,
                                                       y=NMDS_region_pub_ESBT$MDS2)) +
  geom_point(aes(fill = factor(NMDS_region_pub_ESBT$wetland),
                 shape = NMDS_region_pub_ESBT$depth,
                 color = NMDS_region_pub_ESBT$wetland),
             size=5) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(legend.direction="vertical",
        legend.position="right") +
  scale_shape_manual(name="Depth",
                     labels=c("Top", "Middle", "Bottom"),
                     values = c(21,22,24)) +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Site") +
  scale_fill_brewer(type = "qual", palette = "Set1", guide = FALSE) +
  geom_text(label="Stress = 0.076", x = -0.2, y = -0.4, size = 4) +
  coord_equal(ylim = 0, xlim = 0) +
  guides(fill = FALSE, color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  geom_text(label="b. ESBT", x = -0.36, y = 0.45, size = 4) +
  theme_classic()# +  
#guides(fill = guide_legend(override.aes= list(colour = c("white","grey","black")))) +
#theme(legend.key = element_rect(fill = "grey90"))
#glNMDSplotReg <- glNMDSplotReg + geom_text_repel(aes(label = NMDS_region_pub$depth))
glNMDSplotReg_ESBT <- glNMDSplotReg_ESBT + 
  geom_segment(data = envfit_ESBT2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=envfit_ESBT2, aes(x=NMDS1*1.05, y=NMDS2*1.05,
                                   label=rownames(envfit_ESBT2)),size=4)
glNMDSplotReg_ESBT

##########
### LE ###
##########

###Create a distance matrix 
region_Dist_LE <- phyloseq::distance(vst.otu.region[,36:47], method= 'bray')

##NMDS of all samples
NMDS_region_LE <- ordinate(vst.otu.region[,36:47], method= 'NMDS', 
                             distance = region_Dist_LE, formula = NULL)

## plot quick NMDS for comparisons
NMDS_region_plot_LE <- plot_ordination(ALL, NMDS_region_LE, 
                                         color="region", label='depth')
NMDS_region_plot_LE + geom_point(size=1) + theme_bw()

##### envfit depth
ef_LE <- envfit(NMDS_region_LE, NMDS_META2[36:47,c(3:6,8,10)], permu=999, 
                  na.rm = TRUE)
ef_LE

vectors_reg<-as.data.frame(ef_LE$vectors$arrows*sqrt(ef_LE$vectors$r))
pvals_reg=ef_LE$vectors$pvals
r_reg=ef_LE$vectors$r
envfit_LE<-cbind(vectors_reg,pvals_reg, r_reg)
envfit_LE
#keep only those vectors which significantly correlated to structure
envfit_LE <- envfit_LE[c(3,6),]
#vectors are too long to fit on plot, transform with /2
envfit_LE2 <- cbind(envfit_LE[,1:2]/3, envfit_LE[,3:4]/3)
envfit_LE2

NMDS_region_pub_LE = data.frame(MDS1 = NMDS_region_LE$points[,1],
                                  MDS2 = NMDS_region_LE$points[,2])
NMDS_region_pub_LE <- cbind(NMDS_region_pub_LE,NMDS_META[36:47,1:3])
NMDS_region_pub_LE
glNMDSplotReg_LE <- ggplot(NMDS_region_pub_LE, aes(x=NMDS_region_pub_LE$MDS1,
                                                       y=NMDS_region_pub_LE$MDS2)) +
  geom_point(aes(fill = factor(NMDS_region_pub_LE$wetland),
                 shape = NMDS_region_pub_LE$depth,
                 color = NMDS_region_pub_LE$wetland),
             size=5) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(legend.direction="vertical",
        legend.position="right") +
  scale_shape_manual(name="Depth",
                     labels=c("Top", "Middle", "Bottom"),
                     values = c(21,22,24)) +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Site") +
  scale_fill_brewer(type = "qual", palette = "Set1", guide = FALSE) +
  geom_text(label="Stress = 0.040", x = -0.2, y = -0.4, size = 4) +
  coord_equal(ylim = 0, xlim = 0) +
  guides(fill = FALSE, color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  geom_text(label="c. LE", x = -0.4, y = 0.45, size = 4) +
  theme_classic()# +  
#guides(fill = guide_legend(override.aes= list(colour = c("white","grey","black")))) +
#theme(legend.key = element_rect(fill = "grey90"))
#glNMDSplotReg <- glNMDSplotReg + geom_text_repel(aes(label = NMDS_region_pub$depth))
glNMDSplotReg_LE <- glNMDSplotReg_LE + 
  geom_segment(data = envfit_LE2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=envfit_LE2, aes(x=NMDS1*1.05, y=NMDS2*1.05,
                                   label=rownames(envfit_LE2)),size=4)
glNMDSplotReg_LE

###########
### NSB ###
###########

###Create a distance matrix 
region_Dist_NSB <- phyloseq::distance(vst.otu.region[,48:62], method= 'bray')

##NMDS of all samples
NMDS_region_NSB <- ordinate(vst.otu.region[,48:62], method= 'NMDS', 
                             distance = region_Dist_NSB, formula = NULL)

## plot quick NMDS for comparisons
NMDS_region_plot_NSB <- plot_ordination(ALL, NMDS_region_NSB, 
                                         color="region", label='depth')
NMDS_region_plot_NSB + geom_point(size=1) + theme_bw()

##### envfit depth
ef_NSB <- envfit(NMDS_region_NSB, NMDS_META2[48:62,c(3:6,8,10)], permu=999, 
                  na.rm = TRUE)
ef_NSB

vectors_reg<-as.data.frame(ef_NSB$vectors$arrows*sqrt(ef_NSB$vectors$r))
pvals_reg=ef_NSB$vectors$pvals
r_reg=ef_NSB$vectors$r
envfit_NSB<-cbind(vectors_reg,pvals_reg, r_reg)
envfit_NSB
#keep only those vectors which significantly correlated to structure
envfit_NSB <- envfit_NSB[c(1:2,5),]
#vectors are too long to fit on plot, transform with /2
envfit_NSB2 <- cbind(envfit_NSB[,1:2]/3, envfit_NSB[,3:4]/3)
envfit_NSB2

NMDS_region_pub_NSB = data.frame(MDS1 = NMDS_region_NSB$points[,1],
                                  MDS2 = NMDS_region_NSB$points[,2])
NMDS_region_pub_NSB <- cbind(NMDS_region_pub_NSB,NMDS_META[48:62,1:3])
NMDS_region_pub_NSB
glNMDSplotReg_NSB <- ggplot(NMDS_region_pub_NSB, aes(x=NMDS_region_pub_NSB$MDS1,
                                                       y=NMDS_region_pub_NSB$MDS2)) +
  geom_point(aes(fill = factor(NMDS_region_pub_NSB$wetland),
                 shape = NMDS_region_pub_NSB$depth,
                 color = NMDS_region_pub_NSB$wetland),
             size=5) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(legend.direction="vertical",
        legend.position="right") +
  scale_shape_manual(name="Depth",
                     labels=c("Top", "Middle", "Bottom"),
                     values = c(21,22,24)) +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Site") +
  scale_fill_brewer(type = "qual", palette = "Set1", guide = FALSE) +
  geom_text(label="Stress = 0.012", x = -0.2, y = -0.4, size = 4) +
  coord_equal(ylim = 0, xlim = 0) +
  guides(fill = FALSE, color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  geom_text(label="d. NSB", x = -0.38, y = 0.45, size = 4) +
  theme_classic()# +  
#guides(fill = guide_legend(override.aes= list(colour = c("white","grey","black")))) +
#theme(legend.key = element_rect(fill = "grey90"))
#glNMDSplotReg <- glNMDSplotReg + geom_text_repel(aes(label = NMDS_region_pub$depth))
glNMDSplotReg_NSB <- glNMDSplotReg_NSB + 
  geom_segment(data = envfit_NSB2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=envfit_NSB2, aes(x=NMDS1*1.05, y=NMDS2*1.05,
                                   label=rownames(envfit_NSB2)),size=4)
glNMDSplotReg_NSB

############
### WSB ###
############

###Create a distance matrix 
region_Dist_WSB <- phyloseq::distance(vst.otu.region[,63:74], method= 'bray')

##NMDS of all samples
NMDS_region_WSB <- ordinate(vst.otu.region[,63:74], method= 'NMDS', 
                             distance = region_Dist_WSB, formula = NULL)

## plot quick NMDS for comparisons
NMDS_region_plot_WSB <- plot_ordination(ALL, NMDS_region_WSB, 
                                         color="region", label='depth')
NMDS_region_plot_WSB + geom_point(size=1) + theme_bw()

##### envfit depth
ef_WSB <- envfit(NMDS_region_WSB, NMDS_META2[63:74,c(3:6,8,10)], permu=999, 
                  na.rm = TRUE)
ef_WSB

vectors_reg<-as.data.frame(ef_WSB$vectors$arrows*sqrt(ef_WSB$vectors$r))
pvals_reg=ef_WSB$vectors$pvals
r_reg=ef_WSB$vectors$r
envfit_WSB<-cbind(vectors_reg,pvals_reg, r_reg)
envfit_WSB
#keep only those vectors which significantly correlated to structure
envfit_WSB <- envfit_WSB[c(1:2,5),]
#vectors are too long to fit on plot, transform with /2
envfit_WSB2 <- cbind(envfit_WSB[,1:2]/3, envfit_WSB[,3:4]/3)
envfit_WSB2

NMDS_region_pub_WSB = data.frame(MDS1 = NMDS_region_WSB$points[,1],
                                  MDS2 = NMDS_region_WSB$points[,2])
NMDS_region_pub_WSB <- cbind(NMDS_region_pub_WSB,NMDS_META[63:74,1:3])
NMDS_region_pub_WSB
glNMDSplotReg_WSB <- ggplot(NMDS_region_pub_WSB, aes(x=NMDS_region_pub_WSB$MDS1,
                                                       y=NMDS_region_pub_WSB$MDS2)) +
  geom_point(aes(fill = factor(NMDS_region_pub_WSB$wetland),
                 shape = NMDS_region_pub_WSB$depth,
                 color = NMDS_region_pub_WSB$wetland),
             size=5) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(legend.direction="vertical",
        legend.position="right") +
  scale_shape_manual(name="Depth",
                     labels=c("Top", "Middle", "Bottom"),
                     values = c(21,22,24)) +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Site") +
  scale_fill_brewer(type = "qual", palette = "Set1", guide = FALSE) +
  geom_text(label="Stress = 0.062", x = -0.2, y = -0.4, size = 4) +
  coord_equal(ylim = 0, xlim = 0) +
  guides(fill = FALSE, color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  geom_text(label="e. WSB", x = -0.38, y = 0.45, size = 4) +
  theme_classic()# +  
#guides(fill = guide_legend(override.aes= list(colour = c("white","grey","black")))) +
#theme(legend.key = element_rect(fill = "grey90"))
#glNMDSplotReg <- glNMDSplotReg + geom_text_repel(aes(label = NMDS_region_pub$depth))
glNMDSplotReg_WSB <- glNMDSplotReg_WSB + 
  geom_segment(data = envfit_WSB2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=envfit_WSB2, aes(x=NMDS1*1.05, y=NMDS2*1.05,
                                   label=rownames(envfit_WSB2)),size=4)
glNMDSplotReg_WSB

#################
### Comboplot ###
#################

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

RegionNMDSplot <- grid.arrange(arrangeGrob(glNMDSplotReg_BA,
                                       glNMDSplotReg_ESBT,
                                       glNMDSplotReg_LE,
                                       glNMDSplotReg_NSB,
                                       glNMDSplotReg_WSB,
                                       nrow=3))
##############
### ADONIS ###
##############

#BA
#perMANOVA
NMDS_META_BA <- NMDS_META[1:9,]
BA_adonis <- adonis(region_Dist_BA ~ depth, data = NMDS_META_BA)
BA_adonis


#ESBT
#perMANOVA
ESBT_adonis <- adonis(region_Dist_ESBT ~ depth*wetland, data = NMDS_META[10:35,])
ESBT_adonis

#LE
#perMANOVA
LE_adonis <- adonis(region_Dist_LE ~ depth*wetland, data = NMDS_META[36:47,])
LE_adonis

#NSB
#perMANOVA
NSB_adonis <- adonis(region_Dist_NSB ~ depth*wetland, data = NMDS_META[48:62,])
NSB_adonis

#WSB
#perMANOVA
WSB_adonis <- adonis(region_Dist_WSB ~ depth*wetland, data = NMDS_META[63:74,])
WSB_adonis
