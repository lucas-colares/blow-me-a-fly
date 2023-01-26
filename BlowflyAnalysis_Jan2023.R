# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################ Script to conduct the analysis and graphs from: ################
#### Dispersal capacity regulates early arrival during ecological succession ####
############## By Lucas Colares, Anita Herdina and Cristian Dambros #############
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

install.packages(c("geomorph","vegan","FD","ggplot2","stringr","reshape2","ggpubr","RColorBrewer","ggConvexHull","rlist")) #Run only once
library(geomorph)
library(vegan)
library(FD)
library(ggplot2)
library(stringr)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(ggConvexHull)
library(rlist)

##### Create functions #####
aggmorpho<-function(coords,agg,dim1=nrow(coords),dim2=ncol(coords)){
  
  species<-sort(levels(agg))
  
  resu<-lapply(species,function(x){
    test<-coords[,,agg==x]
    apply(test,c(1,2),mean)
  }
  )
  
  resu2<-array(do.call(c,resu),c(dim1,dim2,length(species)))
  
  dimnames(resu2)<-list(NULL,NULL,species)
  
  resu2
}  

#### Read species occurrence data ####

#Importing attribute table (sex, code, etc.)
attribute<-read.csv("specimenID.csv",h=T)
BlowOcc<-read.csv("Blow_Occ.csv", h=T,sep="\t")
head(attribute)
head(BlowOcc)

#### Plot mean abundance variation across phases
data.frame(t(BlowOcc))->BlowOcc
BlowOcc$phase<-substr(rownames(BlowOcc),1,(nchar(rownames(BlowOcc))-1))
aggregate(BlowOcc[,1:(ncol(BlowOcc)-1)],list(BlowOcc$phase),mean)->MeanAbun
melt(MeanAbun)->MeanAbun
spline.d <- as.data.frame(spline(MeanAbun$value))
MeanAbun$variable<-gsub("_"," ",MeanAbun$variable)

AbunPlot<-ggplot(data=MeanAbun, aes(x=factor(Group.1,levels=c("fresh","bloat","active.decay","advance.decay","dry")), y=value, group=variable)) +
  geom_smooth(linetype="dashed",aes(color=variable))+
  geom_point(size=3, aes(color=variable))+
  theme_bw()+
  scale_x_discrete(labels=c("Fresh","Bloat","ActD","AdvD","Dry"), expand=c(0.05, 0))+
  scale_color_manual(values = rev(brewer.pal(n=11,"Spectral")))+
  labs(x="Carrion succession phase",y="Mean abundance",color="Species")+
  theme(legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=14), legend.text = element_text(face="italic",family ="sans", size=12), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); AbunPlot
tiff("Fig1.tiff", units="in", width=7, height=3.5, res=300) #Save plot
plot(AbunPlot)
dev.off()
# Rearrange levels of phase to be in succesion order (not alphabetical) 
attribute$Phase<-factor(attribute$Phase,levels=c("fresh","bloat","active decay","advance decay","dry"))

# Use code ID in the order it appears
attribute$CodeID<-factor(attribute$CodeID,levels=attribute$CodeID)

# Remove males from attribute table
attribute<-droplevels(attribute[attribute$Sex=="F",])

##### Read morphometric data ####
# Load tps and table
# Wing
wing<-readland.tps("wingFixed2.TPS",specID = "image")
# Torax (dorsal)
dorsal<-readland.tps("toraxdorsalFixed2.TPS",specID = "image")
# Torax (lateral)
lateral<-readland.tps("toraxlateralFixed2.TPS",specID = "image")

# Reordering tps tables to match codes in attribute table
wing<-wing[,,match(attribute$CodeID,dimnames(wing)[[3]])]
dorsal<-dorsal[,,match(attribute$CodeID,dimnames(dorsal)[[3]])]
lateral<-lateral[,,match(attribute$CodeID,dimnames(lateral)[[3]])]

### MORPHOMETRIC ANALYSES #####

# Fill missing values
wing<-estimate.missing(wing)
dorsal<-estimate.missing(dorsal)
lateral<-estimate.missing(lateral)

### Individual level ####

# Torax height
# Torax height is the distance between landmarks 1 and 4
# This code replicates the measurement for all specimens
height<-apply(lateral,3,function(x)dist(x[c(1,4),]))

# Shape and Area (Procrustes)

# Wing
wing_gpa<-gpagen(wing,curves = NULL)
wing_matrix<-two.d.array(wing_gpa$coords, sep = ".")
wing_matrix_array<-arrayspecs(wing_matrix, p=17, k=2, sep = NULL)
wing_pcoa<-plotTangentSpace(wing_matrix_array, axis1 = 1, axis2 = 2,groups = attribute$Phase, warpgrids=T)

scores_wing<-wing_pcoa$pc.scores[,1]

gsub("fresh",1,attribute$Phase)->attribute$Phase
gsub("bloat",2,attribute$Phase)->attribute$Phase
gsub("active decay",3,attribute$Phase)->attribute$Phase
gsub("advance decay",4,attribute$Phase)->attribute$Phase
gsub("dry",5,attribute$Phase)->attribute$Phase

rownames(wing_pcoa$pc.scores)<-attribute$Phase

# Dorsal
dorsal_gpa<-gpagen(dorsal,curves=NULL)
dorsal_matrix<-two.d.array(dorsal_gpa$coords, sep = ".")
dorsal_matrix_array<-arrayspecs(dorsal_matrix, p=8, k=2, sep = NULL)
dorsal_pcoa<-plotTangentSpace(dorsal_matrix_array, axis1 = 1, axis2 = 2,groups = attribute$Phase)

scores_dorsal<-dorsal_pcoa$pc.scores[,1]
rownames(dorsal_pcoa$pc.scores)<-attribute$Phase

# Lateral
lateral_gpa<-gpagen(lateral,curves=NULL)
lateral_matrix<-two.d.array(lateral_gpa$coords, sep = ".")
lateral_matrix_array<-arrayspecs(lateral_matrix, p=6, k=2, sep = NULL)
lateral_pcoa<-plotTangentSpace(lateral_matrix_array, axis1 = 1, axis2 = 2,groups = attribute$Phase)

scores_lateral<-lateral_pcoa$pc.scores[,1]
rownames(lateral_pcoa$pc.scores)<-attribute$Phase

## Include shape (pcoa scores) in attribute table
attribute["height"]<-height
attribute["scores_wing"]<-scores_wing
attribute["scores_dorsal"]<-scores_dorsal
attribute["scores_lateral"]<-scores_lateral

## Include Area (Csize) in attribute table
attribute["csize_wing"]<-wing_gpa$Csize
attribute["csize_dorsal"]<-dorsal_gpa$Csize
attribute["csize_lateral"]<-lateral_gpa$Csize

####### Calculate derivated attributes (Volume and Ratio)
tvolume = attribute$csize_dorsal*height
wtratio = attribute$csize_wing/tvolume

attribute["tvolume"]<-tvolume
attribute["wtratio"]<-wtratio

### Species level ####
## Aggregate morphometrics by species

wing_matrix_array_spp<-aggmorpho(wing_matrix_array,factor(attribute$spPhase))
wing_pcoa_spp<-plotTangentSpace(wing_matrix_array_spp, axis1 = 1, axis2 = 2)

scores_wing_spp<-wing_pcoa_spp$pc.scores[,1]
names(scores_wing_spp)<-unlist(str_extract_all(names(scores_wing_spp), '[A-Z]{2,}'))

rownames(wing_pcoa_spp$pc.scores)[grep("fresh", rownames(wing_pcoa_spp$pc.scores))]<-1
rownames(wing_pcoa_spp$pc.scores)[grep("bloat", rownames(wing_pcoa_spp$pc.scores))]<-2
rownames(wing_pcoa_spp$pc.scores)[grep("active", rownames(wing_pcoa_spp$pc.scores))]<-3
rownames(wing_pcoa_spp$pc.scores)[grep("advance", rownames(wing_pcoa_spp$pc.scores))]<-4
rownames(wing_pcoa_spp$pc.scores)[grep("dry", rownames(wing_pcoa_spp$pc.scores))]<-5

dorsal_matrix_array_spp<-aggmorpho(dorsal_matrix_array,factor(attribute$spPhase))
dorsal_pcoa_spp<-plotTangentSpace(dorsal_matrix_array_spp, axis1 = 1, axis2 = 2)

scores_dorsal_spp<-dorsal_pcoa_spp$pc.scores[,1]
names(scores_dorsal_spp)<-unlist(str_extract_all(names(scores_dorsal_spp), '[A-Z]{2,}'))

rownames(dorsal_pcoa_spp$pc.scores)[grep("fresh", rownames(dorsal_pcoa_spp$pc.scores))]<-1
rownames(dorsal_pcoa_spp$pc.scores)[grep("bloat", rownames(dorsal_pcoa_spp$pc.scores))]<-2
rownames(dorsal_pcoa_spp$pc.scores)[grep("active", rownames(dorsal_pcoa_spp$pc.scores))]<-3
rownames(dorsal_pcoa_spp$pc.scores)[grep("advance", rownames(dorsal_pcoa_spp$pc.scores))]<-4
rownames(dorsal_pcoa_spp$pc.scores)[grep("dry", rownames(dorsal_pcoa_spp$pc.scores))]<-5

lateral_matrix_array_spp<-aggmorpho(lateral_matrix_array,factor(attribute$spPhase))
lateral_pcoa_spp<-plotTangentSpace(lateral_matrix_array_spp, axis1 = 1, axis2 = 2)

scores_lateral_spp<-lateral_pcoa_spp$pc.scores[,1]
names(scores_lateral_spp)<-unlist(str_extract_all(names(scores_lateral_spp), '[A-Z]{2,}'))

rownames(lateral_pcoa_spp$pc.scores)[grep("fresh", rownames(lateral_pcoa_spp$pc.scores))]<-1
rownames(lateral_pcoa_spp$pc.scores)[grep("bloat", rownames(lateral_pcoa_spp$pc.scores))]<-2
rownames(lateral_pcoa_spp$pc.scores)[grep("active", rownames(lateral_pcoa_spp$pc.scores))]<-3
rownames(lateral_pcoa_spp$pc.scores)[grep("advance", rownames(lateral_pcoa_spp$pc.scores))]<-4
rownames(lateral_pcoa_spp$pc.scores)[grep("dry", rownames(lateral_pcoa_spp$pc.scores))]<-5

csize_wing_spp = tapply(wing_gpa$Csize,attribute$SpeciesID,mean)
csize_dorsal_spp = tapply(dorsal_gpa$Csize,attribute$SpeciesID,mean)
csize_lateral_spp = tapply(lateral_gpa$Csize,attribute$SpeciesID,mean)

height_spp = tapply(height,attribute$SpeciesID,mean)
tvolume_spp = tapply(tvolume,attribute$SpeciesID,mean)
wtratio_spp = tapply(wtratio,attribute$SpeciesID,mean)

## Include shape (pcoa scores) in attribute table
attribute["height_spp"]<-height_spp[attribute$SpeciesID]
attribute["scores_wing_spp"]<-scores_wing_spp[attribute$SpeciesID]
attribute["scores_dorsal_spp"]<-scores_dorsal_spp[attribute$SpeciesID]
attribute["scores_lateral_spp"]<-scores_lateral_spp[attribute$SpeciesID]

## Include Area (Csize) in attribute table
attribute["csize_wing_spp"]<-csize_wing_spp[attribute$SpeciesID]
attribute["csize_dorsal_spp"]<-csize_dorsal_spp[attribute$SpeciesID]
attribute["csize_lateral_spp"]<-csize_lateral_spp[attribute$SpeciesID]

attribute["tvolume_spp"]<-tvolume_spp[attribute$SpeciesID]
attribute["wtratio_spp"]<-wtratio_spp[attribute$SpeciesID]

##### Aggregate scatterplots for individuals and species #####
rbind(data.frame(wing_pcoa$pc.scores[,1:2],Phase=attribute$Phase,Trait="Wing",Level="Individual level"),
data.frame(dorsal_pcoa$pc.scores[,1:2],Phase=attribute$Phase,Trait="Dorsal",Level="Individual level"),
data.frame(lateral_pcoa$pc.scores[,1:2],Phase=attribute$Phase,Trait="Lateral",Level="Individual level"),
data.frame(wing_pcoa_spp$pc.scores[,1:2],Phase=rownames(wing_pcoa_spp$pc.scores),Trait="Wing",Level="Species level"),
data.frame(dorsal_pcoa_spp$pc.scores[,1:2],Phase=rownames(dorsal_pcoa_spp$pc.scores),Trait="Dorsal",Level="Species level"),
data.frame(lateral_pcoa_spp$pc.scores[,1:2],Phase=rownames(lateral_pcoa_spp$pc.scores),Trait="Lateral",Level="Species level"))->AllPCoAs

reds<-colorRampPalette(c("#a50f15","#fee5d9"))

WingPcoas<-ggplot(AllPCoAs[AllPCoAs$Trait=="Wing",], aes(x=PC1, y=PC2))+ 
  geom_point(aes(fill=Phase),shape=21, color="black", size=3, stroke = 0.25)+
  theme_bw()+
  facet_grid(. ~ Level)+
  scale_fill_manual(values=reds(5),labels=c("Fresh","Bloat","ActD","AdvD","Dry"))+
  labs(x="Axis 1", y="Axis 2", fill='Phase')+ #Mudar o título da legenda
  theme(strip.text = element_text(size=12,family="sans",face = "bold"),strip.background = element_rect(fill="white"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); WingPcoas

DorsalPcoas<-ggplot(AllPCoAs[AllPCoAs$Trait=="Dorsal",], aes(x=PC1, y=PC2))+ 
  geom_point(aes(fill=Phase),shape=21, color="black", size=3, stroke = 0.25)+
  theme_bw()+
  facet_grid(. ~ Level)+
  scale_fill_manual(values=reds(5),labels=c("Fresh","Bloat","ActD","AdvD","Dry"))+
  labs(x="Axis 1", y="Axis 2", fill='Phase')+ #Mudar o título da legenda
  theme(strip.text = element_text(size=12,family="sans",face = "bold"),strip.background = element_rect(fill="white"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); DorsalPcoas

LateralPcoas<-ggplot(AllPCoAs[AllPCoAs$Trait=="Lateral",], aes(x=PC1, y=PC2))+ 
  geom_point(aes(fill=Phase),shape=21, color="black", size=3, stroke = 0.25)+
  theme_bw()+
  facet_grid(. ~ Level)+
  scale_fill_manual(values=reds(5),labels=c("Fresh","Bloat","ActD","AdvD","Dry"))+
  labs(x="Axis 1", y="Axis 2", fill='Phase')+ #Mudar o título da legenda
  theme(strip.text = element_text(size=12,family="sans",face = "bold"),strip.background = element_rect(fill="white"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); LateralPcoas

###### Check differences in the traits between species and phases ####
### Differences between individuals
### Plot traits over succession

vars<-c("scores_wing","scores_dorsal","scores_lateral","scores_wing_spp","scores_dorsal_spp","scores_lateral_spp")
melt(attribute[,c("Phase",vars)])->PCoA_Scores4Shape
PCoA_Scores4Shape$variable<-as.character(PCoA_Scores4Shape$variable)
PCoA_Scores4Shape$Trait<-as.character(PCoA_Scores4Shape$variable)
PCoA_Scores4Shape[grep("wing",PCoA_Scores4Shape$Trait),4]<-"Wing"
PCoA_Scores4Shape[grep("dorsal",PCoA_Scores4Shape$Trait),4]<-"Dorsal"
PCoA_Scores4Shape[grep("lateral",PCoA_Scores4Shape$Trait),4]<-"Lateral"
PCoA_Scores4Shape[-grep("spp",PCoA_Scores4Shape$variable),2]="Individual level"
PCoA_Scores4Shape[grep("spp",PCoA_Scores4Shape$variable),2]="Species level"

WingBoxs<-ggplot(PCoA_Scores4Shape[PCoA_Scores4Shape$Trait=="Wing",], aes(x=Phase, y=value))+ 
  geom_boxplot(aes(fill=Phase))+
  theme_bw()+
  facet_grid(. ~ variable)+
  scale_fill_manual(values=reds(5),labels=c("Fresh","Bloat","ActD","AdvD","Dry"))+
  labs(x="Succession phase", y="Wing shape", fill='Phase')+ #Mudar o título da legenda
  theme(strip.text = element_text(size=12,family="sans",face = "bold"),strip.background = element_rect(fill="white"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); WingBoxs

DorsalBoxs<-ggplot(PCoA_Scores4Shape[PCoA_Scores4Shape$Trait=="Dorsal",], aes(x=Phase, y=value))+ 
  geom_boxplot(aes(fill=Phase))+
  theme_bw()+
  facet_grid(. ~ variable)+
  scale_fill_manual(values=reds(5),labels=c("Fresh","Bloat","ActD","AdvD","Dry"))+
  labs(x="Succession phase", y="Thorax shape (dorsal)", fill='Phase')+ #Mudar o título da legenda
  theme(strip.text = element_text(size=12,family="sans",face = "bold"),strip.background = element_rect(fill="white"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); DorsalBoxs

LateralBoxs<-ggplot(PCoA_Scores4Shape[PCoA_Scores4Shape$Trait=="Lateral",], aes(x=Phase, y=value))+ 
  geom_boxplot(aes(fill=Phase))+
  theme_bw()+
  facet_grid(. ~ variable)+
  scale_fill_manual(values=reds(5),labels=c("Fresh","Bloat","ActD","AdvD","Dry"))+
  labs(x="Succession phase", y="Thorax shape (lateral)", fill='Phase')+ #Mudar o título da legenda
  theme(strip.text = element_text(size=12,family="sans",face = "bold"),strip.background = element_rect(fill="white"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); LateralBoxs

### Aggregate boxplots and scatterplots that represent shape of species and individuals
BoxsAndSc<-ggarrange(WingBoxs,WingPcoas,DorsalBoxs,DorsalPcoas,LateralBoxs,LateralPcoas,nrow=3,ncol=2,labels = c("(a)","(b)","(c)","(d)","(e)","(f)"), legend = "none",widths = c(1,2)); BoxsAndSc
tiff("Fig2.tiff", units="in", width=11, height=3*3, res=300) #Save plot
plot(BoxsAndSc)
dev.off()

## Plot boxplots for the remaining traits
vars2<-c("csize_wing"  ,"csize_dorsal","csize_lateral","tvolume","wtratio","height","csize_wing_spp"  ,"csize_dorsal_spp","csize_lateral_spp","tvolume_spp","wtratio_spp","height_spp")
melt(attribute[,c("Phase",vars2)])->OtherTraits
OtherTraits$variable<-as.character(OtherTraits$variable)
OtherTraits$Trait<-as.character(OtherTraits$variable)
OtherTraits$Trait<-gsub("_spp","",OtherTraits$Trait)
OtherTraits[-grep("spp",OtherTraits$variable),2]="Individual level"
OtherTraits[grep("spp",OtherTraits$variable),2]="Species level"

trait.labs <- c("Wing size", "Thorax size (dorsal)", "Thorax size (lateral)","Thorax volume","Wing-thorax ratio","Thorax height")
names(trait.labs) <- unique(OtherTraits$Trait)

OtherBoxs<-ggplot(OtherTraits, aes(x=Phase, y=value))+ 
  geom_boxplot(aes(fill=Phase))+
  theme_bw()+
  facet_grid(Trait ~ variable,scales="free_y",labeller = labeller(Trait = trait.labs))+
  scale_fill_manual(values=reds(5),labels=c("Fresh","Bloat","ActD","AdvD","Dry"))+
  labs(x="Succession phase", y="Functional trait", fill='Phase')+ #Mudar o título da legenda
  theme(strip.text = element_text(size=12,family="sans",face = "bold"),strip.background = element_rect(fill="white"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"));OtherBoxs
tiff("FigS5.tiff", units="in", width=8, height=12, res=300) #Save plot
plot(OtherBoxs)
dev.off()

### Test differences in the shape of each trait
rbind(PCoA_Scores4Shape,OtherTraits)->AllTraitsMelted
traits.lab<-unique(AllTraitsMelted$Trait)
names(traits.lab)<-c("Wing shape","Thorax shape (dorsal)","Thorax shape (lateral)", "Wing size","Thorax size (dorsal)","Thorax size (lateral)","Thorax volume","Wing-thorax ratio","Thorax height")

ANOVAS<-{}
TUKEYS<-{}
for(x in 1:length(traits.lab)){
  for (y in unique(AllTraitsMelted$variable)) {
    AllTraitsMelted[AllTraitsMelted$Trait==traits.lab[x],]->OnlyTrait
    OnlyTrait[OnlyTrait$variable==y,]->SelTrait
    print(paste0("Level: ",y," // Trait: ",names(traits.lab)[x]))
    aov(value~Phase,data=SelTrait)->AOV1
    as.matrix(summary(AOV1))->SUMMAOV
    print(SUMMAOV)
    ANOVAS[[length(ANOVAS)+1]]<-data.frame(SUMMAOV[[1]],Level=y,Trait=names(traits.lab)[x])
    if(SUMMAOV[[1]][,5][1]<0.05){
      TukeyHSD(AOV1)->TUKEY1
      print(TUKEY1)
      TUKEYS[[length(TUKEYS)+1]]<-data.frame(Comparison=rownames(TUKEY1$Phase),TUKEY1$Phase,Level=y,Trait=names(traits.lab)[x])
    }
  }
}

data.frame(do.call("rbind",ANOVAS))->ANOVAS
data.frame(do.call("rbind",TUKEYS))->TUKEYS
write.csv(ANOVAS,"ANOVAS4Traits.csv") #Save matrix with anova results
write.csv(TUKEYS,"TUKEYS4Traits.csv") #Save matrix with tukey results

### Differences between species
varz<-c("scores_wing","scores_dorsal","scores_lateral","csize_wing"  ,"csize_dorsal","csize_lateral","tvolume","wtratio","height")
melt(attribute[,c("SpeciesID",varz)])->SppTraits

traits.lab<-c("Wing shape","Thorax shape (dorsal)","Thorax shape (lateral)", "Wing size","Thorax size (dorsal)","Thorax size (lateral)","Thorax volume","Wing-thorax ratio","Thorax height")
names(traits.lab)<-unique(SppTraits$variable)

SppBxs<-ggplot(SppTraits, aes(x=factor(SpeciesID,levels=c("CHA","CHM","CHP","CM","SC","LS","CL","PX","HSD","LE","HSG")), y=value))+ 
  geom_boxplot()+
  theme_bw()+
  facet_wrap(variable ~ .,scales="free_y",labeller = labeller(variable = traits.lab))+
  labs(x="Species", y="Functional trait", fill='Phase')+ #Mudar o título da legenda
  theme(strip.text = element_text(size=12,family="sans",face = "bold"),strip.background = element_rect(fill="white"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(angle = 90,family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); SppBxs
tiff("FigS4.tiff", units="in", width=10, height=6, res=300) #Save plot
plot(SppBxs)
dev.off()

ANOVAS_spp<-{}
TUKEYS_spp<-{}
for(x in 1:length(traits.lab)){
  SppTraits[SppTraits$variable==names(traits.lab)[x],]->OnlyTrait
  print(paste0("Trait: ",traits.lab[x]))
  aov(value~SpeciesID,data=OnlyTrait)->AOV2
  as.matrix(summary(AOV2))->SUMMAOV
  print(SUMMAOV)
  ANOVAS_spp[[length(ANOVAS_spp)+1]]<-data.frame(SUMMAOV[[1]],Trait=names(traits.lab)[x])
  if(SUMMAOV[[1]][,5][1]<0.05){
    TukeyHSD(AOV2)->TUKEY2
    print(TUKEY2)
    TUKEYS_spp[[length(TUKEYS_spp)+1]]<-data.frame(Comparison=rownames(TUKEY2$SpeciesID),TUKEY2$SpeciesID,Trait=names(traits.lab)[x])
  }
}

data.frame(do.call("rbind",ANOVAS_spp))->ANOVAS_spp
data.frame(do.call("rbind",TUKEYS_spp))->TUKEYS_spp
write.csv(ANOVAS_spp,"ANOVAS4Spp.csv") #Save matrix with anova results
write.csv(TUKEYS_spp,"TUKEYS4Spp.csv") #Save matrix with tukey results

###### Plot functional space #####
# Select individuals traits
vars<-c("scores_wing","scores_dorsal","scores_lateral","csize_wing"  ,"csize_dorsal","csize_lateral","tvolume","wtratio","height")

traits=attribute[,vars]
rownames(traits)<-attribute$CodeID

#####
### Arrange table with each Sampling unit as a row

# grouping data for each sampling unit
community<-data.frame(row.names = unique(attribute$Sunit)[order(unique(attribute$Sunit))])

community$Phase<-aggregate(attribute$Phase,list(attribute$Sunit),FUN=function(x)x[1])[,-1]

community[,c("TempMax","TempMin","Precipitacao")]<-aggregate(cbind(attribute$TempMax,attribute$TempMin,attribute$Precipitation),list(attribute$Sunit),FUN=mean)[,-1]

########

## Occurrence data for each sampling unit
attribute$N<-1
attach(attribute)

# Individual level
occ<-tapply(N,list(Sunit,CodeID),sum)
occ[is.na(occ)]<-0

detach(attribute)

#### Aggregate occurrence and trait data per species
occ_spp<-aggregate(t(occ),list(attribute$SpeciesID),sum)
traits_spp<-aggregate(traits,list(attribute$SpeciesID),mean)
rownames(occ_spp)<-occ_spp[,1]
rownames(traits_spp)<-traits_spp[,1]

occ_spp<-t(occ_spp[,-1])
traits_spp<-traits_spp[,-1]
###

#### Principal Component Analysis on trait matrix 
PCA=prcomp(scale(traits)) #PCA on trait distance matrix
head(PCA)

PCA_spp=prcomp(scale(traits_spp)) #PCA on trait distance matrix
head(PCA_spp)

##Calculate variance explained
eigenvalues<-eigenvals(PCA)
variance<-eigenvalues/sum(eigenvalues); variance

traits<-PCA$x[,1:2]
comm<-rbind(colSums(occ))

P<-{}
for(i in 1:nrow(comm)){
  A<-traits[comm[i,]>0,]
  P[[i]]<-A[chull(A),]
}

for (x in 1:length(P)) {
  data.frame(P[[x]],Level="Full")->P[[x]]  
}

list.rbind(P)->Full
factor(community$Phase)->community$Phase

AllPhases<-{}
for(y in 1:length(levels(community$Phase))){
  rbind(colSums(occ[community$Phase==levels(community$Phase)[y],]))->PhaseMean
  
  P<-{}
  for(i in 1:nrow(PhaseMean)){
    A<-traits[PhaseMean[i,]>0,]
    P[[i]]<-A[chull(A),]
  }
  
  for (x in 1:length(P)) {
    data.frame(P[[x]],Level=paste0("Phase",y))->P[[x]]  
  }
  
  list.rbind(P)->Phase
  comm<-occ[community$Phase==levels(community$Phase)[y],][rowSums(occ[community$Phase==levels(community$Phase)[y],])>1,]
  P<-{}
  for(i in 1:nrow(comm)){
    A<-traits[comm[i,]>0,]
    P[[i]]<-A[chull(A),]
  }
  
  for (x in 1:length(P)) {
    data.frame(P[[x]],Level=rownames(comm)[x])->P[[x]]  
  }
  
  list.rbind(P)->P
  
  rbind(Phase, P)->AllPhases[[y]]
  
}

rbind(Full,list.rbind(AllPhases))->AllPhases
data.frame(Phase=paste0("Phase",seq(1,5,1)), Label=c("fresh","bloat","active","advance","dry"))->Labels

# Species
traits<-PCA_spp$x[,1:2]
comm<-rbind(colSums(occ_spp))

P<-{}
for(i in 1:nrow(comm)){
  A<-traits[comm[i,]>0,]
  P[[i]]<-A[chull(A),]
}

for (x in 1:length(P)) {
  data.frame(P[[x]],Level="Full")->P[[x]]  
}

list.rbind(P)->Full_sp

AllPhasesSpp<-{}
for(y in 1:5){
  rbind(colSums(occ_spp[community$Phase==levels(community$Phase)[y],]))->comm
  
  P<-{}
  for(i in 1:nrow(comm)){
    A<-traits[comm[i,]>0,]
    P[[i]]<-A[chull(A),]
  }
  
  for (x in 1:length(P)) {
    data.frame(P[[x]],Level=paste0("Phase",y))->P[[x]]  
  }
  
  list.rbind(P)->Phase
  comm<-occ_spp[community$Phase==levels(community$Phase)[y],][rowSums(decostand(occ_spp[community$Phase==levels(community$Phase)[y],],"pa"))>1,]
  
  if(nrow(comm)==0){
    Phase->AllPhasesSpp[[y]]
    next
  } else {
    P<-{}
    for(i in 1:nrow(comm)){
      A<-traits[comm[i,]>0,]
      try(A[chull(A),])->TRY
      if(word(TRY)=="Error"){
        next
      } else {
        P[[i]]<-A[chull(A),]
      }}
    
    P<-P[!sapply(P,is.null)]
    
    for (x in 1:length(P)) {
      data.frame(P[[x]],Level=rownames(comm)[x])->P[[x]]  
    }
    
    list.rbind(P)->P
    
    rbind(Phase, P)->AllPhasesSpp[[y]]
    
  }
  
}

AllPhasesSpp
rbind(Full_sp,list.rbind(AllPhasesSpp))->AllPhasesSpp

rbind(data.frame(AllPhases[grep("Full|1|fresh",AllPhases$Level),],Phase="Fresh"),
      data.frame(AllPhases[grep("Full|2|bloat",AllPhases$Level),],Phase="Bloat"),
      data.frame(AllPhases[grep("Full|3|active",AllPhases$Level),],Phase="ActD"),
      data.frame(AllPhases[grep("Full|4|advance",AllPhases$Level),],Phase="AdvD"),
      data.frame(AllPhases[grep("Full|5|dry",AllPhases$Level),],Phase="Dry"))->FinalInd

rbind(data.frame(AllPhasesSpp[grep("Full|1|fresh",AllPhasesSpp$Level),],Phase="Fresh"),
      data.frame(AllPhasesSpp[grep("Full|2|bloat",AllPhasesSpp$Level),],Phase="Bloat"),
      data.frame(AllPhasesSpp[grep("Full|3|active",AllPhasesSpp$Level),],Phase="ActD"),
      data.frame(AllPhasesSpp[grep("Full|4|advance",AllPhasesSpp$Level),],Phase="AdvD"),
      data.frame(AllPhasesSpp[grep("Full|5|dry",AllPhasesSpp$Level),],Phase="Dry"))->FinalSpp

rbind(data.frame(FinalInd,TaxLevel="Ind"),data.frame(FinalSpp,TaxLevel="Spp"))->Ultimate

factor(Ultimate$Level,levels=c("Full",rev(levels(factor(Ultimate$Level)))[-6]))->Ultimate$Level

unique(levels(factor(Ultimate$Level)))->cols
cols[-grep("Full|Phase",cols)]<-"transparent"
cols[grep("Full",cols)]<-"black"
cols[grep("Phase",cols)]<-"#a50f15"

unique(levels(factor(Ultimate$Level)))->fillz
fillz[-grep("Full|Phase",fillz)]<-"#a50f15"
fillz[grep("Full",fillz)]<-"transparent"
fillz[grep("Phase",fillz)]<-"transparent"

as.character(Ultimate$Level)->alp
alp[-grep("Full|Phase",alp)]<-0.3
alp[grep("Full",alp)]<-0.1
alp[grep("Phase",alp)]<-0.4

level_lab<-c("Individual level", "Species level")
names(level_lab)=c("Ind","Spp")

FuncSpace<-ggplot(Ultimate, aes(x=PC1, y=PC2, fill=Level)) + 
  theme_bw()+
  facet_grid(TaxLevel ~ factor(Phase,levels=c("Fresh","Bloat","ActD","AdvD","Dry")), labeller = labeller(TaxLevel=level_lab))+
  geom_polygon(data = Ultimate, aes(fill = factor(Level), colour = factor(Level), alpha = factor(Level)),show.legend = FALSE, size=0.5,alpha=as.numeric(alp))+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=fillz)+
  guides(fill = guide_legend(override.aes = list(size = 5)))+
  labs(x="Axis 1", y="Axis 2", fill='Phase')+ #Mudar o título da legenda
  theme(strip.text = element_text(size=12,family="sans",face = "bold"),strip.background = element_rect(fill="white"),legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=14), legend.text = element_text(family ="sans", size=12), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); FuncSpace
tiff("Fig3.tiff", units="in", width=7, height=5, res=300) #Save plot
FuncSpace
dev.off()

##### Investigate the influence of succession phases and environmental variables in the traits ####
allParts<-{}
for(i in vars){
  varpart(attribute[,i],as.factor(attribute$Phase),attribute[,c("TempMax","TempMin","Season","Precipitation")])->varp
  as.data.frame(as.numeric(varp$part$indfract$Adj.R.square))->partTemp
  rownames(partTemp)=rownames(varp$part$indfract)
  colnames(partTemp)="Adj. R²"
  partTemp->allParts[[i]]
}

for (x in 1:length(allParts)) {
  cbind(allParts[[x]],names(allParts)[x])->allParts[[x]]
}

list.rbind(allParts)->allParts
cbind(allParts,rep(c("[a]", "[b]", "[a+b]", "Residuals"), length(vars)))->allParts
colnames(allParts)=c("Adj. R²", "Trait", "Part")
allParts[allParts$`Adj. R²`>0,]->allParts2

PartsPlot<-ggplot(allParts2, aes(fill=factor(Part, levels=c("Residuals","[b]", "[a+b]", "[a]")), y=as.numeric(`Adj. R²`), x=factor(Trait, levels=c("height", "wtratio", "tvolume","csize_lateral", "csize_dorsal", "csize_wing", "scores_lateral", "scores_dorsal", "scores_wing")))) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(labels=c("Residuals","[b] Environment", "[a]+[b] Intersect", "[a] - Succession phases"),values = c("grey33","#f0eb65","#14c99c","#00a4e6"))+
  scale_x_discrete(labels=c("Thorax height", "Wing-thorax ratio", "Thorax volume", 
                            "Thorax size (LT)", "Thorax size (DS)", "Wing size", 
                            "Thorax shape (LT)", "Thorax shape (DS)", "Wing shape"))+
  labs(x="Functional trait", y="Proportion of trait variance explained", fill="Predictors") + #Mudar o título da legenda
  theme_bw() + 
  coord_flip()+
  guides(fill=guide_legend(reverse=TRUE), color="none")+
  theme(legend.position = "right", legend.justification = "top",legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=14), legend.text = element_text(family ="sans", size=12), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); PartsPlot
tiff("Fig4.tiff", units="in", width=7, height=3, res=300)
plot(PartsPlot)
dev.off()

### Test exclusive proportions
matrix(NA, 4, 9)->AA
colnames(AA)=c("SHW", "SHTD", "SHTL", "SIW", "SITD", "SITL", "TV", "WTR", "TH")
rownames(AA)=c("[a]", "[a+b]", "[b]", "Residuals")

allParts[allParts$Part=="[a]",][,1]->AA[1,]
allParts[allParts$Part=="[a+b]",][,1]->AA[2,]
allParts[allParts$Part=="[b]",][,1]->AA[3,]
allParts[allParts$Part=="Residuals",][,1]->AA[4,]

data.frame(t(AA))->AA
AA$sig_a<-NA
AA$sig_b<-NA

for(i in 1:length(vars)){
  op <- par(mar = rep(0, 4))
  par("plt")  
  Phase<-as.factor(attribute$Phase)
  Environment<-cbind(attribute[,c("TempMax","TempMin","Season","Precipitation")])
  siga <- rda(attribute[,vars[i]],Phase,Environment)
  siga2 <- rda(attribute[,vars[i]],Environment,Phase)
  print(rev(c("Thorax height", "Wing-thorax ratio", "Thorax volume", 
          "Thorax size (LT)", "Thorax size (DS)", "Wing size", 
          "Thorax shape (LT)", "Thorax shape (DS)", "Wing shape"))[i])
  print(anova(siga, step=10000, perm.max=10000))
  print(anova(siga2, step=10000, perm.max=10000))
  as.matrix(anova(siga, step=10000, perm.max=10000)[4])[1,]->AA$sig_a[i]
  as.matrix(anova(siga2, step=10000, perm.max=10000)[4])[1,]->AA$sig_b[i]
}

write.table(AA, "AdjR&Sig4Parts.csv", sep="\t")