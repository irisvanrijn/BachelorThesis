#----------------------------------Clean memory---------------------------------
rm(list=ls(all=TRUE))

#---------------------------Setting working directory---------------------------
setwd("C:/Users/iriss/Documents/studiegerelateerd/BachelorScriptie_2021")
getwd()


#-------------------------------Loading libraries-------------------------------

#-----------Reading population quantities on which PCA is performed-------------

PopQs<-read_excel("PCA_LOW.xlsx")
# Reading in an excel file with in the first column the species names: make sure that 
# in your Excel file tthere is an underscore _ between the genus and species name (e.g. Cephalopholis_fulva)
# make sure that the names match those in the tree. You can check this by comparing
# smalltree$tip.label with PopQs$Species . This is also done below in the code.
# the other colums are the life history traits (DEB-IPM pars) 
# and derived life history traits (lambda, R0, generation time)

################################################################################
############################ PART A - BUILD THE TREE ###########################
################################################################################

#---------------------------Loading taxa to create tree-------------------------

taxa <- tnrs_match_names(names = c("Caretta caretta",
                                    "Chelydra serpentina",
                                    "Cottus gobio",
                                    "Crocodylus johnsoni",
                                    "Dermochelys coriacea",
                                    "Gasterosteus aculeatus",
                                    "Gehyra variegata",
                                    "Kinosternon flavescens",
                                    "Lamna nasus",
                                    "Natrix natrix",
                                    "Prionace glauca",
                                    "Sceloporus undulatus",
                                    "Sphenodon punctatus",
                                    "Thamnophis elegans",
                                    "Thunnus orientalis",
                                    "Tiliqua rugosa"))

#---------------------------------Small tree------------------------------------

smalltree <- tol_induced_subtree(ott_ids = ott_id(taxa), label_format = "name")
plot(smalltree, cex = .8, label.offset = .1, no.margin = TRUE)

# set labels
smalltree$node.label<-as.character(1:length(smalltree$node.label))

# Drop all the tips in the tree that are not in the data
drop<-smalltree$tip.label[which(!smalltree$tip.label%in%PopQs$Species)]

#Should be TRUE
length(drop)+length(unique(PopQs$Species))==length(smalltree$tip.label)

# Trim tree of not needed species, if any at all
smalltree<-drop.tip(smalltree, drop)

# sorting required for PCA
all(sort(PopQs$Species)==sort(smalltree$tip.label))

#-----------------------------Editing the small tree----------------------------

tree=smalltree
tree$node.label<-as.character(1:length(tree$node.label))

#The tree must have branch lengths to be able to build the phylogenetic 
#variance-covariance matrix into the function phyl.pca for the phylogenetically-informed PCA
tree <- compute.brlen(tree)

# To use the tree in an MCMC analysis, it needs to be rooted
# By default the tree should be unrooted, but we can root it
# using the command root or root.tree
tree <- root(tree, tree$tip.label[1])

tree <- multi2di(tree, random=FALSE)  
is.rooted(tree)   #This is fo checking if the tree became rooted 

any(duplicated(tree$node.label))   
tree$node.label <- unique(tree$node.label)

#Another possible issue could be that the distance between taxa
#Could be 0, and this is not enabled for some analyses (including MCMCglmm)
#Thus we can solve this by adding an 0.001 to all tips
tree$edge.length <-  tree$edge.length + 0.001

#-------------------Transform the tree into an ultametric one-------------------

is.ultrametric(tree) 
tree <- chronos(tree, lambda=0, model="correlated") 

# JOSJE: At this point in the Paniw et al. 2018 code you would add populations as extra branches 

# continue
is.ultrametric(tree)
class(tree) <- "phylo"

any(duplicated(tree$node.label))
tree$node.label <- unique(tree$node.label)

#-----------------------Plot the final, ultametric tree-------------------------

plot(tree, cex = .8, label.offset = .1, no.margin = TRUE)


################################################################################
#################### PART B PHYLOGENETICALLY - INFORMED PCA ####################
################################################################################


#-----------------------------Phylo PCA-----------------------------------------

row.names(PopQs)=PopQs$Species #This shows where the row names can be found

#Standardizing the data
pcaData1=log(PopQs[,c(2:4)])  # identify columns that will be used in PCA (lambda, R0, GT, mu etc.)
pcaData=scale(pcaData1) # here non-standardised data are standardised.

summary(pcaData)

#Checking the means and standard dev
colMeans(pcaData) # should be (almost) zero
apply(pcaData, 2, sd)  # should be (almost) one


#PCA
pca=phyl.pca(tree,pcaData,method="lambda",mode="corr") 

#Some results of the PCA
summary(pca)
pca$lambda # Pagel's lambda: informs on role of phylogeny in data output
pca$logL

#-------------------PCA is ready at this point----------------------------------

diag(pca$Eval)^2 # Apply Kaiser's criterion to determine how many PCA axes to keep

ncomp=2 # keep the first 2 axes

# Varimax correction on the axes
rawLoadings <- pca$L[,1:ncomp]

# find the variance maximizing rotation of loadings      
rotatedLoadings <- varimax(rawLoadings)$loadings

name1=0.604   # PC1 axis: fill this in for your pca output: summary(pca)! fraction variance explained by first axis; this output should be copied from summary(pca) 
name2=0.3238   # PC2 axis: fill this in for your pca output: summary(pca)! fraction variance explained by second axis;this output should be copied from summary(pca)

# Create the inverse of the loading to calculate new scores: data multiplied by rotation matrix 
invLoadings     <- t(pracma::pinv(rotatedLoadings))
scores          <- pcaData %*% invLoadings


x <- list() 
x$scores <- scores[,1:2]
colnames(x$scores)=c("PC1", "PC2")
x$scores[,1:2] <- (-1)*x$scores[,1:2]
x$loadings <- rotatedLoadings[,1:2]
colnames(x$loadings)=c("PC1", "PC2")
x$loadings[,1:2] <- (-1)*x$loadings[,1:2]

#-----------------Save the PCA scores for regression analyses-------------------

write.csv(cbind(x$scores,PopQs$Species),"PCAscoresLOWnodeb.csv",row.names=F)

#----------------------------Pagel's Lambda-------------------------------------  

#Get standard deviation of Pagel's lambda by simple bootstrap
nsim=100
Plambda=rep(NA,nsim)

for(i in 1:nsim){
  
  pcaData1=pcaData[sample(1:nrow(pcaData),nrow(pcaData),replace=T),]
  # pcaData1=pcaData1[!duplicated(pcaData1), ] #switched off because it changes size dataframe
  pca1=phyl.pca(tree,pcaData1,method="lambda",mode="corr")
  Plambda[i]=pca1$lambda   
  
}


################################################################################
#################### PART C Plotting the result of the PCA #####################
################################################################################



### With point labels (species ID)
data1 <- data.frame(x$scores,point.lab=1:nrow(x$scores))

options(ggrepel.max.overlaps = Inf)


PCbiplot <- function(PC, x="PC1", y="PC2") {
  
  #position on the plot
  plot <- ggplot(data1, aes_string(x=x, y=y)) 
  
  #editing the data points
  plot <- plot + geom_point(aes(color = PopQs$Species), size=4.5, alpha=0.5)
  
  #adding lines to the plot
  plot <- plot + geom_hline(aes(yintercept=0), size=.2) + geom_vline(aes(xintercept=0), size=.2)
  
  plot <- plot + scale_size(range = c(0, 6))    
  
  datapc <- data.frame(varnames=colnames(PopQs[2:4]), PC$loadings) #    WARNING: adapt number of columns do your dataset                    
  
  mult <- min(
    (max(data1[,y]) - min(data1[,y])/(max(datapc[,y])-min(datapc[,y]))), 
    (max(data1[,x]) - min(data1[,x])/(max(datapc[,x])-min(datapc[,x])))  
  )
  
  datapc <- transform(datapc,
                      v1 = .50 * mult* (get(x)),    
                      v2 = .50 * mult* (get(y))      
  )
  
  #Adding the numbers and the title
  plot <- plot + geom_text_repel(aes(label=point.lab)) + labs(title = "Low feeding level") + theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"))
  
  #adding the arrows: you have more than three arrows, try adding them yourself or ask Josje!
  plot <- plot + geom_segment(data=datapc[1:1,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=1.5, color="purple")  
  plot <- plot + geom_segment(data=datapc[2:2,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=1.5, color="deepskyblue")  
  plot <- plot + geom_segment(data=datapc[3:3,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=1.5, color="orangered3")
 
  
  #Axis settings
  plot <- plot+ylab(paste("PCA 2 (",round(name2,2)*100,"%)",sep=""))+xlab(paste("PCA 1 (",round(name1,2)*100,"%)",sep=""))
  plot<- plot+theme(axis.text = element_text(size=14))+theme(axis.title = element_text(size=16))
  
  #Legend settings
  plot<- plot + labs(color="Species") + theme(legend.title = element_text(size=14, face="bold"), legend.text = element_text(size=14),legend.key.size = unit(1, "lines"))
  plot<- plot + scale_color_hue(labels = c("Caretta caretta","Chelydra serpentina","Cottus gobio","Crocodylus johnsoni","Dermochelys coriacea","Gasterosteus aculeatus","Gehyra variegata","Kinosternon flavescens","Lamna nasus","Natrix natrix","Prionace glauca","Sceloporus undulatus","Sphenodon punctatus","Thamnophis elegans","Thunnus orientalis","Tiliqua rugosa"))
  
  plot
}

PCbiplot(x)

