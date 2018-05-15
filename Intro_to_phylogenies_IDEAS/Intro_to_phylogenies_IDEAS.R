library(ape)

ftree <- read.tree("sample_newick_tree.txt")

ftree2 <- read.nexus("sample_NEXUS_tree.txt ")

plot(ftree2)

str(ftree2)

#drop.tip

plot(ftree)
ftree3 <- drop.tip(ftree, 4)
plot(ftree3)

plot(ftree)
ftree3 <- drop.tip(ftree, "A_boreas")
plot(ftree3)

plot(ftree)
ftree3 <- drop.tip(ftree, c(21,22))
plot(ftree3)

#extract clade
plot(ftree)
ftree3 <- extract.clade(ftree, 37)
plot(ftree3)

plot(ftree)
nodelabels()

#lining up data and a tree

#loading sample tree
ctree <- read.tree("canid_tree.txt")
str(ctree)
plot(ctree)

#tree is not fully resolved
ctree <- multi2di(ctree)
str(ctree)

#loading sample data
cdata <- read.csv("Canid_traits.csv")
str(cdata)

#using treedata to get body mass data lined up with tree

#extracting and cleaning body mass data mass
cmass <- cdata$AdultBodyMass_g

#assigning names
names(cmass) <- cdata$Binomial

#eliminating missing data
hist(cmass)
range(cmass)
cmass <- cmass[cmass > 0]

#log transforming (optional)
hist(cmass)
hist(log10(cmass))
cmass <- log10(cmass)

#using treedata
library(geiger)

CMout <- treedata(ctree, cmass)
str(CMout)
plot(CMout[[1]])
hist(CMout[[2]])

#final cleanup
CMtree <- CMout[[1]]
plot(CMtree)

#returns dataframe
CMass2 <- CMout[[2]]
head(CMass2)

#returns vector
CMass2 <- CMout[[2]][,1] 

head(CMass2)
str(CMtree)

#note that data is not in same order as tip labels
#to fix
CMass2 <- CMass2[CMtree$tip.label]

head(CMass2)
str(CMtree)

#testing for phylogenetic signal 
library(picante)

phylosignal(CMass2, CMtree)

######################################
#Tree plotting basics

#adding a scale bar

plot(ftree)
add.scale.bar()

plot(ftree)
axisPhylo()

#################################

#continuous character reconstruction

#creating a color pallete
mypallete <- c("red","orange","yellow","green","blue","purple")

#creating bins for species traits
divider <- (max(CMass2+0.0001)-min(CMass2))/6
index <- floor((CMass2-min(CMass2))/divider)

#ACE
CMtree$edge.length <- CMtree$edge.length+0.001
recon <- ace(CMass2, CMtree, CI = F)
recon.values <- recon$ace
recon.index <- floor((recon.values-min(CMass2))/divider)

#creating a plot
plot(CMtree, label.offset = 0.1, edge.width = 2)
nodelabels(pch = 21, cex = 2, bg = mypallete[recon.index+1])
tiplabels(pch = 21, cex = 2, bg = mypallete[index+1])

#add a legend
pts = c(0,1,2,3,4,5)
legend("topright", pch = 21, pt.cex = 2, c("0.00-0.25","0.25-0.50","0.50-0.75","0.75-1.00","1.00-1.26","1.26-1.51"), title = "Log 10 Mass", pt.bg = mypallete[pts+1])

##################################

#Discrete reconstruction

#getting a continuous character

AC <- cdata$ActivityCycle
names(AC) <- cdata$Binomial
AC <- AC[AC > 0]

#aligning tree and data
out2 <- treedata(ctree, AC)

CAtree <- out2[[1]]
CAtree <- multi2di(CAtree)
CAtree$edge.length <- CAtree$edge.length+0.001

CAC <- out2[[2]][,1]
CAC <- CAC[CAtree$tip.label]

#peforming reconstruction
anC <- ace(CAC, CAtree, type = "d")
states <- anC$lik.anc

#making the plot
plot(CAtree, edge.width = 2, label.offset = 1)
co <- c("white", "gray", "black")
tiplabels(pch = 22, bg = co[as.numeric(CAC)], cex = 2, adj = 1)
nodelabels(pie = states, piecol = c("white", "gray", "black"), cex = 0.5)
axisPhylo()

#adding a legend
pts = c(1,2,3)
legend("topright", pch = 22, pt.cex = 2, c("nocturnal","crepuscular","diurnal"), title = "Activity Cycle", pt.bg = co[pts])

