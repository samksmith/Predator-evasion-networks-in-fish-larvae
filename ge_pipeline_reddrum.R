# assembling data, running outlier detection, and fitting models

library(DESeq2)
library(arrayQualityMetrics)

#read in counts
counts = read.table("allcounts_redDrum.txt")
table(is.na(counts))
# how many genes do we have total?
nrow(counts)
ncol(counts)
names(counts)
# remove nonsamples
counts=counts[,-c(83:88)]
# remove data for samples that have no associated behavior data
notraits=c("X140rd","X142rd","X193rd","X194rd","X195rd","X196rd")
vector=c()
for(item in c(1:6)){
  vector=c(vector,(which(colnames(counts)==notraits[item])))
}
counts=counts[-vector]

# how many genes have mean count >=5?
# apply(table, 1=rows 2=columns, what apply to it)
means=apply(counts,1,mean)
#how many time above is satisfied
table(means>=5)
hist(log(means,10)) # histogram on log base 10

# only include genes with mean count >= 5
countData=counts[means>=5,]
nrow(countData) #17503 genes retained

# making a fake table of experimental conditions:
# label everything as "1"
ind=rep(1,ncol(countData))
a.con=data.frame(cbind(ind))
table(is.na(countData))
# because dataset is highly multifactorial and not fully crossed we can't use linear modeling
dds = DESeqDataSetFromMatrix(countData=countData, colData=a.con, design=~1)
# make big dataframe, getting normalized data
str(dds)
rl=rlog(dds) # normalizes and takes log
vsd=assay(rl)
save(dds,rl,countData,vsd,a.con,file="vsd_redDrum_final.RData")

#--------------------------------- WGCNA -------------------------------- #
# installing WGCNA:
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
#install.packages("flashClust",repos="http://cran.us.r-project.org")
#install.packages("WGCNA",dependencies=TRUE,repos="http://cran.us.r-project.org")
#,repos="http://cran.us.r-project.org"

library(WGCNA)
library(flashClust)
library(ape)
options(stringsAsFactors=FALSE)
#allowWGCNAThreads()

#--------------------------------
# run this if you are in /KOG or /GO or /heatmaps 
lnames=load("vsd_redDrum_final.RData")
lnames # "dds" "rl" "countData" vsd"  "a.con" # log-transformed variance-stabilized gene expression, and table or experimental conditions
datt=t(vsd)

# load behavioral datasets
traits1=read.csv("fuimanfish_metadata.csv")
traits2=read.csv("performance_reorg_042120.csv")
# combine them to contain just the data we want
traits=cbind.data.frame(traits1,traits2[,5:27])
row.names(traits)=traits$Matz_Number
# extract sample names out of the gene expression data in the same format
snames=gsub("rd|X","",colnames(vsd))
# use snames to address columns in traits
table(snames %in% traits$Matz_Number) # all true
table(traits$Matz_Number %in% traits$Matz_Number)
traits=traits[snames,]
traits$Matz_Number=snames
save(vsd,datt,traits,file="wgcnaData_fuimanfish.RData")

#-------------- PCA,heatmaps prior to covariate removal
load("wgcnaData_fuimanfish.RData")
library(pheatmap)
# similarity among samples
vsdn=vsd
# make heatmap
colnames(vsdn)=paste(traits[,1])
annotation_col = data.frame(
                    Group = traits$Replicate)
rownames(annotation_col) = colnames(vsdn)
ann_colors = list(
  Group = c("1" = "#FFB266", "2" = "#99CCFF"))
quartz()
pheatmap(cor(vsdn), annotation_col = annotation_col, 
         annotation_legend = FALSE, 
         annotation_colors = ann_colors,
         fontsize = 8)

# Principle coordinate analysis (pcoa)
library(vegan)
library(ape)

# import file that contains the relatedness group designations from ibs matrix pca (genetic groupings)
grp_id=read.csv("grp_geno.csv")
row.names(grp_id)=grp_id$X
table(grp_id$X %in% traits$Matz_Number)
# make a table of relevant covariates = genetic group, replicate group, and size
covariates=data.frame(cbind("group_id"=grp_id$grp_id,"rep"=traits$Replicate,"size"=traits$Size_mm))

adonis(datt ~ covariates[,1] + covariates[,3], data = covariates)

# pcoa based on manhattan distances
# divide by 1000 because manhattan dist are sum of all fold changes among all genes - large number
dds.pcoa=pcoa(dist(t(vsd),method="manhattan")/1000)
scores=dds.pcoa$vectors

# how many good PC's we have? Compared to random ("broken stick") model
quartz()
plot(dds.pcoa$values$Relative_eig)
points(dds.pcoa$values$Broken_stick,col="red",pch=3)

#fit covariates onto PCA
vectors=envfit(scores,covariates,choices=c(1,2)) # how covariates map onto 1st and 2nd PCo
vectors2=envfit(scores,covariates,choices=c(2,3)) # how covariates map onto 2nd and 3rd PCo
# plotting PCoA (change numbers to plot PCo1 and 2 or 2 and 3)
# find variance explained by each pco by looking at dds.pcoa object
quartz()
plot(scores[,2], scores[,3],
     col=ifelse(traits$Assay == "visual", "#FF66B3", ifelse(traits$Assay == "acoustic", "#66FFB2", NA)),
     xlab = "PCo 2 (8.6%)", ylab = "PCo 3 (3.9%)", pch=16) 
legend("bottomright", title = "Assay", legend = c("visual", "acoustic"), col = c('#FF66B3', '#66FFB2'), 
    pch = c(19,19), bty = "n", pt.cex = 0.7, cex = 1, text.col = "black", horiz = F, inset = c(0, 1))
plot(vectors2) # make sure put appropriate vector based on which PCo's displaying

# record the scores from PC1 and 2 into trait table
traits$gpc1=scores[,1]
traits$gpc2=scores[,2]
traits$grp_id=grp_id$grp_id

save(vsd,traits,covariates,file="vsd_traits_pcs_final.RData")

# --------------- removing effects of grp_id and size on data
library(limma)
load("vsd_traits_pcs_final.RData")
# remove grp_id and size as covariates
vsd2=removeBatchEffect(vsd,batch=traits$grp_id,covariates=traits$Size_mm)
datt=t(vsd2)
vsd=vsd2
save(vsd,datt,traits,file="wgcnaData_fuimanfish_allrem_final.RData")

#-------------- PCA after covariate removal -- check that removed effect of replicate, size, and grp_id
load("wgcnaData_fuimanfish_allrem_final.RData")
library(vegan)
library(ape)

# pcoa based on manhattan distances
dds.pcoa=capscale(dist(t(vsd),method="manhattan")~1)
quartz()
plot(dds.pcoa$CA$eig) # how many PCs are useful?

covariates=data.frame(cbind("grp_id"=traits$grp_id,"rep"=traits$Replicate,"size"=traits$Size_mm))
row.names(covariates)=traits$Matz_Number
vectors=envfit(dds.pcoa,covariates)
# if you print vectors you see that there is no relationship b/t vectors and pco values anymore
# Still see large vectors on plot bc they are automatically rescaled to be visible
quartz()
plot(dds.pcoa$CA$u[,1:2],
     col=ifelse(traits$Replicate == "1", "#FFB266", ifelse(traits$Replicate == "2", "#99CCFF", NA)),
     xlab = "PCo 1 (11.2%)", ylab = "PCo 2 (5.2%)", pch=16)
legend("bottomright", title = "Assay", legend = c("visual", "acoustic"), col = c('#FF66B3', '#66FFB2'), 
       pch = c(19,19), bty = "n", pt.cex = 0.7, cex = 1, text.col = "black", horiz = F, inset = c(0, 1))
# plotting PCoA
quartz()
ordihull(scores,conditions[,1],label=T,draw="polygon",col="grey90",cex=2)

# record scores from PC1 as trait variable
traits$gpc1b=dds.pcoa$CA$u[,1]

save(vsd,datt,traits,file="vsd_traits_pcs_allrem_final.RData")

# change names in traits table
pc1vis=traits$PC1Group
pc1vis[traits$Assay=="acoustic"]=NA
pc1ac=traits$PC1Group
pc1ac[traits$Assay=="visual"]=NA
pc2vis=traits$PC2Group
pc2vis[traits$Assay=="acoustic"]=NA
pc2ac=traits$PC2Group
pc2ac[traits$Assay=="visual"]=NA
pc3ac=traits$PC3Group

# recreate traits table with new names
traits1=data.frame(cbind(pc1vis,pc1ac,pc2vis,pc2ac,pc3ac,"gPC1b"=traits$gpc1b))
traits2=cbind.data.frame(traits1,traits[,10:30])
traits=traits2

save(vsd,datt,traits,file="wgcnaData_fuimanfish_allrem_final1.RData")

#################### WGCNA
# Try different betas ("soft threshold") - power factor for calling connections between genes
library(WGCNA)
load("wgcnaData_fuimanfish_allrem_final1.RData")
powers = c(seq(from = 1, to=10, by=0.5))
# Call the network topology analysis function
sft = pickSoftThreshold(datt, powerVector = powers, verbose = 5,networkType="signed")

# Plot the results:
quartz()
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#dev.off()

#####################
# making modules
library(flashClust)
library(ape)

s.th=2.5 # re-specify according to previous section; chose value that corresponds with R^2 cut-off
# calculating correlation or distance network adjacency from given expression data
adjacency = adjacency(datt, power = s.th,type="signed");
# calculates topological overlap matrix and corresponding dissimilarity from given adjacency matrix
TOM = TOMsimilarity(adjacency,TOMType="signed");
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30; 
# adaptive branch pruning of hierarchical clustering dendrograms
dynamicMods = cutreeDynamic(dendro = geneTree,distM = dissTOM,
                            deepSplit = 2,pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
# convert vector or array of numerical labels into corresponding vector or array of colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Calculate eigengenes (1st principal component) of modules in dataset
MEList = moduleEigengenes(datt, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# hierarchical clustering of dissimilarity
METree = flashClust(as.dist(MEDiss), method = "average");

save(dynamicMods,dynamicColors,MEs,METree,geneTree,file="1stPassModules_final.RData")
######################### merging modules
load("wgcnaData_fuimanfish_allrem_final1.RData")
mm=load("1stPassModules_final.RData")
mm
library(WGCNA)

quartz()
MEDissThres = 0 # in the first pass, set this to 0 - no merging 
# (we want to see the module-traits heatmap first, then decide which modules are telling us the same story and better be merged)
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")  # on 2nd pass: does this cut height meet your merging goals? If not, reset MEDissThres and replot

# Call an automatic merging function
merge = mergeCloseModules(datt, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# plotting the dendrogram
quartz()
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,
                    addGuide = FALSE, guideHang = 0.05,lwd=0.3)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Calculate dissimilarity of module eigengenes
quartz()
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# how many genes in each module?
table(moduleColors)
write.csv(moduleColors, file = "modules_varincl_042120.csv" )
save(MEs, geneTree, moduleLabels, moduleColors, file = "networkdata_signed_final.RData")

################### plotting correlations with traits
load(file = "networkdata_signed_final.RData")
load("wgcnaData_fuimanfish_allrem_final1.RData")

# Define numbers of genes and samples
nGenes = ncol(datt);
nSamples = nrow(datt);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datt)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);

moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# module-trait correlations
quartz()
# put in symbols for significance on heatmap
textMatrix = signif(moduleTraitPvalue, 1)
textMatPsym = matrix(,nrow = 30, ncol = 27)
for (row in 1:nrow(textMatrix)){
  for (col in 1:ncol(textMatrix)){
    if(textMatrix[row,col] <= 0.01){
      textMatPsym[row,col] = "***"
    }
    else if(textMatrix[row,col] <= 0.05){
      textMatPsym[row,col] = "**"
    }
    else if(textMatrix[row,col] <= 0.1){
      textMatPsym[row,col] = "*"
    }
    else {
      textMatPsym[row,col] = ""
    }
  }
}

dim(textMatPsym) = dim(moduleTraitCor)
par(mar = c(6, 7.5, 2, 2));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatPsym,
               setStdMargins = FALSE,
               cex.text = 0.8,
               cex.lab = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

print(data.frame(table(moduleColors))) # gives numbers of genes in each module

save(MEs, geneTree, moduleLabels, moduleColors,traits,file="networkdata_signed_fuimanfish_1_final.RData")
# if it was first pass with no module merging, this is where you examine your heatmap
# and dendrogram of module eigengenes to see where you would like to set cut height (MEDissThres parameter) in the previous section
# to merge modules that are telling the same story for your trait data

# good way to do it is to find a group of similar modules in the heat map and then see 
# at which tree height they connect in the dendrogram.

############# scatterplots of gene significance (correlation-based) vs kME
## these plots tell you whether genes that are more central (well-connected) in the network are
# also more highly correlated with a particular trait
library(WGCNA)
load("wgcnaData_fuimanfish_allrem_final1.RData");
load(file ="networkdata_signed_fuimanfish_1_final.RData")

table(moduleColors)
whichTrait="gPC1b"

nGenes = ncol(datt);
nSamples = nrow(datt);
selTrait = as.data.frame(traits[,whichTrait]);
names(selTrait) = whichTrait
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(signedKME(datt, MEs));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(selTrait), sep="");
names(GSPvalue) = paste("p.GS.", names(selTrait), sep="");
quartz()
par(mfrow=c(3,3))
counter=0
for(module in modNames[1:length(modNames)]){
  counter=counter+1
  if (counter>9) {
    quartz()
    par(mfrow=c(3,3))
    counter=1
  }
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste(module,"module membership"),
                     ylab = paste("GS for", whichTrait),
                     col = "grey50",mgp=c(2.3,1,0))
}

################ eigengene-heatmap plot (sanity check - is the whole module driven by just one crazy sample?)
# note: this part does not make much sense for unsigned modules

which.module="blue" 
datME=MEs
datExpr=datt
quartz()
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="sample")

length(datExpr[1,moduleColors==which.module ]) # number of genes in chosen module

################# saving selected modules for GO and KOG analysis (two-parts: Fisher test, MWU test within-module)

library(WGCNA)
load("wgcnaData_fuimanfish_allrem_final1.RData");
load(file = "networkdata_signed_final.RData")

# calculating module memberships for all genes for all modules
allkME =as.data.frame(signedKME(datt, MEs)) 
names(allkME)=gsub("kME","",names(allkME))

whichModule="blue"
table(moduleColors==whichModule) # how many genes are in it? need at least 100 for GOMWU
vsd=t(datt)
# Saving data for Fisher-MWU combo test (GO_MWU)
inModuleBinary=as.numeric(moduleColors==whichModule)
combo=data.frame("gene"=row.names(vsd),"Fish_kME"=allkME[,whichModule]*inModuleBinary)
write.csv(combo,file=paste(whichModule,"_varincl_042120.csv",sep=""),row.names=F,quote=F)

################ plotting heatmap for named top-kME genes
load("wgcnaData_fuimanfish_allrem_final.RData");
load(file = "networkdata_signed_final.RData")
allkME =as.data.frame(signedKME(datt, MEs))
gg=read.table("redDrum_iso2gene.tab",sep="\t")
library(pheatmap)
vsd=t(datt)
whichModule="turquoise"
top=30 # number of named top-kME genes to plot

datME=MEs
datExpr=datt
modcol=paste("kME",whichModule,sep="")
sorted=vsd[order(allkME[,modcol],decreasing=T),]
head(sorted)
# selection top N names genes, attaching gene names
gnames=c();counts=0;hubs=c()
for(i in 1:length(sorted[,1])) {
  if (row.names(sorted)[i] %in% gg$V1) { 
    counts=counts+1
    gn=gg[gg$V1==row.names(sorted)[i],2]
    gn=paste(gn,row.names(sorted)[i],sep=".")
    if (gn %in% gnames) {
      gn=paste(gn,counts,sep=".")
    }
    gnames=append(gnames,gn) 
    hubs=data.frame(rbind(hubs,sorted[i,]))
    if (counts==top) {break}
  }
} 
row.names(hubs)=gnames

contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
contrasting2 = colorRampPalette(rev(c("chocolate1","chocolate1","#FEE090","grey10", "cyan3","cyan","cyan")))(100)
contrasting3 = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan","cyan")))(100)
quartz()
pheatmap(hubs,scale="row",col=contrasting2,border_color=NA,treeheight_col=0,cex=0.7,cluster_rows=F)

# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriches with either up- or down-regulated genes. The advantage - no need to impose arbitrary significance cutoff.
# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 
# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.
# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.
# The fraction next to GO category name indicates the fracton of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.
# Stretch the plot manually to match tree to text
# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu
################################################################
# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.
setwd("")
input="blue_varincl_042120.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="redDrum_iso2go_emapperNR.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF - molecular function, or BP - biological process, or CC - cellular compartment
source("gomwu.functions.R")


# Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.

# Plotting results
quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=0.001,  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=1,
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=0.75,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=2.0, # height of the hierarchical clustering tree
                  # colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results
# -------------------------- ANALYSIS WITH REMOVAL OF LATENT VARIABLE -------------------
setwd("~/Desktop/Reddrum_final")
library(WGCNA)
library(flashClust)
library(ape)
options(stringsAsFactors=FALSE)

# --------------- removing effects of gpc1b from data (replicate/grp_id and size already removed)
library(limma)
load("vsd_traits_pcs_allrem_final.RData")
# identify effect size of latent variable
dds.pcoa=capscale(dist(t(vsd),method="manhattan")~1)
vectors = envfit(dds.pcoa,traits$gpc1b)
vectors #identify affect size
vsd2=removeBatchEffect(vsd,covariates=traits$gpc1b)
datt=t(vsd2)
vsd=vsd2
save(vsd,datt,traits,file="wgcnaData_fuimanfish_fin.RData")

#-------------- PCA after covariate removal -- check that removed effect of replicate, size, grp_id, and gpc1
load("wgcnaData_fuimanfish_fin.RData")

# Principle coordiante analysis
library(vegan)
library(ape)

# pcoa based on manhattan distances
dds.pcoa=pcoa(dist(t(vsd),method="manhattan")/1000)
scores=dds.pcoa$vectors

# how many good PC's we have? Compared to random ("broken stick") model
quartz()
plot(dds.pcoa$values$Relative_eig) # find where eigengene info 
points(dds.pcoa$values$Broken_stick,col="red",pch=3)

covariates=data.frame(cbind("grp_id"=traits$grp_id,"rep"=traits$Replicate,"size"=traits$Size_mm,"gpc1"=traits$gpc1))
row.names(covariates)=traits$Matz_Number
vectors=envfit(scores,covariates,choices=c(1,2))
# plotting PCoA
quartz()
plot(scores[,1], scores[,2],
     col=ifelse(traits$Assay == "visual", "#FF66B3", ifelse(traits$Assay == "acoustic", "#66FFB2", NA)),
     xlab = "PCo 1 (5.8%)", ylab = "PCo 2 (4.7%)", pch=16)
legend("bottomright", title = "Assay", legend = c("visual", "acoustic"), col = c('#FF66B3', '#66FFB2'), 
       pch = c(19,19), bty = "n", pt.cex = 0.7, cex = 1, text.col = "black", horiz = F, inset = c(0, 1))

save(vsd,datt,traits,file="vsd_traits_pcs_final.RData")

# change names for traits
pc1vis=traits$PC1Group
pc1vis[traits$Assay=="acoustic"]=NA
pc1ac=traits$PC1Group
pc1ac[traits$Assay=="visual"]=NA
pc2vis=traits$PC2Group
pc2vis[traits$Assay=="acoustic"]=NA
pc2ac=traits$PC2Group
pc2ac[traits$Assay=="visual"]=NA
pc3ac=traits$PC3Group

# recreate traits table with new names
traits1=data.frame(cbind(pc1vis,pc1ac,pc2vis,pc2ac,pc3ac))
traits2=cbind.data.frame(traits1,traits[,10:30])
traits=traits2

save(vsd,datt,traits,file="wgcnaData_fuimanfish_final.RData")

#################### WGCNA
# Try different betas ("soft threshold") - power factor for calling connections between genes
library(WGCNA)
load("wgcnaData_fuimanfish_final.RData")
powers = c(seq(from = 1, to=10, by=0.5))
# Call the network topology analysis function
sft = pickSoftThreshold(datt, powerVector = powers, verbose = 5,networkType="signed")

quartz()
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#dev.off()

##################### making modules

s.th=5 # re-specify according to previous section
adjacency = adjacency(datt, power = s.th,type="signed");
TOM = TOMsimilarity(adjacency,TOMType="signed");
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30; 
dynamicMods = cutreeDynamic(dendro = geneTree,distM = dissTOM,
                            deepSplit = 2,pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Calculate eigengenes
MEList = moduleEigengenes(datt, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
METree = flashClust(as.dist(MEDiss), method = "average");

save(dynamicMods,dynamicColors,MEs,METree,geneTree,file="1stPassModules_1_final.RData")
######################### merging modules

load("wgcnaData_fuimanfish_final.RData")
mm=load("1stPassModules_1_final.RData")
mm
library(WGCNA)

quartz()
MEDissThres = 0.4 # in the first pass, set this to 0 - no merging (we want to see the module-traits heatmap first, then decide which modules are telling us the same story and better be merged)
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")  # on 2nd pass: does this cut height meet your merging goals? If not, reset MEDissThres and replot

# Call an automatic merging function
merge = mergeCloseModules(datt, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# plotting the dendrogram
quartz()
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,
                    addGuide = FALSE, guideHang = 0.05,lwd=0.3)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Calculate dissimilarity of module eigengenes
quartz()
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# how many genes in each module?
modcol = table(moduleColors)
write.csv(modcol, "allrem_module_list.csv")
save(MEs, geneTree, moduleLabels, moduleColors, file = "networkdata_signed_1_final.RData")

################### plotting correlations with traits:
load(file = "networkdata_signed_1_final.RData")
load("wgcnaData_fuimanfish_final.RData")

# Define numbers of genes and samples
nGenes = ncol(datt);
nSamples = nrow(datt);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datt)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);

moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
modtraitpval = signif(moduleTraitPvalue, 1)

# module-trait correlations
quartz()
textMatrix = signif(moduleTraitPvalue, 1)
# make a matrix that turns p-vals into symbols
textMatPsym = matrix(,nrow = 25, ncol = 26) # you need to know how big your matrix is
for (row in 1:nrow(textMatrix)){
  for (col in 1:ncol(textMatrix)){
    if(textMatrix[row,col] < 0.01){
      textMatPsym[row,col] = "***"
    }
    else if(textMatrix[row,col] < 0.05){
      textMatPsym[row,col] = "**"
    }
    else if(textMatrix[row,col] < 0.1){
      textMatPsym[row,col] = "*"
    }
    else {
      textMatPsym[row,col] = ""
    }
  }
}
dim(textMatPsym) = dim(moduleTraitCor)
par(mar = c(6, 7.5, 2, 2));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatPsym,
               setStdMargins = FALSE,
               cex.text = 0.8,
               cex.lab = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

print(data.frame(table(moduleColors))) # gives numbers of genes in each module

save(MEs, geneTree, moduleLabels, moduleColors,traits,file="networkdata_signed_fuimanfish_2_final.RData")
# if it was first pass with no module merging, this is where you examine your heatmap
# and dendrogram of module eigengenes to see where you would like to see 
# where you woudl like to set cut height (MEDissThres parameter) in the previous section
# to merge modules that are talling the same story for your trait data

# good way to do it is to find a group of similar modules in the heat map and then see 
# at which tree height they connect in the dendrogram.

############# scatterplots of gene significance (correlation-based) vs kME
library(WGCNA)
load("wgcnaData_fuimanfish_final.RData");
load(file ="networkdata_signed_fuimanfish_2_final.RData")

table(moduleColors)
whichTrait="RMeanMobility"

nGenes = ncol(datt);
nSamples = nrow(datt);
selTrait = as.data.frame(traits[,whichTrait]);
names(selTrait) = whichTrait
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(signedKME(datt, MEs));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(selTrait), sep="");
names(GSPvalue) = paste("p.GS.", names(selTrait), sep="");
quartz()
par(mfrow=c(3,3))
counter=0
for(module in modNames[1:length(modNames)]){
  counter=counter+1
  if (counter>9) {
    quartz()
    par(mfrow=c(3,3))
    counter=1
  }
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste(module,"module membership"),
                     ylab = paste("GS for", whichTrait),
                     col = "grey50",mgp=c(2.3,1,0))
}

################
# eigengene-heatmap plot (sanity check - is the whole module driven by just one crazy sample?)
# note: this part does not make much sense for unsigned modules

which.module="plum1" 
datME=MEs
datExpr=datt
quartz()
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="sample")

length(datExpr[1,moduleColors==which.module ]) # number of genes in chosen module

################# saving selected modules for GO and KOG analysis (two-parts: Fisher test, MWU test within-module)

library(WGCNA)
load("wgcnaData_fuimanfish_final.RData");
load(file = "networkdata_signed_fuimanfish_2_final.RData")

# calculating module memberships for all genes for all modules
allkME =as.data.frame(signedKME(datt, MEs)) 
names(allkME)=gsub("kME","",names(allkME))

whichModule="saddlebrown"
table(moduleColors==whichModule) # how many genes are in it? You want at least 100 for GOMWU
vsd=t(datt)
# Saving data for Fisher-MWU combo test (GO_MWU)
inModuleBinary=as.numeric(moduleColors==whichModule)
combo=data.frame("gene"=row.names(vsd),"Fish_kME"=allkME[,whichModule]*inModuleBinary)
write.csv(combo,file=paste(whichModule,"_allrem.csv",sep=""),row.names=F,quote=F)

################ plotting heatmap for named top-kME genes

load("wgcnaData_fuimanfish_final.RData");
load(file = "networkdata_signed_fuimanfish_2_final.RData")
allkME =as.data.frame(signedKME(datt, MEs))
gg=read.table("redDrum_iso2gene_emapperNR.tab",sep="\t",quote="",fill=FALSE)
library(pheatmap)
vsd=t(datt)
whichModule="darkmagenta"
top=30 # number of named top-kME genes to plot

datME=MEs
datExpr=datt
modcol=paste("kME",whichModule,sep="")
sorted=vsd[order(allkME[,whichModule],decreasing=T),]
head(sorted)
# selection top N names genes, attaching gene names
gnames=c();counts=0;hubs=c()
for(i in 1:length(sorted[,1])) {
  if (row.names(sorted)[i] %in% gg$V1) { 
    counts=counts+1
    gn=gg[gg$V1==row.names(sorted)[i],2]
    gn=paste(gn,row.names(sorted)[i],sep=".")
    if (gn %in% gnames) {
      gn=paste(gn,counts,sep=".")
    }
    gnames=append(gnames,gn) 
    hubs=data.frame(rbind(hubs,sorted[i,]))
    if (counts==top) {break}
  }
} 
row.names(hubs)=gnames

contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
contrasting2 = colorRampPalette(rev(c("chocolate1","chocolate1","#FEE090","grey10", "cyan3","cyan","cyan")))(100)
contrasting3 = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan","cyan")))(100)
quartz()
pheatmap(hubs,scale="row",show_rownames=FALSE,col=contrasting2,border_color=NA,treeheight_col=0,cex=0.7,cluster_rows=F)

# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriches with either up- or down-regulated genes. The advantage - no need to impose arbitrary significance cutoff.
# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 
# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.
# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.
# The fraction next to GO category name indicates the fracton of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.
# Stretch the plot manually to match tree to text
# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu
################################################################
# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.
setwd("~/Desktop/RedDrum_final/GO_MWU_050620")
input="darkmagenta_allrem.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="redDrum_iso2go_emapperNR.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="CC" # either MF - molecular function, or BP - biological process, or CC - cellular compartment
source("gomwu.functions.R")


# Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.

# Plotting results
quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=0.001,  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=1,
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=0.75,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=2.0, # height of the hierarchical clustering tree
                  # colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results
