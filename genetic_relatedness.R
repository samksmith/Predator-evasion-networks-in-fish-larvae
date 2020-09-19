setwd("~/Desktop/Reddrum_final")
bams=read.table("bams_onlytraits")[,1] # list of bam files for samples that had behavioral data
#--------------------
# covariance / PCA 
library(vegan)
library(stringr)
library(ggplot2)
co = as.matrix(read.table("mito_traits.covMat")) # covariance based on single-read resampling for mitochondrial genes
co.geno = as.matrix(read.table("myresult.covMat")) # covariance based on single-read resampling for all genes
ids = read.table("bams")[,1] 
nam = gsub(".fq.trim.bam","",ids)
names=gsub(".fq.trim.bam","",bams)
dimnames(co)=list(names,names)
dimnames(co.geno)=list(nam,nam)
# remove samples that don't have traits associated with them from co.geno (covariance matrix w all genes included)
notraits=c("140","142","193","194","195","196")
vector=c()
for(item in c(1:6)){
  vector=c(vector,(which(colnames(co.geno)==notraits[item])))
}
co.geno=co.geno[-vector,-vector]

traits=read.csv("fuimanfish_metadata.csv")
row.names(traits)=traits$Matz_Number

pp0=capscale(as.dist(1-cov2cor(co))~1) # PCoA maternal genes
pp0a=capscale(as.dist(1-cov2cor(co.geno))~1) # PCoA all genes

# eigenvectors: how many are interesting?
plot(pp0a$CA$eig) 

axes2plot=c(1,2)   # plotting PC 1 and 2
quartz()
cc=pp0 # using unconstrained analysis data (mitochondrial genes)
cc1=pp0a # all genes
# cluster into four groups
k=kmeans(cc1$CA$u[,1:3],4)
plot(cc,choices=axes2plot, type = 'n', xlab = "PCo 1 (57%)", ylab = "PCo 2 (0.87%)") # choices - axes to display
abline(h=0, col="white", lty = 1, lwd = 5)
abline(v=0, col="white", lty = 1, lwd = 5)
box()
points(cc,choices=axes2plot,pch=19,
       col=ifelse(traits$Replicate == 1, "#FFB266", ifelse(traits$Replicate == 2, "#99CCFF", NA)))
legend("right", legend = c("Group 1", "Group 2"), col = c('#FFB266', '#99CCFF'), 
  pch = c(19,19), bty = "n", pt.cex = 0.7, cex = 1, text.col = "black", horiz = F, inset = c(0, 1))

# plot ellipses to show clustering
ordiellipse(cc1,choices= axes2plot,groups=k$cluster,draw="polygon",label=T)
n2identify=4
identify(cc1$CA$u[,axes2plot],labels=colnames(co.geno),n=n2identify,cex=0.7)

# record which group individuals in
grp_geno = data.frame(k$cluster)
colnames(grp_geno) = "grp_id"
write.csv(grp_geno, file = "grp_geno.csv")
cov_geno_scores=cc1$CA$u[,1:2]
write.csv(cov_geno_scores, file = "cov_geno_scores.csv")

#-------------
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)

ma = as.matrix(read.table("mito_traits.ibsMat"))
ma.geno = as.matrix(read.table("myresult.ibsMat"))
nam = gsub(".fq.trim.bam","",ids)
names=gsub(".fq.trim.bam","",bams)
dimnames(ma)=list(names,names)
dimnames(ma.geno)=list(nam,nam)
notraits=c("140","142","193","194","195","196")
vec=c()
for(item in c(1:6)){
  vec=c(vec,(which(colnames(ma.geno)==notraits[item])))
}
ma.geno = ma.geno[-vec,-vec]
hc=hclust(as.dist(ma),"ave")
quartz()
plot(hc,cex=0.5)
hc.geno = hclust(as.dist(ma.geno), "ave")
quartz()
plot(hc.geno,cex=0.5)


pp1=capscale(ma~1)

# eigenvectors
quartz()
plot(pp1$CA$eig) 

axes2plot=c(1,2)  
quartz()
library(adegenet) # for transp()
cmd=pp1
k3=kmeans(cmd$CA$u[,1:3],4)
plot(cmd,choices=axes2plot,display="sites",type="n", xlab = "PCoA 1 (44%)", ylab = "PCoA 2 (8.4%)") # choices - axes to display
abline(h=0, col="white", lty = 1, lwd = 5)
abline(v=0, col="white", lty = 1, lwd = 5)
box()
points(cmd,choices=axes2plot,pch=19,
       col=ifelse(traits$Replicate == 1, "#FFB266", ifelse(traits$Replicate == 2, "#99CCFF", NA)))
legend("right", legend = c("Group 1", "Group 2"), col = c('#FFB266', '#99CCFF'), 
       pch = c(19,19), bty = "n", pt.cex = 0.7, cex = 1, text.col = "black", horiz = F, inset = c(0, 1))
ordiellipse(cmd,choices= axes2plot,groups=k3$cluster,draw="polygon",label=T)
cl=kmeans(ma,4)