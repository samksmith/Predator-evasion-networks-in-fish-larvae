options(stringsAsFactors = FALSE)
library(stringr)
library(dplyr)
##################### Making top50 kME gene module files ################
# load all csv files of gene module membership for modules enriched for GO terms
# an example below for the saddlebrown module
saddlebrown = read.csv("saddlebrown.csv")
# keep only genes that belong to the module
g = c()
kme = c()
for(gene in 1:nrow(saddlebrown)){
  if(saddlebrown[gene,2] > 0){
    g = append(g,saddlebrown[gene,1])
    kme = append(kme,saddlebrown[gene,2])
  }
}
saddlebrown = as.data.frame(cbind(g,kme,m = rep("saddlebrown",length(g))))
# sort in order of highest kME to lowest kME
saddlebrown=saddlebrown[order(-kme),]
# keep only the top 50% in kME value
t50saddlebrown = saddlebrown[-c((0.5*nrow(saddlebrown))+1:nrow(saddlebrown)),]
write.csv(t50saddlebrown,"t50saddlebrown.csv",row.names=FALSE)

####################### Gene Annotation Generation #############################
options(stringsAsFactors = FALSE)
# load necessary data
library(stringr)
library(dplyr)
load("wgcnaData_fuimanfish_final.RData");
load(file = "networkdata_signed_fuimanfish_2_final.RData")

# load top 50 kME gene files for each module
t50green = read.csv("t50green.csv")
t50blue = read.csv("t50blue.csv")
t50black = read.csv("t50black.csv")
t50brown = read.csv("t50brown.csv")
t50darkgreen = read.csv("t50darkgreen.csv")
t50skyblue = read.csv("t50skyblue.csv")
t50red = read.csv("t50red.csv")
t50greenyellow = read.csv("t50greenyellow.csv")
t50darkred = read.csv("t50darkred.csv")
t50darkorange = read.csv("t50darkorange.csv")
t50darkmagenta = read.csv("t50darkmagenta.csv")
t50lightcyan1 = read.csv("t50lightcyan1.csv")
t50darkgrey = read.csv("t50darkgrey.csv")
t50lightsteelblue1 = read.csv("t50lightsteelblue1.csv")
t50cyan = read.csv("t50cyan.csv")
t50orangered4 = read.csv("t50orangered4.csv")
t50floralwhite = read.csv("t50floralwhite.csv")
t50ivory = read.csv("t50ivory.csv")
t50plum1 = read.csv("t50plum1.csv")
t50midnightblue = read.csv("t50midnightblue.csv")
t50brown4 = read.csv("t50brown4.csv")
t50saddlebrown = read.csv("t50saddlebrown.csv")

# import table that associates module with significantly associated traits
modTraitCor = read.csv("module-trait_NoPCs.csv")
# import table that has human readable names for each trait variable
trait2name = read.csv("trait_to_name.csv")
modTraitCor$mod = str_replace(modTraitCor$mod,"ME","") # Remove "ME" from module name
# add human readable names to module trait association table
modTrName = merge(modTraitCor,trait2name,by = "tr", all.x = TRUE, all.y = FALSE, sort = FALSE)
modTrName$cor = round(as.numeric(modTrName$cor),digits=2)
# make dataframe with trait association description for each module
# vector for module that is upreg in association w a particular trait
upmodfortrait = c()
# vector for  trait that is upregulated
upreg = c()
# vector for module that is downreg in association w particular trait
downmodfortrait = c()
# vector for trait that is downregulated
downreg = c()
for(row in 1:nrow(modTrName)){
  # if shows positive correlation b/t module and trait, record as upregulated
  if(modTrName[row,3] > 0){
    upmodfortrait = append(upmodfortrait,modTrName[row,2])
    upreg = append(upreg,paste("(R2=+",modTrName[row,3],") ",modTrName[row,5],sep=""))
  }
  # otherwise, shows a negative correlation, record as downreg
  else{
    downmodfortrait = append(downmodfortrait,modTrName[row,2])
    downreg = append(downreg,paste("(R2=",modTrName[row,3],") ",modTrName[row,5],sep=""))
  }
}
# create dataframe with modules and associated upregulated traits
upregulatedtraits = as.data.frame(cbind(upmodfortrait,upreg))
# create a dataframe with modules and associated downregulated traits
downregulatedtraits = as.data.frame(cbind(downmodfortrait,downreg))

# combine all module kme files into a single file
g_mod_kme = as.data.frame(rbind(t50black,t50blue,t50brown,t50darkgreen,t50darkorange,t50darkred,t50green,t50greenyellow,
                t50red,t50skyblue,t50darkmagenta,t50lightcyan1,t50darkgrey,
                t50lightsteelblue1,t50ivory,t50cyan,t50orangered4,t50floralwhite,
                t50plum1,t50midnightblue,t50brown4,t50saddlebrown))
# round kME values to 2 digits
g_mod_kme$kme = round(as.numeric(g_mod_kme$kme),digits=2)
# create dataframe with annotations 
annots = c()
# This for loop creates a vector containing all annotations
for(v in 1:nrow(g_mod_kme)){
  # find traits upregulated in assoc. with given module
  upts = filter(upregulatedtraits,upmodfortrait == g_mod_kme[v,3])
  # combine all traits (upreg) into given string
  upts = paste(upts[,2],sep="",collapse="; ")
  # find traits downregulated in assoc. with given module
  dnts = filter(downregulatedtraits,downmodfortrait == g_mod_kme[v,3])
  # combine all traits (downreg) into given string
  dnts = paste(dnts[,2],sep="",collapse="; ")
  # if there are upreg and downreg traits significantly associated write them
  if(upts != "" & dnts != ""){
    annots = append(annots,paste(g_mod_kme[v,1],"|",g_mod_kme[v,3]," kME=",g_mod_kme[v,2],"|",
                                 "Associations:"
                                 ,upts,"; ",dnts,".",sep=""))
  }
  # if there are only upregulated traits significantly associated with module
  else if(upts != "" & dnts == ""){
    annots = append(annots,paste(g_mod_kme[v,1],"|",g_mod_kme[v,3]," kME=",g_mod_kme[v,2],"|",
                                 "Associations:",upts,".",sep=""))
  }
  # if there are only downregulated traits significantly associated with module
  else{
    annots = append(annots,paste(g_mod_kme[v,1],"|",g_mod_kme[v,3]," kME=",g_mod_kme[v,2],"|",
                                 "Associations:",dnts,".",sep=""))
  }
}
# make the dataframe with the sequence ID and gene annotation
annotations = as.data.frame(cbind(g_mod_kme[,1],annots))

# write dataframe to csv
write.csv(annotations,"gene_annotations.csv",row.names = FALSE)
save(annotations,g_mod_kme,modTrName,file="gene_annotations.RData")
########################### Testing gene annotation scheme against EggNOG annotations ##############################
# load gene_annotation tables generated above
load("gene_annotations.RData")
# gene annotations from somewhere....figure out
gg=read.table("redDrum_iso2gene_emapperNR.tab",sep="\t",quote="",fill=FALSE)

# only keep genes from gg table that we annotated above (top 50% kME in modules w GO terms)
table(gg[,1] %in% g_mod_kme[,1]) # FALSE 8694; TRUE 4602

g_annots = as.data.frame(NULL)
for(row in 1:nrow(gg)){
  if(gg$V1[row] %in% g_mod_kme[,1] == TRUE){
    g_annots = rbind(g_annots,gg[row,])
  }
}
# merge the annotations into a single dataframe
annotation_test = merge(annotations,g_annots, by.x = "V1", by.y = "V1",
                        all.x = TRUE, all.y = TRUE, sort = TRUE)
colnames(annotation_test) = c("sequence","NI_annotation","EggNOG_annotation")
write.csv(annotation_test,"gene_annotation_test.csv",row.names=FALSE)
save(annotations,annotation_test,g_annots,g_mod_kme,modTrName,file="gene_annotations_test.RData")
