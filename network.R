# data preparation
counts <- read.table("counts.txt", sep="\t", header=T, row.names = 1)[1:31]
traits <- rep(1,31)
traits[c(6,18,19,21,25)] <- 0 #Modern corn do not branch
subsp <- c(1,1,1,0,1,NA,1,1,1,1,1,1,1,1,1,1,1,NA,NA,0,NA,0,0,1,NA,1,1,0,1,1,1) # "1" for parviglumis, "0" for mexicana, "NA" for modern maize
require(edgeR)
y <- DGEList(counts = counts)
y <- calcNormFactors(y)
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
#Filter out low read genes
filter_low_counts <- function (x, cutoff = 1) {
  sum(x < cutoff)
}
counts <- data.frame(cpm(y)[-which(apply(cpm(y), MARGIN = 1, FUN = filter_low_counts) > dim(counts)[2]*.9) ,])
counts <- counts[-c(7,8)]
target <- traits[-c(7,8)]

## ---Network analysis---
library(WGCNA)
library(flashClust)
datExpr <- data.frame(t(log(counts+1, base=2)))
# powers <- c(1:30)
# choose power based on SFT criterion
# disableWGCNAThreads() #To get rid of the warning in the notes of the code next line
# sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed") #If not run the code above, Warning message: executing %dopar% sequentially: no parallel backend registered. I think it just means it is running sequentially if no parallel backend has been registered. It will issue this warning only once. 
A <- adjacency(datExpr, power=15, type="signed", distFnc="bicor")
dissTOM <- TOMdist(A,TOMType="signed")
geneTree <- flashClust(as.dist(dissTOM), method="average")
# plotDendroAndColors(geneTree,colors="red",dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
TOM <- TOMsimilarityFromExpr(datExpr, power=15, TOMType = "signed") #Time consuming
probes <- names(datExpr)

# read in candidates
lasso <- read.table("./cutoff-1-cpm/reduced-dimension/lasso.txt")
a <- lasso.genes <- as.character(read.table("./cutoff-1-cpm/reduced-dimension/lasso.txt")[,1])
elnet <- read.table("./cutoff-1-cpm/original-dimension/elnet.txt")
b <- elnet.genes <- as.character(read.table("./cutoff-1-cpm/original-dimension/elnet.txt")[,1])
rf <- read.table("./cutoff-1-cpm/original-dimension/rfimp.txt")
c <- rf.genes <- as.character(rownames(read.table("./cutoff-1-cpm/original-dimension/rfimp.txt")))
xgb <- read.table("./cutoff-1-cpm/reduced-dimension/xgboost.txt")
d <- xgb.genes <- as.character(read.table("./cutoff-1-cpm/reduced-dimension/xgboost.txt")[,1])
nb <- read.table("./cutoff-1-cpm/original-dimension/nbcimp.txt")
e <- nb.genes <- as.character(read.table("./cutoff-1-cpm/original-dimension/nbcimp.txt")[,1])

require(VennDiagram)
require(venn)
pdf("./venn.pdf",width = 7,height= 7)
draw.quintuple.venn(124 ,41, 2368, 118, 913  ,  2 ,  29 , 118 ,  26  , 32  ,  2 ,  40 ,  29  ,495,   26,
                    1  ,  2 ,   2 ,  29 ,  17,   26  ,  1  , 32    ,2,   17 ,   1 ,   1,    2,   17 ,   1,  1,
                    category = c("LASSO", "Elastic Net", "Random Forest", "Boosted Tree", "Naive Bayes"),
                    cat.dist = 0.09)
dev.off()

# # to test if the non-overlapping features are highly correlated
# moduleLabelsManual <- cutreeDynamic(dendro=geneTree, distM = dissTOM, method = "hybrid",
#                                      deepSplit = 2, pamRespectsDendro = F, minClusterSize = 30)
# moduleColorsManual <- labels2colors(moduleLabelsManual)
# lasso.module <- moduleColorsManual[which(probes %in% lasso.genes)]
# elnet.module <- moduleColorsManual[which(probes %in% elnet.genes)]
# rf.module <- moduleColorsManual[which(probes %in% rf.genes)]
# xgb.module <- moduleColorsManual[which(probes %in% xgb.genes)]
# nb.module <- moduleColorsManual[which(probes %in% nb.genes)]
# pdf("./venn-modules.pdf")
# draw.quintuple.venn(124 ,41, 2368, 118, 913  ,  2 ,  29 , 118 ,  26  , 32  ,  2 ,  40 ,  29  ,495,   26,
#                     1  ,  2 ,   2 ,  29 ,  17,   26  ,  1  , 32    ,2,   17 ,   1 ,   1,    2,   17 ,   1,  1,
#                     category = c("LASSO", "Elastic Net", "Random Forest", "Boosted Tree", "Naive Bayes Classifier"))
# dev.off()

x <- c()
x[1] <- length(a)
x[2] <- length(b)
x[3] <- length(c)
x[4] <- length(d)
x[5] <- length(e)
x[6] <- length(intersect(a,b))
x[7] <- length(intersect(a,c))
x[8] <- length(intersect(a,d))
x[9] <- length(intersect(a,e))
x[10] <- length(intersect(b,c))
x[11] <- length(intersect(b,d))
x[12] <- length(intersect(b,e))
x[13] <- length(intersect(c,d))
x[14] <- length(intersect(c,e))
x[15] <- length(intersect(d,e))
x[16] <- length(intersect(a,intersect(b,c)))
x[17] <- length(intersect(a,intersect(b,d)))
x[18] <- length(intersect(a,intersect(b,e)))
x[19] <- length(intersect(a,intersect(c,d)))
x[20] <- length(intersect(a,intersect(c,e)))
x[21] <- length(intersect(d,intersect(a,e)))
x[22] <- length(intersect(b,intersect(c,d)))
x[23] <- length(intersect(b,intersect(c,e)))
x[24] <- length(intersect(d,intersect(b,e)))
x[25] <- length(intersect(d,intersect(c,e)))
x[26] <- length(intersect(intersect(a,b),intersect(c,d)))
x[27] <- length(intersect(intersect(a,b),intersect(c,e)))
x[28] <- length(intersect(intersect(a,b),intersect(d,e)))
x[29] <- length(intersect(intersect(a,c),intersect(d,e)))
x[30] <- length(intersect(intersect(b,c),intersect(d,e)))
x[31] <- length(intersect(a,intersect(intersect(b,c),intersect(d,e))))


# rf.genes <- rf.genes[1:(0.1*length(rf.genes))]
# nb.genes <- nb.genes[1:(0.1*length(nb.genes))]

# get union of candidate genes
candidates <- union(union(union(lasso.genes,elnet.genes),union(rf.genes,xgb.genes)),nb.genes)

# require(ltm)
# counts <- data.frame(t(counts))
# output <- data.frame(candidates)
# for (i in 1:length(candidates)) {
#   output$cor[i] <- cor(counts[,which(colnames(counts) %in% candidates[i])],target)
# }
# 
# write.table(candidates,"candidates.txt",quote=F,sep="\t")
# 
# modProbes=probes[is.finite(match(probes,candidates))]
# # Select the corresponding Topological Overlap
# modTOM = TOM[is.finite(match(probes,candidates)), is.finite(match(probes,candidates))]
# dimnames(modTOM) = list(modProbes, modProbes)
# # Export the network into edge and node list files for Cytoscape
# cyt = exportNetworkToCytoscape(modTOM,
#                                edgeFile="CytoEdge.txt",
#                                nodeFile="CytoNode.txt",
#                                weighted = TRUE, threshold = 0.2,nodeNames=modProbes) #threshold was 0.02

library(GGally)
library(network)
library(sna)
library(ggplot2)
require(RColorBrewer)
require(intergraph)

index <- which(probes %in% candidates)
net <- TOM[index,index]
net[net < 0.02] = 0
diag(net) <- 0
discardindx <- which(apply(net,1,sum) == 0)
net <- net[-discardindx,-discardindx]
degree <- apply(net,1,sum)
# hist(degree)
net <- network(net, directed = FALSE)
degree.int <- c()
for (i in 1:length(degree)) {degree.int[i] <- sum(net[i,])}
col <- colorRampPalette(c('grey', 'black'))(length(degree))[rank(degree)]
set.seed(123)
pdf("./gcn.pdf")
ggnet2(net, color = col,node.size = 1.5,edge.size = 0.1) + 
  guides(color = FALSE, size = FALSE)
dev.off()
gcn <- data.frame(genes=probes[index][-discardindx],degree=degree,degree.int = degree.int)
gcn <- gcn[order(gcn$degree,decreasing=TRUE),]
gcn$coeff <- rep(NA,dim(gcn)[1])
for (i in 1:dim(gcn)[1]) {
  if (sum(gcn$genes[i] == lasso.genes) == 1) {
    gcn$coeff[i] <- lasso$lasso.coeffs[which(lasso$lasso.genes %in% gcn$genes[i])]
  }
  else if (sum(gcn$genes[i] == elnet.genes) == 1) {
    gcn$coeff[i] <- elnet$coefficient[which(elnet$name %in% gcn$genes[i])]
  }
}
write.table(gcn,"GCN.txt",sep="\t",quote=FALSE,col.names=NA)

## ---map decision boundaries for gene candidates---

# LASSO candidate
plot(as.numeric(counts[which(probes %in% "GRMZM2G070649"),]),target) # very linear decision boundary

# Elastic net candidates
pdf("./elnet-candidates-interaction.pdf")
plot(as.numeric(counts[which(probes %in% elnet.genes)[3],]), as.numeric(counts[which(probes %in% elnet.genes)[6],]),type="n", xlab=elnet.genes[3], ylab=elnet.genes[6]) 
text(as.numeric(counts[which(probes %in% elnet.genes)[3],]), as.numeric(counts[which(probes %in% elnet.genes)[6],]),target) # very linear decision boundary
abline(a=2.8, b = -0.5, col="red")
dev.off()

# Random forest candidates
plot(as.numeric(counts[which(probes %in% rf.genes)[1],]), as.numeric(counts[which(probes %in% rf.genes)[2],]),type="n") 
text(as.numeric(counts[which(probes %in% rf.genes)[1],]), as.numeric(counts[which(probes %in% rf.genes)[2],]),target) # very linear decision boundary


