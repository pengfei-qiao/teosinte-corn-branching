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

tb1 = "AC233950.1_FG002"
counts[match(tb1,rownames(counts)),]
probes <- rownames(counts)


# hclust for outlier detection
require(WGCNA)
counts.log <- log(counts+1, base=2)
A <- adjacency(counts.log, type="distance") #Cluster samples usually use Euclidean distance, clustering genes use correlation
k <- as.numeric(apply(A,2,sum)) - 1
Z.k <- scale(k)
#Designate samples as outlying based on Z.k
thresholdZ.k <- -2.5
outlierColor <- ifelse(Z.k < thresholdZ.k, "red", "black")
#Sample 7, 8 is the outlier, discard
counts <- counts[,-c(7,8)]
#Calculate the cluster tree using flashClust
library(flashClust)
sampleTree <- flashClust(as.dist(1-A), method = "average")
plot(sampleTree, main = "Sample dendrogram and trait heatmap")

### ML

## --- PCA for dimension reduction
# data(USArrests)
# pr.out <- prcomp(USArrests,sclae=T)
# scale(USArrests) %*% pr.out$rotation == pr.out$x
result <- prcomp(t(counts), scale = TRUE)
pdf("./figures/PC-cumsum.pdf")
plot(cumsum(result$sdev^2/sum(result$sdev^2)*100),type="o",lwd=2,xlab = "PCs",ylab="Cumulative Variance Expained (%)")
abline(a = 85.22416, b=0, col="red")
# identical(scale(t(counts)) %*% result$rotation, result$x)
dev.off()
pcacounts <- result$x[1:29,1:20]

# get cosine of angle between each transcript's loading and PC
PCs <- matrix(0, nrow=29,ncol=29)
diag(PCs) <- 1
pc.cor <- result$rotation %*% PCs
for (i in 1:nrow(pc.cor)) {
  pc.cor[i,] <- pc.cor[i,]/norm(result$rotation[i,],type="2")
}


## --- LASSO ---
data <- pcacounts
target <- traits[-c(7,8)]
# split into train (validation) and test sets
set.seed(123)
test <- sample(1:nrow(data),1/5*nrow(data)) # leave out test set
data.test <- data[test,] # test set
data <- data[-test,] # training and validation set
ptm <- proc.time()
require(glmnet)
lasso.fit = cv.glmnet(data,as.factor(target[-test]),alpha = 1, lambda = 10^seq(-10,0,length=100), nfolds = nrow(data),family = "binomial",type.measure = "class")
print(lasso.time <- proc.time() - ptm)
plot(lasso.fit, xvar="lambda",main="LASSO")
lasso.pred <- predict(lasso.fit,data.test,s=lasso.fit$lambda.min,"class")
print(lasso.acc <- mean(lasso.pred == as.numeric(target[test])))

tmp_coeffs <- coef(lasso.fit,s=lasso.fit$lambda.min)
print(lasso.pcs <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)) # https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame, the reason for +1 is that the @i method indexes from 0 for the intercept but @Dimnames[[1]] starts at 1
print(lasso.Rsq <- lasso.fit$glmnet.fit$dev.ratio[which(lasso.fit$glmnet.fit$lambda==lasso.fit$lambda.min)])
lasso.pcs <- lasso.pcs[order(abs(lasso.pcs$coefficient),decreasing = TRUE),]

lasso.genes <- c()
lasso.coeffs <- c()
for (i in 2:nrow(lasso.pcs)) {
  pc <- as.character(lasso.pcs$name[i])
  pcindx <- as.numeric(strsplit(pc,"PC")[[1]][2])
  genes <- probes[which(abs(pc.cor[,pcindx]) > 0.7)]
  lasso.coeffs <- c(lasso.coeffs,lasso.pcs$coefficient[i]*pc.cor[which(probes %in% genes),pcindx])
  lasso.genes <- c(lasso.genes,genes)
}
lasso.candidates <- data.frame(lasso.genes, lasso.coeffs)[order(abs(lasso.coeffs),decreasing = TRUE),]

hardboundary <- c()
for (i in 1:length(lasso.genes)){
  teosinte.lower = min(counts[which(probes %in% lasso.genes[i]),][which(target != 0)])
  teosinte.upper = max(counts[which(probes %in% lasso.genes[i]),][which(target != 0)])
  maize.lower = min(counts[which(probes %in% lasso.genes[i]),][which(target == 0)])
  maize.upper = max(counts[which(probes %in% lasso.genes[i]),][which(target == 0)])
  if (teosinte.lower > maize.upper) {hardboundary <-c(hardboundary,i)}
  if (teosinte.upper < maize.lower) {hardboundary <- c(hardboundary,i)}
}
lasso.candidates$hardboundary <- rep(0,length(lasso.genes))
lasso.candidates$hardboundary[hardboundary] <- 1

write.table(lasso.candidates,"lasso.txt",quote=F,sep="\t")
plot(as.numeric(counts[which(probes %in% "GRMZM2G070649"),]),target)
max(counts[which(probes %in% "GRMZM2G070649"),][which(target != 0)])
min(counts[which(probes %in% "GRMZM2G070649"),][which(target == 0)])



# ## --- Elastic net --- features are no longer correlated, so LASSO outperforms elastic net
# data <- pcacounts
# target <- traits[-c(7,8)]
# # split into train (validation) and test sets
# set.seed(123)
# test <- sample(1:nrow(data),1/5*nrow(data)) # leave out test set
# data.test <- data[test,] # test set
# data <- data[-test,] # training and validation set
# ptm <- proc.time()
# require(glmnet)
# elnet.fit = cv.glmnet(data,as.factor(target[-test]), alpha = 0.5, lambda = 10^seq(0,-10,length=100), nfolds = nrow(data),family = "binomial",type.measure = "class")
# print(elnet.time <- proc.time() - ptm)
# plot(elnet.fit,xvar="lambda",main="Elastic net")
# elnet.pred <- predict(elnet.fit,data.test,s=elnet.fit$lambda.min,"class")
# print(elnet.acc <- mean(elnet.pred == as.numeric(target[test])))
# 
# tmp_coeffs <- coef(elnet.fit,s=elnet.fit$lambda.min)
# print(elnet.pcs <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)) # https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame, the reason for +1 is that the @i method indexes from 0 for the intercept but @Dimnames[[1]] starts at 1
# elnet.pcs <- elnet.pcs[order(abs(elnet.pcs[,2]),decreasing = TRUE),]
# print(elnet.Rsq <- elnet.fit$glmnet.fit$dev.ratio[which(elnet.fit$glmnet.fit$lambda==elnet.fit$lambda.min)])
# 
# elnet.genes <- c()
# for (i in 2:nrow(elnet.pcs)) {
#   pc <- as.character(elnet.pcs$name[i])
#   elnet.genes <- c(elnet.genes,rownames(pc.cor)[which(abs(pc.cor[,which(colnames(pc.cor) %in% pc)]) > 0.9)])
# }
# 
# write.table(elnet.genes,"elnet.txt",quote=F,sep="\t")

# ## --- Random Forest --- features are no longer correlated, so no rf
# require(randomForest)
# data <- data.frame(t(counts))
# data$target <- as.factor(traits[-c(7,8)])
# # set.seed(123)
# # test <- sample(1:nrow(data),1/5*nrow(data)) # leave out test set
# # data.test <- data[test,] # test set
# # data <- data[-test,] # training and validation set
# ptm <- proc.time()
# set.seed(123)
# rf.cv <- rfcv(data[,-dim(data)[2]],data$target,cv.fold=nrow(data)) # to select mtry for randomForest
# with(rf.cv, plot(n.var, error.cv, log="x", type="o", lwd=2))
# pdf("./figures/rf.loocv.pdf")
# par(mfrow=c(1,1))
# plot(rf.cv$n.var, rf.cv$error.cv, log="x", type="o", lwd=2, xlab = "Number of Variables", ylab = "LOOCV")
# dev.off()
# data <- t(counts)
# target <- traits[-c(7,8)]
# # split into train (validation) and test sets
# set.seed(123)
# test <- sample(1:nrow(data),1/5*nrow(data)) # leave out test set
# data.test <- data[test,] # test set
# data <- data[-test,] # training and validation set
# set.seed(123)
# rf = randomForest(data, as.factor(target[-test]), mtry = 39, importance = TRUE, ntree = 2000)
# print(rf.time <- proc.time() - ptm)
# rf.imp <- importance(rf,scale = FALSE) # why scaled=F? https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-307
# rf.imp <- rf.imp[which(apply(rf.imp[,3:4],1,sum) != 0),]
# rf.imp <- rf.imp[order(rf.imp[,3],decreasing = TRUE),] # why sort by permutation decrease in accuracy? http://explained.ai/rf-importance/index.html
# write.table(rf.imp,"rfimp.txt",quote = F, sep = "\t")
# rf.genes <- rownames(rf.imp)
# pdf("./figures/rf.varimpplot.pdf")
# par(mfrow=c(1,1))
# varImpPlot(rf, scale = FALSE)
# dev.off()
# rf.pred <- predict(rf,newdata=data.test,type="response")
# print(rf.acc <- mean(rf.pred == as.numeric(target[test])))
# rf # OOB estimate of error rate: 17.24%
# #getTree(rf,1)


## --- XGBoost ---
require(xgboost)
data <- pcacounts
target <- traits[-c(7,8)]
# split into train (validation) and test sets
set.seed(123)
test <- sample(1:nrow(data),1/5*nrow(data)) # leave out test set
data.test <- data[test,] # test set
data <- data[-test,]
ptm <- proc.time()
dtrain <- xgb.DMatrix(data=data,label=target[-test])
numberOfClasses <- 2

# algorithm parameters for XGBoost
param <- list("objective" = "multi:softprob",
              "eval_metric" = "mlogloss",
              "num_class" = numberOfClasses)

# tune parameters by cv
xgbtune <- function(a,b,c,d) {
  xgbcv <- xgb.cv(data = dtrain, 
                params = param,
                nfold=nrow(data),
                nrounds = a, 
                max.depth = b,     # tree depth (not the same as interaction.depth in gbm!)
                eta = c,
                gamma = d)
  e <- xgbcv$evaluation_log
  pdf("./figures/xgboostcv.pdf")
  plot(e$iter, e$train_mlogloss_mean, type = "o",col="blue",ylim=c(0,max(e$test_mlogloss_mean)),lwd = 2,xlab="Number of Trees", ylab="LOOCV Error Rate")
  lines(e$iter,e$test_mlogloss_mean, type = "o",col="red")  
  legend("topright",legend=c("training error","test error"),col=c("blue","red"),lty=1,cex=1.2)
  dev.off()
  print(e$iter[which(e$test_mlogloss_mean == min(e$test_mlogloss_mean))])
  print(min(e$test_mlogloss_mean))
}

xgbtune(100,4,0.5,0)

nround <- 10 # number of rounds/trees
ptm <- proc.time()
set.seed(123)
xgbtree <- xgb.train(params=param, 
                   data = dtrain, 
                   nrounds = nround, 
                   max.depth = 4,     # tree depth (not the same as interaction.depth in gbm!)
                   eta = 0.5)           # shrinkage parameter
print(xgbtree.time <- proc.time() - ptm)    

xgbtree.prob <- predict(xgbtree, data.test)
xgbtree.prob <- t(matrix(xgbtree.prob, nrow=numberOfClasses, ncol=nrow(data.test)) ) # need to convert it to a matrix
xgbtree.pred <-apply(xgbtree.prob, 1, which.max) - 1 # use the class with maximum probability as prediction
table(xgbtree.pred, target[test])
print(xgbtree.acc <- mean(xgbtree.pred == as.numeric(target[test])))  # classification accuracy on test set
print(xgbtree.imp <- xgb.importance(colnames(data), model=xgbtree))
xgb.plot.importance(xgbtree.imp)


xgb.genes <- c()
for (i in 1:nrow(xgbtree.imp)) {
  pc <- as.character(xgbtree.imp$Feature[i])
  pcindx <- as.numeric(strsplit(pc,"PC")[[1]][2])
  genes <- probes[which(abs(pc.cor[,pcindx]) > 0.7)]
  xgb.genes <- c(xgb.genes,genes)
}

write.table(xgb.genes,"xgboost.txt",quote=F,sep="\t")

## ---Naive Bayes---
require(e1071)
data <- pcacounts
target <- traits[-c(7,8)]
# split into train (validation) and test sets
set.seed(123)
test <- sample(1:nrow(data),1/5*nrow(data)) # leave out test set
data.test <- data[test,] # test set
data <- data[-test,] # training and validation set

# having problem installing caret, and implementing LOOCV myself
nb.acc.loocv <- c()
for (i in 1:nrow(data)) {
  nb.fit <- naiveBayes(data[-i,], as.factor(target[-test])[-i])
  nb.pred <- predict(nb.fit,t(data[i,]))
  nb.acc.loocv[i] <- mean(nb.pred == target[i])
}
print(sum(nb.acc.loocv)/length(nb.acc.loocv)) # 0.8333333

# calculate information gain for each feature - stackoverflow when ran in RStudio. So run in terminal:Rscript --max-ppsize=500000 infogain.R
require(FSelector)
data <- data.frame(pcacounts)
data$target <- as.factor(traits[-c(7,8)])
set.seed(123)
nb.imp <- information.gain(target ~., data = data)
nb.imp <- data.frame(gene=rownames(nb.imp),nb.imp)
nb.imp <- nb.imp[nb.imp[,2] != 0,]
nb.imp <- nb.imp[order(nb.imp[,2],decreasing = TRUE),]
write.table(nb.imp,"nbcimp.txt",quote=F,sep = "\t")
nb.imp <- read.table("nbcimp.txt")
nb.genes <- as.character(nb.imp[,1]) # only PC4 with non-zero info gain, 0.3664735. so subset of LASSO's genes.

# # feature selection using permutation
# nb.p.accu <- c()
# set.seed(123)
# for (i in 1:(ncol(data)-1)) {
#   data2 <- data
#   data2[,i] <- sample(data[,i],replace = FALSE)
#   nb.fit <- naiveBayes(target ~., data=data2)
#   nb.pred <- predict(nb.fit,data.test)
#   print(i)
#   nb.p.accu[i] <- mean(nb.pred == data$target[test])
# }

# prediction accuracy on test set
nb.pred <- predict(nb.fit,data.test)
print(nb.acc <- mean(nb.pred == target[test]))


# ## ---SVM--- not suitable to find an optimal decision boundary
# require(e1071)
# data <- data.frame(t(counts))
# data$target <- as.factor(traits[-c(7,8)])
# # split into train (validation) and test sets
# set.seed(123)
# test <- sample(1:nrow(data),1/10*nrow(data)) # leave out test set
# data.test <- data[test,] # test set
# data <- data[-test,]
# 
# # # PCA to determine decision boundary for kernel choice
# counts.scaled <- scale(counts)
# result <- prcomp(t(counts.scaled))
# result$sdev^2/sum(result$sdev^2)*100
# pdf("./figures/pca.pdf")
# plot(result$x[,1], result$x[,2],type="n",xlab = "PC1 - 33.34%", ylab = "PC2 - 7.99%")
# text(result$x[,1], result$x[,2],traits)
# dev.off()
#  
# svm.fit <- svm(target ~., data=data, kernel = "radial", gamma = 0.1, cost = 10, scale=TRUE)
# predict(svm.fit, data.test)
# 
# tune.out = tune(svm, target~., data=data, kernel="radial",
#                 ranges=list(cost=c(0.001,0.01,0.1,1,5,10,100)),
#                 gamma=c(0.5,1,2,3,4), scale = T )
# bestmod = tune.out$best.model
# bestmod


# look at overlaps between features selected from different algorithms
# tmp_coeffs <- coef(lasso.fit,s=lasso.fit$lambda.min)
# lasso.genes <- as.character(data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)[,1])[2]
# tmp_coeffs <- coef(elnet.fit,s=elnet.fit$lambda.min)
# elnet.genes <- as.character(data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)[,1])[2:122]
# rf.genes <- rownames(rf.imp)
# xgb.genes <- "GRMZM2G028594"
# nb.genes <- as.character(nb.imp[,1])

# require(VennDiagram)
# require(venn)
# pdf("./figures/venn.pdf")
# draw.quintuple.venn( 1,  121 ,4148  ,  1 , 428  ,  1 ,   1  ,  0  ,  1 ,  98,    1 ,  91 ,   1 , 299 ,   1 ,   1 ,   0 ,   1 ,   0   , 1   , 0   , 1 ,  76,    1,    1, 0  ,  1  ,  0  ,  0   , 1 ,   0,
#                      category = c("LASSO", "Elastic Net", "Random Forest", "Boosted Tree", "Naive Bayes Classifier"))
# dev.off()
# 
# # to test if the non-overlapping features are highly correlated
# moduleLabelsManual <- cutreeDynamic(dendro=geneTree, distM = dissTOM, method = "hybrid",
#                                      deepSplit = 2, pamRespectsDendro = F, minClusterSize = 30)
# moduleColorsManual <- labels2colors(moduleLabelsManual)
# lasso.module <- moduleColorsManual[which(probes %in% lasso.genes)]
# elnet.module <- moduleColorsManual[which(probes %in% elnet.genes)]
# rf.module <- moduleColorsManual[which(probes %in% rf.genes)]
# xgb.module <- moduleColorsManual[which(probes %in% xgb.genes)]
# nb.module <- moduleColorsManual[which(probes %in% nb.genes)]
# pdf("./figures/venn-modules.pdf")
# draw.quintuple.venn(1 , 121, 4148 ,   1  ,428,    1   , 1    ,0  ,  1   ,14 ,   1   ,14 ,   1 ,  18  ,  1  ,  1,    0  ,  1,    0 ,   1,    0  ,  1 ,  14  ,  1  ,  1, 0  ,  1 ,   0,    0,    1   , 0,
#                     category = c("LASSO", "Elastic Net", "Random Forest", "Boosted Tree", "Naive Bayes Classifier"))
# dev.off()
# 
# x <- c()
# x[1] <- length(a)
# x[2] <- length(b)
# x[3] <- length(c)
# x[4] <- length(d)
# x[5] <- length(e)
# x[6] <- length(intersect(a,b))
# x[7] <- length(intersect(a,c))
# x[8] <- length(intersect(a,d))
# x[9] <- length(intersect(a,e))
# x[10] <- length(intersect(b,c))
# x[11] <- length(intersect(b,d))
# x[12] <- length(intersect(b,e))
# x[13] <- length(intersect(c,d))
# x[14] <- length(intersect(c,e))
# x[15] <- length(intersect(d,e))
# x[16] <- length(intersect(a,intersect(b,c)))
# x[17] <- length(intersect(a,intersect(b,d)))
# x[18] <- length(intersect(a,intersect(b,e)))
# x[19] <- length(intersect(a,intersect(c,d)))
# x[20] <- length(intersect(a,intersect(c,e)))
# x[21] <- length(intersect(d,intersect(a,e)))
# x[22] <- length(intersect(b,intersect(c,d)))
# x[23] <- length(intersect(b,intersect(c,e)))
# x[24] <- length(intersect(d,intersect(b,e)))
# x[25] <- length(intersect(d,intersect(c,e)))
# x[26] <- length(intersect(intersect(a,b),intersect(c,d)))
# x[27] <- length(intersect(intersect(a,b),intersect(c,e)))
# x[28] <- length(intersect(intersect(a,b),intersect(d,e)))
# x[29] <- length(intersect(intersect(a,c),intersect(d,e)))
# x[30] <- length(intersect(intersect(b,c),intersect(d,e)))
# x[31] <- length(intersect(a,intersect(intersect(b,c),intersect(d,e))))
# 
# ## ---map decision boundaries for gene candidates---
# 
# # LASSO candidate
# plot(as.numeric(counts[which(probes %in% lasso.genes),]),target) # very linear decision boundary
# 
# # Elastic net candidates
# plot(as.numeric(counts[which(probes %in% elnet.genes)[1],]), as.numeric(counts[which(probes %in% elnet.genes)[2],]),type="n") 
# text(as.numeric(counts[which(probes %in% elnet.genes)[1],]), as.numeric(counts[which(probes %in% elnet.genes)[2],]),target) # very linear decision boundary
# 
# # Random forest candidates
# plot(as.numeric(counts[which(probes %in% rf.genes)[1],]), as.numeric(counts[which(probes %in% rf.genes)[2],]),type="n") 
# text(as.numeric(counts[which(probes %in% rf.genes)[1],]), as.numeric(counts[which(probes %in% rf.genes)[2],]),target) # very linear decision boundary

