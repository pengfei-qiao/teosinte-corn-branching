counts <- read.table("counts.txt", sep="\t", header=T, row.names = 1)[1:31]
traits <- rep(1,31)
traits[c(6,18,19,21,25)] <- 0 #Modern corn do not branch
require(edgeR)
y <- DGEList(counts = counts)
y <- calcNormFactors(y)
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
#Filter out low read genes
filter_low_counts <- function (x, cutoff = 10) {
  sum(x < cutoff)
}
counts <- data.frame(cpm(y)[-which(apply(cpm(y), MARGIN = 1, FUN = filter_low_counts) > dim(counts)[2]*.9) ,])

tb1 = "AC233950.1_FG002"
counts[match(tb1,rownames(counts)),]

# PCA
counts.scaled <- scale(counts)
result <- prcomp(t(counts.scaled))
result$sdev^2/sum(result$sdev^2)*100
plot(result$x[,1], result$x[,2],type="n")
text(result$x[,1], result$x[,2],traits)
dev.off()

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

## --- LASSO ---
data <- t(counts)
target <- traits[-c(7,8)]
# split into train (validation) and test sets
set.seed(123)
test <- sample(1:nrow(data),1/10*nrow(data)) # leave out test set
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
print(lasso.genes <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)) # https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame, the reason for +1 is that the @i method indexes from 0 for the intercept but @Dimnames[[1]] starts at 1
print(lasso.Rsq <- lasso.fit$glmnet.fit$dev.ratio[which(lasso.fit$glmnet.fit$lambda==lasso.fit$lambda.min)])


## --- Elastic net ---
data <- t(counts)
target <- traits[-c(7,8)]
# split into train (validation) and test sets
set.seed(123)
test <- sample(1:nrow(data),1/10*nrow(data)) # leave out test set
data.test <- data[test,] # test set
data <- data[-test,] # training and validation set
ptm <- proc.time()
require(glmnet)
elnet.fit = cv.glmnet(data,as.factor(target[-test]), alpha = 0.5, lambda = 10^seq(-10,0,length=100), nfolds = nrow(data),family = "binomial",type.measure = "class")
print(elnet.time <- proc.time() - ptm)
plot(elnet.fit,xvar="lambda",main="Elastic net")
elnet.pred <- predict(elnet.fit,data.test,s=elnet.fit$lambda.min,"class")
print(elnet.acc <- mean(elnet.pred == as.numeric(target[test])))

tmp_coeffs <- coef(elnet.fit,s=elnet.fit$lambda.min)
print(elnet.genes <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)) # https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame, the reason for +1 is that the @i method indexes from 0 for the intercept but @Dimnames[[1]] starts at 1
print(elnet.Rsq <- elnet.fit$glmnet.fit$dev.ratio[which(elnet.fit$glmnet.fit$lambda==elnet.fit$lambda.min)])


## --- Random Forest ---
require(randomForest)
data <- data.frame(t(counts))
data$target <- as.factor(traits[-c(7,8)])
# set.seed(123)
# test <- sample(1:nrow(data),1/5*nrow(data)) # leave out test set
# data.test <- data[test,] # test set
# data <- data[-test,] # training and validation set
ptm <- proc.time()
set.seed(123)
rf.cv <- rfcv(data[,-dim(data)[2]],data$target,cv.fold=nrow(data)) # to select mtry for randomForest
with(rf.cv, plot(n.var, error.cv, log="x", type="o", lwd=2))
set.seed(123)
rf = randomForest(target ~., data = data, mtry = 6, importance = TRUE, ntree = 2000)
print(rf.time <- proc.time() - ptm)
rf.imp <- importance(rf,scale = FALSE) # why scaled=F? https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-307
rf.imp <- rf.imp[which(apply(rf.imp[,3:4],1,sum) != 0),]
rf.imp <- rf.imp[order(rf.imp[,3],decreasing = TRUE),] # why sort by permutation decrease in accuracy? http://explained.ai/rf-importance/index.html
varImpPlot(rf, scale = FALSE)
# rf.pred <- predict(rf,newdata=data.test,type="response")
# print(rf.acc <- mean(rf.pred == as.numeric(target[test])))
rf # OOB estimate of error rate: 17.24%
#getTree(rf,1)


## --- XGBoost ---
require(xgboost)
data <- t(counts)
target <- traits[-c(7,8)]
# split into train (validation) and test sets
set.seed(123)
test <- sample(1:nrow(data),1/10*nrow(data)) # leave out test set
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
  plot(e$iter, e$train_mlogloss_mean, type = "o",col="blue",ylim=c(0,max(e$test_mlogloss_mean)))
  lines(e$iter,e$test_mlogloss_mean, type = "o",col="red")
  print(e$iter[which(e$test_mlogloss_mean == min(e$test_mlogloss_mean))])
}

xgbtune(20,4,0.1,0)
# xgbtune(200,4,0.01,0)
# xgbtune(50,4,0.1,100)

nround <- 8 # number of rounds/trees
ptm <- proc.time()
set.seed(123)
xgbtree <- xgb.train(params=param, 
                   data = dtrain, 
                   nrounds = nround, 
                   max.depth = 4,     # tree depth (not the same as interaction.depth in gbm!)
                   eta = 0.1)           # shrinkage parameter
print(xgbtree.time <- proc.time() - ptm)    

xgbtree.prob <- predict(xgbtree, data.test)
xgbtree.prob <- t(matrix(xgbtree.prob, nrow=numberOfClasses, ncol=nrow(data.test)) ) # need to convert it to a matrix
xgbtree.pred <-apply(xgbtree.prob, 1, which.max) - 1 # use the class with maximum probability as prediction
table(xgbtree.pred, target[test])
print(xgbtree.acc <- mean(xgbtree.pred == as.numeric(target[test])))  # classification accuracy on test set
print(xgbtree.imp <- xgb.importance(colnames(data), model=xgbtree))
xgb.plot.importance(xgbtree.imp)
