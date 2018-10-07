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
counts <- counts[,-c(7,8)]
require(e1071)
require(FSelector)
data <- data.frame(t(counts))
data$target <- as.factor(traits[-c(7,8)])
set.seed(123)
nb.imp <- information.gain(target ~., data = data)
nb.imp <- data.frame(gene=rownames(nb.imp),nb.imp)
nb.imp <- nb.imp[nb.imp[,2] != 0,]
nb.imp <- nb.imp[order(nb.imp[,2],decreasing = TRUE),]
write.table(nb.imp,"nbcimp.txt",quote=F,sep = "\t")

