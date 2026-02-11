## Package loading

library("limma")
library("stringr")
library("edgeR")
library(readxl)
library(plotly)
library(heatmaply)
library(dplyr)
library(data.table)


# Data import & preprocessing

## Import Salmon quant data with catchSalmon
### Import count data frame and specify first column as row names

paths <- list.dirs("../SALMON/output",full.names = TRUE, recursive = FALSE)
paths <- paste0(paths, "/")
raw <- catchSalmon(paths)
colnames(raw$counts) <- basename(colnames(raw$counts))

## Import sample metadata (phenodata) for the dataset

phenodata <- read.delim("analytes.txt")
phenodata

SID.Group <- phenodata[ ,c("SampleID","group")]
SID.raw <- unlist(colnames(raw$counts))
#this command combines the data from two data structures (SID.raw and SID.Group), merges them based on their first columns, removes duplicate rows, and stores the result in the variable translat.o
translat.o <- unique(merge(x = as.matrix(SID.raw), y = SID.Group, by.x = 1, by.y = 1, all.x = TRUE, all.y = FALSE, sort = FALSE, suffixes = c(".x",".y"), no.dups = TRUE, incomparables = NA))

sampleID <- translat.o$V1
group <- translat.o$group

## Fix raw$counts colnames

colnames.old <- colnames(raw$counts)
check.order.colnames <- cbind(colnames.old, translat.o)
check.order.colnames

colnames(raw$counts) <- sampleID

head(raw$counts)

raw.leaves <- raw
colnames(raw$counts[ , 1:48])
raw.leaves$counts <- raw$counts[ , 1:48]
group.leaves <- group[1:48]


## Scaling counts

#Calculate scaled counts and create DGEList objects of raw and scaled Salmon counts

scaled_counts <- raw.leaves$counts/raw.leaves$annotation$Overdispersion

scaled <- DGEList(counts=scaled_counts, group=group.leaves)
unscaled <- DGEList(counts=raw.leaves$counts, group=group.leaves)

levels(as.factor(group.leaves))


### Color palette for graphs


col = c("darkred","orangered","darkgreen", "#B0FEA0","orange", "olivedrab", "tomato2", "cyan", "black","brown", "green", "gold", "#ADD8E6", "cyan","ivory4", "brown")

col <- col[as.numeric(as.factor(group.leaves))]

col

### Plot non-scaled and scaled count densities

opar <- par()
par(mfrow=c(1,2), cex = 0.6)
nsamples <- ncol(scaled_counts)

lcpm <- log(as.matrix(unscaled),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,0.25), las=2, main="", xlab="")
title(main="A. Raw Salmon counts", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

lcpm <- log(as.matrix(scaled),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,0.25), las=2, main="", xlab="")
title(main="B. Scaled Salmon counts", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}


par(opar)


keep_scaled <- filterByExpr(scaled, group=group, min.count=20)
scaled_f <- scaled[keep_scaled, ,keep.lib.sizes=FALSE]

cat("All transcripts in dataset: ", dim(scaled$counts)[1], "\n",
    "Transcripts in dataset after removing low expressed ones: ", dim(scaled_f$counts)[1], 
    " (", round(dim(scaled_f$counts)[1]/dim(scaled$counts)[1]*100,2), "%)", "\n", sep = "")

## Calculate normalisation factors


scaled_f <- calcNormFactors(scaled_f, method = "TMM")

### QC boxplots of scaled counts after filtering and normalisation

opar <- par()
par(mfrow=c(1,3), cex = 0.6)
boxplot(log(scaled$counts+1,10), las=2, ylab="log10(scaled counts)", col=col)
boxplot(log(scaled_f$counts+1,10), las=2, ylab="log10(filtered scaled counts)", col=col)
lcpm <- cpm(scaled_f, log=TRUE)
boxplot(lcpm, ylab="filtered and normalised CPMs", col=col, las=2)
par(opar)


### Density plots of scaled counts before and after filtering

opar <- par()
par(mfrow=c(1,2), cex = 0.6)
nsamples <- ncol(scaled_counts)

lcpm <- log(as.matrix(scaled),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,0.3), las=2, main="", xlab="")
title(main="A. Before filtering", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

lcpm <- log(as.matrix(scaled_f),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,0.6), las=2, main="", xlab="")
title(main="B. After filtering", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
par(opar)

### MDS plots


op <- par(pty = "s")
plotMDS(scaled_f, labels=colnames(scaled_f), col = col, cex = 0.4, main = "MDS plot of scaled normalised CPMs")
par(op)



op <- par(pty = "s")
plotMDS(scaled_f, labels=group.leaves, col = col, cex = 0.6, main = "MDS plot of scaled normalised CPMs")
par(op)


  ### MD plots

  for (i in 1:length(colnames(scaled_counts))) {
    plotMD(scaled_f, column = i, main = paste(colnames(scaled_counts)[i], "scaled"))
  }

  

## Statistical model
  
### Create model matrix
  

design <- model.matrix(~0+group.leaves, data = scaled_f$samples)
  
cbind(colnames(design), levels(as.factor(group.leaves)))
  
colnames(design) <- levels(as.factor(group.leaves))
  
### Estimate dispersion from scaled counts

scaled_f <- estimateDisp(scaled_f,design)
scaled_f$common.dispersion

#### Plot of biological coefficient of variation

plotBCV(scaled_f, main = "BCV of scaled counts")

## Contrast matrix for differential expression

contrasts=makeContrasts("Rywal_PVY_DEX_A-Rywal_mock_DEX",
                        "Rywal_PVY_DEX_B-Rywal_mock_DEX",
                         "Rywal_PVY_DEX_B-Rywal_PVY_DEX_A",
                          "L4_mock_DEX-Rywal_mock_DEX",
                          "L4_PVY_DEX_A-Rywal_PVY_DEX_A",
                        "L4_PVY_DEX_A-Rywal_PVY_DEX_A"
                        "L4_PVY_DEX_B-Rywal_PVY_DEX_B",
                        "L4_PVY_DEX_A-L4_mock_DEX",
                        "L4_PVY_DEX_B-L4_mock_DEX",
                        "L4_PVY_DEX_B-L4_PVY_DEX_A",
                        levels=design)



### limma-voom LM fit
#- using the function from edgeR that is better handling 0 counts
#Here using the default normalize.method="none" as the data is already TMM-normalised.

#older: fitL <- edgeR::voomLmFit(scaled_f, design, block = group.leaves, sample.weights = TRUE, plot= TRUE, normalize.method="none")
fitL <- edgeR::voomLmFit(scaled_f, design, sample.weights = TRUE,normalize.method="none",plot = TRUE)
fitL2 = contrasts.fit(fitL, contrasts)
fitL2 <- eBayes(fitL2)

plotSA(fitL2, main = "SA plot after limma")


#Quick check for DEGs:
logFCcut <- 1
dt <- decideTests(fitL2, lfc=logFCcut)

op <- par(mar = c(7,3,3,3))
barplot(summary(dt)[c(1,3),], beside= TRUE, las= 2, cex.axis = 0.5, cex.names = 0.48, legend = TRUE, col= c("darkgreen", "darkred"))
abline(v = 24.5, col = "grey", lwd = 2)
abline(v = 48.5, col = "grey", lwd = 2, lty = "dashed")

par(op)


# Initialize a list to store the results
results_list <- list()


# Loop over all coefficients
for (i in 1:length(colnames(fitL2$contrasts))) {
  results_list[[i]] <- topTable(fitL2, coef=i, number=1000000, sort.by="none")[, c("logFC", "adj.P.Val")]
  colnames(results_list[[i]]) <- c(paste0(colnames(fitL2$contrasts)[i], "_logFC"), paste0(colnames(fitL2$contrasts)[i], "_adj.P.Val"))
}

results_matrix <- do.call(cbind, results_list)
colnames(results_matrix)



### Export DEG results into file

#add CPMs to stats
results_matrix.addcpms <- merge(results_matrix, scaled_f$counts, by.x="row.names", by.y="row.names", all.x= TRUE, all.y= FALSE, sort= FALSE)

head(results_matrix.addcpms)


#Export results to file

write.table(results_matrix.addcpms, file="RNASeq_5dpi.txt", sep="\t", quote=TRUE, row.names=FALSE)








