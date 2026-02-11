library(limma)
library(edgeR)

setwd("../A2_srna-seq/sRNA_counting")

x <- read.table("sRNAs_table.txt",header=TRUE,row.names="sRNAs",sep="\t", stringsAsFactors=FALSE)
group <- factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6))
y <- DGEList(counts=x,group=group,lib.size=c(21485098,20539094,17608510,20494169,20307240,11788921,20934214,12424045,11090105,16954680,16504589,23833313,10783005,10641559,9292463,10913374,10318003,7094812))
keep.exprs <- filterByExpr(y, group=group, min.count=20)
y1 <- y[keep.exprs, , keep.lib.sizes=TRUE]
y1 <- calcNormFactors(y1, method = "TMM")
col=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6)
plotMDS(y1, col=col, cex=0.6)


lcpm <- cpm(y1, log=TRUE)
boxplot(lcpm, ylab="filtered and normalised CPMs", col=col, las=2)

design = cbind(Rywal_mock=c(1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), Rywal_A=c(0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0), Rywal_B=c(0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0), NahG_mock=c(0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0), NahG_A=c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0), NahG_B=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1))


v <- voom(y1,design,plot=TRUE)
fit <- lmFit(v,design)

contrastMatrix = makeContrasts("Rywal_A-Rywal_mock","Rywal_B-Rywal_mock","Rywal_B-Rywal_A","NahG_A-NahG_mock","NahG_B-NahG_mock","NahG_B-NahG_A","NahG_mock-Rywal_mock", "NahG_A-Rywal_A",
                               "NahG_B-Rywal_B",levels=design)
fit2 = contrasts.fit(fit, contrastMatrix)
fit2 <- eBayes(fit2)



results1 <- topTable(fit2, coef=1, number=1000000, sort.by="none")
results2 <- topTable(fit2, coef=2, number=1000000, sort.by="none")
results3 <- topTable(fit2, coef=3, number=1000000, sort.by="none")
results4 <- topTable(fit2, coef=4, number=1000000, sort.by="none")
results5 <- topTable(fit2, coef=5, number=1000000, sort.by="none")
results6 <- topTable(fit2, coef=6, number=1000000, sort.by="none")
results7 <- topTable(fit2, coef=7, number=1000000, sort.by="none")
results8 <- topTable(fit2, coef=8, number=1000000, sort.by="none")
results9 <- topTable(fit2, coef=9, number=1000000, sort.by="none")

all(rownames(results1)==rownames(results2))
all(rownames(results2)==rownames(results3))
all(rownames(results3)==rownames(results4))
all(rownames(results4)==rownames(results5))
all(rownames(results5)==rownames(results6))
all(rownames(results6)==rownames(results7))
all(rownames(results7)==rownames(results8))
all(rownames(results8)==rownames(results9))


#make one expression matrix with logFCs and adj.P.vals
results <- cbind(results1[,1], results1[,5], results2[,1], results2[,5], results3[,1], results3[,5],
                 results4[,1], results4[,5], results5[,1], results5[,5], results6[,1], results6[,5],
                 results7[,1], results7[,5], results8[,1], results8[,5], results9[,1], results9[,5])




colnames(results) <- c(paste(colnames(contrastMatrix)[1], " ", colnames(results1[1])),					
                       paste(colnames(contrastMatrix)[1], " ", colnames(results1[5])),					
                       paste(colnames(contrastMatrix)[2], " ", colnames(results2[1])),					
                       paste(colnames(contrastMatrix)[2], " ", colnames(results2[5])),					
                       paste(colnames(contrastMatrix)[3], " ", colnames(results3[1])),					
                       paste(colnames(contrastMatrix)[3], " ", colnames(results3[5])),					
                       paste(colnames(contrastMatrix)[4], " ", colnames(results4[1])),					
                       paste(colnames(contrastMatrix)[4], " ", colnames(results4[5])),					
                       paste(colnames(contrastMatrix)[5], " ", colnames(results5[1])),					
                       paste(colnames(contrastMatrix)[5], " ", colnames(results5[5])),					
                       paste(colnames(contrastMatrix)[6], " ", colnames(results6[1])),					
                       paste(colnames(contrastMatrix)[6], " ", colnames(results6[5])),					
                       paste(colnames(contrastMatrix)[7], " ", colnames(results7[1])),					
                       paste(colnames(contrastMatrix)[7], " ", colnames(results7[5])),					
                       paste(colnames(contrastMatrix)[8], " ", colnames(results8[1])),					
                       paste(colnames(contrastMatrix)[8], " ", colnames(results8[5])),					
                       paste(colnames(contrastMatrix)[9], " ", colnames(results9[1])),					
                       paste(colnames(contrastMatrix)[9], " ", colnames(results9[5])))

rownames(results) <- rownames(results1)

length(rownames(y1))==length(results[,1])
all(rownames(y1) == results[,1]) # if FALSE, have to do merge, not cbind!
results.raw <- merge(results, y1$counts, by.x="row.names", by.y="row.names", all.x= TRUE, all.y= FALSE, sort= FALSE)
head(results.raw)
write.table(results.raw, file="sRNAseq_logFC_padj_rawReads.txt", sep="\t", quote=TRUE, row.names=FALSE)
#instead of raw counts here I export TMM normalized cpm values with prior.counts 0.5 
results.tmm <- merge(results, cpm(y1), by.x="row.names", by.y="row.names", all.x= TRUE, all.y= FALSE, sort= FALSE)
head(results.tmm)
write.table(results.tmm, file="sRNAseq_logFC_padj_TMMcpm.txt", sep="\t", quote=TRUE, row.names=FALSE)

