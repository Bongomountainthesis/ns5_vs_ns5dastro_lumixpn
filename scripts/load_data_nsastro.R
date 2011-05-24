#!/usr/local/bin/Rscript


library(beadarray)

#We only have bead-summary level data:

dataFile<-"data/raw_data.csv"

BSData<-readBeadSummaryData(dataFile,
		   sep=",",
                   skip=1,
                   columns=list(
		       exprs="AVG_Signal",
		       se.exprs="BEAD_STDERR",
		       NoBeads="Avg_NBEADS",
		       Detection="Detection"
		   )
)	



fix.col<-function(x){
	x[x<=0]<-1 
	return(x)
}

exprs(BSData) <- apply(exprs(BSData),2,fix.col) 

pd <- data.frame(cell.line=rep(c("NS","Astro"),4))
rownames(pd) = c("NS5.1","Astro.1","NS5.2", "Astro.2", "NS5.3", "Astro.3","NS5.4", "Astro.4")
metadata<-data.frame(labelDescription=c("Cell Line (Astrocyte or NS5"))
pd <- new("AnnotatedDataFrame", data = pd, varMetadata=metadata)

phenoData(BSData)<-pd

save(BSData, file="results/BSData.RData")

postscript(file="results/Boxplotprenorm.ps", horizontal=FALSE)
boxplot(as.data.frame(log2(exprs(BSData))),  las = 2, outline = FALSE, ylab = "log2(intensity)", col=c(2,4,2,4,2,4,2,4), main="Log2 Expression")
dev.off()

postscript(file="results/plotMAXYprenorm.ps", horizontal=FALSE)
plotMAXY(exprs(BSData), arrays = 1:8,  pch = 16)
dev.off()

postscript(file="results/plotDENSITYprenorm.ps", horizontal=FALSE)
E <- log2(exprs(BSData))
plot(density(E[,1]))
for(i in 1:8){
  lines(density(E[,i]),col=i)
}
dev.off()



BSData.quantile = normaliseIllumina(BSData, method = "quantile", transform = "log2")


postscript(file="results/Boxplotpostnorm.ps", horizontal=FALSE)
boxplot(as.data.frame((exprs(BSData.quantile))),  las = 2, outline = FALSE, ylab = "log2(intensity)", col=c(2,4,2,4,2,4,2,4), main="Log2 Expression")
dev.off()

postscript(file="results/plotMAXYpostnorm.ps", horizontal=FALSE)
plotMAXY(exprs(BSData.quantile), arrays = 1:8, log = FALSE, pch = 16)
dev.off()

save(BSData.quantile, file="results/BSData.quantile.RData")

postscript(file="results/plotDENSITYpostnorm.ps", horizontal=FALSE)
E <- exprs(BSData.quantile)
plot(density(E[,1]))
for(i in 1:8){
  lines(density(E[,i]),col=i)
}
dev.off()


