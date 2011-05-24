#!/usr/bin/Rscript

# because I can't be bothered to type all the quotes
qw <- function(...) {
  as.character(sys.call()[-1])
}

options(stringsAsFactors = FALSE);
options(scipen=10)

args <- commandArgs(trailingOnly=TRUE)
filename = args[1]
outfile = args[2]
data <- read.csv(filename)


# set values <10 to 10
# This is a kludge and could possibly have horrible repercussions later on when
# we're determining differential expression as there will be some cases in which there
# is 0 variation. Not really sure what else to do. Limma fudges the variance anyway, so
# hopefully it should be reasonably robust to this sort of thing, but even so...
# Should be fine for just getting an idea of the things that are *really* changing, but
# be careful using this dataset for anything more complicated without maybe getting hold of
# the raw data.
fix10<-function(x){
  x[x<10]<-10
  return(x)
}

library(beadarray)
	
BSData<-readBeadSummaryData(filename,
		   sep=",",
                   skip=1,
                   columns=list(
		       exprs="AVG_Signal",
		       se.exprs="BEAD_STDERR",
		       NoBeads="Avg_NBEADS",
		       Detection="Detection"
		   )
)	

exprs(BSData) <- fix10(exprs(BSData))
E = normaliseIllumina(BSData, method="quantile", transform="log2")
data <- exprs(E)

library(limma)

#grab expression value cols
cols <- rep(qw(ns, astro),4)
ns<-which(cols=="ns")
astro<-which(cols=="astro")

design<-matrix(0,nrow=(ncol(data)), ncol=2)
colnames(design)<-c("ns","astro")
design[ns,"ns"]<-1
design[astro,"astro"]<-1


fit<-lmFit(data, design)
cont.matrix<-makeContrasts(nsvsastro=astro-ns, levels=design)
fit<-contrasts.fit(fit, cont.matrix)
ebFit<-eBayes(fit)

write.fit(ebFit, file=outfile , adjust="BH")
data<-read.table(outfile, sep="\t", header=T)

data<- topTable(ebFit, number=nrow(data))
write.csv(data,outfile)




