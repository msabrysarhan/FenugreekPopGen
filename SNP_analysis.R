#----Fenugreek_SNP_Analysis----#
setwd("D:/Sarhan_Research_&_Work/Publications/Feenugreek_SNPs/_trial")

install.packages("digest")
install.packages("Rcpp")
install.packages("adegenet")
install.packages("sf")
install.packages("deldir")
BiocManager::install("SNPRelate")
install.packages("PopGenReport")


BiocManager::install("qvalue")
install.packages("dartR")

library("SNPRelate")
library("qvalue")
library("sf")
library("PopGenReport")
library("Rcpp")
library("adegenet")
library("digest")
library("deldir")
library("dartR")
library("ggplot2")


##import data with headers


trial = SNP_matrix
x = t(trial)
ind = colnames(SNP_matrix)
x = cbind(ind, x)
pop = character(length(trial))
for (i in 1:length(trial)){
  pop[i] = ind[i]
}
pop
loc.all = character(22937)
position = numeric(117)

x = cbind(as.data.frame(pop),x)



write.table(x, row.names = FALSE, sep = ",", file = "Fenugreek_SNPs_STRUCTURE.csv")
platy <- read.genetable("Fenugreek_SNPs_STRUCTURE.csv",
                        oneColPerAll = FALSE, ind = 1, pop = 1 ,sep = "")

platy.gl <- gi2gl(platy)

#write.table(platy.gl, row.names = TRUE, sep = "\t", file = "Fenugreek_SNPs_STRUCTURE.txt")
soda_SNP = as.matrix(soda_SNP)
dim(soda_SNP)
write.table(soda_SNP, sep = "\t", file = "transformed_soda_snps_filtered.txt")


start = Sys.time()
#write.table(gl.report.ld(soda_SNP, ncores = 8, save = FALSE, nchunks = 1), sep = "\t", file = "Fenugreek_SNPs_LD.txt")
LD = gl.report.ld(soda_SNP)
end = Sys.time()

end - start
gl.repo
#??gl.report.ld

write.table(LD, sep = "\t", file = "Fenugreek_SNPs_LD.txt")



gl2gds(platy.gl, outfile="test.gds")

#______________________________________________________#

x <- gl.report.pa(platy.gl, id= "X80", nmin=0, t=0)


length(platy.gl@ind.names)
platy.gl
soda_SNP
??dartr
dim(as.matrix(soda_SNP))
View(as.matrix(soda_SNP))
write.table(as.matrix(soda_SNP), sep = "\t", file = "str.txt")
n = gl.dist.pop(soda_SNP, method="manhattan")
plot(hclust(n))
pcoa_s = pcoa(n)
plot(pcoa_s$vectors)
gl.dist.heatmap(n)
gl.diversity(platy.gl)

gl.report.callrate(platy.gl)
gl2 <- gl.filter.repavg(platy.gl, t=0.5)
gl2
platy.gl@position <- as.integer(runif(nLoc(platy.gl),2,49))
platy.gl@loc.all <- testset.gl@loc.all[1:6]
??gl.dist.pop
n=platy.gl@gen
n
barplot(table(pop(platy.gl)), las=2)
b=gl.dist.pop(testset.gl, method="euclidean", diag=TRUE)
plot(hclust(b))
disttt=gl.dist.pop(platy.gl)
disttt
pc <- gl.pcoa(platy.gl, nfactors=5)
gl.pcoa.plot(pc, platy.gl, xaxis=1, yaxis=2)

library(ape)
install.packages("labdsv")
library("labdsv")
pco(gl.dist.pop(platy.gl))
biplot(gl.dist.pop(platy.gl))
pcoa(gl.dist.pop(platy.gl))
pc
testset.gl[1:180,]
gl.filter.repavg(platy.gl, t=0.5)
gl <- gl.ibd(platy.gl[1:52,])
barplot(table(pop(platy.gl)), las=2)
write.table(platy@tab, sep = "\t", "sample.txt")



x = read.csv("soda")
platy <- read.genetable(soda, ind=38, oneColPerAll=FALSE, sep="" )


read.csv(paste(.libPaths()[1],"/dartR/extdata/platy.csv"))

names
??write.csv
??read.genetable
soda

x = t(SNP_data)

read.genetable(SNP_data)

dim(soda)



platy <- read.genetable("D:/Sarhan_Research_&_Work/Publications/Feenugreek_SNPs/SNP_data.csv", oneColPerAll = 1)
platy <- read.genetable( paste(.libPaths()[1],"/dartR/extdata/platy.csv",sep="" ),
                         ind=1, oneColPerAll=FALSE,sep="/")
platy
read.csv(paste(.libPaths()[1],"/dartR/extdata/platy.csv"))

platy.gl <- gi2gl(platy)

platy.gl
??gi2gl

gl.report.callrate(platy.gl)
gl2 <- gl.filter.repavg(platy.gl, t=0.5)

platy.gl@position <- as.integer(runif(nLoc(platy.gl),2,49))
platy.gl@loc.all <- testset.gl@loc.all[1:6]


barplot(table(pop(gl)), las=2)
gl.dist.heatmap(platy.gl)
gl.dis
pc <- gl.pcoa(platy.gl, nfactors=5)
gl.pcoa.plot(pc, platy.gl, xaxis=1, yaxis=2)
library(ggplot2)
??gl.pcoa.plot
pc
testset.gl[1:180,]
gl.filter.repavg(platy.gl, t=0.5)
gl <- gl.ibd(platy.gl[1:52,])
barplot(table(pop(platy.gl)), las=2)
write.table(platy@tab, sep = "\t", "sample.txt")
#....................................................#
soda = as.matrix(soda)
ss = numeric(length(soda[,1]))
for (i in 1:length(soda[,1])){
  print(length(levels(factor(soda[i,]))))
  ss[i] = length(levels(factor(soda[i,])))
}
plot(ss)
soda[1,1]
levels(factor(soda[3,]))
summary(ss)
summary.factor(ss)

gl2gds(testset.gl)
gl2gds(platy.gl)


class(testset.gl)

class(platy.gl)

pruned_snp_data = gl2gds(platy.gl)
??gl2gds
snpgdsLDpruning()



x  = hapmap_geno
hapmap_geno = platy.gl

as.matrix(x)

library(gdsfmt)

data()

snpgdsCreateGeno("test.gds", 
                 genmat = as.matrix(x),
                 sample.id = x$ind.names,
                 snp.id = x$loc.names,
                 snpfirstdim=F)
genofile <- snpgdsOpen("test.gds")

snpset <- snpgdsLDpruning(genofile, ld.threshold=0.4)

snpset
