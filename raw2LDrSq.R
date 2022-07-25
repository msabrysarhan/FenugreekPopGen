#Genome-Wide SNP Analysis


setwd("D:/Sarhan_Research_&_Work/Publications/Feenugreek_SNPs/trial")
#tassel
#gappit

#=======================================================#
#remove heterogenous pairs and N bases


rm_hetero <- function(matrix){
  errors = c()
  for (i in 1:length(matrix[1,])){
    if (matrix[1,i] == "NN"|matrix[2,i]=="NN"|
        substr(matrix[1,i], 1, 1) != substr(matrix[1,i], 2, 2)|
        substr(matrix[2,i], 1, 1) != substr(matrix[2,i], 2, 2)){
      errors = append(errors, i)
    }
  }
  for (n in errors){
    matrix = matrix[,-n]
  }
  matrix
}



##2X2 table
haplotyping <- function(matrix){
  haplotypes = character(length(matrix[,1]))
  for (i in 1:length(matrix[1,])){
    haplotypes[i] = paste(as.character(matrix[1,i]), 
                          as.character(matrix[2,i]))
  }
  unit = as.data.frame(summary(as.factor(haplotypes)))
  unit
}



LD_D <- function(df){
  total = sum(df[,1])
  D = ((df[1,1]/total)*(df[4,1]/total))-((df[2,1]/total)*(df[3,1]/total))
  D
}



LD_r_sq <- function(df){
  total = sum(df[,1])
  D = ((df[1,1]/total)*(df[4,1]/total))-((df[2,1]/total)*(df[3,1]/total))
  r_sq = (D^2)/(((df[1,1]/total)+(df[3,1]/total))*
                  ((df[2,1]/total)+(df[4,1]/total))*
                  ((df[1,1]/total)+(df[2,1]/total))*
                  ((df[3,1]/total)+(df[4,1]/total)))
  r_sq
}


##matrix construction
LD_r_sq_matrix <- function(SNP_matrix){
  matrix_length = length(SNP_matrix[,1])
  print(matrix_length)
  result = matrix(0, matrix_length, matrix_length)
  for (i in 1:matrix_length){
    for (n in 1:i){
      pair = rm_hetero(SNP_matrix[c(i,n),])
      haplo = haplotyping(pair)
      result[i,n] = LD_r_sq(haplo)
      #print(LD_r_sq(haplo))
    }
    print(i)
  }
  result
}


rrr = LD_r_sq_matrix(soda[c(10000:10100),])

library(pheatmap)
pheatmap(rrr, cluster_rows = FALSE, cluster_cols = FALSE, 
         color = farben(256))
View(rrr)
library(ggplot2)

pheatmap(matrix(0,38000,38000), cluster_rows = FALSE, cluster_cols = FALSE, 
         color = farben(256))


