install.packages("mixOmics") #install mixOmics package
require(mixOmics) #import mixOmics package
data('breast.TCGA') #import TCGAdata
data = list(mRNA = breast.TCGA$data.train$mrna, 
            miRNA = breast.TCGA$data.train$mirna, 
            proteomics = breast.TCGA$data.train$protein)
install.packages("rTensor") #install rTensor package
require(rTensor) #import rTensor package
#generate tensor
Z <- array(NA,c(150,200,184,142))
for (i in c(1:150))
{cat(i," ")
Z[i,,,] <-data.matrix(outer(outer(data$mRNA[i,],data$miRNA[i,],"*"),data$proteomics[i,],"*"))}
#tensor decomposition
HOSVD <- hosvd(as.tensor(Z))
U1 <- HOSVD$U[[1]] #x_{\ell_4,J}
U2 <- HOSVD$U[[2]] #x_{\ell_1,i_1}
U3 <- HOSVD$U[[3]] #x_{\ell_2,i_2}
U4 <- HOSVD$U[[4]] #x_{\ell_3,i_3}
#150 samples embeded (x_{\ell_4,j})
plot(U1[,c(1,4)],col=breast.TCGA$data.train$subtype)
legend(0.05,0.25,names(summary(breast.TCGA$data.train$subtype)),col=1:3,pch=1)
require(MASS) #import MASS package
LD <- lda(U1[,c(1,4)],breast.TCGA$data.train$subtype,CV=T,prior=rep(1/3,3)) #linear discrminant analysys
table(LD$class,breast.TCGA$data.train$subtype) #confusion matrix
ZZ <-order(-abs(HOSVD$Z@data[c(1,4),,,])) 
[1:20];data.frame(arrayInd(ZZ,dim(HOSVD$Z@data[c(1,4),,,])),HOSVD$Z@data[c(1,4),,,][ZZ]) #rank due to absolute G values
P2 <- pchisq(rowSums(scale(U2[,1:2])^2),2,lower.tail=F) #adding P-value to i_1 assuming multi gaussian to x_{\ell_1,i_1},\ell=1,2
P3 <- pchisq(rowSums(scale(U3[,1:2])^2),2,lower.tail=F) #adding P-value to i_2 assuming multi gaussian to x_{\ell_2,i_2},\ell=1,2
P4 <- pchisq(rowSums(scale(U4[,1:4])^2),4,lower.tail=F) #adding P-value to i_3 assuming multi gaussian to x_{\ell_3,i_3},\ell=1,2,3,4

U <- cbind(data$mRNA[,order(P2)[1:10]],cbind(data$miRNA[,order(P3)[1:10]],data$proteomics[,order(P4)[1:10]])) #merging top tens
install.packages("gplots") #install gplots package
require(gplots) #import gplot package
#draw heatmap
heatmap.2(scale(U),col=rgb(seq(0,1,by=0.1),seq(0,1,by=0.1),seq(1,0,by=-0.1)),RowSideColors=c(rep("black",45),rep("red",30),rep("green",75)),hclustfun=function(x){hclust(x,method="average")},trace="none",ColSideColors=c(rep("blue",10),rep("pink",10),rep("cyan",10)))
