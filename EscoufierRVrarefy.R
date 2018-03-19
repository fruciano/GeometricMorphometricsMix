### If you use this function, please cite:
### Fruciano et al. 2013 - Plos One (for using the rarefaction procedure)
### Escoufier 1973 - Biometrics (because the procedure is based on Escoufier RV)

### This function computes rarefied version of Escoufier RV coefficient
### as suggested by Fruciano et al 2013 - Plos One
### This can be useful to compare RV among groups with the same variables
### but different sample sizes (as RV depends on sample size, see Fruciano et al 2013 - Plos One)

### Usage: RVrarefied(Block1, Block2, rep, samplesize)
### Block1, Block2: data frames with the variables in the first block and second block, respectively
### rep: number of resamplings to obtain the rarefied estimate
### samplesize: sample size to which the rarefaction procedure is carried out

### The function outputs a list with:
### mean rarefied RV (Rarefied_RV)
### 2.5%, 50% (median) and 97.5% percentiles of the rarefied RV
### All RV values obtained using the rarefaction procedure

### Notice: the function does NOT perform GPA on each rarefied sample
### this may or may not make a difference in estimates.
### In many cases, it will probably not make much difference
### (e.g., Fig. 2 in Fruciano et al 2013 - Plos One)

### Example: RVrarefied(A,B,1000,40)
### Obtain 1000 rarefied estimates of RV for blocks A and B
### (same individuals, different blocks of variables)
### rarefying at sample size 40


RVrarefied=function(Block1,Block2,rep,samplesize) {
BothBlocks=cbind(Block1,Block2)
endB1=ncol(Block1)
startB2=endB1+1
sizeboth=ncol(BothBlocks)
RV=vector(length=rep)
for (i in 1:rep) {

NewSamp=cbind(Block1,Block2)[sample(1:nrow(Block1),samplesize,replace=TRUE),]
RV[i]=EscoufierRV(NewSamp[,1:endB1],NewSamp[,startB2:sizeboth])
}
Results=list(
Rarefied_RV=mean(RV),
Quantiles=quantile(RV,c(0.025,0.5,0.975)),
AllRarefiedSamples=RV
)
return(Results)
}

EscoufierRV=function(Block1,Block2) {
BothBlocks=cbind(Block1,Block2)
nvarB1=ncol(Block1)
COV=cov(BothBlocks)
RV=sum(diag(COV[1:nvarB1,(nvarB1+1):ncol(COV)] %*%COV[(nvarB1+1):nrow(COV),1:nvarB1]))/(sqrt(sum(diag(COV[1:nvarB1,1:nvarB1]%*%COV[1:nvarB1,1:nvarB1]))*sum(diag(COV[(nvarB1+1):ncol(COV),(nvarB1+1):ncol(COV)]%*%COV[(nvarB1+1):ncol(COV),(nvarB1+1):ncol(COV)]))))
return(RV)
}
