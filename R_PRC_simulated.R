####samples of p SNPs and n individuals from coalescent-simulated data
####runs boostrap power calculations for IBD (widehatR) 
####with balanced sampling over locations (number of individuals = n for all locations) 
####and power calculated for 100 randomly sampled iterations of p SNPs and n samples

####load required packages and set settings #####
library(optparse);library(doParallel);library(stringr);library(reshape2);library(msr);library(readr);library(dplyr)
options(scipen=999) #turn off scientific notation
`%nin%` = Negate(`%in%`)
th <- c(0.5,0.8) #specify threshold values (tau) for widehatR

####get command line input####
option_list = list(
  make_option(c("-a", "--stem"), type="character", default=NULL, 
              help="stem for multipopulation ms file (without .txt)", metavar="character"),
  make_option(c("-n", "--ncores"), type="integer", default=NULL, 
              help="number of cores to use", metavar="integer"),
  make_option(c("-r", "--rep"), type="integer", default=NULL,
              help="which simulation rep to get", metavar="integer"),
  make_option(c("-b", "--snps"), type="integer", default=NULL,
              help="number of snps to use", metavar="integer"),
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path to ms file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

stem <- opt$stem
ncores <- opt$ncores
path <- opt$path
rep <- opt$rep
snps <- opt$snps

####set up cluster####
try(stopCluster(cls))
cls <- makeCluster(ncores,outfile="")
registerDoParallel(cls)

####lists and bins for selecting location pairs and simulation replicates####
cb <- t(combn(1:10,2)) #list of all pairwise comparisons between locations 1:10
tf <- t(combn(1:45,2)) #list of all pairwise comparisons between pairs of locations
bbins <- split(1:1000,ceiling(seq_along(1:1000)/100)) #bins for each ms replicate simulation
dex <- expand.grid(1:45,1:20) #indices for rows/columns of result matrix

####load and parse SNP data####
file <- paste(path,stem,'.txt',sep='')

####select SNPs with MAF > 0.35, get loci number for hmmIBD####
afsub <- function(x){
  xm <- do.call(rbind,x$gametes)
  afs <- apply(xm,2,sum)/length(xm[,1])
  sites <-  which(afs > 0.35 & afs <  0.65)
  print(length(sites))
  sub <- x$gametes[[1]][,sites]
  loci <- x$positions[[1]][sites]*100000
  return(rbind(loci,sub))
}

####get specified ms simulation replicate###
print('loading data from ms replicate')
print(rep)
sk<-3+((rep-1)*1004)
repl <- c(read_lines(file,n_max=3),read_lines(file,skip=sk,n_max=1003)) %>% parse_ms()
out <- afsub(repl)
ms_r <- out[-1,]
ms_rL <- out[1,]
wd <- dim(out)[2]
print(wd)

if (wd < 200){quit(save='no')}


#####function ip: takes bootstrap sample of n individuals and calculates widehatR for tau=0.5,0.8
#####for one of 45 possible edges (specified by the two nodes on the edge, a and b)
#####using two balanced samples (with replacement) from the two populations

#ssub is the sequence, lsub is loci names, n is number of samples
#b is index of cb, specifing the two nodes/edge
#d is dummy variable for apply
ip <- function(ssub,lsub,n,b,d){
  s1 <- sample(bbins[[cb[b,1]]],n,replace=T)
  s2 <- sample(bbins[[cb[b,2]]],n,replace=T)
  seqn1 <- ssub[s1,] 
  seqn2 <- ssub[s2,]
  
  afs1<-apply(ssub[bbins[[cb[b,1]]],],2,sum)/length(ssub[bbins[[cb[b,1]]],][,1]) #get allele frequencies for hmmIBD
  afs1 <- cbind(lsub,afs1,1-afs1)
  afs2<-apply(ssub[bbins[[cb[b,2]]],],2,sum)/length(ssub[bbins[[cb[b,2]]],][,1])
  afs2 <- cbind(lsub,afs2,1-afs2)
  
  n1 <- c('chrom','pos',s1)
  se1<-rbind(n1,cbind(lsub,t(seqn1)))
  
  n2 <- c('chrom','pos',s2)
  se2<-rbind(n2,cbind(lsub,t(seqn2)))
  
  r<-sample(1:100000000,1) #tag all files written out with random number to avoid overwriting while parallel processing
  write.table(afs1,paste(r,'_afs1.tab',sep=''),col.names = FALSE, row.names=FALSE,sep='\t')
  write.table(se1,paste('s1_',r,'.tab',sep=''),col.names=F,row.names=F,quote=F,sep='\t')
  write.table(afs2,paste(r,'_afs2.tab',sep=''),col.names = FALSE, row.names=FALSE,sep='\t')
  write.table(se2,paste('s2_',r,'.tab',sep=''),col.names=F,row.names=F,quote=F,sep='\t')
  
  cmd1 <- paste('./hmmIBD  -i s1_',r,'.tab -I s2_',r,'.tab -o t_',r,' -f ',paste(r,'_afs1.tab',sep=''),' -F ',paste(r,'_afs2.tab',sep=''),sep='')
  cmd3 <- paste('rm s1_',r,'.tab ','s2_',r,'.tab t_',r,'.hmm.txt t_',r,'.hmm_fract.txt ',r,'_afs2.tab ',r,'_afs1.tab',sep='')
  system(cmd1,ignore.stdout = TRUE)
  op<- data.frame(read.table(paste('t_',r,'.hmm_fract.txt',sep=''),header=T))$fract_sites_IBD
  system(cmd3)
 
 
  return(sapply(th,function(a)length(which(op >a))/length(op)))
}

####functions for power calculations###
#where d is the dataset, x is number of samples, k is index for tf
funcp1 <- function(d,k,x,tf) {
  a <- d[x,tf[k,1],];a <- a[which(!is.na(a))]
  b <- d[x,tf[k,2],];b <- b[which(!is.na(b))]
  cc <- length(which(outer(a,b,"-") < 0))
  dd <- length(which(outer(a,b,"-") > 0))
  return(max(c(cc,dd))/(length(a)*length(b)))
}

funcp<- function(d,k,tf){
  return(sapply(1:21,funcp1,d=d,k=k,tf=tf))  
}

####sapply wrapper for ip####
ipp <- function(i,pp,ll){
      res <- sapply(1:100,ip,ssub=pp,lsub=ll,b=dex[i,1],n=seq(5,100,by=5)[dex[i,2]])
      return(res)
}

####unwrapper for after parLapply####
unr <- function(i,x){
i50[dex[i,2],dex[i,1],1:length(x[[i]][1,])] <<- x[[i]][1,]
i80[dex[i,2],dex[i,1],1:length(x[[i]][2,])] <<- x[[i]][2,] 
}

####master function to sample snps, sample individuals, calculate widehatR, and calculate power
i50 <- i80 <- array(NA,c(21,45,100))

master <- function(n) {
  sub <- sort(sample(1:dim(ms_r)[2],n,replace=F))
  l<-cbind(rep(1,length(ms_rL)),ms_rL)
  lsub <- l[sub,];ssub <- ms_r[,sub]

  clusterExport(cl=cls,c('ip','ipp','dex','bbins','cb','th','funcp','funcp1','tf','lsub','ssub'),envir=environment())
  
  print("Calculating IBD values")
  ch<-parLapply(cl=cls,1:900,ipp,pp=ssub,ll=lsub)
  invisible(sapply(1:900,unr,x=ch))
  
  print("Calculating power")
  curves.ib50 <- sapply(1:990,funcp,d=i50,tf=tf)
  curves.ib80 <- sapply(1:990,funcp,d=i80,tf=tf)

  out1<-list()
  out1[[1]] <- curves.ib50
  out1[[2]] <- curves.ib80

  return(out1)
} 

res <- master(snps)
r<-sample(1:100000,1)

res50 <- res[[1]]
res50 <- res50[-21,]
res50 <- cbind(seq(5,100,by=5),res50)
colnames(res50) <- c('snps',1:990)

res80 <- res[[2]]
res80 <- res80[-21,]
res80 <- cbind(seq(5,100,by=5),res80)
colnames(res80) <- c('snps',1:990)

write.table(res50,paste(stem,'_',snps,'_ib50_',rep,'_',r,'.tab',sep=''),row.names=FALSE,quote=FALSE)
write.table(res80,paste(stem,'_',snps,'_ib80_',rep,'_',r,'.tab',sep=''),row.names=FALSE,quote=FALSE)

stopCluster(cls)


