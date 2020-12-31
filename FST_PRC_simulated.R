
####(0) packages and prelims####
library(optparse);library(dplyr);library(doParallel);library(stringr);library(reshape2);library(msr);library(readr)
options(scipen=999) #turn off scientific notation

option_list = list(
  make_option(c("-a", "--stem"), type="character", default=NULL,
              help="stem for multipopulation ms file (without .txt)", metavar="character"),
  make_option(c("-n", "--ncores"), type="integer", default=NULL,
              help="number of cores to use", metavar="integer"),
  make_option(c("-r", "--rep"), type="integer", default=NULL,
              help="which simulation replicate", metavar="integer"),
  make_option(c("-b", "--snps"), type="integer", default=NULL,
              help="number of snps to use", metavar="integer"),
  make_option(c("-p", "--path"), type="character", default=NULL,
              help="output stem", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

stem <- opt$stem
ncores <- opt$ncores
path <- opt$path
snps <- opt$snps
rep <- opt$rep

try(stopCluster(cls))
cls <- makeCluster(ncores,outfile="")
registerDoParallel(cls)

`%nin%` = Negate(`%in%`)
tf <- t(combn(1:45,2))

####(1) load and parse SNP data####
file <- paste(path,stem,'.txt',sep='')

afsub <- function(x){
  xm <- do.call(rbind,x$gametes)
  afs <- apply(xm,2,sum)/length(xm[,1])
  sites <-  which(afs > 0.35 & afs <  0.65)
  sub <- x$gametes[[1]][,sites]
  loci <- x$positions[[1]][sites]*100000
  return(rbind(loci,sub))
}

bbins <- split(1:1000,ceiling(seq_along(1:1000)/100))

cb <- t(combn(1:10,2))

#sk=3+((rep-1)*1004)
#repl <- c(read_lines(file,n_max=3),read_lines(file,skip=sk,n_max=1003)) %>% parse_ms()
#out <- afsub(repl)
#ms_r <- out[-1,]

####(2) functions ####
##(A) FST and its permuter

fst.wc <- function(x,s1,s2) {
  if(max(x[which(!is.na(x))]) == min(x[which(!is.na(x))])) {return(NA)}
  f <- apply(x,2,function(x){sum(x[which(!is.na(x))])/length(x[which(!is.na(x))])})
  f1 <- apply(x[1:s1,],2,function(x){sum(x[which(!is.na(x))])/length(x[which(!is.na(x))])})
  f2 <- apply(x[(s1+1):(s1+s2),],2,function(x){sum(x[which(!is.na(x))])/length(x[which(!is.na(x))])})
  msp <- sum(s1*(f1-f)^2 + s2*(f2-f)^2)
  msg <- sum((1/((s1-1)+(s2-1)))*(s1*f1*(1-f1)+s2*f2*(1-f2)))
  n <- (s1+s2) - ((s1^2+s2^2)/(s1+s2))
  fst <- (msp-msg)/(msp + (n-1)*msg)
  return(fst)
}

perm.f <- function(j,b,d){
  fs <- vector()
  f <- NA
  t <- 0
  for (k in 4:100){
  while (is.na(f) & t <= 100){
    s1 <- sample(bbins[[cb[b,1]]],j,replace=T)
    s2 <- sample(bbins[[cb[b,2]]],j,replace=T)
    seq1 <- seq_s[s1,]
    seq2 <- seq_s[s2,]
    s1 <- dim(seq1)[1]
    s2 <- dim(seq2)[1]
    aln<-apply(rbind(seq1,seq2),2,as.numeric)
    f <- fst.wc(aln,s1,s2)
    t <- t+1}
  return(f)
}
}

##(B) for power calculations

funcp1 <- function(d,k,x,tf) {
  a <- d[x,tf[k,1],];a <- a[which(!is.na(a))]
  b <- d[x,tf[k,2],];b <- b[which(!is.na(b))]
  cc <- length(which(outer(a,b,"-") < 0))
  dd <- length(which(outer(a,b,"-") > 0))
  return(max(c(cc,dd))/(length(a)*length(b)))
}
#where d is the dataset (f.val, ibs1.val, etc), x is number of samples, k is index for tf)

funcp<- function(d,k,tf){
  return(sapply(1:21,funcp1,d=d,k=k,tf=tf))
}

##(C) master function that samples nsnps, gets power for each population pair
master <- function(n,r) {
  sk=3+((r-1)*1004)
  repl <- c(read_lines(file,n_max=3),read_lines(file,skip=sk,n_max=1003)) %>% parse_ms()

  out <- afsub(repl)
  ms_r <- out[-1,]
  wd <- dim(ms_r)[2]
  print(wd)
  
  if (wd < n){quit(save='no')}

  sub <- sort(sample(1:dim(ms_r)[2],n,replace=F))
  seq_s <- ms_r[,sub]
  out <- array(NA,c(21,45,100))
  clusterExport(cl=cls,c('seq_s','cb','bbins','fst.wc','perm.f','funcp','funcp1','tf'),envir=environment())
 for (e in 1:45){
    for (jj in 1:20){
      #print(c(e,seq(5,100,by=5)[jj]))
      res <- foreach(1:100,.errorhandling='remove', .combine=cbind) %dopar% perm.f(b=e,j=seq(5,100,by=5)[jj])
      out[jj,e,] <- res
    }
  }
 curves <- parSapply(cl=cls,1:990,funcp,d=out,tf=tf)
 curves <- cbind(seq(5,100,by=5),curves[-20,])
 colnames(curves) <- c('n_ind',1:990)
 write.table(curves,paste(snps,'SNPFST_',stem,'_',r,'.tab',sep=''),row.names=F,quote=F)
}


master(snps,rep)
