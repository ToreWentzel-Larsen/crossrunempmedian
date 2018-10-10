library(Rmpfr)
library(crossrun)

# finction for generating all subsets of size m1 within n1 objects:
gensubset <- function(n1, m1) {
  chooseall <- choose(n1,m1)
  if (m1==n1) res <- matrix(1, ncol=1, nrow=m1) else
    if (m1==1) res <- diag(rep(1,n1)) else {
      res <- matrix(0,ncol=n1, nrow=choose(n1,m1))
      chooses <- choose(c((n1-1):(m1-1)),m1-1)
      ends <- cumsum(chooses)
      nch <- length(chooses)
      starts <- c(1,1+ends[-nch])
      for (part in 1:nch) {
        startp <- starts[part]
        endp <- ends[part]
        res[startp:endp,part] <- 1
        res[startp:endp,(part+1):n1] <- gensubset(n1-part,m1-1)
      } # end parts
    } # end if not simple case
    return(res)
} # end function gensubset

gensubset(3,3)
gensubset(3,1)
gensubset(4,3)
g9.4 <- gensubset(9,4)
g9.4[1:100,]
g9.4[101:126,]
gensubset(7,3)

# C, number of crossings and L, longest run,
# empirical median, functions:
cf <- function(subs) length(rle(subs)$lengths)-1
lf <- function(subs) max(rle(subs)$lengths)
clem <- function(m1) {
  n1 <- 2*m1
  # only included subsets including the first element:
  chooseall <- choose(n1-1,m1-1)
  subsm <- matrix(0, ncol=n1, nrow=chooseall)
  subsm[,1] <- 1
  subsm[,-1] <- gensubset(n1-1,m1-1)
  c1 <- apply(subsm,1,cf)
  l1 <- apply(subsm,1,lf)
  return(table(c1,l1))
} # end function clem

# computations as long as possible due
# to storage requirements:
cl2 <- clem(m1=2)
cl3 <- clem(m1=3)
cl4 <- clem(m1=4)
cl5 <- clem(m1=5)
cl6 <- clem(m1=6)
cl7 <- clem(m1=7)
cl8 <- clem(m1=8)
cl9 <- clem(m1=9)
cl10 <- clem(m1=10)
cl11 <- clem(m1=11)
Sys.time()
cl12 <- clem(m1=12)
Sys.time() # a little more than 2 minutes
Sys.time()
cl13 <- clem(m1=13)
Sys.time() # a little more than 8 minutes
Sys.time()
cl14 <- clem(m1=14)
Sys.time()
# higher m requires too much storage

# specific checks in small cases (for comparison
# with resutls of brute force calculations):
check4 <- cbind(rep(1,35),gensubset(7,3))
data.frame(c=apply(check4,1,cf),
           l=apply(check4,1,lf))
check5 <- cbind(rep(1,126),gensubset(9,4))
check5 <- data.frame(check5,c=apply(check5,1,cf),
                     l=apply(check5,1,lf))
check5[1:80,]
check5[81:126,]
data.frame(c=apply(check5,1,cf),
           l=apply(check5,1,lf))
table(c=apply(check5[,1:10],1,cf),
      l=apply(check5[,1:10],1,lf))

# look at the results:
cl2
cl3
cl4
cl5
cl6
cl7
cl8
cl9
cl10
cl11
cl12
cl13
cl14


# check with simulations:
simclem <- function(m1=14, nsim = 1e+05) {
  n1 <- 2*m1
  cs <- rep(NA, nsim)
  ls <- rep(NA, nsim)
  for (sim in 1:nsim) {
    series <- rnorm(n1)
    med <- median(series)
    above <- as.numeric(series>med)
    rleabove <- rle(above)$lengths
    cs[sim] <- length(rleabove) - 1
    ls[sim] <- max(rleabove)
  } # end for sim
  return(data.frame(cs=cs,ls=ls))
} #end function simclem

#simulate for m=14 (n=2*14=28):
simcl14 <- simclem(m1=14)
# check means of C and L:
sum(apply(cl14,1,sum)*c(1:27))/sum(cl14)
mean(simcl14$cs)
sum(apply(cl14,2,sum)*c(1:14))/sum(cl14)
mean(simcl14$ls)
# check standard deviations of C and L:
sqrt(sum(apply(cl14,1,sum)*c(1:27)^2)/sum(cl14) -
       (sum(apply(cl14,1,sum)*c(1:27))/sum(cl14))^2)
sd(simcl14$cs)
sqrt(sum(apply(cl14,2,sum)*c(1:14)^2)/sum(cl14) -
       (sum(apply(cl14,2,sum)*c(1:14))/sum(cl14))^2)
sd(simcl14$ls)
# check means of C*L:
sum(diag(1:27) %*% cl14 %*% diag(1:14))/sum(cl14)
mean(simcl14$cs*simcl14$ls)
# check theoretical and empirical cdf:
plot(x=as.numeric(names(table(simcl14$cs))),
     y=(cumsum(cumsumm(cl14)[,14])/(sum(cl14)))[
       as.numeric(names(table(simcl14$cs)))],
     type="l", xlab="Number of crossings", ylab="CDF", las=1)
points(x=as.numeric(names(table(simcl14$cs))),
       y=cumsum(table(simcl14$cs))/sum(table(simcl14$cs)),
       type="l", col="red",lty="dotted")
lines(x=c(4,7), y=c(.9,.9), col="red")
text(x=7, y=.9, pos=4, labels="red: simulations", col="red")
plot(x=as.numeric(names(table(simcl14$ls))),
     y=(cumsum(cumsummcol(cl14)[27,])/(sum(cl14)))[
       as.numeric(names(table(simcl14$ls)))],
     type="l", xlab="Longest run", ylab="CDF", las=1)
points(x=as.numeric(names(table(simcl14$ls))),
       y=cumsum(table(simcl14$ls))/sum(table(simcl14$ls)),
       type="l", col="red",lty="dotted")
lines(x=c(4,6), y=c(.4,.4), col="red")
text(x=6, y=.4, pos=4, labels="red: simulations", col="red")




# function for computing the number of subsets of each size m
# from n subsequent observations, for each value of the number of
# crossings and the longest run. subsets are classified by
# whether they contain the first observation and the computations
# are partition on the length of the initial run.
crossrunem <- function(nmax = 100, prec = 120,
                       printn = FALSE) {
  # conditioning of S= first value included in the subset (S=1)
  # or not (S=0).
  # nfi: number of subsets, first value included,
  # nfn: number of subsets, first value not included.
  # separate code by brute force for low n (n=1,2,3,4):
  nfi <- list(Rmpfr::mpfrArray(0, prec, dim = c(1, 1, 2)))
  nfn <- list(Rmpfr::mpfrArray(0, prec, dim = c(1, 1, 2)))
  dimnames(nfi[[1]]) <- list("c=0","l=1",c("m=0","m=1"))
  dimnames(nfn[[1]]) <- list("c=0","l=1",c("m=0","m=1"))
  nfn[[1]][1,1,1] <- 1 # m=0
  nfi[[1]][1,1,2] <- 1 # m=1
  nfi[[2]] <- Rmpfr::mpfrArray(0, prec, dim = c(2, 2, 3))
  nfn[[2]] <- Rmpfr::mpfrArray(0, prec, dim = c(2, 2, 3))
  dimnames(nfi[[2]]) <- list(paste0("c=",0:1),paste0("l=",1:2),
                             paste0("m=",0:2))
  dimnames(nfn[[2]]) <- list(paste0("c=",0:1),paste0("l=",1:2),
                             paste0("m=",0:2))
  nfn[[2]][1,2,1] <- 1 # m=0
  nfi[[2]][2,1,2] <- 1 # m=1
  nfn[[2]][2,1,2] <- 1
  nfi[[2]][1,2,3] <- 1 # m=2
  nfi[[3]] <- Rmpfr::mpfrArray(0, prec, dim = c(3, 3, 4))
  nfn[[3]] <- Rmpfr::mpfrArray(0, prec, dim = c(3, 3, 4))
  dimnames(nfi[[3]]) <- list(paste0("c=",0:2),paste0("l=",1:3),
                             paste0("m=",0:3))
  dimnames(nfn[[3]]) <- list(paste0("c=",0:2),paste0("l=",1:3),
                             paste0("m=",0:3))
  nfn[[3]][1,3,1] <- 1 # m=0
  nfi[[3]][2,2,2] <- 1 # m=1
  nfn[[3]][2,2,2] <- 1
  nfn[[3]][3,1,2] <- 1
  nfi[[3]][2,2,3] <- 1 # m=2
  nfi[[3]][3,1,3] <- 1
  nfn[[3]][2,2,3] <- 1
  nfi[[3]][1,3,4] <- 1 # m=3
  nfi[[4]] <- Rmpfr::mpfrArray(0, prec, dim = c(4, 4, 5))
  nfn[[4]] <- Rmpfr::mpfrArray(0, prec, dim = c(4, 4, 5))
  dimnames(nfi[[4]]) <- list(paste0("c=",0:3),paste0("l=",1:4),
                             paste0("m=",0:4))
  dimnames(nfn[[4]]) <- list(paste0("c=",0:3),paste0("l=",1:4),
                             paste0("m=",0:4))
  nfn[[4]][1,4,1] <- 1 # m=0
  nfi[[4]][2,3,2] <- 1 # m=1
  nfn[[4]][2,3,2] <- 1
  nfn[[4]][3,2,2] <- 2
  nfi[[4]][2:3,2,3] <- 1 # m=2
  nfi[[4]][4,1,3] <- 1
  nfn[[4]][2:3,2,3] <- 1
  nfn[[4]][4,1,3] <- 1
  nfi[[4]][2,3,4] <- 1 # m=3
  nfi[[4]][3,2,4] <- 2
  nfn[[4]][2,3,4] <- 1 
  nfi[[4]][1,4,5] <- 1 # m=4
  # iterative procedure for higher n (>=5):
  if (nmax>4) for (nn in 5:nmax) {
    nfi[[nn]] <- Rmpfr::mpfrArray(0, prec, dim = c(nn, nn, nn+1))
    nfn[[nn]] <- Rmpfr::mpfrArray(0, prec, dim = c(nn, nn, nn+1))
    dimnames(nfi[[nn]]) <- list(paste0("c=",0:(nn-1)),paste0("l=",1:nn),
                                paste0("m=",0:nn))
    dimnames(nfn[[nn]]) <- list(paste0("c=",0:(nn-1)),paste0("l=",1:nn),
                                paste0("m=",0:nn))
    # separate computation for m=0,n:
    nfn[[nn]][1,nn,1] <- 1 # m=0
    nfi[[nn]][1,nn,nn+1] <- 1 # m=n
    # iterative procedure nfi, for 1 <= m <= n-1:
    for (mm in 1:(nn-1)) for (gg in 1:mm) {
      if (gg>=nn-gg) {
        if (nn-gg==1) 
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <- 
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] + nfn[[1]][1,1,1]
        else
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <- 
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] + 
            crossrun::cumsumm(nfn[[nn-gg]][1:(nn-gg),,mm-gg+1])[,nn-gg]
      }
      if (gg<nn-gg) {
        if (gg==1) {
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <- 
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] + 
            nfn[[nn-gg]][1:(nn-gg),1:gg,mm-gg+1]
        }
        else {
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <- 
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] + 
            crossrun::cumsumm(nfn[[nn-gg]][1:(nn-gg),1:gg,mm-gg+1])[,gg]
        }
        nfi[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] <- 
          nfi[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] + 
          nfn[[nn-gg]][1:(nn-gg),(gg+1):(nn-gg),mm-gg+1]
      } # end low g
    } # end iterative procedure nfi
    # iterative procedure nfn, for 1 <= m <= n-1:
    for (mm in 1:(nn-1)) for (gg in 1:(nn-mm)) {
      if (gg>=nn-gg) {
        if (nn-gg==1) 
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <- 
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] + nfi[[1]][1,1,2]
        else {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <- 
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] + 
            crossrun::cumsumm(nfi[[nn-gg]][1:(nn-gg),,mm+1])[,nn-gg]
        }
      } # end high g
      if (gg<nn-gg) {
        if (gg==1) {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <- 
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] + 
            nfi[[nn-gg]][1:(nn-gg),1:gg,mm+1]
        }
        else {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] + 
            crossrun::cumsumm(nfi[[nn-gg]][1:(nn-gg),1:gg,mm+1])[,gg]
        }
        nfn[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] <- 
          nfn[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] + 
          nfi[[nn-gg]][1:(nn-gg),(gg+1):(nn-gg),mm+1]
      } # end low g
    } # end iterative procedure nfi
    if (printn==TRUE) print(nn)
    if (printn==TRUE) print(Sys.time())
  } # end for nn
  return(list(nfi = nfi, nfn = nfn))
} # end function crossrunem

Sys.time()
em10 <- crossrunem(nmax=10)
Sys.time() # a few seconds
asNumeric(em10$nfn[[5]])
# agreement with brute force calculation.

# reproducing the results from cl (up to n=28):
Sys.time()
em28 <- crossrunem(nmax=28, printn=TRUE)
Sys.time() # about 3 minutes
asNumeric(em28$nfi[[28]][,,15]) 
# check limited nonzero range for m=n/2:
asNumeric(em28$nfi[[28]][-1,1:14,15])
cl14
asNumeric(em28$nfi[[28]][-1,1:14,15]) - cl14 # consistent 
# check for all lower n:
asNumeric(em28$nfi[[4]][-1,1:2,3]) - cl2
asNumeric(em28$nfi[[6]][-1,1:3,4]) - cl3
asNumeric(em28$nfi[[8]][-1,1:4,5]) - cl4
asNumeric(em28$nfi[[10]][-1,1:5,6]) - cl5
asNumeric(em28$nfi[[12]][-1,1:6,7]) - cl6
asNumeric(em28$nfi[[14]][-1,1:7,8]) - cl7
asNumeric(em28$nfi[[16]][-1,1:8,9]) - cl8
asNumeric(em28$nfi[[18]][-1,1:9,10]) - cl9
asNumeric(em28$nfi[[20]][-1,1:10,11]) - cl10
asNumeric(em28$nfi[[22]][-1,1:11,12]) - cl11
asNumeric(em28$nfi[[24]][-1,1:12,13]) - cl12
asNumeric(em28$nfi[[26]][-1,1:13,14]) - cl13
# no differences
# til 60:
Sys.time()
em60 <- crossrunem(nmax=60, printn=TRUE)
Sys.time() # several hours, the last step (n=60)
# about 20 minutes

# konsistency with cl for low n:
sum(abs(asNumeric(em60$nfi[[4]][-1,1:2,3]) - cl2))
sum(abs(asNumeric(em60$nfi[[6]][-1,1:3,4]) - cl3))
sum(abs(asNumeric(em60$nfi[[8]][-1,1:4,5]) - cl4))
sum(abs(asNumeric(em60$nfi[[10]][-1,1:5,6]) - cl5))
sum(abs(asNumeric(em60$nfi[[12]][-1,1:6,7]) - cl6))
sum(abs(asNumeric(em60$nfi[[14]][-1,1:7,8]) - cl7))
sum(abs(asNumeric(em60$nfi[[16]][-1,1:8,9]) - cl8))
sum(abs(asNumeric(em60$nfi[[18]][-1,1:9,10]) - cl9))
sum(abs(asNumeric(em60$nfi[[20]][-1,1:10,11]) - cl10))
sum(abs(asNumeric(em60$nfi[[22]][-1,1:11,12]) - cl11))
sum(abs(asNumeric(em60$nfi[[24]][-1,1:12,13]) - cl12))
sum(abs(asNumeric(em60$nfi[[26]][-1,1:13,14]) - cl13))
sum(abs(asNumeric(em60$nfi[[28]][-1,1:14,15]) - cl14))
# check og equality between nfi andnfn for even n and m=n/2;
for (mm in 1:30) print(c(mm,2*mm,
        sum(abs(asNumeric(em60$nfi[[2*mm]][,,mm+1] - 
                            em60$nfn[[2*mm]][,,mm+1])))))
# sjekk av konsistens mellom em60$nfi[[2*mm]][,,mm+1] og
# em60$nfi[[2*mm]][-1,1:mm,mm+1]:
for (mm in 1:30) print(c(mm,2*mm,
        sum(asNumeric(em60$nfi[[2*mm]][,,mm+1])) -
          sum(asNumeric(em60$nfi[[2*mm]][-1,1:mm,mm+1])))) 

# check with simulations for n=60:
#simulate for m=14 (n=2*14=28):
simcl30 <- simclem(m1=30)
jointem60 <- asNumeric(em60$nfi[[60]][-1,1:30,31])
# check means of C and L:
sum(apply(jointem60,1,sum)*c(1:59))/sum(jointem60)
mean(simcl30$cs)
sum(apply(jointem60,2,sum)*c(1:30))/sum(jointem60)
mean(simcl30$ls)
# check standard deviations of C and L:
sqrt(sum(apply(jointem60,1,sum)*c(1:59)^2)/sum(jointem60) -
       (sum(apply(jointem60,1,sum)*c(1:59))/sum(jointem60))^2)
sd(simcl30$cs)
sqrt(sum(apply(jointem60,2,sum)*c(1:30)^2)/sum(jointem60) -
       (sum(apply(jointem60,2,sum)*c(1:30))/sum(jointem60))^2)
sd(simcl30$ls)
# check means of C*L:
sum(diag(1:59) %*% jointem60 %*% diag(1:30))/sum(jointem60)
mean(simcl30$cs*simcl30$ls)
# check theoretical and empirical cdf:
plot(x=as.numeric(names(table(simcl30$cs))),
     y=(cumsum(cumsumm(jointem60)[,30])/(sum(jointem60)))[
       as.numeric(names(table(simcl30$cs)))],
     type="l", xlab="Number of crossings", ylab="CDF", las=1)
points(x=as.numeric(names(table(simcl30$cs))),
       y=cumsum(table(simcl30$cs))/sum(table(simcl30$cs)),
       type="l", col="red",lty="dotted")
lines(x=c(12,16), y=c(.9,.9), col="red")
text(x=16, y=.9, pos=4, labels="red: simulations", col="red")
plot(x=as.numeric(names(table(simcl30$ls))),
     y=(cumsum(cumsummcol(jointem60)[59,])/(sum(jointem60)))[
       as.numeric(names(table(simcl30$ls)))],
     type="l", xlab="Longest run", ylab="CDF", las=1)
points(x=as.numeric(names(table(simcl30$ls))),
       y=cumsum(table(simcl30$ls))/sum(table(simcl30$ls)),
       type="l", col="red",lty="dotted")
lines(x=c(8,11), y=c(.2,.2), col="red")
text(x=11, y=.2, pos=4, labels="red: simulations", col="red")

# crossrunem is quite time consuming for n above 60.
# function for continuation of an already computed
# crossrunem result is therefore useful:
crossrunemcont <- function(emstart, n1=61, nmax = 100, 
                           prec = 120, printn = FALSE) {
  nfi <- emstart$nfi
  nfn <- emstart$nfn
  for (nn in n1:nmax) {
    nfi[[nn]] <- Rmpfr::mpfrArray(0, prec, dim = c(nn, nn, nn+1))
    nfn[[nn]] <- Rmpfr::mpfrArray(0, prec, dim = c(nn, nn, nn+1))
    dimnames(nfi[[nn]]) <- list(paste0("c=",0:(nn-1)),paste0("l=",1:nn),
                                paste0("m=",0:nn))
    dimnames(nfn[[nn]]) <- list(paste0("c=",0:(nn-1)),paste0("l=",1:nn),
                                paste0("m=",0:nn))
    # separate computation for m=0,n:
    nfn[[nn]][1,nn,1] <- 1 # m=0
    nfi[[nn]][1,nn,nn+1] <- 1 # m=n
    # iterative procedure nfi, for 1 <= m <= n-1:
    for (mm in 1:(nn-1)) for (gg in 1:mm) {
      if (gg>=nn-gg) {
        if (nn-gg==1) 
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <- 
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] + nfn[[1]][1,1,1]
        else
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <- 
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] + 
            crossrun::cumsumm(nfn[[nn-gg]][1:(nn-gg),,mm-gg+1])[,nn-gg]
      }
      if (gg<nn-gg) {
        if (gg==1) {
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <- 
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] + 
            nfn[[nn-gg]][1:(nn-gg),1:gg,mm-gg+1]
        }
        else {
          nfi[[nn]][2:(nn-gg+1),gg,mm+1] <- 
            nfi[[nn]][2:(nn-gg+1),gg,mm+1] + 
            crossrun::cumsumm(nfn[[nn-gg]][1:(nn-gg),1:gg,mm-gg+1])[,gg]
        }
        nfi[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] <- 
          nfi[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] + 
          nfn[[nn-gg]][1:(nn-gg),(gg+1):(nn-gg),mm-gg+1]
      } # end low g
    } # end iterative procedure nfi
    # iterative procedure nfn, for 1 <= m <= n-1:
    for (mm in 1:(nn-1)) for (gg in 1:(nn-mm)) {
      if (gg>=nn-gg) {
        if (nn-gg==1) 
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <- 
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] + nfi[[1]][1,1,2]
        else {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <- 
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] + 
            crossrun::cumsumm(nfi[[nn-gg]][1:(nn-gg),,mm+1])[,nn-gg]
        }
      } # end high g
      if (gg<nn-gg) {
        if (gg==1) {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <- 
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] + 
            nfi[[nn-gg]][1:(nn-gg),1:gg,mm+1]
        }
        else {
          nfn[[nn]][2:(nn-gg+1),gg,mm+1] <-
            nfn[[nn]][2:(nn-gg+1),gg,mm+1] + 
            crossrun::cumsumm(nfi[[nn-gg]][1:(nn-gg),1:gg,mm+1])[,gg]
        }
        nfn[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] <- 
          nfn[[nn]][2:(nn-gg+1),(gg+1):(nn-gg),mm+1] + 
          nfi[[nn-gg]][1:(nn-gg),(gg+1):(nn-gg),mm+1]
      } # end low g
    } # end iterative procedure nfi
    if (printn==TRUE) print(nn)
    if (printn==TRUE) print(Sys.time())
  } # end for nn
  return(list(nfi = nfi, nfn = nfn))
} # end function crossrunemcont

# check for low n:
em15cont <- crossrunemcont(em10, n1=11, nmax=15,printn=TRUE)
