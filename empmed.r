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
# tostorage requirements:
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
cl15

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
# auxiliary funtions:
cumsummrow <- function(mtrx) {
  mtrs <- mtrx
  nrw <- nrow(mtrx)
  for (rw in 1:nrw) mtrs[rw, ] <- cumsum(mtrx[rw, ])
  return(mtrs)
} # end function cumsumrow
cumsummcol <- function(mtrx) {
  mtrs <- mtrx
  ncl <- ncol(mtrx)
  for (cl in 1:ncl) mtrs[, cl] <- cumsum(mtrx[, cl])
  return(mtrs)
} # end function cumsummcol

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
     y=(cumsum(cumsummrow(cl14)[,14])/(sum(cl14)))[
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
