library(MCMCpack)

#Use normal instead of t
#Compute the stick-breaking probabilities
#v: current probability #m: number of clusters
make.probs<-function(v,m){
  if(m==1){
    p=1
  }
  if(m>1){
    p<-v
    p[2:m]<-p[2:m]*cumprod(1-v[2:m-1])
    p[m]<-1-sum(p[2:m-1])
  }
  p
}

############ LOG-GPD FUNCTIONS ############ (Stephenson, 2013)
rBGPD <- function(n, mu, sigma, xi, alpha){
  w = matrix(rexp(2*n), ncol=2, nrow=n)
  w = w/rowSums(w)
  k = ifelse(runif(n)<(1-alpha), 1, 2)
  z = rgamma(n, shape=k, scale=1)
  y = 1/(z*w^alpha)

  if (xi[1] == 0)
    x1 <- mu[,1]-sigma[1]*log(1-exp(-1/y[,1]))
  else x1 <- mu[,1] + sigma[1]*((1-exp(-1/y[,1]))^(-xi[1])-1)/xi[1]

  if (xi[2] == 0)
    x2 <- mu[,2]-sigma[2]*log(1-exp(-1/y[,2]))
  else x2 <- mu[,2] + sigma[2]*((1-exp(-1/y[,2]))^(-xi[2])-1)/xi[2]

  x <- cbind(x1, x2)
return(x)}

############ DENSITY FUNCTIONS ##############
loglike <- function(sigma, xi, alpha, U, Y){

  library(evd)
  if(is.null(dim(Y))){
    Y     <- matrix(Y, ncol=2, nrow=1)
    U     <- matrix(U, ncol=2, nrow=1)
  }

  cdfx  <- pgpd(Y[,1], loc=U[,1], scale=sigma[1], shape=xi[1])
  cdfy  <- pgpd(Y[,2], loc=U[,2], scale=sigma[2], shape=xi[2])
  pdfx  <- dgpd(Y[,1], loc=U[,1], scale=sigma[1], shape=xi[1])
  pdfy  <- dgpd(Y[,2], loc=U[,2], scale=sigma[2], shape=xi[2])
  A     <- (1/alpha-1)*log(-log(cdfx))+
           (1/alpha-1)*log(-log(cdfy))+
           log(pdfx/cdfx)+
           log(pdfy/cdfy)
  V     <- ((-log(cdfx))^(1/alpha)+(-log(cdfy))^(1/alpha))^alpha
  
  den1 <- -V
  den2 <- (1-2/alpha)*log(V)
  den3 <- log(V-1+1/alpha)
  den  <- A+den1+den2+den3
  den[is.na(den)] <- -Inf
return(den)}


############ DENSITY FUNCTIONS ##############
dBGPD <- function(s, U, Y){
  library(evd)
  sigma <- s[[1]]
  xi    <- s[[2]]
  alpha <- s[[3]]
  if(is.null(dim(Y))){
    Y     <- matrix(Y, ncol=2, nrow=1)
    U     <- matrix(U, ncol=2, nrow=1)
  }
  
  cdfx  <- pgpd(Y[,1], loc=U[,1], scale=sigma[1], shape=xi[1])
  cdfy  <- pgpd(Y[,2], loc=U[,2], scale=sigma[2], shape=xi[2])
  pdfx  <- dgpd(Y[,1], loc=U[,1], scale=sigma[1], shape=xi[1])
  pdfy  <- dgpd(Y[,2], loc=U[,2], scale=sigma[2], shape=xi[2])
  A     <- (1/alpha-1)*log(-log(cdfx))+
    (1/alpha-1)*log(-log(cdfy))+
    log(pdfx/cdfx)+
    log(pdfy/cdfy)
  V     <- ((-log(cdfx))^(1/alpha)+(-log(cdfy))^(1/alpha))^alpha
  
  den1 <- -V
  den2 <- (1-2/alpha)*log(V)
  den3 <- log(V-1+1/alpha)
  den  <- A+den1+den2+den3
  
  # res  <- exp(den)
  res <- den
  res[is.na(res)] <- 0
  
  return(res)}



library(Rcpp)
library(inline)

code = '
Rcpp::NumericMatrix lp(l);
int nr = lp.nrow();
int nc = lp.ncol();
Rcpp::NumericVector temp(nc);
Rcpp::NumericVector res(nr);

RNGScope scp;

for(int i=0; i<nr; i++){
NumericVector rn = runif(1);
double loc = 0.0;
int j = 0;
while(rn[0]>=loc){
loc = loc + lp(i, j);
j = j + 1;
}
res[i] = j;
}

return res;
'

mycl <- cxxfunction(signature(l="numeric"), body=code, plugin="Rcpp")

cdxBGP <- function(D, mu, sigma, xi, alpha){
  library(evd)
  
  cdfy  <- pgpd(D[1], loc=mu[1], scale=sigma[1], shape=xi[1])
  cdfx  <- pgpd(D[2], loc=mu[2], scale=sigma[2], shape=xi[2])
  pdfy  <- dgpd(D[1], loc=mu[1], scale=sigma[1], shape=xi[1])
  pdfx  <- dgpd(D[2], loc=mu[2], scale=sigma[2], shape=xi[2])
  A     <- (1/alpha-1)*log(-log(cdfx))+log(pdfx/cdfx)
  V     <- ((-log(cdfx))^(1/alpha)+(-log(cdfy))^(1/alpha))^alpha
  
  den1 <- -V
  den2 <- (1-1/alpha)*log(V)
  den  <- A+den1+den2
  
  res  <- exp(den)
  res[is.na(res)] <- 0
  
  return(res)}




# Here is the main MCMC function

MCMC <- function(Y,iters=5000,burn=1000,update=100,K=3){

 library(MCMCpack)

 # Initial values

  n        <- nrow(Y)
  v        <- rep(0, K)
  p        <- rep(.5, K)
  k        <- rmultinom(n,1,p)
  k        <- apply(k, 2, function(x) which(x==1))
  U        <- Y-1
  logsigma <- matrix(rep(log(apply(Y,2,sd)),K), ncol=2, byrow=T)
  xi       <- matrix(0.01,ncol=2,nrow=K)
  mu       <- matrix(0,ncol=2,nrow=K)
  Q        <- solve(cov(U))
  alpha    <- rep(0.3,K)

 # Cache the log likelihood  

  # scur      <- list() # keep location, scale and shape in a list
  # scur[[1]] <- U
  # scur[[2]] <- exp(logsigma)
  # scur[[3]] <- xi 
  log_cur   <- loglike(exp(logsigma[1,]),xi[1,], alpha[1], U, Y)

 # Keep track of the samples 
  temp <- matrix(0,iters,6)
  colnames(temp) <- c("Sigma1","Sigma2","xi1","xi2","alpha","p")
  keepers <- rep(list(temp), K)
  U_mn <- rep(list(matrix(0, iters, 2)), K)
  U_sig <- matrix(0, iters, 4)
  
  rm(temp)
 # Start the MCMC

 for(iter in 1:iters){
   # param = list()
   # for(j in 1:K){
   #   param[[j]] = list(exp(logsigma[j,]),xi[j,],alpha[j])
   # }
   #  
   # # Update K
   # # Have problem? May not be able to update all together
   # l = matrix(p,ncol=K,nrow=n,byrow=T)+sapply(param,dBGPD,U=U,Y=Y)
   # l[is.na(l)] = 1
   # l = l/rowSums(l)
   # k = mycl(l)
   
   
   
   # k = g
   # # Update p
   for(i in 1:n){
     cank<-k;log_can<-log_cur
     cank[i]<-sample(1:K,1,prob=p)
     log_can[i]<-loglike(exp(logsigma[cank[i],]),xi[cank[i],],alpha[cank[i]], U[i,], Y[i,])
     R<-exp(log_can[i]-log_cur[i])
     if(!is.na(R) & !is.na(cank[i])){if(runif(1)<R){
       k<-cank;log_cur<-log_can
     }}
   }
   
   # p = c(0.5,0.5)
   
   # # Update p
   for(j in 1:K){
     v[j]<-rbeta(1,1+sum(k==j),1+sum(k>j))
   }
   # BBB<-1-sum(log(1-v[-K]))
   # if(!is.na(BBB)){if(BBB>0){
   #   D<-rgamma(1,K-1+1,1-sum(log(1-v[-K])))}}
   # if(D==0){D<-0.01}
   prev_probs<-p
   p<-make.probs(v,K)
   p[p<0] <- 0
   if(is.na(sum(p))){p<-prev_probs}
   print(p)
   # k <- g
   
   for(j in 1:K){     
     ind       <- which(k==j)
     nk        <- length(ind)
     if(nk>=1){
     # Update U

     P         <- chol(solve(Q))
     can       <- 1*matrix(rnorm(2*nk),nk,2)%*%P
     can       <- cbind(U[ind,1]+can[,1],U[ind,2]+can[,2])
     log_can   <- loglike(exp(logsigma[j,]),xi[j,], alpha[j], can, Y[ind,])

     R         <- log_can - log_cur[ind] -
                  0.5*(Q[1,1]*(can[,1]-mu[j,1])^2 + 2*Q[1,2]*(can[,1]-mu[j,1])*(can[,2]-mu[j,2])+Q[2,2]*(can[,2]-mu[j,2])^2)+
                  0.5*(Q[1,1]*(U[ind,1]-mu[j,1])^2 + 2*Q[1,2]*(U[ind,1]-mu[j,1])*(U[ind,2]-mu[j,2])+Q[2,2]*(U[ind,2]-mu[j,2])^2)
     
     R[is.na(R)] <- -Inf
     keep      <- log(runif(nk))<R 
     if((sum(keep)>0)&(nk>1)){
       U[ind,][keep,]     <- can[keep,]
       log_cur[ind][keep] <- log_can[keep]
     }
     if((sum(keep)>0)&(nk==1)){
       U[ind,]      <- can
       log_cur[ind] <- log_can
     } 

   # Update log sigma

     can       <- rnorm(2, logsigma[j,], 0.05)
     log_can   <- loglike(exp(can),xi[j,], alpha[j], U[ind,], Y[ind,])
     R         <- sum(log_can)-sum(log_cur[ind]) +
                  sum(dnorm(can,0,10, log=T)) -
                  sum(dnorm(logsigma[j,],0,10, log=T)) 
     if(!is.na(R)){if((log(runif(1))<R)){
        logsigma[j,] <- can
        log_cur[ind]  <- log_can
     }}

   # Update xi

     can       <- rnorm(2,xi[j,],0.025)
     log_can   <- loglike(exp(logsigma[j,]),can, alpha[j], U[ind,], Y[ind,])
     R         <- sum(log_can)-sum(log_cur[ind]) +
                  sum(dnorm(can,0,0.25, log=T)) -
                  sum(dnorm( xi[j,],0,0.25, log=T)) 
     if(!is.na(R)){if((log(runif(1))<R)){
       xi[j,]       <- can
       log_cur[ind]  <- log_can
     }}

   # Update alpha

     can     <- pnorm(rnorm(1,qnorm(alpha[j]), 0.05))
     log_can <- loglike(exp(logsigma[j,]),xi[j,], can, U[ind,], Y[ind,])
     R       <- sum(log_can)-sum(log_cur[ind]) +
                dnorm(qnorm(can),log=T) -
                dnorm(qnorm(alpha[j]),log=T) 
     if(!is.na(R)){if((log(runif(1))<R)){
       alpha[j]    <- can
       log_cur[ind]  <- log_can
     }}

   # Update mu

    VVV <- solve(nk*Q + diag(2)*0.01)
    MMM <- Q%*%c(sum(U[ind,1]),sum(U[ind,2]))
    mu[j,]  <- VVV%*%MMM + t(chol(VVV))%*%rnorm(2)
    
     }
   }

   # Update Q
   A <- diag(2)
   for(j in 1:K){
     ind <- which(k==j)
     S   <- cbind((U[ind,1]-mu[j,1]),(U[ind,2]-mu[j,2]))
     S   <- t(S)%*%S 
     A   <- A + S
   }
   Q <- rwish(n+2,solve(A))
   U_sig[iter, ] <- c(Q)


   for(j in 1:K){
     # Store the output
     # assign(paste("keepers",j,"[",iter,",]",sep=""),c(exp(logsigma[j,]),xi[j,],df[j],alpha[j]))
     keepers[[j]][iter,] <- c(exp(logsigma[j,]),xi[j,],alpha[j], p[j])
     U_mn[[j]][iter, ] <- c(mu[j, ])
     # U_sig[[j]][iter, ] <- c(get(paste("Q",j,sep="")))}
   }


 } # End MCMC

 out <- list(keepers=keepers,U=U_mn,V=U_sig)


return(out)}


