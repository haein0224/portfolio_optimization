###### 사전작업
setwd('~/Desktop/portfolio_selection')
dyn.load('SLOPE_code/cproxSortedL1.so')
source('SLOPE_code/proxSortedL1.R') # 여기서 ㅊcproxSortedL1을 이용해서 함수를 생성하고, 읽어옴
library(mvtnorm) # generation simulation data
library(emdbook) # lseq lambda grid

###### Data Preparing
# lambda grid 생성
# 한 '행'이 하나의 람다 시퀀스를 지칭
t <- 50 ; r<-3 ; k <- 12
jth <- seq(from = 10^-5, to = 10^2, length=100)
jlist <- 1:100
ilist <- 1:12
lambda_mat <- matrix(0, 100, 12)
for (j in jlist) {
  alpha <- jth[j]/qnorm(1-0.01/(2*k))
  for (i in ilist) {
    lambda_mat[j,i] <- alpha*(qnorm(1-(0.01*i)/(2*k)))
  }
}


# 공분산행렬 추정
col1 <- c(0.77,0.64,0)
col2 <- c(0.9,0,-0.42)
col3 <- c(0, 0.31,0.64)

BB <- matrix(0,3,12)
for ( i in (1:4)) {
  BB[,i] = t(t(col1))
  BB[,i+4] = t(t(col2))
  BB[,i+8] = t(t(col3))
}
BB

II <- diag(12)

sigma <- t(BB)%*%BB + 0.05*II    # simulation sigma

cor <- matrix(0,12,12)
for ( i in 1:12) {
  for (j in 1:12) {
    cor[i,j] <- sigma[i,j]/sqrt(sigma[i,i]*sigma[j,j])
  }
}

heatmap(cor, scale='none')

III <- diag(3) # epsilon 만들기 위한것
IIII <- diag(12)
library(mvtnorm)
epsilon <- rmvnorm(50,rep(0,12),0.05*IIII)

FF <- rmvnorm(50,rep(0,3),III)

RR <- FF%*%BB + epsilon
mu <- apply(RR, 2, mean) # simulation mu

sigma_hat <- cov(RR)


###### lambda sequence 생성
k <- length(mu)
ll <- 100 # 생성할 람다 시퀀스 개수
q <- 0.01
qi <- rep(0,k)
for (i in 1:k) qi[i] <- i*(q/(2*k))
aa <- lseq(from = 10^-5, to = 10^2, length=ll)/qnorm(1-qi[1]) # lseq 사용해야함
lambda <- matrix(0,k,ll)
for ( j in 1:length(aa)) {
  for ( i in 1:k) {
    lambda[i,j] <- aa[j]*qnorm(1-qi[i])
  }
}

####### Algorithm
t0 <- proc.time()[3] # time counting
ww <- vv <- matrix(0, nrow = nrow(lambda), ncol = ncol(lambda))
gg <- jj <- rep(0, ncol(lambda))
for ( z in 1:ncol(lambda) ){
  # algorithm
  # 초기값 지정
  alpha <- rep(0, k)
  beta <- 0
  w <- rep(0,k)
  v <- rep(1/k,k)
  e <- rep(1,k)
  II <- diag(k)
  phi <- 40 # 상대위험회피정도 : phi > 0 (여기서는 임의의 값으로 두었음,,)
  eta <-  15 # round(1/ncol(base),2) # eta > 0 (여기서는 임의의 값으로 두었음,,)
  j <- 0
  tau <- 10^(-4)
  lambda_1 <- lambda[,z]
  
  while (TRUE) {
    # update
    w <- solve(phi*sigma + eta*(II+e%*%t(e)))%*%(mu-alpha-as.vector(beta%*%e)+eta*(v+e))
    v <- proxSortedL1(w+(1/eta)*alpha, lambda_1)
    alpha <- alpha + eta*(w-v)
    beta <- beta + eta*(t(e)%*%w - 1)
    
    # dual gap 계산 및 체크
    G <- abs(-t(alpha + as.vector(beta%*%e))%*%w + beta + sum(lambda_1*sort(abs(v), decreasing = TRUE)))
    if (G < tau) break
    
    # 다음 단계 넘어가기 전에 전 단계 값을 저장
    preA <- alpha
    preB <- beta
    preW <- w
    preV <- v
    
    j <- j+1
    if (j>=10000) break
  }
  ww[,z] <- w
  # vv[,z] <- v
  # jj[z] <- j
  # gg[z] <- G
  print(z)
  if (z == ncol(lambda)) t0 <- proc.time()[3] - t0 # time counting
}

##### Visualization
plot(ww[1,], ylim=c(min(ww),max(ww)), type='l', col=1, lwd=2, main="data = simulation, phi = 40, eta = 15", xlab=xlb, ylab=ylb) # round(eta,3)))
for(i in 2:4) lines(ww[i,], type='l', col=1, lwd=2)
for(i in 5:8) lines(ww[i,], type='l', col=2, lwd=2)
for(i in 9:12) lines(ww[i,], type='l', col=3, lwd=2)
