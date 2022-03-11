setwd('~/Desktop/portfolio_selection')
dyn.load('SLOPE_code/cproxSortedL1.so')
source('SLOPE_code/proxSortedL1.R') # 여기서 ㅊcproxSortedL1을 이용해서 함수를 생성하고, 읽어옴
library(mvtnorm) # generation simulation data
library(emdbook) # lseq lambda grid
library(readxl)
library(quadprog)
# 앞으로 데이터 셋 불러올때 일반화 해서 불러오기 위한 작업
data_set <- c("DOW30", "DAX30", "SP100", "FTSE100", "FTSE250", "SP500")

number <- 1
title <- paste("data/",data_set[number],".xlsx", sep="")
print(title)
data1 <- read_excel(title ,sheet = 1, col_names = TRUE)
data1 <- as.matrix(data1[,seq(from=1,to=ncol(data1),2)])

# 결측치 제거
h <- c()
for ( i in 1:ncol(data1)) { if (sum(is.na(data1[,i])>0)) h <- append(h, i) }
data1 <- data1[,-h] # 결측치가 있는 열의 경우 삭제
dim(data1)

# 기본 input 정의 # window : 250
window <- 250
start <- 100
end <- start+window
base <- data1[start:end,1:ncol(data1)-1] # 마지막 인덱스값은 제외하고 사용

daily_return <- matrix(0,window, dim(base)[2])
for ( j in 1:ncol(base)) {
  for ( i in 2:nrow(base)) {
    daily_return[i-1,j] <- (base[i,j] - base[i-1,j])/base[i-1,j]  # 지금까지 틀렸던거 찾아서 수정..
  }
} # i 인덱스가 1~T까지 돌아감
R <- t(daily_return) # dimension : p*T

# r-bar
rb <- apply(R,mean,MARGIN = 1) # 자산별 250일 평균 계산
# rb <- (R%*%rep(1,window))/window # 같은값 계산하는 코드

# simulation data용
#rb <- mu 
#R <- t(daily_return) # dimension : p*T
#rp <- 0.06  #임의로 지정
#window <- nrow(R)

# rp : 목표 수익률 -> 포트폴리오 평균 수익률을 써야하나..? -> 사실 이건 hyperparameter가 맞긴함!
#rp_prep <- rep(0,window)
#for( i in 2:nrow(base)) rp_prep[i-1] <- (data1[start:end,ncol(base)][i]-data1[start:end,ncol(base)][i-1])/data1[start:end,30][i-1]
#rp <- mean(rp_prep) # 기간동안 포트폴리오 평균 수익률을 사용
rp <- 0.02 # 임의로 지정하는 경우

########## lambda sequence 생성
p <- length(rb) # 고려하는 자산의 개수
ll <- 20 # 생성할 람다 시퀀스 개수
q <- 0.01
qi <- rep(0,p)
for (i in 1:p) qi[i] <- i*(q/(2*p))
aa <- lseq(from = 10^-5, to = 10^2, length=ll)/qnorm(1-qi[1]) # lseq 사용해야함 / real data에 대해서는 범위가 다름
lambda <- matrix(0,p,ll)
for ( j in 1:length(aa)) {
  for ( i in 1:p) {
    lambda[i,j] <- aa[j]*qnorm(1-qi[i])
  }
}

#####################
ww <-  matrix(0, nrow = nrow(lambda), ncol = ncol(lambda))
jj <- rep(0, ncol(lambda))
for ( nth in 1:ncol(lambda)){
  print(nth)
  # algorithm
  # 초기값 지정
  alpha <- rep(0, p) # theta는 shortfall에서 나온 새로운 모수
  beta <- theta <- 0
  gamma <- z <- rep(0,window)
  w <- v <- rep(0,p) # v <- rep(0,k)
  #v <- rep(1/p,p) #초기값은 다 0으로 놓긴했었음
  e <- rep(1,p)
  II <- diag(p)
  q <- 0.05 # quantile값
  KK <- floor(q*window) # big K
  # phi <- 3 #15 #40 #k # 상대위험회피정도 : phi > 0 (여기서는 임의의 값으로 두었음,,)
  eta <- 20 # round(1/ncol(base),2) # eta > 0 (여기서는 임의의 값으로 두었음,,)
  j <- 0
  tau <- 10^(-5) # 이거에 따라서도 차이가 심함 check 필요!!
  lambda_1 <- lambda[,nth]
  
  while (TRUE) {
    
    vv <- matrix(0,length(w),2000) ; zz <- matrix(0,window,2000)
    for (i in 1:2000) {
    
    # update
    w <- solve(as.vector(eta)*(II + R%*%t(R) + rb%*%t(rb) + e%*%t(e)))%*%(as.vector(eta)*(v + R%*%z + as.vector(rp)*rb + e) - (rb + alpha + as.vector(beta)*e - R%*%gamma + as.vector(theta)*rb))
    v <- proxSortedL1(w+(1/eta)*alpha, lambda_1/eta)
    z <- proxZ(R,w,gamma,eta,KK)
    alpha <- alpha + eta*(w-v)
    beta <- beta + eta*(t(e)%*%w - 1)
    theta <- theta + eta*(t(w)%*%rb-rp)
    gamma <- gamma + eta*(z - t(R)%*%w)
    
    vv[,i] <- 2000 ; zz[,i] <- 2000
    }
  
    # dual gap 계산 및 체크 ## 아닌것 같음!
    G <- abs(-t(alpha + as.vector(beta%*%e) + as.vector(theta%*%rb))%*%w + beta + sum(lambda_1*sort(abs(v), decreasing = TRUE))) # 계산 필요
    if (G < tau) break
    
    j <- j+1
    if (j>=3000) break
  }
  ww[,nth] <- w
  print(nth)
}

plot(ww[1,1:8], ylim=c(min(ww),max(ww)),type='l', lwd=1)#, main=paste("data = ", Dow, ", phi = ", phi, ", eta = ",eta)) # round(eta,3)))
for(i in 2:nrow(ww)) lines(ww[i,1:8], type='l') #, col=ddd[i])

plot(ww[1,1:24], ylim=c(-0.5,0.5)) #, type='l', lwd=1)#, main=paste("data = ", Dow, ", phi = ", phi, ", eta = ",eta)) # round(eta,3)))
for(i in 2:nrow(ww)) points(ww[i,1:24])#, type='l') #, col=ddd[i])

vv <- ww # 1:17, 50은 채워져있음!
