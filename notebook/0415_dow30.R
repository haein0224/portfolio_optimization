
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

##################3
jth <- 1
window <- 250
Nwindow <- 21
start<-2000
nreturn <- bestl <- rep(0,12)
ww_gain <- matrix(0, 12, ncol(data1)-1)
t0 <- proc.time()[0]

par(mfrow=c(3,1))
for ( jth in 1:12) {
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
  rp <- 0.05
  # simulation data용
  #rb <- mu
  #R <- t(RR)
  #rp <- 0.05  #임의로 지정
  #window <- ncol(R)
  
  # rp : 목표 수익률 -> 포트폴리오 평균 수익률을 써야하나..? -> 사실 이건 hyperparameter가 맞긴함!
  #rp_prep <- rep(0,window)
  #for( i in 2:nrow(base)) rp_prep[i-1] <- (data1[start:end,ncol(base)][i]-data1[start:end,ncol(base)][i-1])/data1[start:end,30][i-1]
  #rp <- mean(rp_prep) # 기간동안 포트폴리오 평균 수익률을 사용
  # 임의로 지정하는 경우
  
  ########## lambda sequence 생성
  
  p <- length(rb) # 고려하는 자산의 개수
  ll <- 30 # 생성할 람다 시퀀스 개수
  q <- 0.01
  qi <- rep(0,p)
  for (i in 1:p) qi[i] <- i*(q/(2*p))
  aa <- lseq(from = 10^-5, to = 10, length=ll)/qnorm(1-qi[1]) # 1까지로 설정해서 봄 !!(원래는 10^2까지인것 까먹으면 안됨)
  lambda <- matrix(0,p,ll)
  for ( j in 1:length(aa)) {
    for ( i in 1:p) {
      lambda[i,j] <- aa[j]*qnorm(1-qi[i])
    }
  }
  
  #####################
  ww <-  matrix(0, nrow = nrow(lambda), ncol = ncol(lambda))
  #t0 <- proc.time()[0]
  for (nth in 1:ncol(lambda)){
    # algorithm
    # 초기값 지정
    alpha <- rep(0, p) # theta는 shortfall에서 나온 새로운 모수
    beta <- theta <- 0
    gamma <- z <- rep(0,window)
    w <- v <- rep(0,p)
    e <- rep(1,p)
    II <- diag(p)
    q <- 0.04 # quantile값
    KK <- floor(q*window) # big K
    kvec <- c(rep(1/KK, KK),rep(0, ncol(R)-KK))
    eta <- 0.03 #70 #100 #5 #30 #0.01 #2 # 10 #0.02 # round(1/ncol(base),2) # eta > 0 (여기서는 임의의 값으로 두었음,,)
    j <- 0
    tau <- 10^(-2) # 이거에 따라서도 차이가 심함 check 필요!!
    lambda_1 <- lambda[,nth]
    
    while (TRUE) {
      w <- solve( c(eta)*(II + e%*%t(e) + rb%*%t(rb) + R%*%t(R) )) %*% ( c(eta)*(v+e+c(rp)*rb +R%*%z) - (rb + alpha + c(beta)*e + c(theta)*rb - R%*%gamma))
      v <- proxSortedL1(w+(1/eta)*alpha, lambda_1/eta)
      z <- proxZ(R,w,gamma,eta,kvec)
      alpha <- alpha + c(eta)*(w-v)
      beta <- beta + eta*(t(e)%*%w - 1)
      theta <- theta + eta*(t(w)%*%rb-rp)
      gamma <- gamma + c(eta)*(z - t(R)%*%w) # ; gamma2 <- gamma2 + c(eta)*(z2 - t(R)%*%w)
      
      # dual gap 계산 및 체크 ## 아닌것 같음!
      G3 <- t(w-v)%*%(w-v) + t(z-t(R)%*%w)%*%(z-t(R)%*%w) + (t(rb)%*%w-rp)^2 + (t(e)%*%w-1)^2 # primal-dual gap으로 계산해야함!
      #G3 <- abs(-t(alpha+c(beta)*e + c(theta)*rb - R%*%gamma)%*%w - t(gamma)%*%z + c(beta) + c(theta)*rp)
      
      
      if (G3 < tau) break
      
      j <- j+1
      if (j>=10000) break
    }
    ww[,nth] <- w
    #print(nth)
    #if (nth == ncol(lambda)) { t0 <- proc.time()[0]-t0}
  }
  
  real_viz(ww,eta)
  
  # 평가지표 생성
  # 직전 한달 간의 수익률 -> 여기서 best인 lambda sequence를 적용해서 다음 한달간의 수익을 계산
  #cwindow <- 21
  #cstart <- end-cwindow
  #cend <- end
  #base <- data1[cstart:cend,1:ncol(data1)-1] # 마지막 인덱스값은 제외하고 사용
  
  #c_daily_return <- matrix(0,cwindow, dim(base)[2])
  #for ( j in 1:ncol(base)) {
  #  for ( i in 2:nrow(base)) {
  #    c_daily_return[i-1,j] <- (base[i,j] - base[i-1,j])/base[i-1,j]  # 지금까지 틀렸던거 찾아서 수정..
  #  }
  #}
  
  
  c_daily_return <- daily_return[250-21:250,] # 가장 최근 21일 데이터로 계산
  CthR <- c_daily_return%*%round(ww,2) # 1%단위로 투자할 것이므로 이를 이용해 평가
  
  #max(colSums(CthR)) # 각 weight sequence에 이번 21간의 수익률
  bestl[jth] <- which(colSums(CthR)==max(colSums(CthR))) # 다음달에 적용할 selected lambda index
  ww_gain[jth,]<-round(ww[,bestl[jth]],2)
  abline(v=bestl[jth])
  #bestl
  
  Nstart <- end+1
  Nend <- Nstart+Nwindow
  base <- data1[Nstart:Nend,1:ncol(data1)-1] # 마지막 인덱스값은 제외하고 사용
  
  N_daily_return <- matrix(0,Nwindow, dim(base)[2])
  for ( j in 1:ncol(base)) {
    for ( i in 2:nrow(base)) {
      N_daily_return[i-1,j] <- (base[i,j] - base[i-1,j])/base[i-1,j]  # 지금까지 틀렸던거 찾아서 수정..
    }
  }
  nreturn[jth] <- sum(N_daily_return%*%round(ww[,bestl[jth]],2))
  print(paste("month=",jth,"&index=",bestl[jth]))
  print(table(round(ww[,bestl[jth]],2)))
  jth <- jth+1
  start <- start+21
  if(jth == 12) { t0 <- proc.time()[0]-t0}
}

