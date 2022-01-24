# 210124
# 성과 : grouping 되느 모습 확인 / 보완필요한 부분 : 논문이랑 비슷한 양상이 되도록,, 노력하기,,ㅎ 처음 수렴되는 속도가 너무 빠름


# 기본 세팅
dyn.load('Desktop/portfolio_selection/SLOPE_code/cproxSortedL1.so')
source('Desktop/portfolio_selection/SLOPE_code/proxSortedL1.R') # 여기서 ㅊcproxSortedL1을 이용해서 함수를 생성하고, 읽어옴

set.seed(123)
t <- 50 ; r<-3 ; k <- 12 

library(mvtnorm)
# risk factors
risk_factors <- rmvnorm(50, rep(0,r), diag(r))

# vectors of error term
epsilon<- rmvnorm(50, rep(0,k), 0.05*diag(k)) ## 논문에서는 diag(r)로 나와있음,,! 해결 필요

# loading matrix
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

# return matrix
RR <- risk_factors%*%BB + epsilon

# mu vector
mu <- apply(RR, 2, mean)
mu

# covariance matrix
sigma <- t(BB)%*%BB + 0.05*diag(k)
sigma_hat <- cov(RR)


# lambda grid 생성
q <- 0.01
qi <- rep(0,12)
for (i in 1:12) qi[i] <- i*(q/(2*k))
aa <- seq(from = 10^-5, to = 10^2, length=100)/qnorm(1-qi[1])
lambda <- matrix(0,12,100)
for ( j in 1:length(aa)) {
  for ( i in 1:12) {
    lambda[i,j] <- aa[j]*qnorm(1-qi[i])
  }
}



# algorithm
# 초기값 지정
alpha <- rep(0, k)
beta <- 0
w <- rep(0,k)
v <- rep(1/12,k)
e <- rep(1,k)
II <- diag(k)
theta <- 10 # 상대위험회피정도 : theta > 0 (여기서는 임의의 값으로 두었음,,)
eta <- 1 # eta > 0 (여기서는 임의의 값으로 두었음,,)
j <- 0
tau <- 10^(-6)
lambda_1 <- lambda[,1]

while (TRUE) {
  # update
  w <- solve(theta*sigma + eta*(II+e%*%t(e)))%*%(mu-alpha-as.vector(beta%*%e)+eta*(v+e))
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
