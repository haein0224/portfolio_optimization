# linear shrinkage method 
# 출처 : 박세영 교수님 github (https://github.com/ishspsy/project/blob/ishspsy-patch-1/Shrinkage_la.m)

library(psych) # tr계산을 위한 패키지

# shrinkage method에 필요한 lambda값 생성 함수
shrinkgae_la <- function(data) { # observation matrix이 input -> output이 la
  n <- dim(data)[1] ; p <- dim(data)[2] ;
  Si <- cov(data) # sample covariance matrix를 먼저 계산
  
  # lambda of Ledoit and Wolf shrinkage
  # Ledoit and Wolf shrinkage formula : Si=(1-la)*Si+la*(tr(Si)/p)*eye(p) => 여기 논문에서는 la를 반환하는 것!
  normalized_diagonal_sum <- tr(Si)/p
  b <- tr((Si - normalized_diagonal_sum*diag(p))^2)/p
  a <- 0
  
  for (k in 1:n) {
    a <- a+tr((data[k,]%*%t(data[k,])-Si)^2)
  }
  
  a <- a/(n^2*p)
  a <- min(a,b)
  la <- a/b
  return(la) 
}

# mu & cov matrix 생성 함수
mean_cov_comp <- function(Daily_Return,la) {
  r <- Daily_Return
  p <- dim(r)[2]
  
  # mean computation
  mu <- t(apply(r,2,mean)) # 자산별 평균 벡터
  
  # covariance matrix computation
  Si <- cov(r)
  #la <- shrinkgae_la(r)
  Si <- (1-la)*Si + la*(tr(Si)/p)*diag(p)
  
  # report output
  result <- c()
  result$Mean <- mu
  result$Covariance_matrix <- Si
  #result$la <- la

  return(result)
}

# 사용 방법..?
mu <- mean_cov_comp(RR)$
cov <- mean_cov_comp(RR)$Covariance_matrix










