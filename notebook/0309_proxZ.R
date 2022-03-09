# z 업데이트를 위한 함수
# 기본 아이디어

proxZ <- function(R, w, gamma, eta, K) {
  # 1. z를 원래 인덱스로 돌려놓기 위해서 w를 데이터 프레임으로 변형해 인덱스 부여 후 sorting
  y <- t(R)%*%w + (1/eta)*gamma
  idx <- order(y) # 크기 순서를 저장
  y <- sort(y)
  A <- diag(1,ncol(R), ncol(R))
  A[1,1] <- 0 ; for ( i in 1:ncol(R)) { A[i, i-1] <- -1 }
  bvec <- rep(0,ncol(R))
  Dmat <- diag(eta, ncol(R), ncol(R))
  kvec <- c(rep(-1/K, K),rep(0, ncol(R)-K))
  dvec <- (eta*y-kvec)
  
  # 계산
  res <- solve.QP(Dmat,dvec, A, bvec)$solution
  
  # z를 원래 인덱스로 되돌리기
  return(res[idx])
}



# 다르게 접근한 것
proxZ <- function(R, w, gamma, eta, K) {
  y <- t(R)%*%w
  idx <- order(y) # original order를 저장
  y <- sort(y)
  n_gamma <- gamma[idx] # 감마도 ordering해서 사용
  A <- diag(1,ncol(R), ncol(R))
  A[1,1] <- 0 ; for ( i in 1:length(w)) { A[i, i-1] <- -1 }
  bvec <- rep(0,ncol(R))
  kvec <- c(rep(-1/K, K),rep(0, ncol(R)-K))
  dvec <- (-n_gamma - kvec + eta*y)
  Dmat <- diag(eta, ncol(R), ncol(R))
  
  # 계산
  res <- solve.QP(Dmat,dvec, A, bvec)$solution
  
  # z를 원래 인덱스로 되돌리기
  return(res[idx])
}




# 실험용
a <- c(-1,5,8,2,10,9,135) # y
b <- as.data.frame(cbind(order(a), sort(a))) # b$V1 : original order를 저장하고있음
b[,2] # ordered value of a # sorted된 value
b[order(b$V1),2] # ==a (original value)
