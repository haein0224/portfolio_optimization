proxZ <- function(R, w, gamma, eta, K) {
  dt <- data.frame(cbind(t(R)%*%w, gamma))
  colnames(dt) <- c('y','gamma')
  idx <- order(dt$y) # original order를 저장
  y <- sort(dt$y)
  
  n_gamma <- dt[idx,]$gamma #rep(0, length(gamma))
  #n_gamma <- gamma # 감마도 index따라 정렬해서 사용. -> 근데 지금 제대로 안먹음... 슈방
  
  #A <- matrix(0, ncol(R), ncol(R))
  A <- diag(1,ncol(R), ncol(R))
  for ( i in 1:ncol(R)) { A[i, i-1] <- -1 }
  A <- A[2:ncol(R),]
  #A <- rbind(diag(1,ncol(R),ncol(R)), A) # 모든게 맞다고 값을 constraint로 넣은 경우 -> 근데 이러면 안되는거 아닌감?.?
  #bvec <- c(y,rep(0,ncol(R)-1))
  bvec <- rep(0,ncol(R)-1)
  kvec <- c(rep(1/K, K),rep(0, ncol(R)-K))
  dvec <- ( kvec + as.vector(eta)*y - n_gamma)
  Dmat <- diag(eta, ncol(R), ncol(R))
  
  # 계산
  res <- solve.QP(Dmat,dvec, t(A), bvec)$solution #, meq = ncol(R))$solution
  z[idx] <- res
  # z를 원래 인덱스로 되돌리기
  return(z)
}
