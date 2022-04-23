# fast proxZ
# QP 방식으로 업데이트하는 방식은 속도가 너무 오래 걸리기 때문에
# fast proximal 방식을 사용해 업데이트를 진행하고자함

fastproxZ <- function(R, w, gamma, eta, kvec) {
  # 처음에 작 -> 큰 순서로 정렬된
  dtt <- data.frame(cbind(t(R)%*%w, gamma))
  colnames(dtt) <- c('y','gamma')
  idx <- order(dtt$y) # original order를 저장
  y <- sort(dtt$y)
  
  n_gamma <- dtt[idx,]$gamma
  
  #kvec <- c(rep(1/KK, KK),rep(0, ncol(R)-KK))
  
  #prox <- new_prox <- as.vector(1/eta) * (k + as.vector(eta)*y - n_gamma)
  prox <- new_prox <- (kvec + as.vector(eta)*y - n_gamma)/eta # nondecreasing하게 만들고 싶은 대상
  #prox[1:KK] <- prox[1:KK] + 1/KK
  
  # decr <- TRUE # TRUE로 넣어 놓고 시작
  #set.seed(123)
  #prox <- runif(20) # i # 생성한 데이터 (여기선 삭제할 예정)
  n <- 0
  # 반복할 구간
  #for ( n in 1:length(prox)) {
  while (TRUE) {

    # pt <- which(diff(prox) < 0 ) 
    pt <- which(round(diff(prox),3)<0) # strictly decreasing이 일어나는 지점의 index
    ptt <- c(0,(pt+1)[1:length(pt)-1])
    
    #ept <- which(diff(c(prox,prox[length(prox)])) >= 0 ) 
    ept <- which(round(diff(c(prox,prox[length(prox)])),3) >= 0 ) 
    eptt <- c(0, ept[1:length(ept)-1] + 1)
    
    # 구간의 시작 포인트, 마지막 포인트 지정용 벡터 생성
    startpt <- pt[pt!=ptt] ; if (length(startpt)==0) {break} # 더이상 바꿀게 없으면 멈춤
    endpt <- ept[ept!=eptt]
    dt <- data.frame(startpt, endpt = endpt[which(endpt>min(startpt))])
    
    # 구간별로 업데이트
    for ( i in 1:nrow(dt)) {
      start <- dt[i,'startpt'] ; end <- dt[i,'endpt']
      new_prox[start:end] <- sum(prox[start:end])/(end-start+1)
    }
    
    # 다음 단계 필요한 경우
    prox <- new_prox # 새로운 데이터로 업데이트
    n <- n+1
  }
  
  z[idx] <- prox # z를 원래 인덱스로 되돌리기
  
  return(z)
}


proxZ(R,w,gamma,eta,kvec)
fastproxZ(R,w,gamma,eta,kvec)

# 차분 
diff(c(1,2,3))
which(diff(prox)<0)
