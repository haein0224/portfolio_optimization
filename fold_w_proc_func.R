fold_w_proc <- function(R, valid){
  rb <<- apply(R,mean,MARGIN = 1) # 자산별 평균 계산
  window <- ncol(R)
  
  kvec <- c(rep(1/KK, KK), rep(0, window-KK))
  ww<<-  matrix(0, nrow = nrow(lambda), ncol = ncol(lambda)) # reset ww & rww
  for (nth in 1:ncol(lambda)){
    # 초기값으로 reset
    alpha <- rep(0, p) # theta는 shortfall에서 나온 새로운 모수
    beta <- theta <- 0
    gamma <- z <<- rep(0,window) # ""window""가 fold내 데이터수 뺀 값이 되어야함
    w <- v <- rep(0,p)
    j <- 0
    lambda_1 <- lambda[,nth]
    
    while (TRUE) {
      w <- solve( c(eta)*(II + e%*%t(e) + rb%*%t(rb) + R%*%t(R) )) %*% ( c(eta)*(v+e+c(rp)*rb +R%*%z) - (rb + alpha + c(beta)*e + c(theta)*rb - R%*%gamma))
      v <- proxSortedL1(w+(1/eta)*alpha, lambda_1/eta)
      z <<- proxZ(R,w,gamma,eta,kvec)
      alpha <- alpha + c(eta)*(w-v)
      beta <- beta + eta*(t(e)%*%w - 1)
      theta <- theta + eta*(t(w)%*%rb-rp)
      gamma <- gamma + c(eta)*(z - t(R)%*%w) # ; gamma2 <- gamma2 + c(eta)*(z2 - t(R)%*%w)
      
      # dual gap 계산 및 체크 ## 아닌것 같음!
      G3 <- t(w-v)%*%(w-v) + t(z-t(R)%*%w)%*%(z-t(R)%*%w) + (t(rb)%*%w-rp)^2 + (t(e)%*%w-1)^2 # primal-dual gap으로 계산해야함!
      
      if (G3 < tau) break
      j <- j+1
      if (j>=10000) break
    }
    ww[,nth] <<- w
   # print(nth)
  }
}
