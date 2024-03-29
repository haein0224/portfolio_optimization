###### 0504 : SLOPE & LASSO


real_viz <- function(ww, month, w_idx){
  plot(ww[1,], ylim=c(min(ww),max(ww)),type='l', lwd=1, main=paste(" month = ",month," selected_idx = ", w_idx,"th")) # round(eta,3)))
  for(i in 2:nrow(ww)) lines(ww[i,], type='l')
}


#################
par(mfrow=c(4,1))

#################
library(doParallel)
n_core = detectCores() # 사용가능한 코어가 몇개인지 확인하는 용
cl = makeCluster(n_core-3) # 사용가능한 코어중 하나를 빼놓고 할당(코어 5개 사용)
cl
registerDoParallel(cl)

results_lasso <- data.frame("month"=0, "lambda"=0, "portfolio"=0, "lasso"=0,"index"=0, "EW"=0) # 최종 결과를 저장할 데이터 프레임 생성


######### 데이터 준비
data_prep("DOW30", 1)

########## lambda sequence 생성
lamb_seq(30) # 람다 시퀀스 생성 함수 (몇개를 생성할 것인지 지정)

lambda <- cbind(lambda, rep(lambda[1,17],29))
##################
# 기본 input 정의
window <- 250 # 11개월치 사용
start <- 1

t0 <- proc.time()[3] # 시간 측정
for (month in 1:12) {
  end <- start+window
  base <- data1[start:end,1:ncol(data1)-1] # 마지막 인덱스값은 제외하고 사용
  Nwindow <- Twindow <- 21
  
  daily_return <- matrix(0,window, dim(base)[2])
  for ( j in 1:ncol(base)) {
    for ( i in 2:nrow(base)) {
      daily_return[i-1,j] <- (base[i,j] - base[i-1,j])/base[i-1,j]  # 지금까지 틀렸던거 찾아서 수정..
    }
  } # i 인덱스가 1~T까지 돌아감
  
  
  R <- t(daily_return) # dimension : p*T
  
  # r-bar
  rb <- apply(R,mean,MARGIN = 1) # 자산별 250일 평균 계산
  
  
  ww <- foreach(lambda=iter(lambda, by='column'), .combine=cbind,.packages = "quadprog") %dopar% { #, .export=c("proxSortedL1", "proxZ")) %dopar% { #, .combine=cbind) {
    dyn.load('~/Desktop/portfolio_selection/SLOPE_code/cproxSortedL1.so')
    alpha <- rep(0, p)
    beta <- theta <- 0
    gamma <- z <- rep(0,window)
    w <- v <- rep(0,p)
    e <- rep(1,p)
    II <- diag(p)
    q <- 0.1 # quantile값..
    KK <- floor(q*window) # big K
    kvec <- c(rep(1/KK, KK),rep(0, ncol(R)-KK))
    eta <- 4 # 실험용(빠르게 작동하는것 알고있는 값으로,,)
    j <- 0
    rp <- 0.005
    tau <- 10^(-5) # 실험용(빠르게 작동하는것 알고있는 값으로,,)
    
    while (TRUE) {
      w <- solve( c(eta)*(II + e%*%t(e) + rb%*%t(rb) + R%*%t(R) )) %*% ( c(eta)*(v+e+c(rp)*rb +R%*%z) - (rb + alpha + c(beta)*e + c(theta)*rb - R%*%gamma))
      v <- proxSortedL1(w+(1/eta)*alpha, lambda/eta)
      z <- proxZ(R,w,gamma,eta,kvec)
      alpha <- alpha + c(eta)*(w-v)
      beta <- beta + eta*(t(e)%*%w - 1)
      theta <- theta + eta*(t(w)%*%rb-rp)
      gamma <- gamma + c(eta)*(z - t(R)%*%w) 
      
      # dual gap 계산 및 체크 ## 아닌것 같음!
      Gap <- t(w-v)%*%(w-v) + t(z-t(R)%*%w)%*%(z-t(R)%*%w) + (t(rb)%*%w-rp)^2 + (t(e)%*%w-1)^2 # primal-dual gap으로 계산해야함!
      
      if (Gap < tau) break
      
      j <- j+1
      if (j>=5000) break
    }
    return(w) # 결과 내보내기
  }
  
  
  # 평가
  # 결정된 람다를 적용해 다음 21간의 수익률을 계산 및 저장
  Nstart <- end+1 #Tend+1
  Nend <- Nstart+Nwindow
  base <- data1[Nstart:Nend,1:ncol(data1)]
  
  N_daily_return <- matrix(0,Nwindow, dim(base)[2])
  for ( j in 1:ncol(base)) {
    for ( i in 2:nrow(base)) {
      N_daily_return[i-1,j] <- (base[i,j] - base[i-1,j])/base[i-1,j]  # 지금까지 틀렸던거 찾아서 수정..
    }
  }
  
  idx_return <- N_daily_return[,ncol(N_daily_return)]
  N_daily_return <- N_daily_return[,-ncol(N_daily_return)]
  
  ###### 최종적으로 내보낼때는 평가 대상 : 인덱스 수익률 & Equally Weighted 데이터(데이터 프레임에 순서대로 저장)
  #nreturn <- mean(N_daily_return%*%www) # 평균 수익률
  return_17 <- mean(N_daily_return%*%round(ww[,17],2))
  return_lasso <- mean(N_daily_return%*%round(ww[,31],2))
  idx <- mean(idx_return) # index mean return
  ew <- mean(N_daily_return%*%rep(1/p, p)) # equally weighted portfolio mean return
  
  #######  진행상황 및 grouping 확인용 
  #print(paste("month=",month,"& index=",w_idx, "(",sum(www),")"))
  #print(c(nreturn, idx, ew))
  print(c(return_17, return_lasso, idx, ew))
  print(table(round(ww[,17],2)))
  print(table(round(ww[,31],2)))
  
  #results3 <- rbind(results3, c(month, w_idx, nreturn, idx, ew))
  results_lasso <- rbind(results_lasso, c(month, 17, return_17, return_lasso, idx, ew))
  start <- start+21
  
  # 시각화
  #real_viz(ww, month, w_idx)
  #abline(v=w_idx, col="red", lty=2)
  #abline(v=17,col='blue')
  if (month == 12 ) { t0 <- proc.time()[3] - t0 }
}

stopImplicitCluster()

View(results_lasso[-1,-2])
