simul <- function(number, eta) {
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
  base <- data1[500:750,1:ncol(data1)-1] # 마지막 인덱스값은 제외하고 사용
  
  daily_return <- matrix(0,dim(base)[1], dim(base)[2])
  for ( j in 1:ncol(base)) {
    for ( i in 2:nrow(base)) {
      daily_return[i-1,j] <- (base[i,j] - base[i-1,j])/base[i,j]
    }
  }
  
  mu <- as.vector(mean_cov_comp(daily_return[,1:ncol(base)])$Mean)
  sigma <- mean_cov_comp(daily_return[,1:ncol(base)])$Covariance_matrix
  la <- mean_cov_comp(daily_return[,1:ncol(base)])$la
  
  ########## lambda sequence 생성
  k <- length(mu)
  ll <- 30 # 생성할 람다 시퀀스 개수
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
  
  #####################
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
    pi <- k # 상대위험회피정도 : pi > 0 (여기서는 임의의 값으로 두었음,,)
    # eta <-  15 # round(1/ncol(base),2) # eta > 0 (여기서는 임의의 값으로 두었음,,)
    j <- 0
    tau <- 10^(-4)
    lambda_1 <- lambda[,z]
    
    while (TRUE) {
      # update
      w <- solve(pi*sigma + eta*(II+e%*%t(e)))%*%(mu-alpha-as.vector(beta%*%e)+eta*(v+e))
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
    print(z)
  }
  write.csv(ww, file=paste("result_",data_set[number],'.csv',sep=""), row.names=FALSE)
}


simul(6, )

for ( number in 1:5 ) {#length(data_set) ) {
  simul(number, 3)
}


### 시각화
# write.csv(ww, file='weights/FTSE250.csv', row.names = FALSE)
ww <- as.matrix(read.csv("result(k,0.5)/result_FTSE250.csv"))
plot(ww[1,], ylim=c(min(ww),max(ww)), type='l', col=1, lwd=1)#, main=paste("data = ", Dow, ", pi = ", pi, ", eta = ",eta)) # round(eta,3)))
for(i in 2:nrow(ww)) lines(ww[i,], type='l', col=1)

ncol(ww)
