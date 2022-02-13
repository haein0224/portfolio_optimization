# grid..

www <- array(dim = c(5,8,12,30))
pi_list <- seq(3,12)
eta_list <- 5 #c(0.05, 0.1, 0.5, 1, 3, 5, 10, 15)

wwww <-array(dim=c(10,12,30))

for ( x in 1:length(eta_list) ) {
  for ( y in 1:length(pi_list) ) {
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
      pi <- pi_list[y] #k # 상대위험회피정도 : pi > 0 (여기서는 임의의 값으로 두었음,,)
      eta <- eta_list[x] # 15 # round(1/ncol(base),2) # eta > 0 (여기서는 임의의 값으로 두었음,,)
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
      wwww[y,,z] <- w
      #www[y,x, ,z] <- w
      #vv[,z] <- v
      #jj[z] <- j
      #gg[z] <- G
      print(z)
      #if (z == ncol(lambda)) t0 <- proc.time()[3] - t0
    }
  }
}




ylb <- "Weights" ; xlb <- "Grid Lambda"
par(mfrow=c(4,2))
for (c in 1:length(pi_list)) {
   for(j in 1:length(eta_list)){
    title <- paste("pi = ", pi_list[c], ", eta = ", eta_list[j])
    plot(www[c,j,1,], ylim=c(min(www[1,j,,]), max(www[c,j,,])), type='l', col=1, lwd=1, main=title, xlab = xlb, ylab=ylb) # round(eta,3)))
    for(i in 2:4) lines(www[c,j,i,], type='l', col=1, lwd=1)
    for(i in 5:8) lines(www[c,j,i,], type='l', col=2,lwd=1)
    for(i in 9:12) lines(www[c,j,i,], type='l', col=3, lwd=1)
  }
}

c <- eta_list[1]
for (j in 1:length(pi_list)) {
  title <- paste("pi = ", pi_list[j], ", eta = ", c)
  plot(wwww[j,1,], ylim=c(min(wwww[j,,]), max(wwww[j,,])), type='l', col=1, lwd=1, main=title, xlab = xlb, ylab=ylb) # round(eta,3)))
  for(i in 2:4) lines(wwww[j,i,], type='l', col=1, lwd=1)
  for(i in 5:8) lines(wwww[j,i,], type='l', col=2,lwd=1)
  for(i in 9:12) lines(wwww[j,i,], type='l', col=3, lwd=1)
}

par(mfrow=c(1,1))

title <- paste("pi = ", pi_list[j], ", eta = ", c)
plot(wwww[j,1,], ylim=c(min(wwww[j,,]), max(wwww[j,,])), type='l', col=1, lwd=1, main=title, xlab = xlb, ylab=ylb) # round(eta,3)))
for(i in 2:4) lines(wwww[j,i,], type='l', col=1, lwd=1)
for(i in 5:8) lines(wwww[j,i,], type='l', col=2,lwd=1)
for(i in 9:12) lines(wwww[j,i,], type='l', col=3, lwd=1)
