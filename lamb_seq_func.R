###### 람다 시퀀스 생성 함수
lamb_seq <- function(ll) { # ll : 생성하려는 람다 시퀀스의 개수
  p <<- dim(data1)[2]-1 # 고려하는 자산의 개수
  q <- 0.01
  qi <- rep(0,p)
  for (i in 1:p) qi[i] <- i*(q/(2*p))
  aa <- lseq(from = 10^-5, to = 10, length=ll)/qnorm(1-qi[1]) # 1까지로 설정해서 봄 !!(원래는 10^2까지인것 까먹으면 안됨)
  lambda <<- matrix(0,p,ll)
  for ( j in 1:length(aa)) {
    for ( i in 1:p) {
      lambda[i,j] <<- aa[j]*qnorm(1-qi[i])
    }
  }
}
