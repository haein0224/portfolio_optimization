######## 데이터셋 준비 함수
# 아직 다른 시트에 있는 데이터 등은 어떻게 처리할 것인지 결정 및 코드 수정해야함
# data_set <- c("DOW30", "DAX30", "SP100", "FTSE100", "FTSE250", "SP500")

data_prep <- function(data_name, sheet) { # number : 몇번째 데이터셋을 사용할 것인지 지정 
  title <- paste("data/",data_name,".xlsx", sep="")
  print(title)
  data1 <<- read_excel(title ,sheet = sheet, col_names = TRUE) # 두번째 시트사용
  data1 <<- as.matrix(data1[,seq(from=1,to=ncol(data1),2)])
  
  # 결측치 제거
  h <- c()
  for ( i in 1:ncol(data1)) { if (sum(is.na(data1[,i])>0)) h <- append(h, i) }
  data1 <- data1[,-h] # 결측치가 있는 열의 경우 삭제
  data1 <<- data1[order(1:nrow(data1), decreasing=TRUE),] # 역순으로 정렬해야 시간의 순서대로 됨
  print(paste("nrow : ", dim(data1)[1], "ncol : ", dim(data1)[2]))
}
