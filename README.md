# let-s_go_graduation
졸업을 위한 나의 노력,,


### 0124 시뮬결과 이미지
![image](https://user-images.githubusercontent.com/84130229/150819872-f112401c-a96d-47a9-93be-f21ff93f141f.png)

문제점 : 100번까지있는데 시뮬레이션에서는 첫번째 수렴이 (40,50) 사이, 전체 수렴이  (60,70) 사이에서 발생해야함 뭐가 문제일까..!


### 0126 Coefficient of Relative Risk Aversion
사전적 정의 :
![image](https://user-images.githubusercontent.com/84130229/151123289-50287abb-93b6-4bb0-8e1f-6f5b0fd95473.png)

#### 오늘의 시뮬결과 이미지
![image](https://user-images.githubusercontent.com/84130229/151176984-db147671-5043-4c00-8f6a-00bcd9ffbad7.png)

0124의 문제점을 해결한 부분 ) lambda sequence를 짤때 그냥 seq function을 사용하는게 아니라 lseq를 사용했어야함!! log-spaced lambda를 얻는것이 필요했으므로
새로운 문제가 발견된 부분 ) 현재 목표가 논문의 시뮬레이션과 같은 결과를 도출하는 것이기때문에 초기값이 0124처럼 고른 분포를 보여야하며,람다값에는 차이가 없는데 초기값이 다르게 나왔다는 것은 다른 파라미터의 문제인지 유의할 필요가 있음
