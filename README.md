# let-s_go_graduation
연구보고 archive : https://www.notion.so/c26577b9d19942e8a7f3ca868d50b629?v=b8de0316fced46a5aac8e2efb3bf051d 💻✔️


### 0124 시뮬결과 이미지
![image](https://user-images.githubusercontent.com/84130229/150819872-f112401c-a96d-47a9-93be-f21ff93f141f.png)

문제점 : 100번까지있는데 시뮬레이션에서는 첫번째 수렴이 (40,50) 사이, 전체 수렴이  (60,70) 사이에서 발생해야함 뭐가 문제일까..! -> 1026 해결완료


### 0126 Coefficient of Relative Risk Aversion
사전적 정의 :
![image](https://user-images.githubusercontent.com/84130229/151123289-50287abb-93b6-4bb0-8e1f-6f5b0fd95473.png)

#### 오늘의 시뮬결과 이미지
![image](https://user-images.githubusercontent.com/84130229/151176984-db147671-5043-4c00-8f6a-00bcd9ffbad7.png)

0124의 문제점을 해결한 부분 ) lambda sequence를 짤때 그냥 seq function을 사용하는게 아니라 lseq를 사용했어야함!! log-spaced lambda를 얻는것이 필요했으므로

새로운 문제가 발견된 부분 ) 현재 목표가 논문의 시뮬레이션과 같은 결과를 도출하는 것이기때문에 초기값이 0124처럼 고른 분포를 보여야하며,람다값에는 차이가 없는데 초기값이 다르게 나왔다는 것은 다른 파라미터의 문제인지 유의할 필요가 있음


### 0127 - 시뮬결과 이미지
![image](https://user-images.githubusercontent.com/84130229/151315965-a05871b2-7a01-48e0-8220-7962fd4f33bd.png)

논문에서의 결과와 거의 유사(tuning parameter 값이 있기 때문에 완벽히 같은 결과는 불가를 예상..) !   
PHI : 1번째 iteration의 각 weight values' range를 결정 - 상대위험회피계수 높을 수록 위험을 기피한다는 뜻이므로, 너무 큰 값의 매수나 매도가 발생하지 않음..?   
ETA : 수렴하는 속도를 조절 : eta값이 클수록 더 큰 람다가 주어져야 그룹화 되는 양상ㅇ   
오늘 목표 :    
1) 실제 데이터에는 eta = 1.5로 고정하고, phi value를 계산할 수 있는 수식이 있을지 더 develope
2) 실제로 섹터가 나뉘어있던 데이터를 기반으로 해보기,,(SGLasso 논문 참고)


#### 0129 
https://bioinformaticsandme.tistory.com/145   
연구보고할때는 예시로 그룹화가 어떻게 이뤄지는지 특정람다에서 보여드리고자함   
이때 hierchical clustering으로 거리에 따라 그룹화하고 적당한 개수에서 나뉘는것을 보여드리면 될듯함

+) 추가로 0127 파일을 다시 돌린 결과 논문상의 시뮬레이션과 더 유사한 그래프를 얻음!  
![image](https://user-images.githubusercontent.com/84130229/151700679-ba5f806c-ea9a-48ea-9edb-9ccb7149239b.png).   
시뮬레이션 데이터는 동일, phi : 40 , eta : 10 / 수직 선은 long only 구역을 말함

#### 0214
threshold : 0.005%로 언급하고 있어 5번째 자리까지 사용(반올림)   
연구보고용은 phi = 40, eta = 15로 고정 -> 추후에 교수님께 조언 구하고자함   
mission : 논문 다시 읽어보고 빠진 부분있는지 catch up하기!

#### 0304
what's going on now : mean-shortfall + SLOPE combining
writing optimization problem is my mission,,

#### 0317
I solved the optimization problem applying ADMM Algorithm on mean-shortfall + SLOPE
But still there are so many stuffs to be solved
e.g. grouping effect is not that obvious than my expection,,
I'm gonna recalculate the stopping rule and upgrade my algorithm's converging speed!


#### 0320
result from applying simulation data on the mean-shortfall + SLOPE !
![image](https://user-images.githubusercontent.com/84130229/159166378-8d56c676-4b54-4cb9-a379-eff45c5a8dcf.png)

#### 0506
Recently, working on applying SLOPE on DOW30 real data.  
Prameter setting is the difficulitiest part.....!
