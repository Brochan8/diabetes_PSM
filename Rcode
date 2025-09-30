# 두 집단의 데이터 불균형 해결: 성향점수분석
# Propensity score matching ----
# 2024-08-12 By Jinseo Kim 

# I. 데이터 확인 ----

## 1. R 버전 확인하기 ----
version$version.string    # R version 4.4.1 (2024-06-14 ucrt)
RStudio.Version()$version # 2024.4.2.764


## 2. 패키지 불러오기 ----
# install.packages("dplyr")
library(dplyr) # pipe  

## 3. 데이터 불러오기 ----
setwd("C:/Users/LG/Desktop/na/PSM")
getwd()
data <- read.csv("Data04.csv",header = TRUE)

head(data)
tail(data)
dim(data)
View(data)


# II. 분포 확인 ----

## 1. 변수 타입지정 ----

### 1) 연속형 자료 ----
con_var <- c("Age", "CCI")
con <- lapply(data[con_var],  as.numeric) %>% as.data.frame()

### 2) 범주형 자료 ----
var <- names(data)
fac_var = var [!var %in% con_var]
fac <- lapply(data[fac_var], as.factor) %>% as.data.frame()

data2 <- cbind(con, fac)



## 2. 기술통계량 ----

### 1) 연속형 자료 ----
# "Age", "CCI"

# 히스토그램(Histogram)
par(mfrow=c(1,2))

hist(data2$Age, main = "Age", freq = FALSE)

hist(data2$CCI, main = "CCI", freq = FALSE)


# 상자그림(Boxplot) 
boxplot(data2$Age, main = "Age")

boxplot(data2$CCI,main = "CCI")


# 상자그림과 히스토그램 한 화면에 보이도록
par(mfrow=c(2,2))

boxplot(data2$Age,main = "Age")
boxplot(data2$CCI,main = "CCI")
hist(data2$Age, main = "Age")
hist(data2$CCI, main = "CCI")



### 2) 범주형 자료 ----
# 33개 
ncol(data2) # 35 


# III. Characteristics with SMD----

# install.packages("tableone")
library(tableone)

var <- names(data2)

# 1. 그룹 구분 
table1 <- CreateTableOne(var = var,
                         strata = "Treatment",
                         data = data2, 
                         factorVars = fac_var,
                         smd = TRUE)


# (1) 연속형변수가 정규분포 따르는 경우: mean, SD
table1_normal <- print(table1, 
                       smd=TRUE)

View(table1_normal)

# 결과표 다운로드 
write.csv(table1_normal,"table1_normal.csv")



# (2) 연속형변수가 정규분포 따르지 않는 경우: median, IQR
table1_nonnormal <- print(table1, 
                          smd=TRUE,
                          nonnormal=c("Age","CCI"))

View(table1_nonnormal)

# 결과표 다운로드 
write.csv(table1_nonnormal,"table1_nonnormal.csv")



# IV. Matching ----

## 1. 4 to 1 matching ----

install.packages("MatchIt")
library(MatchIt)

data2$Treatment %>% table

match <- matchit(Treatment ~ Age + Sex + CCI,
                 data=data2, 
                 caliper=0.2, 
                 ratio=2)
summary(match)


## 2. 매칭 후 분포 확인 ----

### 1) Jitter for PS ----
par(mfrow=c(1,1))

plot(match, type = "jitter")


### 2) Histogram for PS ----
plot(match, type = "hist") # (참고) 코드 수행 장시간 소요


### 3) 매칭 후 데이터셋 ----
# matched data set으로 향후 분석 수행
m_data <- match.data(match)  

head(m_data)
View(m_data)


# 일부 샘플 50개 무작위로 추출
set.seed(0812) 
sampled_data <- m_data[sample(nrow(m_data), 50), ]


# Histogram 중첩
with(sampled_data, 
     hist(match$distance[Treatment=="Molnupiravir"], 
          col="orange",
          main="Propensity score distribution",
          xlab="ps"))

with(sampled_data, 
     hist(match$distance[Treatment=="Nirmatrelvir-Ritonavir"], 
          col="skyblue", 
          add=T))



# V. 매칭 후 SMD ----

var2 <- names(m_data[c(1:35)]) # table에 들어갈 변수만 선택

# 1. 그룹 구분 
m_table1 <- CreateTableOne(var = var2,
                           strata = "Treatment",
                           data = m_data, 
                           factorVars = fac_var,
                           smd = TRUE)


# 연속형변수 정규성 확인
par(mfrow=c(1,2))
hist(m_data$Age, main = "Age")
hist(m_data$CCI, main = "CCI")


# (1) 정규분포 따르는 경우: mean, SD
m_table1_normal <- print(m_table1, 
                         smd=TRUE)
View(m_table1_normal)

# 결과값 저장 
write.csv(m_table1_normal, "m_table1_normal.csv")


## (2) 정규분포 따르지 않는 경우: median, IQR
m_table1_nonnormal <- print(m_table1, 
                          smd=TRUE,
                          nonnormal=c("Age","CCI"))

View(m_table1_nonnormal)

# 결과값 저장 
# write.csv(m_table1_nonnormal, "m_table1_nonnormal.csv")




# 성향점수 이용한 생존분석 예시
m_cox_model <- coxph(Surv(Time, Event) ~ Treatment + cluster(distance), 
                     data = m_data)




# VI. IPTW ----
# 필요한 패키지 로드
library(survey)

# 성향 점수 추정 (로지스틱 회귀 사용)
ps_model <- glm(Treatment ~ Age + Sex + CCI, data = data2, family = binomial)

# 성향 점수 계산
data2$ps <- predict(ps_model, type = "response")

# IPTW 가중치 계산
data2$weight <- ifelse(data2$Treatment == 1, 1/data2$ps, 1/(1 - data2$ps))


# 가중치를 이용한 생존분석 예시
w_cox_model <- coxph(Surv(Time, Event) ~ Treatment, 
                     data = data2, 
                     weights = weight)
