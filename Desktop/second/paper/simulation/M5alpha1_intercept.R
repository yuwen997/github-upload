#' @param 

# 三个alpha:alpha1 c(1,0,0,0): randomaziton; alpha2: 系数小一点,就是现在这个； 第三个不除以4； 
#* linear(1: beta相同 2: beta不同) 3: nonlinear;
# fit 1: 错的propensity regression 2: right propensity score 3*3*2种 算relatively bias和mse 
#n=5000 
# time=1000
library(Matching)
n=5000
p=4
M=5
time=1000
simulation<- function(M,time,type,beta_1,beta_2,alpha,pro){
  error1=rep(0,time)
  error4=rep(0,time)
  pb <- txtProgressBar(min = 0, max = time, style = 3)
  for(i in 1:time){
    x=matrix(runif(p*n,1-sqrt(3),1+sqrt(3)),n,p,T)
    colnames(x)<- paste0("x",1:p)
    x<- cbind(rep(1,n),x)
    e1<- rnorm(n,0,1)
    e2<- rnorm(n,0,2)
    if(type=='linear'){
      trueY1<- x%*%beta_1+e1
      trueY2<- x%*%beta_2+e2
    }
    if(type=='non linear'){
      trueY1<- (x^2)%*%beta_1+e1
      trueY2<- (x^2)%*%beta_2+e2
    }
    true=mean(trueY1)-mean(trueY2)
    propen<- 1/(1+exp(-(x%*%alpha)))
    A<- rbinom(n,1,propen)
    B<- 1-A
    Y<- ifelse(A==1,trueY1,trueY2)
    data<- data.frame(x,A,Y)
    X<- data.matrix(data)
    if(pro=='right'){
      prophat<- glm(A~1,family=binomial('logit'),data=data,control = list(maxit = 50))
      phat<- predict.glm(prophat) 
    }
    if(pro=='wrong'){
      prophat<- glm(A~x1+x2+x3,family=binomial('logit'),data=data,control = list(maxit = 50))
      phat<- predict.glm(prophat)
    }
    
    out1 <- Matching::Match(Y = Y, Tr = A, X =phat ,M=M,
                            distance.tolerance = 0, ties = TRUE, Weight = 2)
    out<- Match(Y=Y,Tr=B,X =phat ,M=M,
                distance.tolerance = 0, ties = TRUE, Weight = 2)
    # 对于A=0的部分的match
    ymc0<- Y[out$index.treated]#control自己
    ymt0<- Y[out$index.control]#control的match
    xmt0<- x[out$index.control,]
    xmc0<- x[out$index.treated,]
    # 对于A=1的部分的match
    xmc<- x[out1$index.control,]
    ymc<- Y[out1$index.control]#treat的match
    xmt<- x[out1$index.treated,]
    ymt<- Y[out1$index.treated]#treat自己
    
    
    y0<- c(ymc0,ymc)
    y1<- c(ymt0,ymt)
    x0<- rbind(xmc0,xmc)
    x1<- rbind(xmt0,xmt)
    trt<- c(rep(1,length(y1)),rep(0,length(y0)))
    
    #4
    newdata<- data.frame(c(y1,y0),trt,rbind(x1,x0))
    names(newdata)[1:2]<- c('Y',"A")
    f<- lm(Y~., data=newdata)
    error4[i]=coef(f)[2]-true
    #1
    error1[i]=mean(y1)-mean(y0)-true
    
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  return(list(error1=error1,error4=error4))
}
mse<- function(data){
  return((mean(data))^2+var(data))
}

alpha1<- matrix(c(1,0,0,0,0),p+1,1)
alpha2<-  matrix(c(1,1,1,-1,-1)/4,p+1,1)
alpha3<- matrix(c(1,1,1,-1,-1),p+1,1)

beta_1<- matrix(c(1,2,1,-1,-1)/2,p+1,1)
beta_21<- matrix(c(1,2,1,-1,-1)/2,p+1,1)
beta_22<- matrix(c(1,3.1,4.2,1,1)/2,p+1,1)

#------
case111<- simulation(M,time,'linear',beta_1,beta_21,alpha1,'wrong')
mean(case111$error1);var(case111$error1);
mean(case111$error4);var(case111$error4);
mse(case111$error1)
mse(case111$error4)
(var(case111$error1)-var(case111$error4))/var(case111$error1)
#alpha1 linear1 right propensity score
case112<- simulation(M,time,'linear',beta_1,beta_21,alpha1,'right')
mean(case112$error1);var(case112$error1);
mean(case112$error4);var(case112$error4);
mse(case112$error1)
mse(case112$error4)
(var(case112$error1)-var(case112$error4))/var(case112$error1)
#alpha1 linear2 wrong propensity score
case121<- simulation(M,time,'linear',beta_1,beta_22,alpha1,'wrong')
mean(case121$error1);var(case121$error1);
mean(case121$error4);var(case121$error4);
mse(case121$error1)
mse(case121$error4)
(var(case121$error1)-var(case121$error4))/var(case121$error1)
#alpha1 linear2 right propensity score
case122<- simulation(M,time,'linear',beta_1,beta_22,alpha1,'right')
mean(case122$error1);var(case122$error1);
mean(case122$error4);var(case122$error4);
mse(case122$error1)
mse(case122$error4)
(var(case122$error1)-var(case122$error4))/var(case122$error1)
# alpha1 non linear wrong propensity score
case131<- simulation(M,time,'non linear',beta_1,beta_21,alpha1,'wrong')
mean(case131$error1);var(case131$error1);
mean(case131$error4);var(case131$error4);
mse(case131$error1)
mse(case131$error4)
(var(case131$error1)-var(case131$error4))/var(case131$error1)
# alpha1 non linear right propensity score
case132<- simulation(M,time,'non linear',beta_1,beta_21,alpha1,'right')
mean(case132$error1);var(case132$error1);
mean(case132$error4);var(case132$error4);
mse(case132$error1)
mse(case132$error4)
(var(case132$error1)-var(case132$error4))/var(case132$error1)

save(case111,case112,case121,case122,case131,case132,file='M5alpha1_intercept.RData')











