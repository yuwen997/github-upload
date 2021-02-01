library(Matching)
library(mgcv)
n=5000
m=n*10
p=4
M=5
time=1000
#error1 traditional error2 过去的 error3 new
simulation<- function(M,time,type,beta_1,beta_2,alpha,pro){
  error1=rep(0,time)
  error2=rep(0,time)
  error3=rep(0,time)
  pb <- txtProgressBar(min = 0, max = time, style = 3)
  
  for(i in 1:time){
    # x=matrix(runif(p*n,1-sqrt(3),1+sqrt(3)),n,p,T)
    # colnames(x)<- paste0("x",1:p)
    # x<- cbind(rep(1,n),x)
    # e1<- rnorm(n,0,1)
    # e2<- rnorm(n,0,2)
    x=matrix(runif(p*m,1-sqrt(3),1+sqrt(3)),m,p,T)
    colnames(x)<- paste0("x",1:p)
    x<- cbind(rep(1,m),x)
    e1<- rnorm(m,0,0.1)
    e2<- rnorm(m,0,0.2)
    if(type=='linear'){
      trueY1<- x%*%beta_1+e1
      trueY2<- x%*%beta_2+e2
    }
    if(type=='non linear'){
      trueY1<- (x^2)%*%beta_1+e1
      trueY2<- (x^2)%*%beta_2+e2
    }
    propen<- 1/(1+exp(-(x%*%alpha)))# 在这里加一个intercept x%*%alpha+1 or +c
    A<- rbinom(m,1,propen)
    B<- 1-A
    Y<- ifelse(A==1,trueY1,trueY2)
    if(sum(A)<n/2 | sum(B)<n/2){
      next
    }
    index1<- sample(which(A==1),n/2)
    index2<- sample(which(B==1),n/2)
    x<- x[c(index1,index2),]
    Y<- Y[c(index1,index2)]
    A<- A[c(index1,index2)]
    B<- B[c(index1,index2)]
    # true=mean(trueY1[c(index1,index2)])-mean(trueY2[c(index1,index2)])
    true=sum(beta_1)-sum(beta_2)
    # 改true
    
    data<- data.frame(x,A,Y)
    X<- data.matrix(data)
    if(pro=='right'){
      prophat<- glm(A~x1+x2+x3+x4,family=binomial('logit'),data=data,control = list(maxit = 50))
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
    
    ymc0<- Y[out$index.treated]#control自己
    ymt0<- Y[out$index.control]#control的match
    ymt<- Y[out1$index.treated]#真正的treat
    ymc<- Y[out1$index.control]#treat的match
    if(length(ymc0)!=M*sum(B) | length(ymt) != M*sum(A)){
      next
    }
    
    #-------
    
    trt<- rep(1,M*sum(A))
    col<- rep(0,M*sum(A))
    xmt<- cbind(x[out1$index.treated,],col)
    xmc<- cbind(x[out1$index.control,],col)

 

    col<- rep(0,M*sum(B))
    trt<- rep(1,M*sum(B))
    xmc0<- cbind(x[out$index.control,],trt)
    xmt0<- cbind(x[out$index.treated,],trt)

    Y1<- Y[out1$index.treated]
    data1<- data.frame(Y1,x[out1$index.treated,])
    Y0<- Y[out$index.treated]
    data0<- data.frame(Y0,x[out$index.treated,])

    model1<- gam(Y1~x1+x2+x3+x4,data=data1)
    muhat1i<- predict.gam(object=model1,newdata=data0)
    data1j<- data.frame(Y[out$index.control],x[out$index.control,])
    names(data1j)<- c("Y1","inter",paste0("x",1:p))
    muhat1j<- predict.gam(model1,data1j)

    model0<- gam(Y0~x1+x2+x3+x4,data=data0)
    muhat0i<- predict.gam(object=model0,newdata=data1)
    data0j<- data.frame(Y[out1$index.control],x[out1$index.control,])
    names(data0j)<- c("Y0",'inter',paste0("x",1:p))
    muhat0j<- predict.gam(model0,data0j)

    new_ymt0<- ymt0+muhat1i-muhat1j
    new_ymt0<- matrix(new_ymt0,ncol=M,byrow=T,nrow=length(new_ymt0)/M)
    new_ymt0<- apply(new_ymt0,1,mean)
    new_ymc<- unique(ymc0)
    new_ymt<- unique(ymt)
    new_ymc0<- ymc+muhat0i-muhat0j
    new_ymc0<- matrix(new_ymc0,ncol=M,byrow=T,nrow=length(new_ymc0)/M)
    new_ymc0<- apply(new_ymc0,1,mean)

    y1<- c(new_ymt,new_ymt0)
    y0<- c(new_ymc,new_ymc0)
    x1<- unique(rbind(cbind(x[out1$index.treated,],rep(1,sum(A))),xmt0))
    x0<- unique(rbind(cbind(x[out$index.treated,],rep(0,sum(B))),xmt))


    newdata<- data.frame(c(y1,y0),rbind(x1,x0))
    names(newdata)<- c("Y","inter","x1","x2","x3","x4","A")
    f<- lm(Y~x1+x2+x3+x4+A, data=newdata)
    error3[i]=coef(f)[6]-true


    # error3[i]=coef(f)[6]

    # #------
    # model1<- gam(Y1~x1+x2+x3+x4+A,data=data)
    # newd<- data.frame(rbind(xmt0,xmt,xmc0,xmc))
    # names(newd)<- c("x1","x2","x3","x4","A")
    # muhat1<- predict.gam(model1,newd)
    # # ymt0 1:M*(sum(B)); ymt M*(sum(B))+1: M*n; ymc0 M*n+1: M*n+M*sum(B) ymc M*(n+sum(B))+1: 2*M*n
    # # new_ymt0<- ymt0+muhat1[1:sum(B)]-muhat1[(n+1):(n+sum(B))]
    # new_ymt0<- ymt0+muhat1[1:(M*sum(B))]-muhat1[(M*n+1):(M*n+M*sum(B))]
    # new_ymt0<- matrix(new_ymt0,ncol=M,byrow=T,nrow=length(new_ymt0)/M)
    # new_ymt0<- apply(new_ymt0,1,mean)
    # new_ymc<- unique(ymc0)
    # new_ymt<- unique(ymt)
    # # new_ymc0<- ymc+muhat1[(sum(B)+1):n]-muhat1[(n+sum(B)+1):(2*n)]
    # new_ymc0<- ymc+muhat1[(M*sum(B)+1):(M*n)]-muhat1[(M*n+M*sum(B)+1):(M*2*n)]
    # new_ymc0<- matrix(new_ymc0,ncol=M,byrow=T,nrow=length(new_ymc0)/M)
    # new_ymc0<- apply(new_ymc0,1,mean)
    #
    # y1<- c(new_ymt,new_ymt0)
    # y0<- c(new_ymc,new_ymc0)
    # x1<- unique(rbind(cbind(x[out1$index.treated,],rep(1,sum(A))),xmt0))
    # x0<- unique(rbind(cbind(x[out$index.treated,],rep(0,sum(B))),xmt))
    #
    #
    #
    #
    # newdata<- data.frame(c(y1,y0),rbind(x1,x0))
    # names(newdata)<- c("Y","x1","x2","x3","x4","A")
    # f<- lm(Y~., data=newdata)
    # error3[i]=coef(f)[6]-true
    #

    #-------
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

    #过去的方法
    newdata<- data.frame(c(y1,y0),trt,rbind(x1,x0))
    names(newdata)[1:2]<- c('Y',"A")
    f<- lm(Y~., data=newdata)
    error2[i]=coef(f)[2]-true
    # error2[i]=coef(f)[2]
    #1

    #------
    error1[i]=mean(y1)-mean(y0)-true
    # error1[i]=mean(y1)-mean(y0)
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  return(list(error1=error1,error2=error2,error3=error3))
}

mse<- function(data){
  return((mean(data))^2+var(data))
}



alpha2<-  matrix(c(1,c(1,1,1,-1)/4),p+1,1)
alpha3<- matrix(c(1,c(1,1,1,-1)),p+1,1)

beta_1<- matrix(c(1,2,1,-1,-1)/2,p+1,1)
beta_21<- matrix(c(1,2,1,-1,-1)/2,p+1,1)
beta_22<- matrix(c(1,3.1,4.2,1,1)/2,p+1,1)

# alpha2<-  matrix(c(1,1,1,-1)/4,p,1)
# alpha3<- matrix(c(1,1,1,-1),p,1)
# #
# beta_1<- matrix(c(2,1,-1,-1)/2,p,1)
# beta_21<- matrix(c(2,1,-1,-1)/2,p,1)
# beta_22<- matrix(c(3.1,4.2,1,1)/2,p,1)


#alpha2 linear1 wrong propensity score
case211<- simulation(M,time,'linear',beta_1,beta_21,alpha2,'wrong')
mean(case211$error1);var(case211$error1);
mean(case211$error2);var(case211$error2);
mean(case211$error3);var(case211$error3);
mse(case211$error1)
mse(case211$error2)
mse(case211$error3)
(var(case211$error1)-var(case211$error2))/var(case211$error1)
(var(case211$error1)-var(case211$error3))/var(case211$error1)

#alpha2 linear1 right propensity score
case212<- simulation(M,time,'linear',beta_1,beta_21,alpha2,'right')
mean(case212$error1);var(case212$error1);
mean(case212$error2);var(case212$error2);
mean(case212$error3);var(case212$error3);
mse(case212$error1)
mse(case212$error2)
mse(case212$error3)
(var(case212$error1)-var(case212$error2))/var(case212$error1)
(var(case212$error1)-var(case212$error3))/var(case212$error1)

#alpha2 linear2 wrong propensity score
case221<- simulation(M,time,'linear',beta_1,beta_22,alpha2,'wrong')
mean(case221$error1);var(case221$error1);
mean(case221$error2);var(case221$error2);
mean(case221$error3);var(case221$error3);
mse(case221$error1)
mse(case221$error2)
mse(case221$error3)
(var(case221$error1)-var(case221$error2))/var(case221$error1)
(var(case221$error1)-var(case221$error3))/var(case221$error1)

#alpha2 linear2 right propensity score
case222<- simulation(M,time,'linear',beta_1,beta_22,alpha2,'right')
mean(case222$error1);var(case222$error1);
mean(case222$error2);var(case222$error2);
mean(case222$error3);var(case222$error3);
mse(case222$error1)
mse(case222$error2)
mse(case222$error3)
(var(case222$error1)-var(case222$error2))/var(case222$error1)
(var(case222$error1)-var(case222$error3))/var(case222$error1)

# alpha2 non linear wrong propensity score
case231<- simulation(M,time,'non linear',beta_1,beta_21,alpha2,'wrong')
mean(case231$error1);var(case231$error1);
mean(case231$error2);var(case231$error2);
mean(case231$error3);var(case231$error3);
mse(case231$error1)
mse(case231$error2)
mse(case231$error3)
(var(case231$error1)-var(case231$error2))/var(case231$error1)
(var(case231$error1)-var(case231$error3))/var(case231$error1)

# alpha2 non linear right propensity score
case232<- simulation(M,time,'non linear',beta_1,beta_21,alpha2,'right')
mean(case232$error1);var(case232$error1);
mean(case232$error2);var(case232$error2);
mean(case232$error3);var(case232$error3);
mse(case232$error1)
mse(case232$error2)
mse(case232$error3)
(var(case232$error1)-var(case232$error2))/var(case232$error1)
(var(case232$error1)-var(case232$error3))/var(case232$error1)

save(case211,case212,case221,case222,case231,case232,file='M5alpha2_var.RData')












