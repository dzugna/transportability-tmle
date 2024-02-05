
rm(list=ls())
# General packages 
pkgs <- c("boot")
vapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  require(pkg, character.only = TRUE, quietly = TRUE)
}, FUN.VALUE = logical(length = 1L))

#TMLE function
eattmle<-function(x,y,site,w,nsitemodel,nxmodel,noutmodel){
	datw<-w
	n.dat<-nrow(w)
	#calculate components of clever covariate: P(S=1|W) and P(Z=1|S=1,W)
	cps<-predict(glm(formula=nsitemodel,data=data.frame(cbind(site=site,datw)),family="binomial"),type="response")
	cpx<-predict(glm(formula=nxmodel,data=data.frame(cbind(x=x,datw)),family="binomial"),type="response")  #probability
	#calculate clever covariate
	g0w<-((1-cpx)*cps)/(1-cps)
	g1w<-(cpx*cps)/(1-cps)
	h0w<-((1-x)*I(site==1))/g0w
	h1w<-(x*I(site==1))/g1w
	ymodel<-glm(formula=noutmodel,family="binomial",data=data.frame(cbind(datw,x=x,site=site, y=y)),subset=site==1)
	data_new0<-data_new1<-datw
	data_new0$x<-0
	data_new1$x<-1
	#initial prediction
	q<-cbind(predict(ymodel,type="link",newdata=data.frame(cbind(datw,x=x))),predict(ymodel,type="link",newdata=data_new0),predict(ymodel,type="link",newdata=data_new1))   #log(odds)
	epsilon<-coef(glm(y~-1+offset(q[,1])+h0w+h1w,family="binomial",subset=site==1))  #epsilon calculated on real data
	#update initial prediction
	q1<-q+c((epsilon[1]*h0w+epsilon[2]*h1w),epsilon[1]/g0w,epsilon[2]/g1w) 
	tmlest<-mean(plogis(q1[,3][site==0]))-mean(plogis(q1[,2][site==0]))  #q1[,2]: z=1 then h1w but mean on [site=0] and then g1w 
	#get efficient influence curve values for everyone
	ps0<-mean(I(site==0))
	eic<-(((x*h1w/ps0)-((1-x)*h0w/ps0))*(y-plogis(q[,1])))+(I(site==0)/ps0*plogis(q1[,3]))-(I(site==0)/ps0*plogis(q1[,2]))-(tmlest/ps0)
	return(list("est"=tmlest,"var"=var(eic)/n.dat,"eic"=eic[site==0])) 
}

# Load datasets
final.data <- readRDS("Data/final.data.rds")

dim(final.data)
names(final.data)
summary(final.data)
head(final.data)
#parity: 0,1,2+
#medu (maternal education): 1=high, 2=medium, 3=low
#maternal age: continuous
#source: 0 (target population), 1 (study population)

#continuous variables centered at their mean value; generation of dummy variables from categorical varriables
final.data$c.mage<-final.data$mage-mean(final.data$mage)
final.data$parity1<-0
final.data$parity1[final.data$parity==1]<-1
final.data$parity2<-0
final.data$parity2[final.data$parity==2]<-1
final.data$medu2<-0
final.data$medu2[final.data$medu==2]<-1
final.data$medu3<-0
final.data$medu3[final.data$medu==3]<-1
S<-final.data$source
W<-data.matrix(as.data.frame(cbind(rep(1,nrow(final.data)),final.data$c.mage,final.data$parity1,final.data$parity2,final.data$medu2,final.data$medu3),rownames=T))

#exposure's simulation 
delta <- c(0,0.1, -0.3, -0.5, 0,0)  #Z~B(pigreco), pigreco=expit(W*delta)
set.seed(7773)
X<- rbinom(nrow(final.data), 1, plogis(W %*% delta))

#outcome's simulation
Z<-data.matrix(as.data.frame(cbind(rep(1,nrow(final.data)),X,final.data$c.mage,final.data$parity1,final.data$parity2,final.data$medu2,final.data$medu3),rownames=T))
beta <- c(0, 0.6, 0.2, -0.1,-0.3,0.4,0.6)
set.seed(4321)
Y<-rbinom(nrow(final.data), 1, plogis(Z %*% beta+0.3*X*Z[,6]+0.5*X*Z[,7]))

prop.table(table(Y,X),2)
prop.table(table(Y,S),2)
prop.table(table(X,S),2)

dat<-data.frame(y=Y,x=X,w1=final.data$c.mage,w2=final.data$parity1,w3=final.data$parity2,w4=final.data$medu2,w5=final.data$medu3,s=S)
wmat<-data.frame(w1=final.data$c.mage,w2=final.data$parity1,w3=final.data$parity2,w4=final.data$medu2,w5=final.data$medu3)

sitemodel<-"site ~ w1 + w2 + w3 + w4 + w5" 
xmodel<-"x~w1+w2+w3" 
outmodel<-"y ~ x + w1 + w2 + w3 + w4 + w5 + x:w4 + x:w5" 

#ATE in the study and target populations (confidence intervals calculated by bootstrap)
library(boot)
boot_ci= function(mydata,indices) {
  d=mydata[indices,]
  d0<-d1<-d
  d0$x<-0
  d1$x<-1
  model.y<-glm(formula=outmodel,data=data.frame(d),family="binomial")
  pred<-cbind(predict(model.y,type="link",newdata=data.frame(d0)),
  predict(model.y,type="link",newdata=data.frame(d1)))
  ate0<-	mean(plogis(pred[,2][d$s==0]))-mean(plogis(pred[,1][d$s==0]))
  ate1<-	mean(plogis(pred[,2][d$s==1]))-mean(plogis(pred[,1][d$s==1]))
  return(c(ate0,ate1))
  }
               
replicates=100
results <- boot(data=dat, statistic=boot_ci, R=replicates)
ate0<-results$t0[1]
ate1<-results$t0[2]
conf0<-boot.ci(results, type="perc",index=1)
ate0 #ate in the target population
conf0 #confidence interval
conf1<-boot.ci(results, type="perc",index=2)
ate1 #ate in the study population
conf1 #confidence interval

#tmle to transport estimates from the study population to the target population
tmle<-eattmle(x=dat$x, y =dat$y, site=dat$s, w =wmat, nsitemodel=sitemodel , nxmodel=xmodel , noutmodel=outmodel)
est.tmle<-eattmle(x=dat$x, y =dat$y, site=dat$s, w =wmat, nsitemodel=sitemodel , nxmodel=xmodel , noutmodel=outmodel)$est 
var.tmle<-tmle$var
est<-c(est.tmle,est.tmle-1.96*sqrt(var.tmle),est.tmle+1.96*sqrt(var.tmle))
est

#Note that the results differ slightly from the results reported in the text because the covariate data are generated by randomly sampling with replacement from the original datasets. 