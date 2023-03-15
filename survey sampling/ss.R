#########分层随机抽样均值估计1####################################################################
ss.mean1=function(Nh, nh, yh, s2h, alpha)
{
  stra.num=length(Nh)
  Wh=Nh/sum(Nh)
  fh=nh/Nh
  
  yh.est  =rep(0, stra.num)
  yh.var  =rep(0, stra.num)
  yh.sd   =rep(0, stra.num)
  yh.left =rep(0, stra.num)
  yh.right=rep(0, stra.num)
  
  for (h in 1:stra.num)
  {
    yh.est[h]=yh[h]
    yh.var[h]=((1-fh[h])/nh[h])*s2h[h]
    yh.sd[h]=sqrt(yh.var[h])
    
    ci.stra=ci(yh.est[h], yh.sd[h], alpha)
    yh.left[h] =ci.stra$left
    yh.right[h]=ci.stra$right
  }
  
  stra.result=cbind(Nh, nh, Wh, yh.est, yh.var, yh.sd, yh.left, yh.right)
  
  mean.est=sum(Wh*yh.est)
  mean.var=sum(Wh^2*yh.var)
  mean.sd =sqrt(mean.var)
  
  ci.result=ci(mean.est, mean.sd, alpha)
  mean.left =ci.result$left
  mean.right=ci.result$right
  
  mean.result=matrix(c(mean.est, mean.var, mean.sd, mean.left, mean.right), nrow=1)
  colnames(mean.result)=c("mean.est", "mean.var", "mean.sd", "mean.left", "mean.right")
  rownames(mean.result)="SS"
  return(list(stra.result=as.data.frame(stra.result), mean.result=as.data.frame(mean.result)))
}

#########分层随机抽样均值估计2######################################################################
ss.mean2=function(Nh, mysample, stra.index, alpha)
{
  stra.num=length(Nh)
  Wh=Nh/sum(Nh)
  
  yh.est  =rep(0, stra.num)
  yh.var  =rep(0, stra.num)
  yh.sd   =rep(0, stra.num)
  yh.left =rep(0, stra.num)
  yh.right=rep(0, stra.num)
  
  for (h in 1:stra.num)
  {
    sample.hth=mysample[stra.index==h]
    stra.result=srs.mean(Nh[h], sample.hth, alpha)
    
    yh.est[h]  =stra.result$ybar
    yh.var[h]  =stra.result$ybar.var
    yh.sd[h]   =stra.result$ybar.sd
    yh.left[h] =stra.result$left
    yh.right[h]=stra.result$right
  }
  
  stra.result=cbind(Nh, Wh, yh.est, yh.var, yh.sd, yh.left, yh.right)
  
  mean.est=sum(Wh*yh.est)
  mean.var=sum(Wh^2*yh.var)
  mean.sd =sqrt(mean.var)
  
  ci.result=ci(mean.est, mean.sd, alpha)
  mean.left =ci.result$left
  mean.right=ci.result$right
  
  mean.result=matrix(c(mean.est, mean.var, mean.sd, mean.left, mean.right), nrow=1)
  colnames(mean.result)=c("mean.est", "mean.var", "mean.sd", "mean.left", "mean.right")
  rownames(mean.result)="SS"
  return(list(stra.result=as.data.frame(stra.result), mean.result=as.data.frame(mean.result)))
}

########分层随机抽样总体比例估计1#################################################################
ss.prop1=function(Nh, nh, ah, alpha)
{
  stra.num=length(Nh)
  Wh=Nh/sum(Nh)
  fh=nh/Nh
  
  ph.est  =rep(0, stra.num)
  ph.var  =rep(0, stra.num)
  ph.sd   =rep(0, stra.num)
  ph.left =rep(0, stra.num)
  ph.right=rep(0, stra.num)
  
  for (h in 1:stra.num)
  {
    ph.est[h]=ah[h]/nh[h]
    ph.var[h]=((1-fh[h])/(nh[h]-1))*ph.est[h]*(1-ph.est[h])
    ph.sd[h] =sqrt(ph.var[h])
    
    ci.stra=ci(ph.est[h], ph.sd[h], alpha)
    ph.left[h] =ci.stra$left
    ph.right[h]=ci.stra$right
  }
  
  stra.result=cbind(Nh, nh, ah, Wh, ph.est, ph.var, ph.sd, ph.left, ph.right)
  
  prop.est=sum(Wh*ph.est)
  prop.var=sum(Wh^2*ph.var)
  prop.sd =sqrt(prop.var)
  
  ci.result=ci(prop.est, prop.sd, alpha)
  prop.left =ci.result$left
  prop.right=ci.result$right
  
  prop.result=matrix(c(prop.est, prop.var, prop.sd, prop.left, prop.right), nrow=1)
  colnames(prop.result)=c("prop.est", "prop.var", "prop.sd", "prop.left", "prop.right")
  rownames(prop.result)="Prop"
  return(list(stra.result=as.data.frame(stra.result), prop.result=as.data.frame(prop.result)))
}

########分层随机抽样总体比例估计2####################################################################
ss.prop2=function(Nh, nh, ah, alpha)
{
  stra.num=length(Nh)
  Wh=Nh/sum(Nh)
  
  ph.est  =rep(0, stra.num)
  ph.var  =rep(0, stra.num)
  ph.sd   =rep(0, stra.num)
  ph.left =rep(0, stra.num)
  ph.right=rep(0, stra.num)
  
  for (h in 1:stra.num)
  {
    stra.result=srs.prop(Nh[h], nh[h], ah[h], alpha)
    ph.est[h]  =stra.result$p.est
    ph.var[h]  =stra.result$p.var
    ph.sd[h]   =stra.result$p.sd
    ph.left[h] =stra.result$left
    ph.right[h]=stra.result$right
  }
  stra.result=cbind(Nh, nh, ah, Wh, ph.est, ph.var, ph.sd, ph.left, ph.right)
  
  prop.est=sum(Wh*ph.est)
  prop.var=sum(Wh^2*ph.var)
  prop.sd =sqrt(prop.var)
  
  ci.result=ci(prop.est, prop.sd, alpha)
  prop.left =ci.result$left
  prop.right=ci.result$right
  
  prop.result=matrix(c(prop.est, prop.var, prop.sd, prop.left, prop.right), nrow=1)
  colnames(prop.result)=c("prop.est", "prop.var", "prop.sd", "prop.left", "prop.right")
  rownames(prop.result)="Prop"
  return(list(stra.result=as.data.frame(stra.result), prop.result=as.data.frame(prop.result)))
}

######各层样本量分配，分配后的比例######################################################################
ss.weight=function(Wh, S2h, Ch=NULL, allocation)
{
  if (allocation=="Prop") 
  {
    wh=Wh
  }
  if (allocation=="Opt") 
  {
    wh=(Wh*sqrt(S2h)/sqrt(Ch))/sum(Wh*sqrt(S2h)/sqrt(Ch))
  }
  if (allocation=="Neyman") 
  {
    wh=(Wh*sqrt(S2h))/sum(Wh*sqrt(S2h))
  }
  return(wh)
}

######各层样本量分配，分配后的样本量################################################################
ss.size=function(n, Wh, S2h, Ch=NULL, allocation)
{
  wh=ss.weight(Wh, S2h, Ch, allocation)
  nh=wh*n 
  return(list(n=n, allocation=allocation, wh=wh, nh=round(nh)))
}

#######精度要求转换为方差V#################################################################
VCdr=function(method, bound, Ybar=NULL, alpha=NULL)
{
  if (method=="V") 
  {
    var.bound=bound
  }
  
  if (method=="C")
  {
    var.bound=(bound*Ybar)^2
  }
  
  if (method=="d") 
  { 
    quan=qnorm(1-alpha/2)
    var.bound=(bound/quan)^2
  }
  
  if (method=="r") 
  {
    quan=qnorm(1-alpha/2)
    var.bound=(bound*Ybar/quan)^2
  }
  return(var.bound)
}

#######估计总体均值时样本量的确定##################################################################
ss.mean.size=function(Nh, S2h, Ch=NULL, allocation, method, bound, Ybar=NULL, alpha=NULL)
{
  N=sum(Nh)
  Wh=Nh/N
  
  wh=ss.weight(Wh, S2h, Ch, allocation)
  var.bound=VCdr(method, bound, Ybar, alpha)
  
  n=sum(Wh^2*S2h/wh)/(var.bound+sum(Wh*S2h)/N)
  n=ceiling(n)
  nh=floor(wh*n)
  nh=fenpei(wh,nh,n)
  return(list(method=method, bound=bound, allocation=allocation, n=n, nh=nh))  
}

########估计总体比例时样本量的确定############################################################
ss.prop.size=function(Nh, Ph, Ch=NULL, allocation, method, bound, Ybar=NULL, alpha=NULL)
{
  S2h=(Nh/(Nh-1))*Ph*(1-Ph)
  size.result=ss.mean.size(Nh, S2h, Ch, allocation, method, bound, Ybar, alpha)
  return(size.result)  
}

#######样本量分配的细节问题##################################################################
fenpei=function(wh,nh,n)
{
  m=n-sum(nh)
  a=order(wh,decreasing = T)
  for(i in 1:m){nh[a[i]]=nh[a[i]]+1}
  return(nh)
}