#######srs均值估计#####################################################################
srs.mean=function(N=NULL, mysample, alpha)
{
  n=length(mysample)
  f=ifelse(is.null(N), 0, n/N)
  
  ybar=mean(mysample)
  ys2 =var(mysample)
  
  ybar.var=((1-f)/n)*ys2
  ybar.sd =sqrt(ybar.var)
  
  ci.result=ci(ybar, ybar.sd, alpha)
  d    =ci.result$d
  r    =ci.result$r
  left =ci.result$left
  right=ci.result$right
  
  return(list(ybar=ybar, ybar.var=ybar.var, ybar.sd=ybar.sd, 
              d=d, r=r, left=left, right=right))
}

#######srs总量估计####################################################################

srs.total=function(N, mysample, alpha)
{
  n=length(mysample)
  f=n/N
  
  ybar=mean(mysample)
  ys2 =var(mysample)
  
  ytot.est=N*ybar
  ytot.var=N^2*((1-f)/n)*ys2
  ytot.sd =sqrt(ytot.var)
  
  ci.result=ci(ytot.est, ytot.sd, alpha)
  d    =ci.result$d
  r    =ci.result$r
  left =ci.result$left
  right=ci.result$right
  
  return(list(ytot.est=ytot.est, ytot.var=ytot.var, ytot.sd=ytot.sd, 
              d=d, r=r, left=left, right=right))
}

########srs总体比例估计#####################################################################

srs.prop=function(N=NULL, n, event.num, alpha)
{
  f=ifelse(is.null(N), 0, n/N)
  
  p.est=event.num/n
  p.var=((1-f)/(n-1))*p.est*(1-p.est)
  p.sd =sqrt(p.var)

  ci.result=ci(p.est, p.sd, alpha)
  d    =ci.result$d
  r    =ci.result$r
  left =ci.result$left
  right=ci.result$right
  
  return(list(p.est=p.est, p.var=p.var, p.sd=p.sd, 
              d=d, r=r, left=left, right=right))
}

#######srs总体特征指标数估计###################################################################
srs.num=function(N=NULL, n, event.num, alpha)
{
  f=ifelse(is.null(N), 0, n/N)
  
  p.est=event.num/n
  p.var=((1-f)/(n-1))*p.est*(1-p.est)
    
  A.est=N*p.est
  A.var=N^2*p.var
  A.sd =sqrt(A.var)
  
  ci.result=ci(A.est, A.sd, alpha)
  d    =ci.result$d
  r    =ci.result$r
  left =ci.result$left
  right=ci.result$right
  
  return(list(A.est=round(A.est), A.var=A.var, A.sd=A.sd, 
              d=d, r=r, left=round(left), right=round(right)))
}

#######估计总体均值时样本量的确定####################################################################
srs.mean.size=function(N=NULL, Mean.his=NULL, Var.his, method, bound, alpha)
{
  quan=qnorm(1-alpha/2)
  
  if (method=="V") {n0=Var.his/bound}
  if (method=="CV"){n0=Var.his/(bound*Mean.his)^2}
  if (method=="d") {n0=(quan/bound)^2*Var.his}
  if (method=="r") {n0=(quan/(bound*Mean.his))^2*Var.his}
  
  size=ifelse(is.null(N), n0, n0/(1+n0/N))
  return(list(method=method, n0=ceiling(n0), size=ceiling(size)))
}

########估计总体比例时样本量的确定###################################################################
srs.prop.size=function(N=NULL, Prop.his, method, bound, alpha)
{
  quan=qnorm(1-alpha/2)
  
  if (method=="V") {n0=Prop.his*(1-Prop.his)/bound}
  if (method=="CV"){n0=(1-Prop.his)/(bound^2*Prop.his)}
  if (method=="d") {n0=(quan^2)*Prop.his*(1-Prop.his)/(bound^2)}
  if (method=="r") {n0=(quan^2)*(1-Prop.his)/(bound^2*Prop.his)}
  
  size=ifelse(is.null(N), n0, n0/(1+(n0-1)/N))
  return(list(method=method, n0=ceiling(n0), size=ceiling(size)))
}
#######n0——修正前样本量，size——修正后样本量##################################

