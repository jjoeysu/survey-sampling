######均值估计1####################################################################
cs.mean1=function(N=NULL,n,M,y.sample,cluster,alpha)
{
  f=ifelse(is.null(N), 0, n/N)
  
  nf=(1-f)/n
  
  yi.bar=rep(0,n)
  
  for(i in 1:n)
  {
    yi.bar[i]=mean(y.sample[cluster==i])
  }
  
  mean.est=sum(y.sample)/(n*M)
  sb2=M/(n-1)*sum((yi.bar-mean.est)^2)
  mean.var=nf/M*sb2
  mean.sd=sqrt(mean.var)
  
  ci.result=ci(mean.est,mean.sd,alpha)
  mean.left =ci.result$left
  mean.right=ci.result$right
  
  mean.result=matrix(c(mean.est, mean.var, mean.sd, mean.left, mean.right), nrow=1)
  colnames(mean.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(mean.result)="Mean_cs"
  return(mean.result=as.data.frame(mean.result))
}

######均值估计2#######################################################################
cs.mean2=function(N=NULL,n,M,yi.bar,sw,alpha)
{
  f=ifelse(is.null(N), 0, n/N)
  
  nf=(1-f)/n
  
  mean.est=sum(yi.bar)/n
  sb2=M/(n-1)*sum((yi.bar-mean.est)^2)
  mean.var=nf/M*sb2
  mean.sd=sqrt(mean.var)
  
  ci.result=ci(mean.est,mean.sd,alpha)
  mean.left =ci.result$left
  mean.right=ci.result$right
  
  mean.result=matrix(c(mean.est, mean.var, mean.sd, mean.left, mean.right), nrow=1)
  colnames(mean.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(mean.result)="Mean_cs"
  return(mean.result=as.data.frame(mean.result))
}

######总量估计1###################################################################
cs.total1=function(N,n,M,y.sample,cluster,alpha)
{
  cs.mean.result=cs.mean1(N,n,M,y.sample,cluster,alpha)
  ybar.cs.est  =cs.mean.result$Est
  ybar.cs.var  =cs.mean.result$Var
  ybar.cs.left =cs.mean.result$Left
  ybar.cs.right=cs.mean.result$Right
  
  ytot.cs.est  =N*ybar.cs.est              
  ytot.cs.var  =N^2*ybar.cs.var
  ytot.cs.sd   =sqrt(ytot.cs.var)
  ytot.cs.left =N*ybar.cs.left
  ytot.cs.right=N*ybar.cs.right
  
  total.result=matrix(c(ytot.cs.est, ytot.cs.var, ytot.cs.sd, ytot.cs.left, ytot.cs.right), nrow=1)
  colnames(total.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(total.result)="Total_cs"
  return(total.result=as.data.frame(total.result))
}

######总量估计2#################################################################
cs.total2=function(N=NULL,n,M,yi.bar,sw,alpha)
{
  cs.mean.result=cs.mean2(N=NULL,n,M,yi.bar,sw,alpha)
  ybar.cs.est  =cs.mean.result$Est
  ybar.cs.var  =cs.mean.result$Var
  ybar.cs.left =cs.mean.result$Left
  ybar.cs.right=cs.mean.result$Right
  
  ytot.cs.est  =N*ybar.cs.est              
  ytot.cs.var  =N^2*ybar.cs.var
  ytot.cs.sd   =sqrt(ytot.cs.var)
  ytot.cs.left =N*ybar.cs.left
  ytot.cs.right=N*ybar.cs.right
  
  total.result=matrix(c(ytot.cs.est, ytot.cs.var, ytot.cs.sd, ytot.cs.left, ytot.cs.right), nrow=1)
  colnames(total.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(total.result)="Total_cs"
  return(total.result=as.data.frame(total.result))
}

######设计效应估计1######################################################################
cs.deff1=function(N=NULL,n,M,y.sample,cluster,alpha)
{
  f=ifelse(is.null(N), 0, n/N)
  
  nf=(1-f)/n
  
  yi.bar=rep(0,n)
  y.bar=rep(0,n*M)
  
  for(i in 1:n)
  {
    yi.bar[i]=mean(y.sample[cluster==i])
    y.bar[((i-1)*M+1):(i*M)]=yi.bar[i]
  }
  
  mean.est=sum(y.sample)/(n*M)
  sb2=M/(n-1)*sum((yi.bar-mean.est)^2)
  sw2=1/(n*(M-1))*sum((y.sample-y.bar)^2)
  
  if(is.null(N))
  {
    s2=1/M*(sb2+(M-1)*sw2)
    
    rho.c=(sb2-sw2)/(sb2+(M-1)*sw2)
    mean.var=nf/M*s2*(1+(M-1)*rho.c)
  }
  else
  {
    s2=1/(N*M-1)*((N-1)*sb2+N*(M-1)*sw2)
    
    rho.c=(M*(N-1)*sb2-(N*M-1)*s2)/((M-1)*(N*M-1)*s2)
    mean.var=nf*(N*M-1)/(M^2*(N-1))*s2*(1+(M-1)*rho.c)
  }
  mean.var.ran=nf/M*s2
  deff=mean.var/mean.var.ran
  
  deff.result=cbind(sb2,sw2,s2,mean.var,mean.var.ran,rho.c,deff)
  
  colnames(deff.result)=c("sb2","sw2","s2","Var","Var.ran","rho.c", "Deff")
  rownames(deff.result)="Deff_cs"
  return(deff.result=as.data.frame(deff.result))
}

######设计效应估计2#####################################################################
cs.deff2=function(N=NULL,n,M,yi.bar,sw,alpha)
{
  f=ifelse(is.null(N), 0, n/N)
  
  nf=(1-f)/n
  
  y.bar=rep(0,n*M)
  
  for(i in 1:n)
  {
    y.bar[((i-1)*M+1):(i*M)]=yi.bar[i]
  }
  
  mean.est=sum(yi.bar)/n
  sb2=M/(n-1)*sum((yi.bar-mean.est)^2)
  sw2=1/n*sum(sw^2)
  
  if(is.null(N))
  {
    s2=1/M*(sb2+(M-1)*sw2)
    
    rho.c=(sb2-sw2)/(sb2+(M-1)*sw2)
    mean.var=nf/M*s2*(1+(M-1)*rho.c)
  }
  else
  {
    s2=1/(N*M-1)*((N-1)*sb2+N*(M-1)*sw2)
    
    rho.c=(M*(N-1)*sb2-(N*M-1)*s2)/((M-1)*(N*M-1)*s2)
    mean.var=nf*(N*M-1)/(M^2*(N-1))*s2*(1+(M-1)*rho.c)
  }
  mean.var.ran=nf/M*s2
  deff=mean.var/mean.var.ran
  deff.result=cbind(sb2,sw2,s2,mean.var,mean.var.ran,rho.c,deff)
  
  colnames(deff.result)=c("sb2","sw2","s2","Var","Var.ran","rho.c", "Deff")
  rownames(deff.result)="Deff_cs"
  return(deff.result=as.data.frame(deff.result))
}
######总体比例估计########################################################################
cs.prop=function(N=NULL,n,M=NULL,m.sample=NULL,a.sample,p.sample=NULL,alpha)
{
  f=ifelse(is.null(N), 0, n/N)
  
  nf=(1-f)/n
  
  if(is.null(M))
  {
    M=1/n*sum(m.sample)
    p=sum(a.sample)/sum(m.sample)
    p.var=nf/(M^2)*1/(n-1)*(sum(a.sample^2)+p^2*sum(m.sample^2)-2*p*sum(a.sample*m.sample))
    p.sd=sqrt(p.var)
    
    srs.prop.re=srs.prop(NULL, n*M, sum(a.sample), alpha)
    p.srs.var=srs.prop.re$p.var
  }
  else
  {
    p=1/(n*M)*sum(a.sample)
    if(is.null(p.sample))
    {p.sample=a.sample/M}
    p.var=nf/(n-1)*sum((p.sample-p)^2)
    p.sd=sqrt(p.var)
    if(is.null(N))
    {
      srs.prop.re=srs.prop(NULL, n*M, sum(a.sample), alpha)
    }
    else
    {
      srs.prop.re=srs.prop(N*M, n*M, sum(a.sample), alpha)
    }
    p.srs.var=srs.prop.re$p.var
  }
  
  deff=p.var/p.srs.var
  
  ci=ci(p,p.sd,alpha)
  left=ci$left
  right=ci$right
  
  prop.result=cbind(p,p.var,p.sd,left,right,p.srs.var,deff)
  
  colnames(prop.result)=c("p","Var","SD","Left","Right","Var.srs","Deff")
  rownames(prop.result)="Prop_cs"
  return(prop.result=as.data.frame(prop.result))
}
