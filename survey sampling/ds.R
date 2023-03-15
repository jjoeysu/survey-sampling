######分层的均值简单估计################################################################
stD.mean=function(N=NULL, nh.1st, nh.2nd, ybarh, s2h, alpha)
{
  N.inv=ifelse(is.null(N), 0, 1/N)
  n.1st=sum(nh.1st)
  wh.1st=nh.1st/n.1st
  
  ybar.stD.est=sum(wh.1st*ybarh)
  
  ybar.stD.var1=(1/n.1st-N.inv)*sum(wh.1st*(ybarh-ybar.stD.est)^2)
  ybar.stD.var2=sum((1/nh.2nd-1/nh.1st)*wh.1st^2*s2h)
  ybar.stD.var =ybar.stD.var1+ybar.stD.var2
  
  ybar.stD.sd =sqrt(ybar.stD.var)
  
  ci.result=ci(ybar.stD.est, ybar.stD.sd, alpha)
  ybar.stD.left =ci.result$left
  ybar.stD.right=ci.result$right
  
  mean.result=matrix(c(ybar.stD.est, ybar.stD.var, ybar.stD.sd, ybar.stD.left, ybar.stD.right), nrow=1)
  colnames(mean.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(mean.result)="stD_Mean"
  return(mean.result=as.data.frame(mean.result))
}

######分层的总量简单估计################################################################
stD.total=function(N=NULL, nh.1st, nh.2nd, ybarh, s2h, alpha)
{
  ybar.stD.result=stD.mean(N, nh.1st, nh.2nd, ybarh, s2h, alpha)
  ybar.stD.est  =ybar.stD.result$Est
  ybar.stD.var  =ybar.stD.result$Var
  ybar.stD.left =ybar.stD.result$Left
  ybar.stD.right=ybar.stD.result$Right
  
  ytot.stD.est  =N*ybar.stD.est              
  ytot.stD.var  =N^2*ybar.stD.var
  ytot.stD.sd   =sqrt(ytot.stD.var)
  ytot.stD.left =N*ybar.stD.left
  ytot.stD.right=N*ybar.stD.right
  
  total.result=matrix(c(ytot.stD.est, ytot.stD.var, ytot.stD.sd, ytot.stD.left, ytot.stD.right), nrow=1)
  colnames(total.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(total.result)="stD_Total"
  return(total.result=as.data.frame(total.result))
}

######均值比估计################################################################
RD.mean=function(N=NULL, n.1st, xbar.1st, y.sample, x.sample, alpha)
{ 
  N.inv=ifelse(is.null(N), 0, 1/N)
  n.2nd=length(y.sample)
  
  ybar=mean(y.sample)
  xbar=mean(x.sample)
  
  sy2 =var(y.sample)
  sx2 =var(x.sample)
  syx =cov(y.sample, x.sample)
  
  ratio.est=ybar/xbar
  
  ybar.RD.est=xbar.1st*ratio.est 
  
  ybar.RD.var1=(1/n.1st-N.inv)*sy2
  ybar.RD.var2=(1/n.2nd-1/n.1st)*(sy2+ratio.est^2*sx2-2*ratio.est*syx)
  ybar.RD.var =ybar.RD.var1+ybar.RD.var2
  
  ybar.RD.sd=sqrt(ybar.RD.var)
  
  ci.result=ci(ybar.RD.est, ybar.RD.sd, alpha)
  ybar.RD.left =ci.result$left
  ybar.RD.right=ci.result$right
  
  mean.result=matrix(c(ybar.RD.est, ybar.RD.var, ybar.RD.sd, ybar.RD.left, ybar.RD.right), nrow=1)
  colnames(mean.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(mean.result)="Mean_RD"
  return(mean.result=as.data.frame(mean.result))
}

######总量比估计################################################################
RD.total=function(N=NULL, n.1st, xbar.1st, y.sample, x.sample, alpha)
{
  ybar.RD.result=RD.mean(N, n.1st, xbar.1st, y.sample, x.sample, alpha)
  ybar.RD.est  =ybar.RD.result$Est
  ybar.RD.var  =ybar.RD.result$Var
  ybar.RD.left =ybar.RD.result$Left
  ybar.RD.right=ybar.RD.result$Right
  
  ytot.RD.est  =N*ybar.RD.est              
  ytot.RD.var  =N^2*ybar.RD.var
  ytot.RD.sd   =sqrt(ytot.RD.var)
  ytot.RD.left =N*ybar.RD.left
  ytot.RD.right=N*ybar.RD.right
  
  total.result=matrix(c(ytot.RD.est, ytot.RD.var, ytot.RD.sd, ytot.RD.left, ytot.RD.right), nrow=1)
  colnames(total.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(total.result)="Total_RD"
  return(total.result=as.data.frame(total.result))
}

######均值回归估计################################################################
lrD.mean=function(N=NULL, n.1st, xbar.1st, y.sample, x.sample, alpha, method="Min", beta0=NULL)
{
  N.inv=ifelse(is.null(N), 0, 1/N)
  n.2nd=length(y.sample)
  
  ybar=mean(y.sample)
  xbar=mean(x.sample)
  
  sy2 =var(y.sample)
  sx2 =var(x.sample)
  syx =cov(y.sample, x.sample)
  
  if (method=="Min")
  {
    beta=syx/sx2
    
    ybar.lrD.est=ybar+beta*(xbar.1st-xbar)
    
    ybar.lrD.var1=(1/n.1st-N.inv)*sy2
    ybar.lrD.var2=(1/n.2nd-1/n.1st)*((n.2nd-1)/(n.2nd-2))*(sy2-syx^2/sx2)
    ybar.lrD.var =ybar.lrD.var1+ybar.lrD.var2
    
    ybar.lrD.sd =sqrt(ybar.lrD.var)
    
    ci.result=ci(ybar.lrD.est, ybar.lrD.sd, alpha)
    ybar.lrD.left =ci.result$left
    ybar.lrD.right=ci.result$right
    
    ybar.lrD.result=matrix(c(ybar.lrD.est, ybar.lrD.var, ybar.lrD.sd, ybar.lrD.left, ybar.lrD.right), nrow=1)
    colnames(ybar.lrD.result)=c("Est", "Var", "SD", "Left", "Right")
    rownames(ybar.lrD.result)=c("Mean_lrD")
  }
  
  if (method=="Constant")
  {
    beta=beta0
    
    ybar.lrD.est=ybar+beta*(xbar.1st-xbar)
    
    ybar.lrD.var1=(1/n.1st-N.inv)*sy2
    ybar.lrD.var2=(1/n.2nd-1/n.1st)*(sy2+beta^2*sx2-2*beta*syx)
    ybar.lrD.var =ybar.lrD.var1+ybar.lrD.var2
    
    ybar.lrD.sd =sqrt(ybar.lrD.var)
    
    ci.result=ci(ybar.lrD.est, ybar.lrD.sd, alpha)
    ybar.lrD.left =ci.result$left
    ybar.lrD.right=ci.result$right
    
    ybar.lrD.result=matrix(c(ybar.lrD.est, ybar.lrD.var, ybar.lrD.sd, ybar.lrD.left, ybar.lrD.right), nrow=1)
    colnames(ybar.lrD.result)=c("Est", "Var", "SD", "Left", "Right")
    rownames(ybar.lrD.result)=c("Mean_lrD")
  }
  return(ybar.lrD.result=as.data.frame(ybar.lrD.result))
}

######总量回归估计################################################################
lrD.total=function(N=NULL, n.1st, xbar.1st, y.sample, x.sample, alpha, method="Min", beta0=NULL)
{
  ybar.lrD.result=lrD.mean(N, n.1st, xbar.1st, y.sample, x.sample, alpha, method, beta0)
  ybar.lrD.est  =ybar.lrD.result$Est
  ybar.lrD.var  =ybar.lrD.result$Var
  ybar.lrD.left =ybar.lrD.result$Left
  ybar.lrD.right=ybar.lrD.result$Right
  
  ytot.lrD.est  =N*ybar.lrD.est              
  ytot.lrD.var  =N^2*ybar.lrD.var
  ytot.lrD.sd   =sqrt(ytot.lrD.var)
  ytot.lrD.left =N*ybar.lrD.left
  ytot.lrD.right=N*ybar.lrD.right
  
  total.result=matrix(c(ytot.lrD.est, ytot.lrD.var, ytot.lrD.sd, ytot.lrD.left, ytot.lrD.right), nrow=1)
  colnames(total.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(total.result)="Total_lrD"
  return(total.result=as.data.frame(total.result))
}
