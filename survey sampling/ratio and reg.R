#######比估计(R的估计)################################################################
ratio=function(y.sample, x.sample, N=NULL, auxiliary=FALSE, Xbar=NULL, alpha)
{
  n=length(y.sample)
  f=ifelse(is.null(N), 0, n/N)
  nf=(1-f)/n
  
  ybar=mean(y.sample)
  xbar=mean(x.sample)
  
  sy2 =var(y.sample)
  sx2 =var(x.sample)
  syx =cov(y.sample, x.sample)
  
  cy2 =nf*sy2/(ybar^2)
  cx2 =nf*sx2/(xbar^2)
  cyx =nf*syx/(ybar*xbar)
  
  ratio.est=ybar/xbar
  
  if(auxiliary==FALSE)
  {
    ratio.var=(nf/xbar^2)*(sy2+ratio.est^2*sx2-2*ratio.est*syx)
    ratio.sd =sqrt(ratio.var)
  }
  else
  {
    ratio.var=(nf/Xbar^2)*(sy2+ratio.est^2*sx2-2*ratio.est*syx)
    ratio.sd =sqrt(ratio.var)
  }  
  
  
  ### CI: Method 1 #######
  ci1=ci(ratio.est, ratio.sd, alpha)
  left1 =ci1$left
  right1=ci1$right
  
  
  ### CI: Method 2 #######
  quan=qnorm(1-alpha/2)
  
  est2=1-quan^2*cyx
  var2=(cy2+cx2-2*cyx)-quan^2*(cy2*cx2-cyx^2)
  
  ci2=ci(est2, sqrt(var2), alpha)
  left2 =ratio.est*(ci2$left)/(1-quan^2*cx2)
  right2=ratio.est*(ci2$right)/(1-quan^2*cx2)
  
  ### CI: Method 3 #######
  var3=ratio.est^2*(cy2+cx2-2*cyx)
  
  ci3=ci(ratio.est, sqrt(var3), alpha)
  left3 =ci3$left
  right3=ci3$right
  
  ratio.ci=t(cbind(c(left1, right1), c(left2, right2), c(left3, right3)))
  colnames(ratio.ci)=c("Left", "Right")
  rownames(ratio.ci)=c("Classic", "Exact", "Exact2")
  return(list(ratio.est=ratio.est, ratio.var=ratio.var, ratio.sd=ratio.sd, ratio.ci=as.data.frame(ratio.ci)))
}

######均值比估计################################################################
ratio.mean=function(y.sample, x.sample, N=NULL, Xbar, alpha)
{ 
  ratio.result=ratio(y.sample, x.sample, N, auxiliary=TRUE, Xbar, alpha)
  ratio.est=ratio.result$ratio.est
  ratio.var=ratio.result$ratio.var
  ratio.ci =ratio.result$ratio.ci
  
  ybarR.est=Xbar*ratio.est              
  ybarR.var=Xbar^2*ratio.var
  ybarR.sd =sqrt(ybarR.var)
  ybarR.ci =Xbar*ratio.ci
  
  return(list(ybarR.est=ybarR.est, ybarR.var=ybarR.var, ybarR.sd=ybarR.sd, 
              ybarR.ci=as.data.frame(ybarR.ci)))
}

######总量比估计#################################################################
ratio.total=function(y.sample, x.sample, N, Xbar, alpha)
{
  ybarR.result=ratio.mean(y.sample, x.sample, N, Xbar, alpha)
  ybarR.est=ybarR.result$ybarR.est
  ybarR.var=ybarR.result$ybarR.var
  ybarR.ci =ybarR.result$ybarR.ci
  
  ytot.est=N*ybarR.est              
  ytot.var=N^2*ybarR.var
  ytot.sd =sqrt(ytot.var)
  ytot.ci =N*ybarR.ci
  
  return(list(ytot.est=ytot.est, ytot.var=ytot.var, ytot.sd=ytot.sd, 
              ytot.ci=as.data.frame(ytot.ci)))
}

######均值回归估计###############################################################
reg.mean=function(y.sample, x.sample, N=NULL, Xbar, alpha, method="Min", beta0=NULL)
{
  n=length(y.sample)
  f=ifelse(is.null(N), 0, n/N)
  nf=(1-f)/n
  
  ybar=mean(y.sample)
  xbar=mean(x.sample)
  
  sy2 =var(y.sample)
  sx2 =var(x.sample)
  syx =cov(y.sample, x.sample)
  
  if (method=="Min")
  {
    beta=syx/sx2
    ybar.reg.est=ybar+beta*(Xbar-xbar)
    ybar.reg.var=nf*(n-1)/(n-2)*(sy2-syx^2/sx2)
    ybar.reg.sd =sqrt(ybar.reg.var)
    
    ci=ci(ybar.reg.est, ybar.reg.sd, alpha)
    left =ci$left
    right=ci$right
    
    ybar.reg.result=matrix(c(ybar.reg.est, ybar.reg.var, ybar.reg.sd, left, right), nrow=1)
    colnames(ybar.reg.result)=c("Est", "Var", "SD", "Left", "Right")
    rownames(ybar.reg.result)=c("Mean_Reg")
  }
  
  if (method=="Constant")
  {
    beta=beta0
    ybar.reg.est=ybar+beta*(Xbar-xbar)
    ybar.reg.var=nf*(sy2+beta^2*sx2-2*beta*syx)
    ybar.reg.sd =sqrt(ybar.reg.var)
    
    ci=ci(ybar.reg.est, ybar.reg.sd, alpha)
    left =ci$left
    right=ci$right
    
    ybar.reg.result=matrix(c(ybar.reg.est, ybar.reg.var, ybar.reg.sd, left, right), nrow=1)
    colnames(ybar.reg.result)=c("Est", "Var", "SD", "Left", "Right")
    rownames(ybar.reg.result)=c("Mean_Reg")
  }
  return(ybar.reg.result=as.data.frame(ybar.reg.result))
}

######总量回归估计####################################################################
reg.total=function(y.sample, x.sample, N=NULL, Xbar, alpha, method="Min", beta0=NULL)
{
  ybar.reg.result=reg.mean(y.sample, x.sample, N, Xbar, alpha, method, beta0)
  ybar.reg.est  =ybar.reg.result$Est
  ybar.reg.var  =ybar.reg.result$Var
  ybar.reg.left =ybar.reg.result$Left
  ybar.reg.right=ybar.reg.result$Right
  
  ytot.reg.est  =N*ybar.reg.est              
  ytot.reg.var  =N^2*ybar.reg.var
  ytot.reg.sd   =sqrt(ytot.reg.var)
  ytot.reg.left =N*ybar.reg.left
  ytot.reg.right=N*ybar.reg.right
  
  ytot.reg.result=matrix(c(ytot.reg.est, ytot.reg.var, ytot.reg.sd, ytot.reg.left, ytot.reg.right), nrow=1)
  colnames(ytot.reg.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(ytot.reg.result)=c("Total_Reg")
  
  return(ytot.reg.result=as.data.frame(ytot.reg.result))
}

######均值分别比估计##################################################################
s.ratio.mean=function(Nh, y.sample, x.sample, stra.index, Xbar, alpha)
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
    y.hth=y.sample[stra.index==h]
    x.hth=x.sample[stra.index==h]
    stra.ratio=ratio.mean(y.hth, x.hth, Nh[h], Xbar[h], alpha)
    
    yh.est[h]  =stra.ratio$ybarR.est
    yh.var[h]  =stra.ratio$ybarR.var
    yh.sd[h]   =stra.ratio$ybarR.sd
    
    yh.ci      =stra.ratio$ybarR.ci
    yh.ci.left =yh.ci$Left
    yh.ci.right=yh.ci$Right
    yh.left[h] =yh.ci.left[1]
    yh.right[h]=yh.ci.right[1]
  }
  stra.result=cbind(Nh, Wh, yh.est, yh.var, yh.sd, yh.left, yh.right)
  
  yRS.est=sum(Wh*yh.est)
  yRS.var=sum(Wh^2*yh.var)
  yRS.sd =sqrt(yRS.var)
  
  yRS.ci=ci(yRS.est, yRS.sd, alpha)
  yRS.left =yRS.ci$left
  yRS.right=yRS.ci$right
  
  yRS.result=matrix(c(yRS.est, yRS.var, yRS.sd, yRS.left, yRS.right), nrow=1)
  colnames(yRS.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(yRS.result)="Mean_RS"
  return(list(stra.result=as.data.frame(stra.result), yRS.result=as.data.frame(yRS.result)))
}

######总量分别比估计####################################################################
s.ratio.total=function(Nh, y.sample, x.sample, stra.index, Xbar, alpha)
{
  N=sum(Nh)
  Wh=Nh/sum(Nh)
  
  mean.RS.result=s.ratio.mean(Nh, y.sample, x.sample, stra.index, Xbar, alpha)
  RS.stra=mean.RS.result$stra.result
  RS.mean=mean.RS.result$yRS.result
  
  yh.totR.est  =N*RS.stra$yh.est
  yh.totR.var  =N^2*RS.stra$yh.var
  yh.totR.sd   =sqrt(yh.totR.var)
  yh.totR.left =N*RS.stra$yh.left
  yh.totR.right=N*RS.stra$yh.right
  
  stra.result=cbind(Nh, Wh, yh.totR.est, yh.totR.var, yh.totR.sd, yh.totR.left, yh.totR.right)
  
  ytot.RS.est  =N*RS.mean$Est             
  ytot.RS.var  =N^2*RS.mean$Var
  ytot.RS.sd   =sqrt(ytot.RS.var)
  ytot.RS.left =N*RS.mean$Left
  ytot.RS.right=N*RS.mean$Right
  
  ytot.RS.result=matrix(c(ytot.RS.est, ytot.RS.var, ytot.RS.sd, ytot.RS.left, ytot.RS.right), nrow=1)
  colnames(ytot.RS.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(ytot.RS.result)="TOTAL_RS"
  return(list(stra.result=as.data.frame(stra.result), ytot.RS.result=as.data.frame(ytot.RS.result)))
}

######均值联合比估计####################################################################
c.ratio.mean=function(Nh, y.sample, x.sample, stra.index, Xbar, alpha)
{
  yst.result=ss.mean2(Nh, y.sample, stra.index, alpha)$mean.result
  xst.result=ss.mean2(Nh, x.sample, stra.index, alpha)$mean.result
  
  ratio.est=yst.result$mean.est/xst.result$mean.est
  yRC.est=Xbar*ratio.est
  
  stra.num=length(Nh)
  Wh=Nh/sum(Nh)
  
  nh  =rep(0, stra.num)
  sy2 =rep(0, stra.num)
  sx2 =rep(0, stra.num)
  syx =rep(0, stra.num)
  
  for (h in 1:stra.num)
  {
    y.hth=y.sample[stra.index==h]
    x.hth=x.sample[stra.index==h]
    
    nh[h]=length(y.hth)
    
    sy2[h] =var(y.hth)
    sx2[h] =var(x.hth)
    syx[h]=cov(y.hth, x.hth)
  }
  
  fh=nh/Nh
  nf=(1-fh)/nh
  
  stra.result=cbind(Nh, Wh, nh, fh)
  
  yRC.var=sum(Wh^2*nf*(sy2+ratio.est^2*sx2-2*ratio.est*syx))
  yRC.sd =sqrt(yRC.var)
  
  yRC.ci=ci(yRC.est, yRC.sd, alpha)
  yRC.left =yRC.ci$left
  yRC.right=yRC.ci$right
  
  yRC.result=matrix(c(yRC.est, yRC.var, yRC.sd, yRC.left, yRC.right), nrow=1)
  colnames(yRC.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(yRC.result)="Mean_RC"
  return(list(stra.result=as.data.frame(stra.result), yRC.result=as.data.frame(yRC.result)))
}

######总量联合比估计####################################################################
c.ratio.total=function(Nh, y.sample, x.sample, stra.index, Xbar, alpha)
{
  N=sum(Nh)
  
  mean.RC.result=c.ratio.mean(Nh, y.sample, x.sample, stra.index, Xbar, alpha)
  RC.stra=mean.RC.result$stra.result
  RC.mean=mean.RC.result$yRC.result
  
  ytot.RC.est  =N*RC.mean$Est             
  ytot.RC.var  =N^2*RC.mean$Var
  ytot.RC.sd   =sqrt(ytot.RC.var)
  ytot.RC.left =N*RC.mean$Left
  ytot.RC.right=N*RC.mean$Right
  
  ytot.RC.result=matrix(c(ytot.RC.est, ytot.RC.var, ytot.RC.sd, ytot.RC.left, ytot.RC.right), nrow=1)
  colnames(ytot.RC.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(ytot.RC.result)="TOTAL_RC"
  return(list(stra.result=as.data.frame(RC.stra), ytot.RC.result=as.data.frame(ytot.RC.result)))
}

######均值分别回归估计####################################################################
s.reg.mean=function(Nh, y.sample, x.sample, stra.index, Xbar, alpha)
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
    y.hth=y.sample[stra.index==h]
    x.hth=x.sample[stra.index==h]
    stra.reg=reg.mean(y.hth, x.hth, Nh[h], Xbar[h], alpha,method="Min", beta0=NULL)
    
    yh.est[h]  =stra.reg$Est
    yh.var[h]  =stra.reg$Var
    yh.sd[h]   =stra.reg$SD
    
    yh.left[h] =stra.reg$Left
    yh.right[h]=stra.reg$Right
  }
  stra.result=cbind(Nh, Wh, yh.est, yh.var, yh.sd, yh.left, yh.right)
  
  yRegS.est=sum(Wh*yh.est)
  yRegS.var=sum(Wh^2*yh.var)
  yRegS.sd =sqrt(yRegS.var)
  
  yRegS.ci=ci(yRegS.est, yRegS.sd, alpha)
  yRegS.left =yRegS.ci$left
  yRegS.right=yRegS.ci$right
  
  yRegS.result=matrix(c(yRegS.est, yRegS.var, yRegS.sd, yRegS.left, yRegS.right), nrow=1)
  colnames(yRegS.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(yRegS.result)="Mean_RegS"
  return(list(stra.result=as.data.frame(stra.result), yRegS.result=as.data.frame(yRegS.result)))
}

######总量分别回归估计###################################################################
s.reg.total=function(Nh, y.sample, x.sample, stra.index, Xbar, alpha)
{
  N=sum(Nh)
  Wh=Nh/sum(Nh)
  
  mean.RegS.result=s.reg.mean(Nh, y.sample, x.sample, stra.index, Xbar, alpha)
  RegS.stra=mean.RegS.result$stra.result
  RegS.mean=mean.RegS.result$yRegS.result
  
  yh.totReg.est  =N*RegS.stra$yh.est
  yh.totReg.var  =N^2*RegS.stra$yh.var
  yh.totReg.sd   =sqrt(yh.totReg.var)
  yh.totReg.left =N*RegS.stra$yh.left
  yh.totReg.right=N*RegS.stra$yh.right
  
  stra.result=cbind(Nh, Wh, yh.totReg.est, yh.totReg.var, yh.totReg.sd, yh.totReg.left, yh.totReg.right)
  
  ytot.RegS.est  =N*RegS.mean$Est             
  ytot.RegS.var  =N^2*RegS.mean$Var
  ytot.RegS.sd   =sqrt(ytot.RegS.var)
  ytot.RegS.left =N*RegS.mean$Left
  ytot.RegS.right=N*RegS.mean$Right
  
  ytot.RegS.result=matrix(c(ytot.RegS.est, ytot.RegS.var, ytot.RegS.sd, ytot.RegS.left, ytot.RegS.right), nrow=1)
  colnames(ytot.RegS.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(ytot.RegS.result)="TOTAL_RegS"
  return(list(stra.result=as.data.frame(stra.result), ytot.RegS.result=as.data.frame(ytot.RegS.result)))
}

######均值联合回归估计####################################################################
c.reg.mean=function(Nh, y.sample, x.sample, stra.index, Xbar, alpha,method="Min", beta0=NULL)
{
  yst.result=ss.mean2(Nh, y.sample, stra.index, alpha)$mean.result
  xst.result=ss.mean2(Nh, x.sample, stra.index, alpha)$mean.result
  stra.num=length(Nh)
  Wh=Nh/sum(Nh)
  
  nh  =rep(0, stra.num)
  sy2 =rep(0, stra.num)
  sx2 =rep(0, stra.num)
  syx =rep(0, stra.num)
  
  for (h in 1:stra.num)
  {
    y.hth=y.sample[stra.index==h]
    x.hth=x.sample[stra.index==h]
    
    nh[h]=length(y.hth)
    
    sy2[h] =var(y.hth)
    sx2[h] =var(x.hth)
    syx[h]=cov(y.hth, x.hth)
  }
  
  fh=nh/Nh
  nf=(1-fh)/nh
  
  if (method=="Min")
  {
    beta=sum(Wh^2*nf*syx)/sum(Wh^2*nf*sx2)
    yRegC.est=yst.result$mean.est+beta*(Xbar-xst.result$mean.est)
  }
  
  if (method=="Constant")
  {
    beta=beta0
    yRegC.est=yst.result$mean.est+beta*(Xbar-xst.result$mean.est)
  }
  
  stra.result=cbind(Nh, Wh, nh, fh)
  
  yRegC.var=sum(Wh^2*nf*(sy2+beta^2*sx2-2*beta*syx))
  yRegC.sd =sqrt(yRegC.var)
  
  yRegC.ci=ci(yRegC.est, yRegC.sd, alpha)
  yRegC.left =yRegC.ci$left
  yRegC.right=yRegC.ci$right
  
  yRegC.result=matrix(c(yRegC.est, yRegC.var, yRegC.sd, yRegC.left, yRegC.right), nrow=1)
  colnames(yRegC.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(yRegC.result)="Mean_RegC"
  return(list(stra.result=as.data.frame(stra.result), yRegC.result=as.data.frame(yRegC.result)))
}

######总量联合回归估计####################################################################
c.reg.total=function(Nh, y.sample, x.sample, stra.index, Xbar, alpha,method="Min", beta0=NULL)
{
  N=sum(Nh)
  
  mean.RegC.result=c.reg.mean(Nh, y.sample, x.sample, stra.index, Xbar, alpha,method="Min", beta0=NULL)
  RegC.stra=mean.RegC.result$stra.result
  RegC.mean=mean.RegC.result$yRegC.result
  
  ytot.RegC.est  =N*RegC.mean$Est             
  ytot.RegC.var  =N^2*RegC.mean$Var
  ytot.RegC.sd   =sqrt(ytot.RegC.var)
  ytot.RegC.left =N*RegC.mean$Left
  ytot.RegC.right=N*RegC.mean$Right
  
  ytot.RegC.result=matrix(c(ytot.RegC.est, ytot.RegC.var, ytot.RegC.sd, ytot.RegC.left, ytot.RegC.right), nrow=1)
  colnames(ytot.RegC.result)=c("Est", "Var", "SD", "Left", "Right")
  rownames(ytot.RegC.result)="TOTAL_RegC"
  return(list(stra.result=as.data.frame(RegC.stra), ytot.RegC.result=as.data.frame(ytot.RegC.result)))
}
