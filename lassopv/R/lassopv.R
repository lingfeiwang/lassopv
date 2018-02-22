# Copyright 2016-2018 Lingfei Wang
# 
# This file is part of lassopv.
# 
# Lassopv is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Lassopv is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with lassopv.  If not, see <http://www.gnu.org/licenses/>.
# 
lassopv<-function (x,y,normalize=TRUE,H0=c("spherical","normal"),log.p=FALSE,max.predictors=NULL,trace = FALSE,Gram,eps = .Machine$double.eps,max.steps,use.Gram=TRUE) 
{
	dx=dim(x)[2]
	ds=dim(x)[1]
	if(ds!=length(y))
		stop("Incorrect dimensions for x or y. The number of rows in x should match the length of y.")
	if(is.null(max.predictors))
		max.predictors=dx
	else if(max.predictors<=0)
		stop("Incorrect max.predictors")

	isdone=rep(F,dx)
	if(dx==1)
	{
		ans=cor(drop(x),y)		
		isdone[1]=T
	}
	else
	{
		#Normalize
		x=t(t(x)-colMeans(x))
		sx=sqrt(colMeans(x**2))
		y=y-mean(y)
		sy=sqrt(mean(y**2))
		if(normalize)
		{
			x=t(t(x)/sx)
			sx=sqrt(colMeans(x**2))
		}
		a=lars(x,y,type='lasso',trace=trace,normalize=F,intercept=F,Gram=Gram,eps=eps,max.steps=max.steps,use.Gram=use.Gram)
		#Crop away ascending lambda values (supposed bug in lars)
		nstep=which(a$lambda[2:length(a$lambda),drop=F]>a$lambda[1:length(a$lambda)-1])
		if(length(nstep)==0)
			nstep=length(a$lambda)
		else
			nstep=nstep[1]
		#lambda values
		lambda=a$lambda[1:nstep,drop=F]/ds
		#Residue
		yres=y-predict.lars(a,x,s=a$lambda[1:nstep,drop=F],mode='lambda')$fit
		#Residue variance
		pv=t(t(yres)-colMeans(yres))
		pv=colMeans(pv**2)
		wcl=lambda/sqrt(pv)
		#ans=lambda_i/(sigma_i*sigma_yres)
		ans=rep(0.,dx)
		ntot=0
		for (i in 1:nstep)
		{
			if(!is.finite(wcl[i]))
				next
			pa=unlist(a$actions[i])
			for(j in 1:length(pa))
			{
				if((pa[j]<=0)||isdone[pa[j]])
					next
				ans[pa[j]]=wcl[i]/sx[pa[j]]
				isdone[pa[j]]=T
				ntot=ntot+1
			}
			if(ntot>=max.predictors)
				break
		}
	}

	if(H0[1]=="normal")
	{
		ans=pchisq((ans**2)*ds,1,lower.tail=F,log.p=log.p)
	}
	else if(H0[1]=="spherical")
	{
		ans=pt(abs(ans)*sqrt((ds-2)/(1-ans**2)),ds-2,lower.tail=F,log.p=log.p)
		if(log.p)
			ans=ans+log(2)
		else
			ans=ans*2
	}
	else
		stop("Unknown null distribution type H0.")
	if(log.p)
	{
		ans[!isdone]=0
		#Account for numerical precision limitations
		ans[ans>0]=0
	}
	else
	{
		ans[!isdone]=1
		#Account for numerical precision limitations
		ans[ans>1]=1
	}
	names(ans)=colnames(x)
	ans
}

require("lars")
require("stats")
























