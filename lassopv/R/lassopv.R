# Copyright 2016,2017 Lingfei Wang
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
lassopv<-function (x,y,trace = FALSE,normalize=TRUE,Gram,eps = .Machine$double.eps,max.steps,use.Gram=TRUE,log.p=FALSE,max.predictors=NULL) 
{
	dx=dim(x)[2]
	ds=dim(x)[1]
	if( ds!=length(y) )
		stop("Incorrect dimensions for x or y.")
	if(is.null(max.predictors))
		max.predictors=dx
	else if(max.predictors<=0)
		stop("Incorrect max.predictors")

	if(dim(x)[2]==1)
		ans=cor(drop(x),y)
	else
	{
		#Normalize
		x=t(t(x)-colMeans(x))
		y=y-mean(y)
		sx=sqrt(colMeans(x**2))
		sy=sqrt(mean(y**2))
		if(normalize==T)
		{
			x=t(t(x)/sx)
			sx=sqrt(colMeans(x**2))
			y=y/sy
			sy=sqrt(mean(y**2))
		}
		a=lars(x,y,type='lasso',trace=trace,normalize=F,intercept=F,Gram=Gram,eps=eps,max.steps=max.steps,use.Gram=use.Gram)
		#Residue variance
		yv=mean(y**2)
		pv=(1-as.double(a$R2)[1:length(a$lambda)])*yv
		#lambda values
		wcl=a$lambda/(ds*sqrt(pv))
		pa=as.integer(unlist(a$actions))
		ans=rep(0.,dx)
		ntot=0
		for (i in 1:length(pa))
		{
			if((pa[i]>0)&&(is.finite(wcl[i]))&&(ans[pa[i]]==0.))
			{
				ans[pa[i]]=wcl[i]
				ntot=ntot+1
			}
			if(ntot>=max.predictors)
				break
		}
	}
	ans=pchisq((ans**2)*ds,1,lower.tail=F,log.p=log.p)
	names(ans)=colnames(x)
	ans
}

require("lars")


























