% Copyright 2016,2017 Lingfei Wang
% 
% This file is part of lassopv.
% 
% Lassopv is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% Lassopv is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with lassopv.  If not, see <http://www.gnu.org/licenses/>.
% 
\name{lassopv}
\alias{lassopv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of Nonparametric P-Value Estimation for Predictors in Lasso
}
\description{
This function estimates p-values for predictors x against target variable y in lasso regression, using the regularization strength when each predictor enters the active set of regularization path for the first time as the statistic. This is based on the assumption that predictors that (first) become active earlier tend to be more significant. Null distribution for each predictor is computed analytically under approximation, which aims at efficiency and accuracy for small p-values.
}
\usage{
lassopv(x,y,trace = FALSE,normalize=TRUE,Gram,eps = .Machine$double.eps,

max.steps,use.Gram=TRUE,log.p=FALSE,max.predictors=NULL) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Input matrix of predictor variables.}
  \item{y}{Input vector of target variable.}
  \item{trace}{Whether traces lasso regression. See lars in package lars.}
  \item{normalize}{Whether every predictor is scaled to unit variance first. Every predictor is forcefully shifted to zero mean regardless of this argument.}
  \item{Gram}{Optional Gram used by lasso regression in lars.}
  \item{eps}{Precision for enet function.}
  \item{max.steps}{The optional maximum number steps for lasso regression. See lars in package lars.}
  \item{use.Gram}{Whether to use Gram in lasso regression. See lars in package lars.}
  \item{log.p}{Whether to output log p-values instead.}
  \item{max.predictors}{The number of top predictors to estimate p-values for. Defaults to all predictors.}
}
\value{
Vector of p-values for predictors. Predictors never entered the active set of regularization path within the given max.steps or not within the top (max.predictors) predictors have p-value=1. If log.p is set, outputs log p-values instead.
}
%\note{}

%\seealso{}
\examples{
library(lars)
library(lassopv)
data(diabetes)
attach(diabetes)
pv=lassopv(x,y)
}
