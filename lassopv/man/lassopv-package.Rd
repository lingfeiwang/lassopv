% Copyright 2016-2018 Lingfei Wang
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
\name{lassopv-package}
\alias{lassopv-package}
\docType{package}
\title{
Nonparametric P-Value Estimation for Predictors in Lasso
}
\description{
Estimate the p-values for predictors x against target variable y in lasso regression, using the regularization strength when each predictor enters the active set of regularization path for the first time as the statistic. This is based on the assumption that predictors (of the same variance) that (first) become active earlier tend to be more significant. Three null distributions are supported: normal and spherical, which are computed separately for each predictor and analytically under approximation, which aims at efficiency and accuracy for small p-values.
}
\details{
This R package provides a simple and efficient method to estimate the p-value of every predictor on a given target variable. The method is based on lasso regression and compares when every predictor enters the active set of the regulatization path against a normally distributed null predictor. The null distribution is computed analytically under approximation, whose errors are small for significant predictors. The whole computation only requires a single lasso regression over the regularization path, and is capable of analyzing high dimensional datasets.
}
\author{
Lingfei Wang <Lingfei.Wang.github@outlook.com>
}
\references{
Lingfei Wang and Tom Michoel, Comparable variable selection with lasso, https://arxiv.org/pdf/1701.07011. 2017, 2018.
}
%\keyword{ package }
\examples{
library(lars)
library(lassopv)
data(diabetes)
attach(diabetes)
pv=lassopv(x,y)
}
