# A new Integrated Mean Variance Correlation and Its Use in High-Dimensional Data Analysis
The goal of package `newIMVC` is to provide an easy way to implement the proposed methods in Xiong et al. (2024), which include a new robust correlation between continuous variables and its use in hypothesis test, feature screening and FDR control.
## Installation
As the dependent package `limma` has been removed from the CRAN repository, one should install the package `limma` from github before the installation of the package `newIMVC`,
```R
install.packages("devtools")
library("devtools")
install_github("gangwug/limma")
```
To install `newIMVC`,
```R
install_github("scottpanhan/newIMVC")
```
## Example
Here is an example showing how to use main functions in package `newIMVC`.
```R
library("newIMVC")
library("MASS")
library("mvtnorm")

n=200
x=rnorm(n)
y=x^2+rt(n,2)
IMVC(y,x,K=10,type="nonlinear")
#> [1] 0.3005832

n=200
p=500
#sigma_z=diag(rep(1,pp))
pho1=0.8
mean_x=rep(0,p)
sigma_x=matrix(NA,nrow = p,ncol = p)
for (i in 1:p) {
  for (j in 1:p) {
    sigma_x[i,j]=pho1^(abs(i-j))
  }
}
x=rmvnorm(n, mean = mean_x, sigma = sigma_x,method = "chol")
x1=x[,1]
x2=x[,2]
x3=x[,12]
x4=x[,22]
y=2*x1+0.5*x2+3*x3*ifelse(x3<0,1,0)+2*x4+rnorm(n)
IMVCS(y,x,K=5,d=round(n/log(n)),type="nonlinear")
#> [1]   1   2  22  12  11  10   3  21  13   9  23  15   4  14   8  20  16   5   6  24  49  18   7 472  19
#> [26]  25  17 294 473 426  47 326 395  26 178 394 461  27
```
