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
