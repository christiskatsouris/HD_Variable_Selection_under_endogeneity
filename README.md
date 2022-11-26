# Variable Selection and Goodness-of-fit Testing in High-Dimensional Linear Models under Endogeneity

## Description 

The R package ['Variable_Selection_GoF_testing_under_endogeneity'](https://github.com/christiskatsouris/Variable_Selection_GoF_testing_under_endogeneity) implements the coding procedure for the research project titled: "Variable Selection and Goodness-of-fit Testing in High-Dimensional Linear Models under Endogeneity" which is joint work with [Dr Panayiota Touloupou](https://www.birmingham.ac.uk/staff/profiles/maths/touloupou-panayiota.aspx), School of Mathematics, University of Birmingham, United Kingdom. 

## Methodology

In this paper, we propose a novel dimension reduction methodology which operates within a sequential variable selection framework. Specifically, a high dimensional vector of available covariates is scanned sequentially in order to identify the best correctly specified model and associated selected covariates, that is, the correct functional form. Furthermore, the statistical model under consideration tackles the aspect of endogeneity in high dimensional linear models. Our empirical application considers an important application from the field of medical sciences. 

## Installation (under development) 

The R package ['Variable_Selection_GoF_testing_under_endogeneity'](https://github.com/christiskatsouris/Variable_Selection_GoF_testing_under_endogeneity) will be able to be installed from Github.

## Usage: 

```R

# After development the package will be able to be installed using
install.packages("Variable_Selection_GoF_testing_under_endogeneity")
library("Variable_Selection_GoF_testing_under_endogeneity")

```

## Estimation Examples:

```R


```


## Simulation Studies:

```R

#############################################################################################
library(MASS)

# function to simulate multivariate normal distribution
# given gene number, sample size and correlation coefficient

multi_norm <- function(gene_num,sample_num,R) 
{# begin of function 

  # initial covariance matrix
  V <- matrix(data=NA, nrow=gene_num, ncol=gene_num)
  
  V <- matrix(data=NA, nrow=gene_num, ncol=gene_num)
  
  # mean for each gene
  meansmodule <- runif(gene_num, min=-3, max=3)
  # variance for each gene
  varsmodule <- runif(gene_num, min=0, max=5)
  
  for (i in 1:gene_num)
  {
    # a two-level nested loop to generate covariance matrix
    for (j in 1:gene_num) 
    {
      if (i == j) 
      {
        # covariances on the diagonal
        V[i,j] <- varsmodule[i]
      } 
      else 
      {
        # covariances
        V[i,j] <- R * sqrt(varsmodule[i]) * sqrt(varsmodule[j])
      }
    }
  }
  
  # simulate multivariate normal distribution
  # given means and covariance matrix
  
  X <- t(mvrnorm(n = sample_num, meansmodule, V))
  
  return(X)
  
}# end of function

multi_norm(5,5,0.4)

```

<p align="center">
  
<img src="https://github.com/christiskatsouris/Variable_Selection_GoF_testing_under_endogeneity/blob/main/data/MCstudy.png" width="780"/>

</p>  


<p align="center">
  
<img src="https://github.com/christiskatsouris/Variable_Selection_GoF_testing_under_endogeneity/blob/main/data/MCdesign1.png" width="700"/>

</p>  


# Bibliography

## Key References:

$\textbf{[1]}$ Chudik, A., Kapetanios, G., & Pesaran, M. H. (2018). A one covariate at a time, multiple testing approach to variable selection in high‐dimensional linear regression models. Econometrica, 86(4), 1479-1512. 

$\textbf{[2]}$ Fan, J., & Liao, Y. (2014). Endogeneity in high dimensions. Annals of statistics, 42(3), 872.

$\textbf{[3]}$ Lockhart, R., Taylor, J., & Tibshirani, R. (2014). A significance test for the lasso. Annals of statistics, 42(2), 413.

$\textbf{[4]}$ Gold, D., Lederer, J., & Tao, J. (2020). Inference for high-dimensional instrumental variables regression. Journal of Econometrics, 217(1), 79-111.

$\textbf{[5]}$ Breunig, C., Mammen, E., & Simoni, A. (2020). Ill-posed estimation in high-dimensional models with instrumental variables. Journal of Econometrics, 219(1), 171-200.

$\textbf{[6]}$  Windmeijer, F., Farbmacher, H., Davies, N., & Davey Smith, G. (2019). On the use of the lasso for instrumental variables estimation with some invalid instruments. Journal of the American Statistical Association, 114(527), 1339-1350.

$\textbf{[7]}$ Gautier, E., & Rose, C. (2011). High-dimensional instrumental variables regression and confidence sets. arXiv preprint arXiv:1105.2454.

$\textbf{[8]}$ Ke, T., Jin, J., & Fan, J. (2014). Covariance assisted screening and estimation. Annals of statistics, 42(6), 2202.

$\textbf{[9]}$ Shah, R. D., & Bühlmann, P. (2018). Goodness‐of‐fit tests for high dimensional linear models. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 80(1), 113-135.

$\textbf{[10]}$ Battey, H., Fan, J., Liu, H., Lu, J., & Zhu, Z. (2018). Distributed testing and estimation under sparse high dimensional models. Annals of statistics, 46(3), 1352.

$\textbf{[11]}$  Yousuf, K. (2018). Variable screening for high dimensional time series. Electronic Journal of Statistics, 12(1), 667-702.

$\textbf{[12]}$  Campbell, F., & Allen, G. I. (2017). Within group variable selection through the exclusive lasso. Electronic Journal of Statistics, 11(2), 4220-4257.

$\textbf{[13]}$ Zhang, C. H. (2010). Nearly unbiased variable selection under minimax concave penalty. The Annals of statistics, 38(2), 894-942.

$\textbf{[14]}$ Hagemann, A. (2012). A simple test for regression specification with non-nested alternatives. Journal of Econometrics, 166(2), 247-254.

$\textbf{[15]}$ Katsouris, C. (2021a). Optimal Portfolio Choice and Stock Centrality for Tail Risk Events. [arXiv:2112.12031](https://arxiv.org/abs/2112.12031).

$\textbf{[16]}$  Katsouris, C. (2021b). Forecast Evaluation in Large Cross-Sections of Realized Volatility. [arXiv:2112.04887](https://arxiv.org/abs/2112.04887).



# Acknowledgments

The second author (Christis G. Katsouris) greatfully acknowledges financial support from the [Department of Economics](http://business-school.exeter.ac.uk/about/departments/economics/) of the [Faculty of Environment, Science and Economy](https://www.exeter.ac.uk/departments/ese/) at the University of Exeter, United Kingdom. 

Christis G. Katsouris is a Lecturer in Economics at the [University of Exeter Business School](http://business-school.exeter.ac.uk/). He is also a member of the [Time Series and Machine Learning Group](https://www.personal.soton.ac.uk/cz1y20/Reading_Group/mlts-group-2022.html) at the [School of Mathematical Sciences](https://www.southampton.ac.uk/about/faculties-schools-departments/school-of-mathematical-sciences) (Statistics Division) of the University of Southampton. 

# Code of Coduct

Please note that the ['Variable_Selection_GoF_testing_under_endogeneity'](https://github.com/christiskatsouris/Variable_Selection_GoF_testing_under_endogeneity) R project will be released with a Contributor Code of Coduct (under construction). By contributing to this project, you agree to abide by its terms.

# Declarations

The author (Christis G. Katsouris) declares no conflicts of interest.

Notice that the academic research presented here is considered to be as open access to the academic and non-academic community. Therefore, we would appreciate it if appropriate acknolwedgement is given to statistical methodologies and econometric procedures developed by academic researchers and made available to the wider applied data scientist community.   

# Historical Background

> Standing on the shoulders of giants.
> 
> ''_If I have been able to see further, it was only because I stood on the shoulders of giants._"
> Isaac Newton, 1676 


#### Charles Stein

Charles Stein  (22 March 1920 – 24 November 2016) was an American mathematical statistician and professor of statistics at Stanford University. 
He received his Ph.D in 1947 at Columbia University with advisor Abraham Wald. He held faculty positions at Berkeley and the University of Chicago before moving permanently to Stanford in 1953. He is known for Stein's paradox in decision theory, which shows that ordinary least squares estimates can be uniformly improved when many parameters are estimated; for Stein's lemma, giving a formula for the covariance of one random variable with the value of a function of another when the two random variables are jointly normally distributed; and for Stein's method, a way of proving theorems such as the Central Limit Theorem that does not require the variables to be independent and identically distributed (Source: Wikepedia).
