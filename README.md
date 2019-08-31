# GGmend: A Mendelian randomization method for finding gene-on-gene regulatory effects in the presence of unobserved confounders


Abstract
Several studies have shown that correlated gene expression levels can be produced by complex biological regulatory networks. Unfortunately, confounding factors, such as shared environmental factors or *trans*-regulators affecting two genes, may spuriously induce correlations between the gene expressions. These false correlations prevent us from being able to make causal statements about the effect of one gene on another: gene-on-gene effects. By leveraging ideas from Mendelian randomization, it is possible to circumvent the confounding effects and estimate gene-on-gene effects. *cis*-eQTLs are genetic variants that regulate gene expression in nearby genes. These variants are subject to Mendelian randomization and can be used as instrumental variables to determine the direction and effect size of gene-on-gene regulating effects. Here we introduce a new method, GGmend (Gene-on-Gene effect estimator using Mendelian randomization) that accurately identifies gene-on-gene regulatory effects utilizing a quartet of two genes and their two *cis*-eQTLs. The GGmend is robust to expression correlations induced by confounding effects and that the family-wise error rate is well-calibrated simulated datasets. Besides, our method outperformed existing causal inference algorithms on the DREAM5 dataset. Further, we applied our method to yeast data and identified several putative causal regulators including previously identified causal genes by establishing gene-on-gene regulatory network. GGmend is publicly available at https://github.com/junghyunJJ/ggmend.

## Installation
> *We currently only support R 3.5+.*
```
install.packages("devtools")
library(devtools)
install_github("junghyunJJ/ggmend")
```

## Example

#### 1. Load exmaple dataset in R.
```
library(ggmend)
data(data_ggmend)
str(data_ggmend)
List of 2
 $ egenes : num [1:14, 1:100] 0.655 0.758 0.479 0.957 0.32 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:14] "G1" "G2" "G3" "G4" ...
  .. ..$ : NULL
 $ cissnps: int [1:14, 1:100] 0 0 0 1 1 1 1 1 1 0 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:14] "G1" "G2" "G3" "G4" ...
  .. ..$ : NULL
```


> **Note:**
> - *ys* (egenes) and *xs* (cissnps) require the same dimension.
> - Parallelization (threads > 1) works only on *nix (Linux, Unix such as macOS) system. please check **.Platform$OS.type** function.


#### 2. Run ggmend
```
res <- ggmend(data_ggmend$egenes, data_ggmend$cissnps,
               effectth = 0, corth = 0, standardize = TRUE,
               maxperms = 1000, threads = 1,seed = 1)
# Loading required package: pbapply
# n = 100 / m = 14
# Data filtering: |effectsize| > 0 and |correlation| > 0
# standardizing
# Number of snp/gene pair: 91
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 04s

```

#### 3. Convert data for visNetwork
```
head(res)
#   g1 g2 row col      Bs1g1      Bs2g2       GGcor      Bg1g2       Bg2g1        Bgg pval       FDR
# 1 G1 G2   1   2  0.6856368  0.4931179 -0.38072544 0.39706536 0.016778401 0.39706536 0.00 0.0000000
# 2 G1 G3   1   3  0.6856368 -0.5087897 -0.03265593 0.04406143 0.106112580 0.10611258 0.48 0.5824000
# 3 G2 G3   2   3  0.4931179 -0.5087897 -0.36814014 0.04331202 0.003397309 0.04331202 0.79 0.8169318
# 4 G1 G4   1   4  0.6856368  0.5628240 -0.08399198 0.08409654 0.042935828 0.08409654 0.42 0.5383099
# 5 G2 G4   2   4  0.4931179  0.5628240  0.25401625 0.12472683 0.122887443 0.12472683 0.47 0.5779730
# 6 G3 G4   3   4 -0.5087897  0.5628240 -0.46654818 0.20976835 0.013751757 0.20976835 0.19 0.4433333
```
filtering results using FDR < 0.05
```
netdata <- convertnet(res, fdrth = 0.05)
# filtering < 0.05...
# nodes: 13 / edges: 15

str(netdata)
# List of 2
#  $ edges:'data.frame':	15 obs. of  3 variables:
#   ..$ from: chr [1:15] "G1" "G6" "G3" "G8" ...
#   ..$ to  : chr [1:15] "G2" "G7" "G10" "G10" ...
#   ..$ Bgg : num [1:15] 0.397 0.453 0.567 0.503 0.461 ...
#  $ nodes:'data.frame':	13 obs. of  2 variables:
#   ..$ id   : Factor w/ 13 levels "G1","G10","G11",..: 1 2 3 4 5 6 7 8 9 10 ...
#   ..$ label: Factor w/ 13 levels "G1","G10","G11",..: 1 2 3 4 5 6 7 8 9 10 ...
```

#### 3. Plot network
```
plotnet(netdata$nodes, netdata$edges)
# Loading required package: visNetwork
```
http://rpubs.com/junghyunjj/523836
