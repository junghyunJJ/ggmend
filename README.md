# ggmend
GGmend: A Mendelian randomization method for finding gene-on-gene regulatory effects in the presence of unobserved confounders


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

#### load exmaple dataset in R.
```
> library(ggmend)
> data(data_ggmend)
> str(data_ggmend)
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


#### run ggmend
```
> res <- ggmend(data_ggmend$egenes, data_ggmend$cissnps,
+               effectth = 0, corth = 0, standardize = TRUE,
+               maxperms = 1000, threads = 1,seed = 1)
Loading required package: pbapply
n = 100 / m = 14
Data filtering: |effectsize| > 0 and |correlation| > 0
standardizing
Number of snp/gene pair: 91
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 04s

```
