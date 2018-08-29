# Generalized Group Lasso for Patient Subgroup Identification
Wenxuan Deng

## Introduction

Prognostic biomarkers and predictive biomarkers.

Why decision trees is not workable? Because the sample size in real clinical datasets is too small, typically no more than 100 patients.

Group lasso [^group]

Elastic net [^elastic] adaptive weights for elastic net [^elasticAdaptive]

Hierarchical Group lasso for interactions [^InteractionGroup]

Overlapping group lasso [^Overlapping1][^Overlapping2][^theoOverlapping]

Sparse Group Lasso [^notesparse] [^sparsegroup]

Structured group lasso [^binyu]

Group lasso for logistic regression [^logistic]

Other variable selection methods:

GUIDE: a regression tree [^GiGs][^GUIDE]

SIS: screening [^SIS][^iSIS]

SIR: [^SIR][^logiSIR]

Stepwise selection: [^stepwise]



## Methods

#### Model

#### Ordinary Linear Model

$$Y=X_0\beta_0 + X_T\beta_\tau + X_1\beta_1 + X_T\otimes X_1 \beta_2+\epsilon$$

Where $$X_0$$ is the baseline variables, $$X_T$$ is the treatment variable, $$X_1$$ is the high dimensional design matrix of genes, i.e. gene expression levels, SNP and mutations, and $$X_T\otimes X_1$$ is the interaction between genes and treatment. $$\beta=(\beta_0, \beta_\tau, \beta_1, \beta_2)$$ is the corresponding coefficients. $$\epsilon$$ is random error.

Let $$X=[X_0,X_T,X_1^{(1)},\dots,X_m^{(1)},X_TX_1^{(1)},\dots,X_TX_1^{(m)}]$$ and $$\beta=[\beta_0,\beta_\tau,\beta_1^{(1)},\dots,\beta_1^{(m)},\beta_2^{(1)},\dots,\beta_2^{(m)}]$$. For each gene $$l$$, its prognostic and predictive design matrix is denoted as $$X^{(l)}=[X_1^{(l)},X_TX_1^{(l)}]$$ and its corresponding coefficients are $$\beta^{(l)}=[\beta_1^{(l)},\beta_2^{(l)}]$$

#### Loss Function

We used group lasso and elastic net for variables selection when $$n\ll p$$, and assumed the hierarchical relationship between prognostic biomarkers and predictive biomarkers, that is the predictive biomarkers should be a prognostic biomarkers. The loss function is 

$$\min_{\theta} f(\beta|Y,X_0,X_T,X_1) + g(\beta)$$

$$g(\beta)=\lambda_1\sum_i\phi_i|\beta_2^{(i)}|+\lambda_1\sum_i\psi_i\sqrt{(\beta_1^{(1)})^2+(\beta_1^{(1)})^2}+\lambda_2(\parallel\beta_1\parallel_2^2+\parallel\beta_2\parallel_2^2)$$

Where $$\beta=(\beta_0,\beta_\tau,\beta_1,\beta_2)$$ is the parameter, and $$f(\beta|Y,X_0,X_T,X_1) $$ is $$L$$-2 loss function. When the model is the ordinary linear model, the $$L$$-2 loss function is $$\parallel Y-(X_0\beta_0 + X_T\beta_\tau + X_1\beta_1+ X_T\otimes X_1 \beta_2) \parallel^2$$. Penalty function $$g(\beta)$$ can construct a complex hierarchical selection of $$\beta_1$$ and $$\beta_2$$, that nonzero $$\beta_2$$ is a sufficient but not necessary condition for nonzero $$\beta_1$$. The contour plot for a pair of $$\beta_1$$ and $$\beta_2$$ is shown in Figure 1. $$\lambda_1$$ and $$\lambda_2$$ are regularization parameters.

![](/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/CodeForFigures/Contour.png)
*Figure 1*

#### Criterion and Adaptive Weights

KKT [^unique]

For group $$\hat{\beta^{(l)}}$$, the KKT condition is 

$${X^{(l)}}^T(Y-X\hat{\beta})=\lambda_1\phi_l$$ 



## Algorithms

Fast iterative shrinkage-thresholding algorithm with backtracking[^FISTA]

Proximal operator for group lasso [^fastGrouplasso]

Adaptive restart for rippling behavior [^restart]

Adaptive stepwise of cyclic Barzilai-Borwein spectral approach [^stepsize]

![](/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/CodeForFigures/Algorithm.png)

## Experiments

Signal to noise ratio: $$SNR=\frac{Var(X\beta)}{Var(\epsilon)}$$


## Reference


[^unique]: Tibshirani, Ryan J. "The lasso problem and uniqueness." Electronic Journal of Statistics 7 (2013): 1456-1490.

[^theoOverlapping]: Percival, Daniel. "Theoretical properties of the overlapping groups lasso." Electronic Journal of Statistics 6 (2012): 269-288.

[^Overlapping1]: Jacob, Laurent, Guillaume Obozinski, and Jean-Philippe Vert. "Group lasso with overlap and graph lasso." Proceedings of the 26th annual international conference on machine learning. ACM, 2009.

[^Overlapping2]: Obozinski, Guillaume, Laurent Jacob, and Jean-Philippe Vert. "Group lasso with overlaps: the latent group lasso approach." arXiv preprint arXiv:1110.0413 (2011).

[^group]: Yuan, Ming, and Yi Lin. "Model selection and estimation in regression with grouped variables." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 68.1 (2006): 49-67.

[^binyu]: Zhao, Peng, Guilherme Rocha, and Bin Yu. "The composite absolute penalties family for grouped and hierarchical variable selection." The Annals of Statistics 37.6A (2009): 3468-3497.

[^InteractionGroup]: Lim, Michael, and Trevor Hastie. "Learning interactions via hierarchical group-lasso regularization." Journal of Computational and Graphical Statistics 24.3 (2015): 627-654.

[^notesparse]: Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. "A note on the group lasso and a sparse group lasso." arXiv preprint arXiv:1001.0736 (2010).

[^sparsegroup]: Simon, Noah, et al. "A sparse-group lasso." Journal of Computational and Graphical Statistics 22.2 (2013): 231-245.

[^elastic]: Zou, Hui, and Trevor Hastie. "Regularization and variable selection via the elastic net." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 67.2 (2005): 301-320.

[^elasticAdaptive]: Zou, Hui, and Hao Helen Zhang. "On the adaptive elastic-net with a diverging number of parameters." Annals of statistics 37.4 (2009): 1733.

[^logistic]: Meier, Lukas, Sara Van De Geer, and Peter Bühlmann. "The group lasso for logistic regression." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 70.1 (2008): 53-71.z

[^GiGs]: Loh, Wei‐Yin, Xu He, and Michael Man. "A regression tree approach to identifying subgroups with differential treatment effects." Statistics in medicine 34.11 (2015): 1818-1833.

[^GUIDE]: Loh, Wei-Yin. "Regression tress with unbiased variable selection and interaction detection." Statistica Sinica (2002): 361-386.

[^SIS]: Fan, Jianqing, and Jinchi Lv. "Sure independence screening for ultrahigh dimensional feature space." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 70.5 (2008): 849-911.

[^iSIS]: Fan, Jianqing, Richard Samworth, and Yichao Wu. "Ultrahigh dimensional feature selection: beyond the linear model." Journal of machine learning research 10.Sep (2009): 2013-2038.

[^SIR]: Jiang, Bo, and Jun S. Liu. "Sliced inverse regression with variable selection and interaction detection." arXiv preprint arXiv:1304.4056 652 (2013).

[^logiSIR]: Li, Yang, and Jun S. Liu. "Robust variable and interaction selection for logistic regression and general index models." Journal of the American Statistical Association (2018): 1-16.

[^stepwise]: Miller, Alan J. "Selection of subsets of regression variables." Journal of the Royal Statistical Society. Series A (General) (1984): 389-425.

[^FISTA]: Beck, Amir, and Marc Teboulle. "A fast iterative shrinkage-thresholding algorithm for linear inverse problems." SIAM journal on imaging sciences 2.1 (2009): 183-202.

[^fastGrouplasso]: Liu, Jun, and Jieping Ye. "Fast overlapping group lasso." arXiv preprint arXiv:1009.0306 (2010).

[^restart]: O’donoghue, Brendan, and Emmanuel Candes. "Adaptive restart for accelerated gradient schemes." Foundations of computational mathematics 15.3 (2015): 715-732.

[^stepsize]: Wright, Stephen J., Robert D. Nowak, and Mário AT Figueiredo. "Sparse reconstruction by separable approximation." IEEE Transactions on Signal Processing 57.7 (2009): 2479-2493.


