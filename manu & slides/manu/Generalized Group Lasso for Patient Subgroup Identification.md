# Generalized Group Lasso for Patient Subgroup Identification
Wenxuan Deng

## Introduction

Prognostic biomarkers and predictive biomarkers.

Why decision trees is not workable? Because the sample size in real clinical datasets is too small, typically no more than 100 patients.

## Methods

#### Model

#### Ordinary Linear Model

$$Y=X\beta + W\tau + G\alpha + W\otimes G \gamma+\epsilon$$

Where $$X$$ is the baseline variables, $$W$$ is the treatment variable, $$G$$ is the high dimensional design matrix of genes, i.e. gene expression levels, SNP and mutations, and $$W\otimes G$$ is the interaction between genes and treatment. $$\theta=(\beta, \tau, \alpha, \gamma)$$ is the corresponding coefficients. $$\epsilon$$ is random error. 

#### Loss Function

We used group lasso for variables selection when $$n\ll p$$, and assumed the hierarchical relationship between prognostic biomarkers and predictive biomarkers, that is the predictive biomarkers should be a prognostic biomarkers. The loss function is 

$$\min_{\theta} f(\theta|Y,X,W,G) + \lambda \sum_i \eta_i^I |\gamma_i| + \lambda \sum_i \sqrt{\alpha_i^2 + \gamma_i^2} $$

Where $$f(\theta|Y,X,W,G) $$ is $$L$$-2 loss function. When the model is the ordinary linear model, the $$L$$-2 loss function is $$\parallel Y-(X\beta + W\tau + G\alpha + W\otimes G \gamma) \parallel^2$$.



