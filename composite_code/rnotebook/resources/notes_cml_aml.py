"""
Notes for CML-AML manuscript

"""

"""
Introduction/Background

A number of prior studies applied machine learning algorithms to AML data.
These studies used a variety of genomics data, including gene expression and genome sequencing, 
to predict various clinical factors including survival, time to relapse, and risk of progression.
Prior related work applied a variety of methods, including multivariate regression, lasso, and others.

Prior work identified shared and unique characteristics of importance for molecular 
subtyping in pediatric and adult AML cases, and this suggests there are age-dependent 
and - independent molecular characteristics in AML as a disease. We leveraged these
prior findings as a kind of validation for characteristics implicated in disease 
progression risk from ensemble machine learning approach.


To our knowlege this is the first study to apply an ensemble machine learning approach
to predict risk of disease progression in pediatric AML cases. Our workflow with each 
prediction method was essentially a two-step process of initial model fitting, followed by
feature implication and consensus across disparate models. Thanks to ongoing 
innovations and improvements to algorithms and code underlying various model fitting
methods programmatically, we were able to run various model fitting algorithms 
efficiently on a budget of time and memory resources. Our investigation paves the
way for wider involvement by the research, data science, and bioinformatics communities
to shed further light on important biological characteristics underlying complex 
diseases, such as AML, with potential clinical relevance.



"""

"""
Preprocessing, Summary Statistics, and Preliminary Analyses

1. DEGs
To reduce dimensionality of the gene expression data, we initially applied a filter
for the most differentially expressed genes (DEGs, t-test, p-adj < 0.05). After 
applying this filter, we were left with a subset of just under 2,000 genes.

2. Gene Expression Covariance and Correlation
We were further interested in evidence of correlated gene expression in the subset
of DEGs identified. We applied gene ontology (GO) tree traversal method as well as
basic bivariate correlation analyses, to identify gene sets likely to be correlated 
(Spearman test, p-adj < 0.05). For subsequent analyses, we included the full set 
of DEGs, as well as the subset of DEGs predicted not to be correlated.

3. Summary Statistics
As mentioned, the classifier of interest was chosen in part for its clinical relevance,
and in part because TARGET AML sample subsets of each classifier group were approximately
balanced on important demographic factors (e.g. gender, age at diagnosis, event 
free survival, etc.). See Table 1 for details.
(Table1)


"""


"""
Machine learning methods

1. Lasso

We used lasso (least absolute shrinkage and selection operator), a type of generalized linear model algorithm using penalized 
maximum likelihood. To apply lasso we used the R glmnet package. Lasso is a 
type of regression that returns coefficient assignments to each factor in a model. 
Factors of lower importance for fitted model prediction are assigned coefficients with absolute 
value closer to zero. It is important to note that, among correlated features, 
lasso will pick one and discard the others. This is by contrast to ridge regression,
which uses a penalty that shrinks the coeffients of correlated predictors togther.

https://cran.r-project.org/web/packages/glmnet/index.html
https://eight2late.wordpress.com/2017/07/11/a-gentle-introduction-to-logistic-regression-and-lasso-regularisation-using-r/

2. SVM

3. Random Forest

4. Random Forest with Boosting

5. AutoML

6. MLseq

...Others...

"""