# RNAseq_Cancer_Biomarkers
Machine Learning to Detect Cancer Biomarkers from RNAseq Data

# Hackathon Project Workflow
!["ML RNA-seq biomarkers, flow chart"](https://github.com/NCBI-Hackathons/RNAseq_Cancer_Biomarkers/blob/master/ml-fhack_day1-flowchart_v2_SeanMaden.jpg "Day 1 Flowchart")

# Background
Leukemia is a cancer of the blood arising in white blood cells of the bone marrow. It poses a substantial population burden as the most common pediatric cancer ([Steliarova-Foucher et al 2017](https://www.ncbi.nlm.nih.gov/pubmed/28410997); [SEER data](https://seer.cancer.gov/statfacts/html/amyl.html)). Acute Myelogenous Leukemeia (AML) is a type of leukemia impacting the myeloblast stem cells that accounts for a substantial population burden, arising at a current rate of approximately 20,000 cases per year with 27.4% 5-year survival 2. AML is molecularly heterogeneous, with several clinically relevant subtypes, including perhaps dozens of subtypes defined by factors ranging from cell differentiation state to cytogenic and sequencing assays ([Yi G. et al 2019](https://www.sciencedirect.com/science/article/pii/S2211124718320588?via%3Dihub); [Tyner J. W. et al 2018](https://www.nature.com/articles/s41586-018-0623-z)). Pediatric AML is characterized at a molecular level by rare somatic mutations, absence of common adult AML mutations, and relatively frequent structural variants ([Bolouri H et al 2018](https://www.nature.com/articles/nm.4439)). Here, we apply several machine learning approaches for feature selection of RNA-seq data from both pediatric and adult AML cases. Our goal was to better understand gene expression-based heterogeneity underlying AML cases, as well as age-related and -unrelated dysregulation patterns. We used clinical and assay data from pediatric cancer patients from the Therapeutically Applicable Research To Generate Effective Treatments (TARGET) initiative (https://ocg.cancer.gov/programs/target/).

# Methods and Analysis Overview
We were interested in applying machine learning principles for feature selection, to identify the most important genes and gene sets for predicting clinically-relevant classifiers in pediatric and adult AML cases. We initially analyzed the TARGET pediatric AML cohort. We combined primary peripheral blood and bone marrow samples, and assessed potential sources of bias and batch effects. Classifiers of main interest included age, stage, and survival. We further pre-filtered TMM-normalized expression based on extent of differential expression between these classifiers of interest, using multiple thresholds.

We then applied a series of machine learning algorithms for feature selection using various R and Python packages: 1. Support Vector Machines (SVM) using [e1071](https://cran.r-project.org/web/packages/e1071/index.html); 2. Random Forest with boosting in Python; 3. Neural Networks with [keras](https://cran.r-project.org/web/packages/keras/index.html); 4. Logistic regression in R; 5. Elastic net with Lasso using [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html); 6. [AutoML](https://pypi.org/project/automl/). 

# Links to Shared Documents

## 1. Manuscript
https://docs.google.com/document/d/1DPAmUFfggAnAjsMIPTs1hV90k25ZKckYLi18b3dBot0/edit#
