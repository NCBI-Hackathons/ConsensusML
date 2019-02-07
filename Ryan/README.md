GLMBoosted Regression on Event Free Survival Time 

Input data was cleaned and prepared as outlined on the main README. 

[WT dataset](https://github.com/NCBI-Hackathons/ConsensusML/blob/master/scripts/WT_assay_clinical.csv.zip)

First distribution of Event Free Survial Time in Days (EFST) was plotted and the median (331 days) calculated. ![Distribution_of_EVST_median](https://github.com/NCBI-Hackathons/ConsensusML/blob/master/Ryan/Distribution_of_EVST_median.pdf) The data was then bisected on the median into two cohorts. Then each normalized gene expression data was averaged for each cohort and the difference between each cohort calculated. ![Log_variance_plot](hhttps://github.com/NCBI-Hackathons/ConsensusML/blob/master/Ryan/log_variance_plot.pdf). Then I removed the genes which had a variance from zero of less than 0.05.

Next, I split the data in half into test and training sets and used [R MLR package](https://mlr.mlr-org.com/) to build a GLM boosted regression model predicting EVST based off expression of the genes that survivied the variance trimming. 

Regression model performance can be seen [here](https://github.com/NCBI-Hackathons/ConsensusML/blob/master/Ryan/Regression_on_training_data.pdf) X - axis is ground truth EFST and y axis is predicted EFST based exclusively on normalized gene expression counts for the genes selected in the earlier steps. Average error is 103, only test set data points are shown. Dotted lines represent the original median of the distribution. This allows me to treat the model as a binary classifier based off if the model predicts the EVST on the same side of the median that the ground truth EVST was. 

TPR = 31/34 = 91%
Precision = 31/46 = 67%
Specificity = 17/32 = 53%
False Discovery Rate = 32%

However, the most interesting data is the genes (or features) which the GLM ended up using to create predictions - these features and their coefiicients are listed [here](https://github.com/NCBI-Hackathons/ConsensusML/blob/master/Ryan/feature_list.csv)

For example some of the genes that were corelated with low EVST are sumarized below:
HPS4,HSPB1 - Heat shock proteins that have been implicated in relapse of [HBV associated hepatocellular cancer](https://doi.org/10.7150/ijms.10735) 

TPCN2 which has been associated with [relapse of prostate cancer after surgery](https://doi.org/10.1007/s12094-018-02029-z)



