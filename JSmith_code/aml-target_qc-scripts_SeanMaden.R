# Script for AML TARGET samples processing and summaries

# author: SKM
# contact: maden@ohsu.edu

#-------------------
# Read in the files
#-------------------
fn = 'AML_dataframe.csv'
aml.cd <- read.csv(fn,stringsAsFactors = F)
fn = 'TARGET_NBL_AML_RT_WT_TMMCPM_log2_Norm_Counts.csv'
pancan.tmm <- read.csv(fn, stringsAsFactors = F)

#----------------
# Data summaries
#----------------
table(aml.cd$Gender)
# Female   Male 
# 96     91

summary(aml.cd$Age.at.Diagnosis.in.Days)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 137     966    3273    3205    5251    8231

summary(aml.cd$Age.at.Diagnosis.in.Days/365)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3753  2.6466  8.9671  8.7819 14.3863 22.5507

summary(as.numeric(aml.cd$Bone.marrow.leukemic.blast.percentage....))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 14.00   58.25   73.60   70.92   89.00  100.00       5

table(aml.cd$CNS.disease)
# No Yes 
# 175  12

# CLASSIFIERS OF INTEREST
# 1. young vs old
hist(aml.cd$Age.at.Diagnosis.in.Days)
med.age <- median(aml.cd$Age.at.Diagnosis.in.Days)
abline(v=med.age,col="blue")

class.age <- ifelse(aml.cd$Age.at.Diagnosis.in.Days<med.age,'young','old')
table(class.age, aml.cd$Gender)

# class.age Female Male
# old       46   48
# young     50   43

# 2. survival time
summary(aml.cd$Overall.Survival.Time.in.Days)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 112.0   501.5  1229.0  1474.1  2364.0  4022.0
summary(aml.cd$Event.Free.Survival.Time.in.Days)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 77.0   269.0   392.0   823.9   861.5  3630.0 

# survival by age grouping
summary(aml.cd[class.age=='young',]$Overall.Survival.Time.in.Days)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 112     600    1464    1555    2364    4022
summary(aml.cd[class.age=='young',]$Event.Free.Survival.Time.in.Days)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 85.0   269.0   373.0   884.3   809.0  3630.0

summary(aml.cd[class.age=='old',]$Overall.Survival.Time.in.Days)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 163.0   447.8   924.0  1393.8  2351.8  3291.0 
summary(aml.cd[class.age=='old',]$Event.Free.Survival.Time.in.Days)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 77.0   262.5   409.0   764.1   887.8  3144.0

# barcode matrix
bcm <- matrix(nrow=0,ncol=5)
for(i in 1:nrow(aml.cd)){
  bcm <- rbind(bcm,unlist(strsplit(aml.cd$entity_submitter_id[i],"-")))
}

x <- colnames(aml.tmm)
bcm2 <- matrix(nrow=0,ncol=5)
for(i in 2:length(x)){bcm2 <- rbind(bcm2,unlist(strsplit(x[i],'\\.')))}

# 01A 02A 03A 04A 06A 09A 11A 40A 
# 262  12  26  41   1 119  11   2 

# 01R 02R 03R 04R 05R 
# 353  52  48  13   8 

# Notes:
# want diagnostic samples only first...
# summarize samples by tissue type(s)


library(ggplot2)

ggdat <- as.data.frame(matrix(ncol=2,nrow=0))
ggdat <- rbind(ggdat,data.frame(group='young.overallsurv',survival.time=aml.cd[class.age=='young',]$Overall.Survival.Time.in.Days))
ggdat <- rbind(ggdat,data.frame(group='young.efsurv',survival.time=aml.cd[class.age=='young',]$Event.Free.Survival.Time.in.Days))
ggdat <- rbind(ggdat,data.frame(group='old.overallsurv',survival.time=aml.cd[class.age=='old',]$Overall.Survival.Time.in.Days))
ggdat <- rbind(ggdat,data.frame(group='old.efsurv',survival.time=aml.cd[class.age=='old',]$Event.Free.Survival.Time.in.Days))

ggplot(ggdat, aes(x=ggdat$survival.time, col=ggdat$group))+geom_density()+
  theme(panel.background = element_rect(fill = 'white',colour = 'black'),
        rect = element_rect(fill = 'white',colour = "white"),
        panel.grid.major = element_line(colour = 'grey75', size=0.2),
        panel.grid.minor = element_line(colour = 'white'),
        legend.position = 'right',
        legend.background = element_rect(fill = "white", 
                                         colour ="white"),
        legend.key = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5))+
  labs(color="Group Survival") + 
  ggtitle("Survival Time by Age Classifier")

#------------------------------
# Gene annotations and GRanges
#------------------------------
# gene annotations
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
edb.genes <- genes(edb)

# gene granges

#----------------------------------
# Compile as Summarized Experiment
#----------------------------------
# https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
