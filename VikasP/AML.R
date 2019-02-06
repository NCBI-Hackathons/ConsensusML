library(MLSeq)
library(DESeq2)
library(xlsx)
library(limma)
library(VennDiagram)

#reading in clinical data files
AML_clinical<-read.csv('AML_assay_clinical.csv')

#raw count expression data
all_expression_data<-read.csv('TARGET_AML_High.Std.Risk_vs_LowRisk_DEGs', header=TRUE)

#setting rownames for MLseq
rownames(all_expression_data)<-all_expression_data$Genes
all_expression_data$Genes<-NULL

#removing duplicate samples
AML_colnames<-c()
for(i in 1:ncol(all_expression_data)){ 
  split_string<-strsplit(colnames(all_expression_data[i]),"[.]")[[1]][3]
  AML_colnames<-append(AML_colnames, split_string)
}

uniquecheck<-duplicated(AML_colnames)
todelete<-c()
for(i in 1:length(uniquecheck)){ 
  if(uniquecheck[i]==TRUE){ 
    todelete<-append(todelete,i)    
  }
}

all_expression_data[,todelete]<-NULL


#deletes any columns with non-unqiue 5 string identifier from AML_clinical
AML_clinical_subset<-c()
for(i in 1:nrow(AML_clinical)){ 
  split_string<-strsplit(as.character(AML_clinical$TARGET.USI[i]),"-|\\s")[[1]][3]
  AML_clinical_subset<-append(AML_clinical_subset, split_string)
}

uniquecheck<-duplicated(AML_clinical_subset)
todelete<-c()
for(i in 1:length(uniquecheck)){ 
  if(uniquecheck[i]==TRUE){ 
    todelete<-append(todelete,i)    
  }
}

AML_clinical<-AML_clinical[-todelete,]



AML_clinical_patients<-(AML_clinical$TARGET.USI)
AML_clinical_patients<-as.character(AML_clinical_patients)
all_expression_data_patients<-colnames(all_expression_data)

count=0
notunique<-c()
all_expression_aml_positives<-c()
delete_aml_clinical_unknown_risk<-c()
delete_expression_unknown_risk<-c()
AML_clinical$Risk.group<-as.character(AML_clinical$Risk.group)
for(i in 1:nrow(AML_clinical)){ 
  for(j in 1:ncol(all_expression_data)) { 
  split_string<-strsplit(as.character(AML_clinical$TARGET.USI[i]), "-|\\s")[[1]][3]
  if(grepl(split_string, all_expression_data_patients[j])){ 
      count = count + 1
      all_expression_aml_positives<-append(all_expression_aml_positives,j)
      if(AML_clinical$Risk.group[i]=='Unknown'){ 
        print(colnames(AML_clinical$Risk.group[i]))
        delete_aml_clinical_unknown_risk<-append(delete_aml_clinical_unknown_risk,i)
        delete_expression_unknown_risk<-append(delete_expression_unknown_risk,j)
        AML_clinical$Risk.group[i]<-('DELETE')
        all_expression_data[2,j]<-NA
        
        }
    }
  }
}


all_expression_aml_subset<-all_expression_data[,c(0,all_expression_aml_positives)]

#removes Unkonwn risk of death samples 
todelete<-c()
for(i in 1:nrow(AML_clinical)){ 
  if(AML_clinical$Risk.group[i]=='DELETE'){ 
      todelete<-append(todelete,i)      
    }
  }

AML_clinical<-AML_clinical[-todelete,]

todelete<-c()
for(i in 1:ncol(all_expression_aml_subset)){ 
  if(is.na(all_expression_aml_subset[2,i])){ 
      todelete<-append(todelete,i)
    }
  }
all_expression_aml_subset<-all_expression_aml_subset[,-todelete]



risk_group<-c()

for(i in 1:(ncol(all_expression_aml_subset))){ 
  for (j in 1:nrow(AML_clinical)){ 
  split_string<-strsplit(colnames(all_expression_aml_subset[i]),"[.]")[[1]][3]
      if(grepl(split_string, AML_clinical$TARGET.USI[j])){
        risk_group<-append(risk_group,as.character(AML_clinical$Risk.group[j]) )
      }
    }
  }


for( i in 1:length(risk_group)){
  if(risk_group[i]=="High"){
    risk_group[i]<-'H'
  }
  if(risk_group[i]=="Standard"){
    risk_group[i]<-'H'
  }
  if(risk_group[i]=='Low'){
    risk_group[i]<-'L'
  }
}

deleterows<-c()
rows<-rownames(all_expression_aml_subset)
for(i in 1:length(rows)){ 
  if(grepl("__",rows[i])){ 
    deleterows<-append(deleterows,i)
  }
}

all_expression_aml_subset<-all_expression_aml_subset[-deleterows,]


class<-DataFrame(condition=factor(risk_group))

set.seed(2128)
vars <- sort(apply(all_expression_aml_subset, 1, var, na.rm = TRUE), decreasing = TRUE)
data <- all_expression_aml_subset[names(vars)[1:100], ]
nTest <- ceiling(ncol(data) * 0.3)
ind <- sample(ncol(data), nTest, FALSE)

data.train <- as.matrix(data[ ,-ind] + 1)
data.test <- as.matrix(data[ ,ind] + 1)
classtr <- DataFrame(condition = class[-ind, ])
classts <- DataFrame(condition = class[ind, ])

data.trainS4 = DESeqDataSetFromMatrix(countData = data.train, colData = classtr,
                                      design = formula(~condition))

data.testS4 = DESeqDataSetFromMatrix(countData = data.test, colData = classts,
                                     design = formula(~condition))


fit <- classify(data = data.trainS4, method = "svmRadial",
                preProcessing = "deseq-rlog", ref = "H",
                control = trainControl(method = "repeatedcv", number = 2,
                                       repeats = 2, classProbs = TRUE))
show(fit)



set.seed(2128)
# Voom based Nearest Shrunken Centroids.
fit <- classify(data = data.trainS4, method = "voomNSC",
                normalize = "deseq", ref = "H",
                control = voomControl(tuneLength = 20))
trained(fit) ## Trained model summary




set.seed(2128)
# Support vector machines with radial basis function kernel
fit.svm <- classify(data = data.trainS4, method = "svmRadial",
                    preProcessing = "deseq-vst", ref = "H", tuneLength = 10,
                    control = trainControl(method = "repeatedcv", number = 5,
                                           repeats = 10, classProbs = TRUE))
show(fit.svm)


trained(fit.svm)


plot(fit.svm)


# Define control list
ctrl.svm <- trainControl(method = "repeatedcv", number = 5, repeats = 1)
ctrl.plda <- discreteControl(method = "repeatedcv", number = 5, repeats = 1,
                             tuneLength = 10)
ctrl.voomDLDA <- voomControl(method = "repeatedcv", number = 5, repeats = 1,
                             tuneLength = 10)
# Support vector machines with radial basis function kernel
fit.svm <- classify(data = data.trainS4, method = "svmRadial",
                    preProcessing = "deseq-vst", ref = "H", tuneLength = 10,
                    control = ctrl.svm)
# Poisson linear discriminant analysis
fit.plda <- classify(data = data.trainS4, method = "PLDA", normalize = "deseq",
                     ref = "H", control = ctrl.plda)
# Voom-based diagonal linear discriminant analysis
fit.voomDLDA <- classify(data = data.trainS4, method = "voomDLDA",
                         normalize = "deseq", ref = "H", control = ctrl.voomDLDA)



trained(fit.voomDLDA)


pred.svm <- predict(fit.svm, data.testS4)
pred.svm

pred.svm <- relevel(pred.svm, ref = "H")
actual <- relevel(classts$condition, ref = "H")
tbl <- table(Predicted = pred.svm, Actual = actual)
confusionMatrix(tbl, positive = "H")



set.seed(2128)
# Define control lists.
ctrl.continuous <- trainControl(method = "repeatedcv", number = 5, repeats = 10)
ctrl.discrete <- discreteControl(method = "repeatedcv", number = 5, repeats = 10,
                                 tuneLength = 10)
ctrl.voom <- voomControl(method = "repeatedcv", number = 5, repeats = 10,
                         tuneLength = 10)
# 1. Continuous classifiers, SVM and NSC
fit.svm <- classify(data = data.trainS4, method = "svmRadial",
                    preProcessing = "deseq-vst", ref = "H", tuneLength = 10,
                    control = ctrl.continuous)
fit.NSC <- classify(data = data.trainS4, method = "pam",
                    preProcessing = "deseq-vst", ref = "H", tuneLength = 10,
                    control = ctrl.continuous)
# 2. Discrete classifiers
fit.plda <- classify(data = data.trainS4, method = "PLDA", normalize = "deseq",
                     ref = "H", control = ctrl.discrete)
fit.plda2 <- classify(data = data.trainS4, method = "PLDA2", normalize = "deseq",
                      ref = "H", control = ctrl.discrete)
fit.nblda <- classify(data = data.trainS4, method = "NBLDA", normalize = "deseq",
                      ref = "H", control = ctrl.discrete)
# 3. voom-based classifiers
fit.voomDLDA <- classify(data = data.trainS4, method = "voomDLDA",
                         normalize = "deseq", ref = "H", control = ctrl.voom)
fit.voomNSC <- classify(data = data.trainS4, method = "voomNSC",
                        normalize = "deseq", ref = "H", control = ctrl.voom)
# 4. Predictions
pred.svm <- predict(fit.svm, data.testS4)
pred.NSC <- predict(fit.NSC, data.testS4)




voomNSC<-selectedGenes(fit.voomNSC)
voomDLDA<-selectedGenes(fit.voomDLDA)
PLDA<-selectedGenes(fit.plda)
PLDA2<-selectedGenes(fit.plda2)
nblda<-selectedGenes(fit.nblda)
svm<-selectedGenes(fit.svm)


write.csv(voomNSC, 'voomnsc.csv')
write.csv(PLDA,'plda.csv')
write.csv(PLDA2,'plda2.csv')
#annotate the genes 
#look at them in ommim 


#genes IDs from biomart search of the csv files created
all_genes<-c("MPO", "MPO", "RPS20", "CD74", "HSPA5", "PABPC1", "ACTB", "ACTB", "ACTB", "ACTB", "ACTB", "ACTB", "ACTB", "FTL", "FTL", "FTL", "FTL", "FTL", "FTL", "LGALS1", "RPL3", "RPS19", "RPS19", "RPS19", "RPL28", "PFN1", "PFN1", "GAPDH", "RPS12", "DUSP1", "SRGN", "ZFP36", "RPS15A", "RPS6", "RPS24", "RPS24", "RPS11", "RPL13A", "SH3BGRL3", "RPS3A", "EEF1A1", "RPL29", "LAPTM5", "RPL9", "RPL27A", "B2M", "B2M", "B2M", "B2M", "EEF2", "EEF2", "FTH1", "FTH1", "CXCL8", "FOS", "AZU1", "CFL1", "EIF1", "RPL4", "RPS27", "RPS27", "RPS17", "RPS17", "ACTG1", "ACTG1", "ACTG1", "ACTG1", "ACTG1", "ACTG1", "ACTG1", "H1FX", "PTMA", "HBA2", "HBA2", "HBA2", "HBA2", "HBA2", "HBA2", "HBA2", "ELANE", "ELANE", "ELANE", "RPS26", "RPS26", "RPL37A", "RPL12", "RPS4X", "MT-ND6", "MT-ND6", "MT-ND6", "MT-ND6", "MT-ND6", "MT-ND6", "MT-CYB", "MT-CYB", "MT-CYB", "MT-CYB", "MT-CYB", "MT-ND2", "MT-ND2", "MT-ND2", "MT-ND2", "MT-ND5", "MT-ND5", "MT-ND5", "MT-ND5", "MT-ND5", "MT-ND5", "MT-CO1", "MT-CO1", "MT-CO1", "MT-CO1", "MT-CO1", "MT-CO1", "MT-CO1", "MT-CO1", "MT-ND3", "MT-ND3", "MT-ND3", "MT-ND4", "MT-ND4", "MT-ND4", "MT-ND4", "MT-ND4", "MT-ND4", "MT-ND1", "MT-ND1", "MT-ND1", "MT-ND1", "MT-ND1", "MT-ND1", "MT-ATP6", "MT-ATP6", "MT-ATP6", "MT-ATP6", "MT-ATP6", "MT-ATP6", "MT-ATP6", "MT-ATP6", "MT-ATP6", "MT-CO3", "MT-CO3", "MT-CO3", "MT-CO3", "MT-CO3", "MT-CO3", "HLA-DRA", "TMSB4X", "MT-RNR2", "MT-ND4L", "MT-ND4L", "MT-ATP8", "MT-ATP8", "MT-ATP8", "RPS18", "HLA-B", "HLA-B", "HLA-B", "HLA-B", "HLA-B", "HLA-B", "HLA-B", "HLA-B", "HLA-B", "HLA-B", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "HBB", "MTATP6P1", "MALAT1")
all_genes<-unique(all_genes)
all_genes

