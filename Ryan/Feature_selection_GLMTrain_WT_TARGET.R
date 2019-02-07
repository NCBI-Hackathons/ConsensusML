# Feature selection, training and basic statistics
# of GLMBoosted regression model for AML WT RNA-Seq Data
library(plyr)
library(mlr)
df <- wt
# Examine distribution of variable of intrest
hist(df$Event.Free.Survival.Time.in.Days,breaks = 30)
abline(v=313, col = 'red')

# dichotomize data set based on if they're above or below the mean 
median(df$Event.Free.Survival.Time.in.Days)
for (x in seq(1,length(df$Event.Free.Survival.Time.in.Days))) {
  if( df$Event.Free.Survival.Time.in.Days[x] <= 313) {
    df$cohort[x] <- "young"
  } else {
    df$cohort[x] <- "old"
  }
}

head(df$cohort)

# split into two cohorts based on EVST relative to the median
old_group <- df[df$cohort== "young",]
young_group <- df[df$cohort== "old",]

count(old_group, c('Gender'))
count(young_group, c('Gender'))
count(df, c('Gender'))
length(old_group$Event.Free.Survival.Time.in.Days)
length(young_group$Event.Free.Survival.Time.in.Days)
hist(old_group$Event.Free.Survival.Time.in.Days)
hist(young_group$Event.Free.Survival.Time.in.Days)

dim(old_group)
dim(young_group)



diff <- data.frame(seq(42,21448))
colnames(diff) <- c('X')
l = vector("list", 21448)

# get the between group variance of all features
for (x in seq(42,21448)) {
  diff$y[x-41] <- mean(old_group[,x]) - mean(young_group[,x])
}
plot(diff$y)
abline(h=0.05, col='red')
abline(h=-0.05, col='red')

# pick out the features that aren't super close to zero vairance between sets ~16000 genes
expression_set_0.05 <- c(9,which(abs(diff$y) > 0.05) + 41)
es_0.05 <- wt[,expression_set_0.05]

# build a glmboosted model to regress for EFST
regr.task_es <- makeRegrTask(id='es_0.05', data = es_0.05, target = 'Event.Free.Survival.Time.in.Days')
regr.lrn = makeLearner('regr.glmboost')
n = getTaskSize(regr.task_es)
# split data in half into test and train
train.set = seq(1,n,by=2)
test.set = seq(1,n,by=2)

# train
mod = train(regr.lrn, regr.task_es, subset = train.set)
# test
task.pred = predict(mod, task = regr.task_es, subset = test.set)

# parse and display results 
results = as.data.frame(task.pred)

r_line_df = seq(0,5000)
regr_ln.pred = predict(mod, task = regr.task_es, )
plt <- ggplot(data = results, aes(x=results$truth, y = results$response)) + geom_point() + 
  theme_classic() + geom_vline(xintercept = 331, linetype = 'dashed') + geom_hline(yintercept = 331, linetype = 'dashed') + 
  xlim(0,5000) + ylim(0,5000) + geom_abline(intercept = 0, slope = 1)
show(plt)
ggsave('Regression_on_training_data.pdf', height = 10, width = 10)
cp= 0
cn = 0
pp = 0
pn = 0
tp = 0
tn = 0
fp = 0
fn = 0
for (x in seq(1, length(results$id))) {
  if (results$truth[x] > 331){
    cp = cp+ 1
    if (results$response[x] > 331) {
      tp = tp + 1
      pp = pp + 1
    } else {
      fn = fn + 1
      pn = pn + 1
    } 
  } else { 
    cn = cn + 1
    if (results$response[x] > 331) {
      fp = fp + 1
      pp = pp + 1
    } else {
      tn = tn + 1
      pn = pn +1 
    }
  }

  
}


