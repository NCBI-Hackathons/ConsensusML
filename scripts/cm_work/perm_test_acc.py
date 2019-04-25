import pandas as pd
import numpy as np

from sklearn.metrics import log_loss, accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.ensemble import  RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import permutation_test_score


df = pd.read_csv('deg.csv')
df.drop(columns=['Unnamed: 0'], inplace=True)

testing = pd.read_csv('../../composite_code/rnotebook/data/TARGET_AML_Testing_Samples.csv')
training = pd.read_csv('../../composite_code/rnotebook/data/TARGET_AML_Training_Samples.csv')

training_id = list(training['x'])
testing_id = list(testing['x'])

train_df = df[df['TARGET.USI'].isin(training_id)]
test_df = df[df['TARGET.USI'].isin(testing_id)]

model_columns = train_df.columns[1:-1]

X_train = train_df[model_columns].astype(float)
y_train = train_df.iloc[:, -1]

X_test = test_df[model_columns].astype(float)
y_test = test_df.iloc[:, -1]

X = df[model_columns]
y = df.iloc[:, -1]

#call models
xgb = XGBClassifier(learning_rate = 0.01, max_depth = 3, n_estimators = 700, random_state=8, n_jobs=-1)
rf = RandomForestClassifier(n_estimators=1000, max_depth=20, random_state=8, n_jobs=-1)

# do xgb permutation test for accuracy
xgb_score, perm_scores, p_value = permutation_test_score(xgb, X, y, n_permutations=5000, cv=3,
                                                        n_jobs=-1, scoring='accuracy')
xgb_perm_df = pd.DataFrame({'xgb_score': [xgb_score] * 5000, 'perm_scores': perm_scores,
                        'p_value': [p_value] * 5000})
pd.DataFrame.to_csv(xgb_perm_df, 'xgb_permutation_acc.csv')

# do rf permutation test for accuracy
rf_score, perm_scores, p_value = permutation_test_score(rf, X, y, n_permutations=5000, cv=3,
                                                        n_jobs=-1, scoring='accuracy')
rf_perm_df = pd.DataFrame({'rf_score': [rf_score] * 5000, 'perm_scores': perm_scores,
                        'p_value': [p_value] * 5000})
pd.DataFrame.to_csv(rf_perm_df, 'rf_permutation_acc.csv')