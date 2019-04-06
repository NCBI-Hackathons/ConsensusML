import pandas as pd
import numpy as np
from tqdm import tqdm
import warnings
from sklearn.ensemble import  RandomForestClassifier
from collections import Counter

from data_prep import create_x_y_data

X_train, y_train, X_test, y_test, X, y = create_x_y_data()
columns = X.columns

warnings.filterwarnings("ignore", category=RuntimeWarning)

def boruta_bootstrap(n_rounds, X, y):
    rf_boruta = RandomForestClassifier(n_jobs=-1)
    
    total_boruta_features = []
    for n in tqdm(range(n_rounds)):
        bootstrap_X = X.sample(n=200, replace=True)
        bootstrap_y = y[bootstrap_X.index]
        
        feat_selector = BorutaPy(rf_boruta, n_estimators='auto', verbose=0, max_iter = 100, random_state=8)
        feat_selector.fit(bootstrap_X.values, bootstrap_y.values)
        
        boruta_mask = feat_selector.support_
        rf_boruta_features = columns[boruta_mask]
        total_boruta_features.append(rf_boruta_features)
    return total_boruta_features

def feature_counter(bf):
    boruta_counter = Counter()
    for boruta_round in bf:
        for feat in boruta_round:
            boruta_counter[feat] +=1
    df = pd.DataFrame(boruta_counter, index=[0]).T
    return df

bf = boruta_bootstrap(1, X, y)
df = feature_counter(bf)
pd.DataFrame.to_csv(df, 'bootstrap_boruta.csv')