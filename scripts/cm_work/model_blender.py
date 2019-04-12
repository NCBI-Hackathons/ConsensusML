from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import log_loss, accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.ensemble import  RandomForestClassifier
from xgboost import XGBClassifier
import pandas as pd 
import numpy as np 
from tqdm import tqdm
from data_prep import cross_val_data




def important_gene_mask(columns, coefs):
    """
    function for finding features that are not zeroed by rf or xgb
    """
    mask = coefs != 0
    important_genes = columns[mask]
    return important_genes


def kfold_gene_intersection(columns, model_weights):
    """
    function for finding union of important feature of kfolds

    inputs
    ------
    columns: the features the model is predicting on
    model_weights: from model fit after kfolds

    outputs:
    intersection aka important features 
    """
    features1 = set(important_gene_mask(columns, model_weights[0]))
    features2 = set(important_gene_mask(columns, model_weights[1]))
    features3 = set(important_gene_mask(columns, model_weights[2]))
    features4 = set(important_gene_mask(columns, model_weights[3]))
    features5 = set(important_gene_mask(columns, model_weights[4]))

    total_feature_union = set.union(features1, features2, features3, features4, features5)
    total_feature_intersection = set.intersection(features1, features2, features3, features4, features5)

    return total_feature_union, total_feature_intersection




def kfold_feature_finder(csv, model_list):
    """
    a function that:
    1) creates x and y
    2) for each model, creates a 5-fold split
    3) returns a dictionary of these features the model found important

    inputs
    -------
    csv: string path to csv
    model_list: list of called models [xgb, rf, lasso]

    outputs
    -------
    """
    X, y, X_val, y_val = cross_val_data(csv)
    columns = X.columns

    kf = KFold(n_splits=5, shuffle=True)

    feature_dict = {}
    for model in tqdm(model_list):
        ll_performance = []
        model_weights = []
        for train_index, test_index in kf.split(X):
            X_train, X_test = X.iloc[train_index], X.iloc[test_index]
            y_train, y_test = y.iloc[train_index], y.iloc[test_index]
            model.fit(X_train, y_train)
            y_pred = model.predict_proba(X_test)
            log_ll = log_loss(y_test, y_pred)

            ll_performance.append(log_ll)
            if str(model).split('(')[0] == 'LogisticRegression':
                model_weights.append(model.coef_[0])
            else:
                model_weights.append(np.array(model.feature_importances_)) #turn into array for xgb and rando forest
        feature_union, feature_intersection = kfold_gene_intersection(columns, model_weights)
        feature_dict[str(model).split('(')[0]] = {'Union': feature_union, 'Intersection': feature_intersection}
    return feature_dict