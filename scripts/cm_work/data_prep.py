import pandas as pd
import numpy as np

from sklearn.metrics import log_loss, accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.ensemble import  RandomForestClassifier
from xgboost import XGBClassifier
import random
from sklearn.model_selection import permutation_test_score
import datetime

def create_x_y_data():
    """
    function for creating pre-chosen x_train and y_train
    for permutation tests
    """
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

    return X_train, y_train, X_test, y_test, X, y

def cross_val_data(csv):
    """
    function for preparing data for kfold validation

    inputs
    ------
    csv: csv path as string

    returns
    -------
    X: dataframe to train on
    y: labels
    final_X: validation dataframe
    final_y: validation labels
    """
    df = pd.read_csv('deg.csv')
    df.drop(columns=['Unnamed: 0'], inplace=True)
    columns = df.columns[1:-1]

    final_validation = df.sample(frac=0.1)

    final_X = final_validation[columns]
    final_y = final_validation.iloc[:, -1]

    data = df.drop(index= final_validation.index)
    X = data[columns]
    y = data.iloc[:, -1]

    return X, y, final_X, final_y