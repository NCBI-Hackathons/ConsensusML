from model_blender import important_gene_mask
from sklearn.metrics import log_loss
import numpy as np

def gene_weight_finder(model, X_train, X_test, y_train, y_test):
    """
    function that returns the most important features, weights and # of features

    inputs
    -------
    model: tree based model
    X_train:
    X_test:
    y_train:
    y_test

    outputs
    -------
    all_important: list of all genes with feature importance > 0
    
    top_20_feature_names: top 20 most important column names (gene)
    based on feature importance

    top_20_weights: weights of the top 20 columns

    num_feats: number of features that are not 0

    number_important: number of features with feature importance > 0

    log_loss: log loss score
    """
    columns = X_train.columns
    model.fit(X_train, y_train)
    y_pred = model.predict_proba(X_test)
    ll = log_loss(y_test, y_pred)
    
    top_20_features = np.argsort(model.feature_importances_)[-20:][::-1]
    top_20_feature_names = columns[top_20_features]
    top_20_weights = model.feature_importances_[top_20_features]

    number_important = len(important_gene_mask(columns, model.feature_importances_))
    all_important = important_gene_mask(columns, model.feature_importances_)

    return all_important, top_20_feature_names, top_20_weights, number_important, ll




