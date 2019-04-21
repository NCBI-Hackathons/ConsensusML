import pandas as pd
import numpy as np
from sklearn.metrics import log_loss
from sklearn.ensemble import  RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split, KFold
from boruta import BorutaPy





def model_fit_score(model, X_train, X_test, y_train, y_test):
    """
    returns model log loss
    """
    model.fit(X_train, y_train)
    y_pred = model.predict_proba(X_test)
    model_ll = log_loss(y_test, y_pred)
    return model_ll

def important_gene_mask(columns, coefs):
    """
    inputs
    ------
    columns: columns of df
    coefs: beta weights of lasso
    
    results
    ------
    important_genes: name of genes with weight != 0
    gene_weights: beta weights of genes
    """
    
    mask = coefs[0] != 0
    
    gene_weights = coefs[0][mask]
    important_genes = columns[mask]

    return dict(zip(important_genes, gene_weights))

def important_gene_mask_tree(columns, coefs):
    """
    gene finder for tree based models since coef_ and feature_importances
    work differently.
    
    inputs
    ------
    columns: columns of df
    coefs: beta weights of lasso
    
    results
    ------
    important_genes: name of genes with weight != 0
    gene_weights: beta weights of genes
    """
    
    mask = coefs != 0
    
    gene_weights = coefs[mask]
    important_genes = columns[mask]

    return dict(zip(important_genes, gene_weights))


class ConsensusML():
    def __init__(self, lasso, xgb, rf, X, y, param=None):
        self.param = param
        self.lasso = lasso
        self.xgb = xgb
        self.rf = rf
        self.X = X
        self.y = y
        self.model_dict = None # dictionary from feat_selection3
        self.columns = None #feature columns aka genes
        self.gene_intersection = None # intersection of Lasso, XGB and RF features
        self.kfold_performances = None
        self.kfold_weights_list = None
        self.lasso_kf = None
        self.xgb_kf = None
        self.rf_kf = None
        self.rf_boruta_features = None
        self.total_consensus = None
        
    
    def feat_selection3(self, X_train, X_test, y_train, y_test):
        """
        fit the three models on X_train and y_train

        return dictionary with:
        keys: Lasso, XGB, RF
        values: {Log Loss': , 'Genes':, 'Total Genes':}
        """
        
        self.columns = X_train.columns #columns are features here

        model_names = ['Lasso', 'XGB', 'RF']
        model_dict = {}       #container of dictionary
        i = 0
        for model in [self.lasso, self.xgb, self.rf]:
            name = model_names[i] #create key for dictionary
            log_loss = model_fit_score(model, X_train, X_test, y_train, y_test)
            if model == self.lasso:
                mask = model.coef_ != 0
                lasso_columns = self.columns[mask[0]]
                model_dict[name] = {'Log Loss': log_loss, 'Genes': lasso_columns, 'Total Genes': len(lasso_columns)}
            elif model == self.xgb or model == self.rf:
                mask = model.feature_importances_ != 0
                tree_columns = self.columns[mask]
                model_dict[name] = {'Log Loss': log_loss, 'Genes': tree_columns, 'Total Genes': len(tree_columns)}
            i += 1
        self.model_dict = model_dict
        return self.model_dict

    def feature_intersection_weights(self):
        """
        function that creates a set that is the intersection fo the three models

        returns intersection and weights from lasso and xgb

        outputs
        ------
        gene_intersection: set of feature importance
        lasso_weights: weights of those featutres
        xgb_feature_importance: weights of those features
        """
        self.gene_intersection = set.intersection(set(self.model_dict['Lasso']['Genes']), set(self.model_dict['XGB']['Genes']),
                                            set(self.model_dict['RF']['Genes'])   )
        intersection_mask = [x in self.gene_intersection for x in self.columns]

        lasso_weights = self.lasso.coef_[0][intersection_mask]
        xgb_feature_importance = self.xgb.feature_importances_[intersection_mask]
        return self.gene_intersection, lasso_weights, xgb_feature_importance

    def kfold_tune(self, param_list):
        """
        input X and y plus param_list
        param_list: list of dictionaries with params to try GridSearch on
        default param_list [ 
            {'C': range(1, 4, 1)},
            {'max_depth': range(3, 10)},
            {'max_depth': range(3, 7), 'n_estimators': range(200, 800, 200)}
        ]
        lasso --> xgb --> rf
        """

        models = [self.lasso, self.xgb, self.rf]
        for i in range(3):
            params = param_list[i]
            model = models[i]
            if model == self.lasso:
                grid_search = GridSearchCV(model, params, cv=5, n_jobs=-1, scoring='neg_log_loss', verbose=True)
                grid_search.fit(self.X, self.y)
                self.lasso = LogisticRegression(penalty='l1', solver='saga', max_iter=10000, **grid_search.best_params_)
            if model == self.rf:
                grid_search = GridSearchCV(model, params, cv=5, n_jobs=-1, scoring='neg_log_loss', verbose=True)
                grid_search.fit(self.X, self.y)
                self.rf = RandomForestClassifier(n_estimators=1000, n_jobs=-1, **grid_search.best_params_)
            if model == self.xgb:
                grid_search = GridSearchCV(model, params,
                                        cv=5, n_jobs=-1, scoring='neg_log_loss', verbose=True)
                grid_search.fit(self.X, self.y)
                self.xgb = XGBClassifier(n_jobs=-1, **grid_search.best_params_)
    
    def kfold_weights(self):
        """
        Runs tuned model on X and y with kfold validation.

        Outputs
        -------
        kfold_weights: average weight of each feature over kfold for each model
        kfold_perofmrnaces: average performance of each model over kfolds
        """

        lasso_performance = []
        rf_ll_performance = []
        xgb_ll_performance = []

        lasso_weights = []
        rf_weights = []
        xgb_weights = []

        kf = KFold(n_splits=5, shuffle=True)
        for train_index, test_index in kf.split(self.X):
            X_train, X_test = self.X.iloc[train_index], self.X.iloc[test_index]
            y_train, y_test = self.y.iloc[train_index], self.y.iloc[test_index]
            
            self.lasso.fit(X_train, y_train)
            self.rf.fit(X_train, y_train)
            self.xgb.fit(X_train, y_train)
            
            p_lr = self.lasso.predict_proba(X_test)
            p_rf = self.rf.predict_proba(X_test)
            p_xgb = self.xgb.predict_proba(X_test)
            
            log_ll = log_loss(y_test, p_lr)
            rf_ll = log_loss(y_test, p_rf)
            xgb_ll = log_loss(y_test, p_xgb)
            
            lasso_performance.append(log_ll)
            rf_ll_performance.append(rf_ll)
            xgb_ll_performance.append(xgb_ll)
            
            lasso_weights.append(self.lasso.coef_)  
            rf_weights.append(self.rf.feature_importances_)
            xgb_weights.append(self.xgb.feature_importances_)
        
        self.kfold_performances = [lasso_performance, rf_ll_performance, xgb_ll_performance] #save log loss performances

        l1, l2, l3, l4, l5 = [important_gene_mask(self.columns, lasso_weights[i]) for i in range(5)]
        rf1, rf2, rf3, rf4, rf5 = [important_gene_mask_tree(self.columns, rf_weights[i]) for i in range(5)]
        xgb1, xgb2, xgb3, xgb4, xgb5 = [important_gene_mask_tree(self.columns, xgb_weights[i]) for i in range(5)]

        lasso_int = set.intersection(set(l1), set(l2), set(l3), set(l4), set(l5))
        rf_int = set.intersection(set(rf1), set(rf2), set(rf3), set(rf4), set(rf5))
        xgb_int = set.intersection(set(xgb1), set(xgb2), set(xgb3), set(xgb4), set(xgb5))

        lasso_weight = {}
        for gene in lasso_int:
            lasso_weight[gene] = l1[gene]
            lasso_weight[gene] += l2[gene]
            lasso_weight[gene] += l3[gene]
            lasso_weight[gene] += l4[gene]
            lasso_weight[gene] += l5[gene]
            lasso_weight[gene] = lasso_weight[gene] / 5
        rf_weight = {}
        for gene in rf_int:
            rf_weight[gene] = rf1[gene]
            rf_weight[gene] += rf2[gene]
            rf_weight[gene] += rf3[gene]
            rf_weight[gene] += rf4[gene]
            rf_weight[gene] += rf5[gene]
            rf_weight[gene] = rf_weight[gene] / 5
        xgb_weight = {}
        for gene in xgb_int:
            xgb_weight[gene] = xgb1[gene]
            xgb_weight[gene] += xgb2[gene]
            xgb_weight[gene] += xgb3[gene]
            xgb_weight[gene] += xgb4[gene]
            xgb_weight[gene] += xgb5[gene]
            xgb_weight[gene] = xgb_weight[gene] / 5
        
        self.kfold_weights_list = [lasso_weight, rf_weight, xgb_weight]
        
        return self.kfold_performances, self.kfold_weights_list

    def boruta_selection(self):
        """feature selection via boruta"""
        rf_boruta = RandomForestClassifier(n_jobs=-1)
        feat_selector = BorutaPy(rf_boruta, n_estimators='auto', verbose=2, max_iter = 100, random_state=8)
        feat_selector.fit(self.X.values, self.y.values)

        selected = self.X.values[:, feat_selector.support_]
        print(selected.shape)

        # get the name of columns that boruta thinks is important
        boruta_mask = feat_selector.support_
        self.rf_boruta_features = self.columns[boruta_mask]

    def feature_consensus(self):
        """
        create intersections and sets with all gene features
        """
        self.gene_intersection # original intersection of ensemble
        self.rf_boruta_features #boruta features
        self.lasso_kf = self.kfold_weights_list[0] #lasso kfold features
        self.rf_kf = self.kfold_weights_list[1] # rf kfold features
        self.xgb_kf = self.kfold_weights_list[2] #xgb kfold features

        self.total_consensus = set.intersection(set(self.gene_intersection), set(self.rf_boruta_features), set(self.lasso_kf),
                                            set(self.rf_kf), set(self.xgb_kf) )
        return self.total_consensus


    def best_combo(self, model, lasso=False):
        """
        compare kfold scores between different combinations of
        1) boruta
        2) boruta + gene_intersection
        3) boruta + kf_lasso
        4) boruta + kf_lasso + xgb_lasso
        5) gene_intersection + kf_lasso + xgb_lasso
        """

        boruta_gene = list(set.union(self.gene_intersection, set(self.rf_boruta_features)))  #2
        boruta_lasso = list(set.union(set(self.rf_boruta_features), set(self.lasso_kf))) #3
        boruta_lasso_xgb = list(set.union(set(self.rf_boruta_features), set(self.lasso_kf), set(self.xgb_kf))) #4
        gene_lasso_xgb = list(set.union(set(self.gene_intersection), set(self.lasso_kf), set(self.xgb_kf))) #5

        boruta_df = self.X[self.rf_boruta_features]     #1
        boruta_gene_df = self.X[boruta_gene]            #2
        boruta_lasso_df = self.X[boruta_lasso]          #3
        boruta_lasso_xgb_df = self.X[boruta_lasso_xgb]  #4
        gene_lasso_xgb_df = self.X[gene_lasso_xgb]      #5

        df_list = [boruta_df, boruta_gene_df, boruta_lasso_df, boruta_lasso_xgb_df, gene_lasso_xgb_df]
        df_names = ['Boruta', 'Boruta + Intersection', 'Boruta + Lasso KF', 'Boruta + Lasso + XGB KF',
                    'Intersection + Lasso + XGB KF']

        i = 0
        feature_performance = {}
        for df in df_list:
            X_train, X_test, y_train, y_test = train_test_split(df, self.y)
            log_loss_score = model_fit_score(model, X_train, X_test, y_train, y_test)
            combo_name = df_names[i]
            if lasso:
                feature_weights = model.coef_
            else:
                feature_weights = model.feature_importances_
            feature_performance[combo_name] = {'Log Loss': log_loss_score,
                                            'Model Weights': feature_weights}

            i +=1
        print([(key, feature_performance[key]['Log Loss']) for key in feature_performance.keys()])
        
        best_feature_pos = np.argmin([feature_performance[key]['Log Loss'] for key in feature_performance.keys()])
        best_df = df_list[best_feature_pos]

        return feature_performance, best_df

