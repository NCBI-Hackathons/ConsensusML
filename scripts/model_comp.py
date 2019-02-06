from sklearn.metrics import log_loss
from sklearn.model_selection import train_test_split
from sklearn.ensemble import  RandomForestClassifier, GradientBoostingClassifier
from xgboost import XGBClassifier

# def model_comp(df):
#     xgboost_model = XGBClassifier(learning_rate = 0.01, max_depth = 3, n_estimators = 300, random_state=8)
#     gradient_boost_model = GradientBoostingClassifier(learning_rate=0.01, max_depth=4, max_features='log2', min_samples_leaf=4, n_estimators=280, subsample=0.25, random_state=8)
#     random_forest_model = RandomForestClassifier(n_estimators=300, max_depth=3, verbose=1, random_state=8)

#     xgboost_model.fit(X_train, y_train)
#     gradient_boost_model.fit(X_train, y_train)
#     random_forest_model.fit(X_train, y_train)

#     p_random_forest = random_forest_model.predict_proba(X_test)
#     p_gradient_boost =  gradient_boost_model.predict_proba(X_test)
#     p_xgboost = xgboost_model.predict_proba(X_test)

#     #need to concat these two arrays somehow...
#     ensemble_p_0 = (p_random_forest[:,0] + p_gradient_boost[:,0] + p_xgboost[:,0])/3
#     ensemble_p_1 = (p_random_forest[:,1] + p_gradient_boost[:,1] + p_xgboost[:,1])/3
#     # ensemble_p = np.concatenate(ensemble_p_0, ensemble_p_1, axis=0)

#     random_forest_ll = log_loss(y_test, p_random_forest)
#     gradient_boost_ll = log_loss(y_test, p_gradient_boost)
#     xgboost_ll = log_loss(y_test, p_xgboost)
#     ensemble_ll = log_loss(y_test, ensemble_p_1)
                                        

#     print("Ensemble Log Loss " + str(ensemble_ll))
#     print("Gradient Boost Log Loss " + str(gradient_boost_ll))
#     print("Random Forest Log Loss " + str(random_forest_ll))


def data_prep_columns(df, var):
    """

    Parameters
    ----------------
    df: dataframe of patients
    var: str either 'Max' or 'Min', sets whether to take max variance or minimum variance

    Returns
    ----------------
    columns of either most/least variance from the two subgroups
    """
    if var =='Max':
        reverse_ = True
    else:
        reverse_ = False

    df1 = df[df['Diagnostic ID'].isin(('09A', '03A', '01A'))]
    df2 = df1[df1['Risk group'] != 'Unknown'].copy()
    df2['Low Risk'] = df2['Risk group'].apply(lambda x: (x == 'Low') * 1)

    low_df = df2[df2['Low Risk'] == 1].copy()
    ldf = low_df[low_df.columns[84:-4]].copy()
    ldf_var = list(enumerate(ldf.var(axis=0)))
    highest_var_low = sorted(ldf_var, key=lambda x:x[1], reverse=reverse_)[:1000]

    low_columns_to_take = [indx for indx, val in highest_var_low]
    low_columns_to_take_full = [item + 84 for item in low_columns_to_take]

    high_df = df2[df2['Low Risk'] == 0].copy()
    hdf = high_df[high_df.columns[84:-4]].copy()

    hdf_var = list(enumerate(hdf.var(axis=0)))
    highest_var_high = sorted(hdf_var, key=lambda x:x[1], reverse=reverse_)[:1000]

    high_columns_to_take = [indx for indx, val in highest_var_high]
    high_columns_to_take_full = [item + 84 for item in high_columns_to_take]

    low_columns_to_take_full.extend(x for x in high_columns_to_take_full if x not in low_columns_to_take_full)
    return low_columns_to_take_full

def model_prep(df, columns_to_take):
    data = df.iloc[:, columns_to_take]
    y = df.loc[data.index, 'Low Risk']
    data['label'] = y.copy()

    holdout = data.sample(frac=.2, random_state=8)
    train = data.drop(holdout.index).copy()
    X = train.iloc[:, :-1]
    y = train['label']

    X_train, X_test, y_train, y_test = train_test_split(X.values, y, test_size=0.33, random_state=8)
    return X_train, X_test, y_train, y_test, holdout