from sklearn.metrics import log_loss
from sklearn.model_selection import train_test_split
from sklearn.ensemble import  RandomForestClassifier, GradientBoostingClassifier
from xgboost import XGBClassifier

def model_comp(X_train, X_test, y_train, y_test):
    xgboost_model = XGBClassifier(learning_rate = 0.01, max_depth = 3, n_estimators = 700, random_state=8)
    gradient_boost_model = GradientBoostingClassifier(learning_rate=0.01, max_depth=4, max_features='log2', min_samples_leaf=4, n_estimators=280, subsample=0.25, random_state=8)
    random_forest_model = RandomForestClassifier(n_estimators=300, max_depth=3, verbose=1, random_state=8)

    xgboost_model.fit(X_train, y_train)
    gradient_boost_model.fit(X_train, y_train)
    random_forest_model.fit(X_train, y_train)

    p_random_forest = random_forest_model.predict_proba(X_test)
    p_gradient_boost =  gradient_boost_model.predict_proba(X_test)
    p_xgboost = xgboost_model.predict_proba(X_test)

    #need to concat these two arrays somehow...
    ensemble_p_0 = (p_random_forest[:,0] + p_gradient_boost[:,0] + p_xgboost[:,0])/3
    ensemble_p_1 = (p_random_forest[:,1] + p_gradient_boost[:,1] + p_xgboost[:,1])/3
    # ensemble_p = np.concatenate(ensemble_p_0, ensemble_p_1, axis=0)

    random_forest_ll = log_loss(y_test, p_random_forest)
    gradient_boost_ll = log_loss(y_test, p_gradient_boost)
    xgboost_ll = log_loss(y_test, p_xgboost)
    ensemble_ll = log_loss(y_test, ensemble_p_1)
                                        

    print("Ensemble Log Loss " + str(ensemble_ll))
    print("Gradient Boost Log Loss " + str(gradient_boost_ll))
    print("Random Forest Log Loss " + str(random_forest_ll))
    print("XGBoost Log Loss " + str(xgboost_ll))
    return xgboost_model, random_forest_model


def data_prep_columns(df, var):
    """

    Parameters
    ----------------
    df: dataframe of patients
    var: str either 'Max' or 'Min', sets whether to take top 1000 columns of
    max variance or minimum variance

    Returns
    ----------------
    df2: transformed dataframe with 'Diagostic ID' and no 'Unknown' values for risk_group
    new column 'Low Risk' with boolean value 0 or 1; THIS COLUMN HAS ALL COLUMNS (clinical data, manifest data)

    low_columns_to_take_full: columns (genes) of either most/least variance from the two subgroups
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
    return df2, low_columns_to_take_full

def model_prep(df, columns_to_take):
    """

    Parameters
    ----------------
    df: - takes transformed dataframe from 'data_prep_columns'

    columns_to_take: - columns you want to input (can be from 'data_prep_columns')
    or chosen in a different way

    Returns
    ----------------
    X_train, X_test, y_train, y_test as expected from train_test_split

    X: is a dataframe with only genes as columns
    """
    data = df.iloc[:, columns_to_take]
    y = df.loc[data.index, 'Low Risk']
    data['label'] = y.copy()

    holdout = data.sample(frac=.2, random_state=8)
    train = data.drop(holdout.index).copy()
    X = train.iloc[:, :-1]
    y = train['label']

    X_train, X_test, y_train, y_test = train_test_split(X.values, y, test_size=0.33, random_state=8)
    return X_train, X_test, y_train, y_test, holdout, X

def model_prep_loc(df, columns_to_take):
    data = df.loc[:, columns_to_take]
    y = df.loc[data.index, 'Low Risk']
    data['label'] = y.copy()

    holdout = data.sample(frac=.2, random_state=8)
    train = data.drop(holdout.index).copy()
    X = train.iloc[:, :-1]
    y = train['label']

    X_train, X_test, y_train, y_test = train_test_split(X.values, y, test_size=0.33, random_state=8)
    return X_train, X_test, y_train, y_test, holdout, X

def xgboost_tuner(X_train, X_test, y_train, y_test, n_estimators_num):
    """

    Parameters
    ----------------
    X_train, X_test, y_train, y_test derived form train_test_split
    n_estimators_num is a list [3, 4, 5, 6, 7] of len(5)

    Returns
    ----------------
    xg boost tuning performance
    """
    trees = [100 * x for x in n_estimators_num]
    a, b, c, d, e = trees

    xgb1 = XGBClassifier(learning_rate=0.01, n_estimators=a, random_state=8)
    xgb2 = XGBClassifier(learning_rate=0.01, n_estimators=b, random_state=8)
    xgb3 = XGBClassifier(learning_rate=0.01, n_estimators=c, random_state=8)
    xgb4 = XGBClassifier(learning_rate=0.01, n_estimators=d, random_state=8)
    xgb5 = XGBClassifier(learning_rate=0.01, n_estimators=e, random_state=8)

    xgb1.fit(X_train, y_train)
    xgb2.fit(X_train, y_train)
    xgb3.fit(X_train, y_train)
    xgb4.fit(X_train, y_train)
    xgb5.fit(X_train, y_train)


    xgb_p1 =  xgb1.predict_proba(X_test)
    xgb_p2 =  xgb2.predict_proba(X_test)
    xgb_p3 =  xgb3.predict_proba(X_test)
    xgb_p4 =  xgb4.predict_proba(X_test)
    xgb_p5 =  xgb5.predict_proba(X_test)



    xgb1_ll = log_loss(y_test, xgb_p1)
    xgb2_ll = log_loss(y_test, xgb_p2)
    xgb3_ll = log_loss(y_test, xgb_p3)
    xgb4_ll = log_loss(y_test, xgb_p4)
    xgb5_ll = log_loss(y_test, xgb_p5)


    print(f"XGB n_estimators: {a} log loss " + str(xgb1_ll))
    print(f"XGB n_estimators: {b} log loss " + str(xgb2_ll))
    print(f"XGB n_estimators: {c} log loss " + str(xgb3_ll))
    print(f"XGB n_estimators: {d} log loss " + str(xgb4_ll))
    print(f"XGB n_estimators: {e} log loss " + str(xgb5_ll))

