from sklearn.metrics import log_loss
from sklearn.model_selection import train_test_split
from sklearn.ensemble import  RandomForestClassifier, GradientBoostingClassifier
from xgboost import XGBClassifier

def model_comp(df):
    xgboost_model = XGBClassifier(learning_rate = 0.01, max_depth = 3, n_estimators = 300, random_state=8)
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