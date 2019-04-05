import pandas as pd
import matplotlib.pyplot as plt

def histogram_plot(csv, model_type, loss='log loss'):
    """
    Function to plot histogram of permutation test

    inputs:
    -------
    csv: string csv address
    model_type: str (Random Forest or XGBoost)
    loss: str, type of loss (log loss or accuracy)

    returns:
    --------
    visualization 
    """

    df = pd.read_csv(csv)
    df.drop(columns=['Unnamed: 0'], inplace=True)
    df.columns = ['real score', 'perm scores', 'p-value']

    if loss == 'log loss':
        loss = 'Log Loss'
        df['perm scores'] = df['perm scores'].apply(lambda x: x * -1.0)
        df['real score'] = df['real score'].apply(lambda x: x * -1.0)
        actual_score = round(df['real score'].mean(), 3)
        actual_score_marker = actual_score + 0.1
        y_ = 50
        x_ = 0.4
        y1 = 100
        x1 = actual_score_marker
        
    if loss == 'accuracy':
        loss = 'Accuracy'
        y_ = 50
        x_ = 0.75
        y1 = 100
        x1 = 0.75
        
    actual_score = round(df['real score'].mean(), 3)
    actual_score_marker = actual_score + 0.1
    p_value = round(df['p-value'].mean(), 4)
    

    plt.figure(figsize=(10, 6))
    plt.hist(df['perm scores'], bins=50)
    plt.axvline(x = actual_score, color='r', )
    plt.title(f'Permutation Test {model_type} {loss} vs Actual {loss}')
    plt.xlabel(f'{loss}')
    plt.ylabel('Occurences')
    plt.text(y = y_, x = x_, s = f'P-value: {p_value}')
    plt.text(y = y1, x = x1, s = f'Actual {loss}: {actual_score}' );