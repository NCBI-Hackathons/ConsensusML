import pandas as pd
import matplotlib.pyplot as plt

def histogram_plot(csv, model_type):
    """
    Function to plot histogram of permutation test

    inputs:
    -------
    csv: string csv address
    model_type: str (Random Forest or XGBoost)

    returns:
    --------
    visualization 
    """
    df = pd.read_csv(csv)
    df.drop(columns=['Unnamed: 0'], inplace=True)
    df.columns = ['real score', 'perm scores', 'p-value']
    df['perm scores'] = df['perm scores'].apply(lambda x: x * -1.0)
    df['real score'] = df['real score'].apply(lambda x: x * -1.0)
    actual_score = round(df['real score'].mean(), 3)
    actual_score_marker = actual_score + 0.1
    p_value = round(df['p-value'].mean(), 4)
    

    plt.figure(figsize=(10, 6))
    plt.hist(df['perm scores'], bins=50)
    plt.axvline(x = actual_score, color='r', )
    plt.title(f'Permutation Test {model_type} Log Loss vs Actual Log Loss')
    plt.xlabel('Log Loss')
    plt.ylabel('Occurences')
    plt.text(y =50, x = 0.4, s = f'P-value: {p_value}')
    plt.text(y = 100, x = actual_score_marker, s = f'Actual Log Loss: {actual_score}' );