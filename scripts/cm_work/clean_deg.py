def assay_transpose(assay_df):
    """
    assay_t_df = assay_transpose(assay_df)

    Parameters
    ----------------
    assay_df: dataframe of assay data with genes as rows and patients as columns

    Returns
    ----------------
    dataframe transposed with patients as rows and genes as columns 
    """
    assay_pt_df = assay_df.T
    assay_pt_df.columns = assay_pt_df.iloc[0]
    assay_pt_df = assay_pt_df.iloc[1:]
    assay_pt_df['entity_submitter_id'] = list(assay_pt_df.index)
    return assay_pt_df

def column_reorder(clinical_df):
    columns = clinical_df.columns
    reorder_columns = [columns[0], columns[-1]]
    for column in columns[1:-1]:
        reorder_columns.append(column)
    return reorder_columns

def assay_clinical_merge(assay_df, clinical_df):
    """
    AML_genes = assay_clinical_merge(assay_t_df, AML_df)

    Parameters
    ----------------
    assay_df: dataframe of assay data with patients as rows and genes as columns
    clinical_df: dataframe of clinical + manifest data with patients as rows
    
    Returns
    ----------------
    merged dataframe on 'entity_submitted_id' with 'Diagnostic ID' as second columns
    """
    assay_df1 = assay_transpose(assay_df)
    df = clinical_df.merge(assay_df1, how='left', on='entity_submitter_id')
    df['Low Risk'] = df['Risk.group'].apply(lambda x: x == 'Low') * 1.0
    columns = ['TARGET.USI']
    columns.extend(assay_df1.columns[:-1])
    columns.append('Low Risk')
    df = df[columns]
    return df
