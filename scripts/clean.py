def manifest_clinical_merge(manifest_df, clinical_df, target):
    """
    AML_df = manifest_clinical_merge(manifest_df, aml_disc_df, 'TARGET-AML')

    Parameters
    ----------------
    manifest_df: dataframe of metadata of study
    clinical_df: dataframe specific for disease with patient as rows
    target: string of target of disease 

    Returns
    ----------------
    dataframe transposed with patients as rows and genes as columns 
    """
    target_df = manifest_df[manifest_df['project.project_id'] == target]
    target_df['TARGET USI'] = target_df.loc[:, 'entity_submitter_id'].apply(lambda x: x[:16])
    final_df = clinical_df.merge(target_df, on='TARGET USI')
    return final_df


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
    # assay_pt_df['TARGET USI'] = list(assay_pt_df.index)
    # assay_pt_df['TARGET USI'] = assay_pt_df['TARGET USI'].apply(lambda x: x[0:16])
    # assay_pt_df['2nd ID'] = list(assay_pt_df.index)
    # assay_pt_df['2nd ID'] = assay_pt_df['2nd ID'].apply(lambda x: x[16:])
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
    clinical_df['entity_submitter_id'] = clinical_df['entity_submitter_id'].apply(lambda x: x.replace('-', '.'))
    # clinical_df['TARGET USI'] = clinical_df['TARGET USI'].apply(lambda x: x.replace('-', '.'))
    # df = clinical_df.merge(assay_df, how='left', on='TARGET USI')
    df = clinical_df.merge(assay_df, how='left', on='entity_submitter_id')
    df['Diagnostic ID'] = df.loc[:, 'entity_submitter_id'].apply(lambda x: x[-7:-4])
    columns = column_reorder(df)
    df = df[columns]

    return df
