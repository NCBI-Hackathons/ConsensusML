def manifest_clinical_merge(manifest_df, clinical_df, target):
    """
    inputs:
    manifest_df - dataframe with all study metadata
    clinical_df - dataframe specific for disease
    target - string ex 'TARGET-WT'

    output:
    returns combined dataframe by merging on TARGET USI (ex TARGET-50-CAAAAA)

    ex: WT_df = manifest_clinical_merge(manifest_df, wt_disc_df, 'TARGET-WT')
    """
    target_df = manifest_df[manifest_df['project.project_id'] == target]
    target_df['TARGET USI'] = target_df.loc[:, 'entity_submitter_id'].apply(lambda x: x[:16])
    final_df = clinical_df.merge(target_df, on='TARGET USI')
    return final_df


def assay_transpose(assay_df):
    """
    input assay_df
    transposes df so that rows are patient assays
    """
    assay_pt_df = assay_df.T
    assay_pt_df.columns = assay_pt_df.iloc[0]
    assay_pt_df = assay_pt_df.iloc[1:]
    assay_pt_df['TARGET USI'] = list(assay_pt_df.index)
    assay_pt_df['TARGET USI'] = assay_pt_df['TARGET USI'].apply(lambda x: x[0:16])
    assay_pt_df['2nd ID'] = list(assay_pt_df.index)
    assay_pt_df['2nd ID'] = assay_pt_df['2nd ID'].apply(lambda x: x[16:])
    return assay_pt_df

def assay_clinical_merge(assay_df, clinical_df):
    """
    input assay df, clinical df, and disease (AML, WT, or NBL)
    return merged data
    """
    clinical_df['TARGET USI'] = clinical_df['TARGET USI'].apply(lambda x: x.replace('-', '.'))
    df = clinical_df.merge(assay_df, how='left', on='TARGET USI')
    return df
