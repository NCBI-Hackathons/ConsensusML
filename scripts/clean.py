def manifest_clinical_merge(manifest_df, clinical_df, target):
    """
    manifest_df - dataframe with all study metadata
    clinical_df - dataframe specific for disease
    target - string ex 'TARGET-WT'
    returns combined dataframe by merging on TARGET USI (ex TARGET-50-CAAAAA)
    """
    target_df = manifest_df[manifest_df['project.project_id'] == target]
    target_df['TARGET USI'] = target_df.loc[:, 'entity_submitter_id'].apply(lambda x: x[:16])
    final_df = clinical_df.merge(target_df, on='TARGET USI')
    return final_df

    