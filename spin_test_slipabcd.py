
# %%
from neuromaps import  nulls, stats
import neuromaps as neuromaps
from nilearn.datasets import fetch_atlas_surf_destrieux
import abagen
import pandas as pd
import nibabel as nib

# %%
def process_beta_dataframe(df, df_labels_lh, df_labels_rh, beta_col, hemi_col='hemisphere', label_col='label'):
    """
    Process a DataFrame of beta values by hemisphere and label mapping,
    returning a 1D numpy array of ordered beta values.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing at least hemisphere, label, and beta columns.
    labels_lh : dict
        Mapping of left hemisphere labels to keys.
    labels_rh : dict
        Mapping of right hemisphere labels to keys.
    beta_col : str
        Name of the column containing beta values. Default 'bw_beta'.
    hemi_col : str
        Name of the column containing hemisphere info. Default 'hemisphere'.
    label_col : str
        Name of the column containing label info. Default 'label'.

    Returns
    -------
    np.ndarray
        Ordered 1D array of beta values combining left and right hemispheres.
    """
    
    # Split by hemisphere
    df_left = df[df[hemi_col] == 'L'].copy()
    df_right = df[df[hemi_col] == 'R'].copy()
    
    
    # Merge label-Key mapping
    df_left = pd.merge(df_left, df_labels_lh, on=label_col)
    df_right = pd.merge(df_right, df_labels_rh, on=label_col)
    beta_col='bw_beta',
    # Sort by Key
    df_left = df_left.sort_values(by='Key')
    df_right = df_right.sort_values(by='Key')
    
    # Combine hemispheres
    df_ordered = pd.concat([df_left, df_right])
    
    # Return 1D vector of beta values
    return df_ordered[beta_col].values

# %%
atlas_files = abagen.fetch_desikan_killiany(surface = True)
atlas_lh = nib.load(atlas_files['image'][0])
atlas_rh = nib.load(atlas_files['image'][1])
labels_lh = atlas_lh.labeltable.get_labels_as_dict()
labels_rh = atlas_rh.labeltable.get_labels_as_dict()
# Convert label dicts into DataFrames
df_labels_lh = pd.DataFrame(list(labels_lh.items()), columns=['Key', label_col])
df_labels_rh = pd.DataFrame(list(labels_rh.items()), columns=['Key', label_col])

# %%
#read in the csvs
abcd_t1_bw_betas = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/ABCD-braincharts/2.0_results/dkt_betas/dsk_birthweight_all_t1.csv")
abcd_t1_ga_betas = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/ABCD-braincharts/2.0_results/dkt_betas/dsk_gestAge_all_t1.csv")
abcd_t2_bw_betas = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/ABCD-braincharts/2.0_results/dkt_betas/dsk_birthweight_all_t2.csv")
abcd_t2_ga_betas = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/ABCD-braincharts/2.0_results/dkt_betas/dsk_gestAge_all_t2.csv")
slip_bw_betas_vol = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-distinct-allgrades-25.11/results_csv/dkt_betas/bw.qc_volume.csv")
slip_ga_betas_vol = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-distinct-allgrades-25.11/results_csv/dkt_betas/ga.qc_volume.csv")
slip_bw_betas_sa = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-distinct-allgrades-25.11/results_csv/dkt_betas/bw.qc_surfArea.csv")
slip_ga_betas_sa = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-distinct-allgrades-25.11/results_csv/dkt_betas/ga.qc_surfArea.csv")
slip_bw_betas_th = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-distinct-allgrades-25.11/results_csv/dkt_betas/bw.qc_thick.csv")
slip_ga_betas_th = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-distinct-allgrades-25.11/results_csv/dkt_betas/ga.qc_thick.csv")
slip_peaks = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/centile_peaks.csv")

# %%
slip_peaks_vol = slip_peaks[slip_peaks['domain' == "vol"]].copy()
slip_peaks_sa = slip_peaks[slip_peaks['domain' == "sa"]].copy()
slip_peaks_th = slip_peaks[slip_peaks['domain' == "th"]].copy()

# %%
abcd_t1_bw_beta_vector = process_beta_dataframe(abcd_t1_bw_betas, df_labels_lh, df_labels_rh, 'bw_beta')
abcd_t2_bw_beta_vector = process_beta_dataframe(abcd_t2_bw_betas, df_labels_lh, df_labels_rh, 'bw_beta')
slip_bw_vol_beta_vector = process_beta_dataframe(slip_bw_betas_vol, df_labels_lh, df_labels_rh, 'bw_beta')
slip_bw_sa_beta_vector = process_beta_dataframe(slip_bw_betas_sa, df_labels_lh, df_labels_rh, 'bw_beta')
slip_bw_th_beta_vector = process_beta_dataframe(slip_bw_betas_th, df_labels_lh, df_labels_rh, 'bw_beta')
slip_peaks_vol_vector = process_beta_dataframe(slip_peaks_vol, df_labels_lh, df_labels_rh, 'rank_peak50')

# %%
#Generating indices of null spins instead of rotating data, so can apply indices to the map being rotated
nulls_idx = nulls.alexander_bloch(data=None, atlas='fsaverage', density='10k', parcellation=[atlas_lh, atlas_rh], n_perm=10000, seed=1234)

# %%
#Using the indices to rotate the abcd maps (will compare to slip)
rotated_abcd_t1_bw = abcd_t1_bw_beta_vector[nulls_idx]
rotated_abcd_t2_bw = abcd_t2_bw_beta_vector[nulls_idx]
rotated_abcd_t1_ga = abcd_t1_ga_beta_vector[nulls_idx]
rotated_abcd_t2_ga = abcd_t2_ga_beta_vector[nulls_idx]
rotated_slip_bw = slip_bw_beta_vector[nulls_idx]

# %%
#Test t1 bw and slip bw
abcd_t1_slip_bw_corr, abcd_t1_slip_bw_pval = stats.compare_images(abcd_t1_bw_beta_vector, slip_bw_beta_vector, nulls=rotated_abcd_t1_bw, metric='spearmanr')
print(abcd_t1_slip_bw_corr, abcd_t1_slip_bw_pval)
#Test t1 ga and slip ga
abcd_t1_slip_ga_corr, abcd_t1_slip_ga_pval = stats.compare_images(abcd_t1_ga_beta_vector, slip_ga_beta_vector, nulls=rotated_abcd_t1_ga, metric='spearmanr')
print(abcd_t1_slip_ga_corr, abcd_t1_slip_ga_pval)
#Test t2 bw and slip bw
abcd_t2_slip_bw_corr, abcd_t2_slip_bw_pval = stats.compare_images(abcd_t2_bw_beta_vector, slip_bw_beta_vector, nulls=rotated_abcd_t2_bw, metric='spearmanr')
print(abcd_t2_slip_bw_corr, abcd_t2_slip_bw_pval)
#Test t2 ga and slip ga
abcd_t2_slip_ga_corr, abcd_t2_slip_ga_pval = stats.compare_images(abcd_t2_ga_beta_vector, slip_ga_beta_vector, nulls=rotated_abcd_t2_ga, metric='spearmanr')
print(abcd_t2_slip_ga_corr, abcd_t2_slip_ga_pval)


#Test slip bw and slip ga
slip_bwga_corr, slip_bwga_pval = stats.compare_images(slip_bw_beta_vector, slip_ga_beta_vector, nulls=rotated_slip_bw, metric='spearmanr')
print(slip_bwga_corr, slip_bwga_pval)

#Test abct t1 bw and abcd t1 ga
abcd_t1_bwga_corr, abcd_t1_bwga_pval = stats.compare_images(abcd_t1_bw_beta_vector, abcd_t1_ga_beta_vector, nulls=rotated_abcd_t1_bw, metric='spearmanr')
print(abcd_t1_bwga_corr, abcd_t1_bwga_pval)

#Test abct t2 bw and abcd t2 ga
abcd_t2_bwga_corr, abcd_t2_bwga_pval = stats.compare_images(abcd_t2_bw_beta_vector, abcd_t2_ga_beta_vector, nulls=rotated_abcd_t2_bw, metric='spearmanr')
print(abcd_t2_bwga_corr, abcd_t2_bwga_pval)
