
# %%
import sys
print(sys.executable)
print(sys.version)

# %%
from neuromaps import  nulls, stats
import neuromaps as neuromaps
from nilearn.datasets import fetch_atlas_surf_destrieux
import abagen
import pandas as pd
import nibabel as nib
import numpy as np

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
    #beta_col='bw_beta',
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
df_labels_lh = pd.DataFrame(list(labels_lh.items()), columns=['Key','label'])
df_labels_rh = pd.DataFrame(list(labels_rh.items()), columns=['Key','label'])

# %%
#read in the csvs
slip_bw_betas_vol = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/analysis_results/median-th-mpr/results_csv/dkt_betas/bw.qc_volume.csv")
slip_bwp_betas_vol = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/analysis_results/median-th-mpr/results_csv/dkt_betas/bwp.qc_volume.csv")
slip_ga_betas_vol = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/analysis_results/median-th-mpr/results_csv/dkt_betas/ga.qc_volume.csv")
slip_bw_betas_sa = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/analysis_results/median-th-mpr/results_csv/dkt_betas/bw.qc_surfArea.csv")
slip_bwp_betas_sa = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/analysis_results/median-th-mpr/results_csv/dkt_betas/bwp.qc_surfArea.csv")
slip_ga_betas_sa = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/analysis_results/median-th-mpr/results_csv/dkt_betas/ga.qc_surfArea.csv")
slip_bw_betas_th = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/analysis_results/median-th-mpr/results_csv/dkt_betas/bw.qc_thick.csv")
slip_bwp_betas_th = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/analysis_results/median-th-mpr/results_csv/dkt_betas/bwp.qc_thick.csv")
slip_ga_betas_th = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/analysis_results/median-th-mpr/results_csv/dkt_betas/ga.qc_thick.csv")

#slip_bw_global_betas_vol = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-allgrades-25.11/results_csv/dkt_betas/bw.qc.global_volume.csv")
#slip_ga_global_betas_vol = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-allgrades-25.11/results_csv/dkt_betas/ga.qc.global_volume.csv")
#slip_bw_global_betas_sa = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-allgrades-25.11/results_csv/dkt_betas/bw.qc.global_surfArea.csv")
#lip_ga_global_betas_sa = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-allgrades-25.11/results_csv/dkt_betas/ga.qc.global_surfArea.csv")
#slip_bw_global_betas_th = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-allgrades-25.11/results_csv/dkt_betas/bw.qc.global_thick.csv")
#slip_ga_global_betas_th = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-allgrades-25.11/results_csv/dkt_betas/ga.qc.global_thick.csv")

# %%
#PREP SA AXIS RANKS FOR COMPARISON
#sa_axis_dsk_full = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/ABCD-braincharts/CSV/sydnor.dk.maps.csv")
#sa_axis_dsk_full['mapname'].unique()

# %%
# rename and recode hemisphere column
sa_axis_dsk_full = sa_axis_dsk_full.rename(columns={"hemi": "hemisphere"})
sa_axis_dsk_full["hemisphere"] = sa_axis_dsk_full["hemisphere"].replace({
    "left": "L",
    "right": "R"
})

# %% SA AXIS
sa_axis_dsk = sa_axis_dsk_full[sa_axis_dsk_full['mapname'] == 'Sensorimotor-Association Axis'].copy()
bloodflow_dsk = sa_axis_dsk_full[sa_axis_dsk_full['mapname'] == 'Cerebral Blood Flow'].copy()
glycolysis_dsk = sa_axis_dsk_full[sa_axis_dsk_full['mapname'] == 'Aerobic Glycosis'].copy()
gene_dsk = sa_axis_dsk_full[sa_axis_dsk_full['mapname'] == 'Gene Expression'].copy()
func_dsk = sa_axis_dsk_full[sa_axis_dsk_full['mapname'] == 'Functional Hierarchy'].copy()
anat_dsk = sa_axis_dsk_full[sa_axis_dsk_full['mapname'] == 'Anatomical Hierarchy'].copy()
allo_dsk = sa_axis_dsk_full[sa_axis_dsk_full['mapname'] == 'Allometric Scaling'].copy()#*****
evo_dsk = sa_axis_dsk_full[sa_axis_dsk_full['mapname'] == 'Evolutionary Expansion'].copy()

sa_axis_rank_vector = process_beta_dataframe(sa_axis_dsk, df_labels_lh, df_labels_rh, 'rank')
bloodflow_rank_vector = process_beta_dataframe(bloodflow_dsk, df_labels_lh, df_labels_rh, 'rank')
glycolysis_rank_vector = process_beta_dataframe(glycolysis_dsk, df_labels_lh, df_labels_rh, 'rank')
gene_rank_vector = process_beta_dataframe(gene_dsk, df_labels_lh, df_labels_rh, 'rank')
func_rank_vector = process_beta_dataframe(func_dsk, df_labels_lh, df_labels_rh, 'rank')
anat_rank_vector = process_beta_dataframe(anat_dsk, df_labels_lh, df_labels_rh, 'rank')
allo_rank_vector = process_beta_dataframe(allo_dsk, df_labels_lh, df_labels_rh, 'rank')
evo_rank_vector = process_beta_dataframe(evo_dsk, df_labels_lh, df_labels_rh, 'rank')

# %% Read SWYC betas
#read in the csvs
#slip_swyc_betas_vol = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/swyc_output/")
slip_swyc_betas_sa = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/swyc_output/swyc_median_surfArea_apr13.csv")
#slip_swyc_betas_th = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/swyc_output/")

# %% process beta frames
slip_bw_vol_beta_vector = process_beta_dataframe(slip_bw_betas_vol, df_labels_lh, df_labels_rh, 'bw_beta')
slip_bw_sa_beta_vector = process_beta_dataframe(slip_bw_betas_sa, df_labels_lh, df_labels_rh, 'bw_beta')
slip_bw_th_beta_vector = process_beta_dataframe(slip_bw_betas_th, df_labels_lh, df_labels_rh, 'bw_beta')
slip_bwp_vol_beta_vector = process_beta_dataframe(slip_bwp_betas_vol, df_labels_lh, df_labels_rh, 'bwp_beta')
slip_bwp_sa_beta_vector = process_beta_dataframe(slip_bwp_betas_sa, df_labels_lh, df_labels_rh, 'bwp_beta')
slip_bwp_th_beta_vector = process_beta_dataframe(slip_bwp_betas_th, df_labels_lh, df_labels_rh, 'bwp_beta')
slip_ga_vol_beta_vector = process_beta_dataframe(slip_ga_betas_vol, df_labels_lh, df_labels_rh, 'ga_beta')
slip_ga_sa_beta_vector = process_beta_dataframe(slip_ga_betas_sa, df_labels_lh, df_labels_rh, 'ga_beta')
slip_ga_th_beta_vector = process_beta_dataframe(slip_ga_betas_th, df_labels_lh, df_labels_rh, 'ga_beta')

# %% controlling for global centiles
slip_bw_global_vol_beta_vector = process_beta_dataframe(slip_bw_global_betas_vol, df_labels_lh, df_labels_rh, 'bw_beta')
slip_bw_global_sa_beta_vector = process_beta_dataframe(slip_bw_global_betas_sa, df_labels_lh, df_labels_rh, 'bw_beta')
slip_bw_global_th_beta_vector = process_beta_dataframe(slip_bw_global_betas_th, df_labels_lh, df_labels_rh, 'bw_beta')
slip_ga_global_vol_beta_vector = process_beta_dataframe(slip_ga_global_betas_vol, df_labels_lh, df_labels_rh, 'ga_beta')
slip_ga_global_sa_beta_vector = process_beta_dataframe(slip_ga_global_betas_sa, df_labels_lh, df_labels_rh, 'ga_beta')
slip_ga_global_th_beta_vector = process_beta_dataframe(slip_ga_global_betas_th, df_labels_lh, df_labels_rh, 'ga_beta')

# %% slip_swyc_vol_beta_vector = process_beta_dataframe(slip_swyc_betas_vol, df_labels_lh, df_labels_rh, 'swyc_beta')
slip_swyc_sa_beta_vector = process_beta_dataframe(slip_swyc_betas_sa, df_labels_lh, df_labels_rh, 'swyc_beta')
#slip_swyc_th_beta_vector = process_beta_dataframe(slip_swyc_betas_th, df_labels_lh, df_labels_rh, 'swyc_beta')

# %% generate and save nulls idx, 10k
#Generating indices of null spins instead of rotating data, so can apply indices to the map being rotated
nulls_idx = nulls.alexander_bloch(data=None, atlas='fsaverage', density='10k', parcellation=[atlas_lh, atlas_rh], n_perm=10000, seed=1234)
np.savetxt("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/dsk_nulls.csv", nulls_idx, delimiter=",")

# %% Rotate Maps
rotated_slip_bw_vol = slip_bw_vol_beta_vector[nulls_idx]
rotated_slip_bw_sa = slip_bw_sa_beta_vector[nulls_idx]
rotated_slip_bw_th = slip_bw_th_beta_vector[nulls_idx]
rotated_slip_bwp_vol = slip_bwp_vol_beta_vector[nulls_idx]
rotated_slip_bwp_sa = slip_bwp_sa_beta_vector[nulls_idx]
rotated_slip_bwp_th = slip_bwp_th_beta_vector[nulls_idx]
rotated_slip_ga_vol = slip_ga_vol_beta_vector[nulls_idx]
rotated_slip_ga_sa = slip_ga_sa_beta_vector[nulls_idx]
rotated_slip_ga_th = slip_ga_th_beta_vector[nulls_idx]

# %% rotated_slip_swyc_vol = slip_swyc_vol_beta_vector[nulls_idx]
rotated_slip_swyc_sa = slip_swyc_sa_beta_vector[nulls_idx]
#rotated_slip_swyc_th = slip_swyc_th_beta_vector[nulls_idx]
rotated_sa = sa_axis_rank_vector[nulls_idx]
rotated_allo = allo_rank_vector[nulls_idx]
rotated_evo = evo_rank_vector[nulls_idx]

# %% Test between BW and GA
slip_bwga_vol_corr, slip_bwga_vol_pval = stats.compare_images(slip_bw_vol_beta_vector, slip_ga_vol_beta_vector, nulls=rotated_slip_bw_vol, metric='spearmanr')
print(f"SLIP BW-GA volume corr: {slip_bwga_vol_corr}, p-val: {slip_bwga_vol_pval}")
slip_bwga_sa_corr, slip_bwga_sa_pval = stats.compare_images(slip_bw_sa_beta_vector, slip_ga_sa_beta_vector, nulls=rotated_slip_bw_sa, metric='spearmanr')
print(f"SLIP BW-GA surfArea corr: {slip_bwga_sa_corr}, p-val: {slip_bwga_sa_pval}")
slip_bwga_th_corr, slip_bwga_th_pval = stats.compare_images(slip_bw_th_beta_vector, slip_ga_th_beta_vector, nulls=rotated_slip_bw_th, metric='spearmanr')
print(f"SLIP BW-GA thickness corr: {slip_bwga_th_corr}, p-val: {slip_bwga_th_pval}")

# %% Test between BWp and GA
slip_bwpga_vol_corr, slip_bwpga_vol_pval = stats.compare_images(slip_bwp_vol_beta_vector, slip_ga_vol_beta_vector, nulls=rotated_slip_bwp_vol, metric='spearmanr')
print(f"SLIP BWP-GA volume corr: {slip_bwpga_vol_corr}, p-val: {slip_bwpga_vol_pval}")
slip_bwpga_sa_corr, slip_bwpga_sa_pval = stats.compare_images(slip_bwp_sa_beta_vector, slip_ga_sa_beta_vector, nulls=rotated_slip_bwp_sa, metric='spearmanr')
print(f"SLIP BWP-GA surfArea corr: {slip_bwpga_sa_corr}, p-val: {slip_bwpga_sa_pval}")
slip_bwpga_th_corr, slip_bwpga_th_pval = stats.compare_images(slip_bwp_th_beta_vector, slip_ga_th_beta_vector, nulls=rotated_slip_bwp_th, metric='spearmanr')
print(f"SLIP BWP-GA thickness corr: {slip_bwpga_th_corr}, p-val: {slip_bwpga_th_pval}")

# %% Test between Vol, SA, and Th
slip_ga_volsa_corr, slip_ga_volsa_pval = stats.compare_images(slip_ga_vol_beta_vector, slip_ga_sa_beta_vector, nulls=rotated_slip_ga_vol, metric='spearmanr')
print(f"SLIP GA Vol-SA corr: {slip_ga_volsa_corr}, p-val: {slip_ga_volsa_pval}")
slip_ga_volth_corr, slip_ga_volth_pval = stats.compare_images(slip_ga_vol_beta_vector, slip_ga_th_beta_vector, nulls=rotated_slip_ga_vol, metric='spearmanr')
print(f"SLIP GA Vol-Th corr: {slip_ga_volth_corr}, p-val: {slip_ga_volth_pval}")
slip_bwp_volth_corr, slip_bwp_volth_pval = stats.compare_images(slip_bwp_th_beta_vector, slip_bwp_th_beta_vector, nulls=rotated_slip_bwp_th, metric='spearmanr')
print(f"SLIP BWP Vol-Th corr: {slip_bwp_volth_corr}, p-val: {slip_bwp_volth_pval}")

# %% Test between birth factors and SWYC
#slip_bwswyc_vol_corr, slip_bwswyc_vol_pval = stats.compare_images(slip_bw_vol_beta_vector, slip_swyc_vol_beta_vector, nulls=rotated_slip_bw_vol, metric='spearmanr')
#print(f"SLIP BW-SWYC volume corr: {slip_bwswyc_vol_corr}, p-val: {slip_bwswyc_vol_pval}")
slip_bwswyc_sa_corr, slip_bwswyc_sa_pval = stats.compare_images(slip_bw_sa_beta_vector, slip_swyc_sa_beta_vector, nulls=rotated_slip_bw_sa, metric='spearmanr')
print(f"SLIP BW-SWYC surfArea corr: {slip_bwswyc_sa_corr}, p-val: {slip_bwswyc_sa_pval}")
#slip_bwswyc_th_corr, slip_bwswyc_th_pval = stats.compare_images(slip_bw_th_beta_vector, slip_swyc_th_beta_vector, nulls=rotated_slip_bw_th, metric='spearmanr')
#print(f"SLIP BW-SWYC thickness corr: {slip_bwswyc_th_corr}, p-val: {slip_bwswyc_th_pval}")

#slip_swycga_vol_corr, slip_swycga_vol_pval = stats.compare_images(slip_swyc_vol_beta_vector, slip_ga_vol_beta_vector, nulls=rotated_slip_swyc_vol, metric='spearmanr')
#print(f"SLIP SWYC-GA volume corr: {slip_swycga_vol_corr}, p-val: {slip_swycga_vol_pval}")
slip_swycga_sa_corr, slip_swycga_sa_pval = stats.compare_images(slip_swyc_sa_beta_vector, slip_ga_sa_beta_vector, nulls=rotated_slip_ga_sa, metric='spearmanr')
print(f"SLIP SWYC-GA surfArea corr: {slip_swycga_sa_corr}, p-val: {slip_swycga_sa_pval}")
#slip_swycga_th_corr, slip_swycga_th_pval = stats.compare_images(slip_swyc_th_beta_vector, slip_ga_th_beta_vector, nulls=rotated_slip_swyc_th, metric='spearmanr')
#print(f"SLIP SWYC-GA thickness corr: {slip_swycga_th_corr}, p-val: {slip_swycga_th_pval}")


# %% Test SA Axis Rank
slip_bwsaaxis_vol_corr, slip_bwsaaxis_vol_pval = stats.compare_images(slip_bw_vol_beta_vector, sa_axis_rank_vector, nulls=rotated_slip_bw_vol, metric='spearmanr')
print(f"SLIP BW-SA-AXIS volume corr: {slip_bwsaaxis_vol_corr}, p-val: {slip_bwsaaxis_vol_pval}")
slip_bwsaaxis_sa_corr, slip_bwsaaxis_sa_pval = stats.compare_images(slip_bw_sa_beta_vector, sa_axis_rank_vector, nulls=rotated_slip_bw_sa, metric='spearmanr')
print(f"SLIP BW-SA-AXIS surfArea corr: {slip_bwsaaxis_sa_corr}, p-val: {slip_bwsaaxis_sa_pval}")
slip_bwsaaxis_th_corr, slip_bwsaaxis_th_pval = stats.compare_images(slip_bw_th_beta_vector, sa_axis_rank_vector, nulls=rotated_slip_bw_th, metric='spearmanr')
print(f"SLIP BW-SA-AXIS thickness corr: {slip_bwsaaxis_th_corr}, p-val: {slip_bwsaaxis_th_pval}")

slip_gasaaxis_vol_corr, slip_gasaaxis_vol_pval = stats.compare_images(slip_ga_vol_beta_vector, sa_axis_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP GA-SA-AXIS volume corr: {slip_gasaaxis_vol_corr}, p-val: {slip_gasaaxis_vol_pval}")
slip_gasaaxis_sa_corr, slip_gasaaxis_sa_pval = stats.compare_images(slip_ga_sa_beta_vector, sa_axis_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP GA-SA-AXIS surfArea corr: {slip_gasaaxis_sa_corr}, p-val: {slip_gasaaxis_sa_pval}")
slip_gasaaxis_th_corr, slip_gasaaxis_th_pval = stats.compare_images(slip_ga_th_beta_vector, sa_axis_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP GA-SA-AXIS thickness corr: {slip_gasaaxis_th_corr}, p-val: {slip_gasaaxis_th_pval}")

#slip_swycsaaxis_vol_corr, slip_swycsaaxis_vol_pval = stats.compare_images(slip_swyc_vol_beta_vector, sa_axis_rank_vector, nulls=rotated_sa, metric='spearmanr')
#print(f"SLIP SWYC-SA-AXIS volume corr: {slip_swycsaaxis_vol_corr}, p-val: {slip_swycsaaxis_vol_pval}")
slip_swycsaaxis_sa_corr, slip_swycsaaxis_sa_pval = stats.compare_images(slip_swyc_sa_beta_vector, sa_axis_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP SWYC-SA-AXIS surfArea corr: {slip_swycsaaxis_sa_corr}, p-val: {slip_swycsaaxis_sa_pval}")
#slip_swycsaaxis_th_corr, slip_swycsaaxis_th_pval = stats.compare_images(slip_swyc_th_beta_vector, sa_axis_rank_vector, nulls=rotated_sa, metric='spearmanr')
#print(f"SLIP SWYC-SA-AXIS thickness corr: {slip_swycsaaxis_th_corr}, p-val: {slip_swycsaaxis_th_pval}")

# %% SA axis for BW and GA when controlling for global centiles
slip_bwsaaxis_global_vol_corr, slip_bwsaaxis_global_vol_pval = stats.compare_images(slip_bw_global_vol_beta_vector, sa_axis_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP BW-SA-AXIS GLOBAL COV volume corr: {slip_bwsaaxis_global_vol_corr}, p-val: {slip_bwsaaxis_global_vol_pval}")
slip_bwsaaxis_global_sa_corr, slip_bwsaaxis_global_sa_pval = stats.compare_images(slip_bw_global_sa_beta_vector, sa_axis_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP BW-SA-AXIS GLOBAL COV surfArea corr: {slip_bwsaaxis_global_sa_corr}, p-val: {slip_bwsaaxis_global_sa_pval}")
slip_bwsaaxis_global_th_corr, slip_bwsaaxis_global_th_pval = stats.compare_images(slip_bw_global_th_beta_vector, sa_axis_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP BW-SA-AXIS GLOBAL COV thickness corr (inc under one): {slip_bwsaaxis_global_th_corr}, p-val: {slip_bwsaaxis_global_th_pval}")

slip_gasaaxis_global_vol_corr, slip_gasaaxis_global_vol_pval = stats.compare_images(slip_ga_global_vol_beta_vector, sa_axis_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP GA-SA-AXIS GLOBAL COV volume corr: {slip_gasaaxis_global_vol_corr}, p-val: {slip_gasaaxis_global_vol_pval}")
slip_gasaaxis_global_sa_corr, slip_gasaaxis_global_sa_pval = stats.compare_images(slip_ga_global_sa_beta_vector, sa_axis_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP GA-SA-AXIS GLOBAL COV surfArea corr: {slip_gasaaxis_global_sa_corr}, p-val: {slip_gasaaxis_global_sa_pval}")
slip_gasaaxis_global_th_corr, slip_gasaaxis_global_th_pval = stats.compare_images(slip_ga_global_th_beta_vector, sa_axis_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP GA-SA-AXIS GLOBAL COV thickness corr: {slip_gasaaxis_global_th_corr}, p-val: {slip_gasaaxis_global_th_pval}")


# %% Test Allometry Rank
slip_bwallo_vol_corr, slip_bwallo_vol_pval = stats.compare_images(slip_bw_vol_beta_vector, allo_rank_vector, nulls=rotated_slip_bw_vol, metric='spearmanr')
print(f"SLIP BW-ALLO volume corr: {slip_bwallo_vol_corr}, p-val: {slip_bwallo_vol_pval}")
slip_bwallo_sa_corr, slip_bwallo_sa_pval = stats.compare_images(slip_bw_sa_beta_vector, allo_rank_vector, nulls=rotated_slip_bw_sa, metric='spearmanr')
print(f"SLIP BW-ALLO surfArea corr: {slip_bwallo_sa_corr}, p-val: {slip_bwallo_sa_pval}")
slip_bwallo_th_corr, slip_bwallo_th_pval = stats.compare_images(slip_bw_th_beta_vector, allo_rank_vector, nulls=rotated_slip_bw_th, metric='spearmanr')
print(f"SLIP BW-ALLO thickness corr: {slip_bwallo_th_corr}, p-val: {slip_bwallo_th_pval}")

slip_gaallo_vol_corr, slip_gaallo_vol_pval = stats.compare_images(slip_ga_vol_beta_vector, allo_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP GA-ALLO volume corr: {slip_gaallo_vol_corr}, p-val: {slip_gaallo_vol_pval}")
slip_gaallo_sa_corr, slip_gaallo_sa_pval = stats.compare_images(slip_ga_sa_beta_vector, allo_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP GA-ALLO surfArea corr: {slip_gaallo_sa_corr}, p-val: {slip_gaallo_sa_pval}")
slip_gaallo_th_corr, slip_gaallo_th_pval = stats.compare_images(slip_ga_th_beta_vector, allo_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP GA-ALLO thickness corr: {slip_gaallo_th_corr}, p-val: {slip_gaallo_th_pval}")

#slip_swycallo_vol_corr, slip_swycallo_vol_pval = stats.compare_images(slip_swyc_vol_beta_vector, allo_rank_vector, nulls=rotated_sa, metric='spearmanr')
#print(f"SLIP SWYC-ALLO volume corr: {slip_swycallo_vol_corr}, p-val: {slip_swycallo_vol_pval}")
slip_swycallo_sa_corr, slip_swycallo_sa_pval = stats.compare_images(slip_swyc_sa_beta_vector, allo_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP SWYC-ALLO surfArea corr: {slip_swycallo_sa_corr}, p-val: {slip_swycallo_sa_pval}")
#slip_swycallo_th_corr, slip_swycallo_th_pval = stats.compare_images(slip_swyc_th_beta_vector, allo_rank_vector, nulls=rotated_sa, metric='spearmanr')
#print(f"SLIP SWYC-ALLO thickness corr: {slip_swycallo_th_corr}, p-val: {slip_swycallo_th_pval}")


# %% Allometry Rank for BW and GA when controlling for global centiles
slip_bwallo_global_vol_corr, slip_bwallo_global_vol_pval = stats.compare_images(slip_bw_global_vol_beta_vector, allo_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP BW-ALLO GLOBAL COV volume corr: {slip_bwallo_global_vol_corr}, p-val: {slip_bwallo_global_vol_pval}")
slip_bwallo_global_sa_corr, slip_bwallo_global_sa_pval = stats.compare_images(slip_bw_global_sa_beta_vector, allo_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP BW-ALLO GLOBAL COV surfArea corr: {slip_bwallo_global_sa_corr}, p-val: {slip_bwallo_global_sa_pval}")
slip_bwallo_global_th_corr, slip_bwallo_global_th_pval = stats.compare_images(slip_bw_global_th_beta_vector, allo_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP BW-ALLO GLOBAL COV thickness corr (inc under one): {slip_bwallo_global_th_corr}, p-val: {slip_bwallo_global_th_pval}")

slip_gaallo_global_vol_corr, slip_gaallo_global_vol_pval = stats.compare_images(slip_ga_global_vol_beta_vector, allo_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP GA-ALLO GLOBAL COV volume corr: {slip_gaallo_global_vol_corr}, p-val: {slip_gaallo_global_vol_pval}")
slip_gaallo_global_sa_corr, slip_gaallo_global_sa_pval = stats.compare_images(slip_ga_global_sa_beta_vector, allo_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP GA-ALLO GLOBAL COV surfArea corr: {slip_gaallo_global_sa_corr}, p-val: {slip_gaallo_global_sa_pval}")
slip_gaallo_global_th_corr, slip_gaallo_global_th_pval = stats.compare_images(slip_ga_global_th_beta_vector, allo_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP GA-ALLO GLOBAL COV thickness corr: {slip_gaallo_global_th_corr}, p-val: {slip_gaallo_global_th_pval}")

# %% Evo rank testing (SWYC only)
slip_bwevo_vol_corr, slip_bwevo_vol_pval = stats.compare_images(slip_bw_vol_beta_vector, evo_rank_vector, nulls=rotated_slip_bw_vol, metric='spearmanr')
print(f"SLIP BW-EVO volume corr: {slip_bwevo_vol_corr}, p-val: {slip_bwevo_vol_pval}")
slip_bwevo_sa_corr, slip_bwevo_sa_pval = stats.compare_images(slip_bw_sa_beta_vector, evo_rank_vector, nulls=rotated_slip_bw_sa, metric='spearmanr')
print(f"SLIP BW-EVO surfArea corr: {slip_bwevo_sa_corr}, p-val: {slip_bwevo_sa_pval}")
slip_bwevo_th_corr, slip_bwevo_th_pval = stats.compare_images(slip_bw_th_beta_vector, evo_rank_vector, nulls=rotated_slip_bw_th, metric='spearmanr')
print(f"SLIP BW-EVO thickness corr: {slip_bwevo_th_corr}, p-val: {slip_bwevo_th_pval}")

slip_gaevo_vol_corr, slip_gaevo_vol_pval = stats.compare_images(slip_ga_vol_beta_vector, evo_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP GA-EVO volume corr: {slip_gaevo_vol_corr}, p-val: {slip_gaevo_vol_pval}")
slip_gaevo_sa_corr, slip_gaevo_sa_pval = stats.compare_images(slip_ga_sa_beta_vector, evo_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP GA-EVO surfArea corr: {slip_gaevo_sa_corr}, p-val: {slip_gaevo_sa_pval}")
slip_gaevo_th_corr, slip_gaevo_th_pval = stats.compare_images(slip_ga_th_beta_vector, evo_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP GA-EVO thickness corr: {slip_gaevo_th_corr}, p-val: {slip_gaevo_th_pval}")

#slip_swycevo_vol_corr, slip_swycevo_vol_pval = stats.compare_images(slip_swyc_vol_beta_vector, evo_rank_vector, nulls=rotated_sa, metric='spearmanr')
#print(f"SLIP SWYC-EVO volume corr: {slip_swycevo_vol_corr}, p-val: {slip_swycevo_vol_pval}")
slip_swycevo_sa_corr, slip_swycevo_sa_pval = stats.compare_images(slip_swyc_sa_beta_vector, evo_rank_vector, nulls=rotated_sa, metric='spearmanr')
print(f"SLIP SWYC-EVO surfArea corr: {slip_swycevo_sa_corr}, p-val: {slip_swycevo_sa_pval}")
#slip_swycevo_th_corr, slip_swycevo_th_pval = stats.compare_images(slip_swyc_th_beta_vector, evo_rank_vector, nulls=rotated_sa, metric='spearmanr')
#print(f"SLIP SWYC-EVO thickness corr: {slip_swycevo_th_corr}, p-val: {slip_swycevo_th_pval}")
# %%
