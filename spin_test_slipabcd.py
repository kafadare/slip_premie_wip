from neuromaps import  nulls, stats
import neuromaps as neuromaps
from nilearn.datasets import fetch_atlas_surf_destrieux
import abagen
import pandas as pd
import nibabel as nib

atlas_files = abagen.fetch_desikan_killiany(surface = True)
atlas_lh = nib.load(atlas_files['image'][0])
atlas_rh = nib.load(atlas_files['image'][1])
labels_lh = atlas_lh.labeltable.get_labels_as_dict()
labels_rh = atlas_rh.labeltable.get_labels_as_dict()

#timepoint 1
#read in the betas
abcd_t1_bw_betas = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/ABCD-braincharts/2.0_results/dkt_betas/dsk_birthweight_all_t1.csv")
abcd_t1_ga_betas = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/ABCD-braincharts/2.0_results/dkt_betas/dsk_gestAge_all_t1.csv")
abcd_t2_bw_betas = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/ABCD-braincharts/2.0_results/dkt_betas/dsk_birthweight_all_t2.csv")
abcd_t2_ga_betas = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/ABCD-braincharts/2.0_results/dkt_betas/dsk_gestAge_all_t2.csv")
slip_bw_betas = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/results_csv/dkt_betas/bw_volume.csv")
slip_ga_betas = pd.read_csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/results_csv/dkt_betas/ga_volume.csv")

#Split by hemisphere, which also removes non-cortical. Should be 34 for each hemisphere
##abcd bw
abcd_t1_bw_left = abcd_t1_bw_betas[abcd_t1_bw_betas['hemisphere'] == 'L'].copy()
abcd_t1_bw_right = abcd_t1_bw_betas[abcd_t1_bw_betas['hemisphere'] == 'R'].copy()
abcd_t2_bw_left = abcd_t2_bw_betas[abcd_t2_bw_betas['hemisphere'] == 'L'].copy()
abcd_t2_bw_right = abcd_t2_bw_betas[abcd_t2_bw_betas['hemisphere'] == 'R'].copy()
##abcd ga
abcd_t1_ga_left = abcd_t1_ga_betas[abcd_t1_ga_betas['hemisphere'] == 'L'].copy()
abcd_t1_ga_right = abcd_t1_ga_betas[abcd_t1_ga_betas['hemisphere'] == 'R'].copy()
abcd_t2_ga_left = abcd_t2_ga_betas[abcd_t2_ga_betas['hemisphere'] == 'L'].copy()
abcd_t2_ga_right = abcd_t2_ga_betas[abcd_t2_ga_betas['hemisphere'] == 'R'].copy()
##slip bw
slip_bw_left = slip_bw_betas[slip_bw_betas['hemisphere'] == 'L'].copy()
slip_bw_right = slip_bw_betas[slip_bw_betas['hemisphere'] == 'R'].copy()
##slip ga
slip_ga_left = slip_ga_betas[slip_ga_betas['hemisphere'] == 'L'].copy()
slip_ga_right = slip_ga_betas[slip_ga_betas['hemisphere'] == 'R'].copy()

# Convert label dicts into DataFrames
df_labels_lh = pd.DataFrame(list(labels_lh.items()), columns=['Key', 'label'])
df_labels_rh = pd.DataFrame(list(labels_rh.items()), columns=['Key', 'label'])

# Merge label-Key mapping with betas
##abcd bw
abcd_t1_bw_left = pd.merge(abcd_t1_bw_left, df_labels_lh, on='label')
abcd_t1_bw_right = pd.merge(abcd_t1_bw_right, df_labels_rh, on='label')
abcd_t2_bw_left = pd.merge(abcd_t2_bw_left, df_labels_lh, on='label')
abcd_t2_bw_right = pd.merge(abcd_t2_bw_right, df_labels_rh, on='label')
##abcd ga
abcd_t1_ga_left = pd.merge(abcd_t1_ga_left, df_labels_lh, on='label')
abcd_t1_ga_right = pd.merge(abcd_t1_ga_right, df_labels_rh, on='label')
abcd_t2_ga_left = pd.merge(abcd_t2_ga_left, df_labels_lh, on='label')
abcd_t2_ga_right = pd.merge(abcd_t2_ga_right, df_labels_rh, on='label')
##slip bw
slip_bw_left = pd.merge(slip_bw_left, df_labels_lh, on='label')
slip_bw_right = pd.merge(slip_bw_right, df_labels_rh, on='label')
##slip ga
slip_ga_left = pd.merge(slip_ga_left, df_labels_lh, on='label')
slip_ga_right = pd.merge(slip_ga_right, df_labels_rh, on='label')

# Sort by Key
##abcd bw
abcd_t1_bw_left = abcd_t1_bw_left.sort_values(by='Key')
abcd_t1_bw_right = abcd_t1_bw_right.sort_values(by='Key')
abcd_t2_bw_left = abcd_t2_bw_left.sort_values(by='Key')
abcd_t2_bw_right = abcd_t2_bw_right.sort_values(by='Key')
##abcd ga
abcd_t1_ga_left = abcd_t1_ga_left.sort_values(by='Key')
abcd_t1_ga_right = abcd_t1_ga_right.sort_values(by='Key')
abcd_t2_ga_left = abcd_t2_ga_left.sort_values(by='Key')
abcd_t2_ga_right = abcd_t2_ga_right.sort_values(by='Key')
##slip bw
slip_bw_left = slip_bw_left.sort_values(by='Key')
slip_bw_right = slip_bw_right.sort_values(by='Key')
##slip ga
slip_ga_left = slip_ga_left.sort_values(by='Key')
slip_ga_right = slip_ga_right.sort_values(by='Key')

# Combine hemispheres
##abcd bw
abcd_t1_bw_betas_ordered = pd.concat([abcd_t1_bw_left, abcd_t1_bw_right])
abcd_t2_bw_betas_ordered = pd.concat([abcd_t2_bw_left, abcd_t2_bw_right])
##abcd ga
abcd_t1_ga_betas_ordered = pd.concat([abcd_t1_ga_left, abcd_t1_ga_right])
abcd_t2_ga_betas_ordered = pd.concat([abcd_t2_ga_left, abcd_t2_ga_right])
##slip bw
slip_bw_betas_ordered = pd.concat([slip_bw_left, slip_bw_right])
##slip ga
slip_ga_betas_ordered = pd.concat([slip_ga_left, slip_ga_right])

# Get 1D vectors of beta values
##abcd bw
abcd_t1_bw_beta_vector = abcd_t1_bw_betas_ordered['bw_beta'].values
abcd_t2_bw_beta_vector = abcd_t2_bw_betas_ordered['bw_beta'].values
##abcd bw
abcd_t1_ga_beta_vector = abcd_t1_ga_betas_ordered['ga_beta'].values
abcd_t2_ga_beta_vector = abcd_t2_ga_betas_ordered['ga_beta'].values
##slip bw
slip_bw_beta_vector = slip_bw_betas_ordered['bw_beta'].values
##slip ga
slip_ga_beta_vector = slip_ga_betas_ordered['ga_beta'].values

#Generating indices of null spins instead of rotating data, so can apply indices to the map being rotated
nulls_idx = nulls.alexander_bloch(data=None, atlas='fsaverage', density='10k', parcellation=[atlas_lh, atlas_rh], n_perm=10000, seed=1234)

#Using the indices to rotate the abcd maps (will compare to slip)
rotated_abcd_t1_bw = abcd_t1_bw_beta_vector[nulls_idx]
rotated_abcd_t2_bw = abcd_t2_bw_beta_vector[nulls_idx]
rotated_abcd_t1_ga = abcd_t1_ga_beta_vector[nulls_idx]
rotated_abcd_t2_ga = abcd_t2_ga_beta_vector[nulls_idx]
rotated_slip_bw = slip_bw_beta_vector[nulls_idx]

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
