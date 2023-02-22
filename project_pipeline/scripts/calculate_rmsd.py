

# Load and select native and pred pdbs.
cmd.delete('all')
cmd.load(gt_fn, 'native')
cmd.load(pred_fn, 'pred')

# select region1 (autoinhibitory domain) and region2 ("functional" domain)

# First define what those regions are.
region1_num = df.loc[pdb, 'region_1']
region2_num = df.loc[pdb, 'region_2']
if ',' in region1_num or region2_num:
    
# Superimpose receptors

# Calculate rmsd

