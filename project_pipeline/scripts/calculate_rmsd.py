from pymol import cmd

# Load and select native and pred pdbs.
def load_and_select(gt_fn, pred_fn, region1, region2):
    cmd.delete('all')
    cmd.load(gt_fn, 'native')
    cmd.load(pred_fn, 'pred')

    for obj in ['native','pred']:
        # select region1 and region2
        cmd.select(f'{obj}_1', f'{obj} and {region1}')
        cmd.select(f'{obj}_2', f'{obj} and {region2}')
        cmd.color('red', f'{obj}_1')
        cmd.color('green',f'{obj}_2')

def superimpose_region(region_num):
    # superimpose given region and calculate rmsd
    super = cmd.super(f'native_{region_num}',f'pred_{region_num}')
    cmd.color('purple','native_2')
    cmd.color('yellow','native_1')
    cmd.color('blue','pred_2')
    cmd.color('orange','pred_1')
    

# First define what those regions are.
# superimpose receptor chains and calculate rmsd for ligand
def calculate_rmsd(gt_pdb_fn, pred_pdb_fn, complex_fn, region1, region2, verbose=False):
        '''Calculate rmsd between gt and pred regions and whole proteins
            Region1 is autoinhibitory region, region2 is domain
        '''
        
        # Rmsds are in order of whole protein, region2, region1
        rmsds = []
        load_and_select \
            (gt_pdb_fn, pred_pdb_fn,
            region1, region2)

        # Superimpose region2 (domains) and calculate rmsd for whole protein and only region2. Save complex based on region2 alignment.
        superimpose_region(2)
        cmd.multisave(complex_fn, 'all', format='pdb')
        # Calculate rmsd of whole protein
        rmsd = cmd.rms_cur('native','pred')
        rmsds.append(rmsd)
        # Calculate rmsd of region2
        rmsd = cmd.rmsd_cur('native_2, pred_2')
        rmsds.append(rmsd)

        # Superimpose region1 (autoinhibitory regions) and calculate rmsd
        superimpose_region(1)
        rmsd = cmd.rms_cur('native_1', 'pred_1')
        rmsds.append(rmsd)

        # save two objects after superimposing receptor chain
        cmd.color('gray','native')
        cmd.color('red','pred')
        for obj in ['native','pred']:
            cmd.color('yellow', f'{obj}_interface_R')
            cmd.color('blue',f'{obj}_interface_L')

        rmsds = [round(rmsd,3) for rmsd in rmsds]
        if verbose: print(rmsds)
        return rmsds

# Superimpose receptors

# Calculate rmsd