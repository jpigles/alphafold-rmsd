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

def superimpose_region(region_num, complex_fn):
    # superimpose given region and calculate rmsd
    super = cmd.super('native_2','pred_2')
    cmd.color('purple','native_2')
    cmd.color('yellow','native_1')
    cmd.color('gray','pred_2')
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

        # Superimpose region2 (domains) and calculate rmsd for whole protein and only region2
        superimpose_region(2, complex_fn)`
        cmd.multisave(complex_fn, 'all', format='pdb')
        rmsd = cmd.rms_cur('native_L','pred_L')
        rmsds.append(rmsd)


        # save two objects after superimposing receptor chain
        cmd.color('gray','native')
        cmd.color('red','pred')
        for obj in ['native','pred']:
            cmd.color('yellow', f'{obj}_interface_R')
            cmd.color('blue',f'{obj}_interface_L')

        # calculate rmsd for interface idr only
        rmsd = cmd.rms_cur('native_interface_L','pred_interface_L')
        rmsds.append(rmsd)

        rmsds = [round(rmsd,3) for rmsd in rmsds]
        if verbose: print(rmsds)
        return rmsds

# Superimpose receptors

# Calculate rmsd