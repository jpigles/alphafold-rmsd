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

        # color 
        cmd.color('yellow', f'{obj}_interface_R')
        cmd.color('blue',f'{obj}_interface_L')ibitory domain) and region2 ("functional" domain)

def superimpose_receptors(complex_fn):
    # superimpose receptor chains and calculate rmsd for idr, assume existence of corresp sels
    super = cmd.super('native_R','pred_R')
    cmd.color('purple','native_R')
    cmd.color('yellow','native_L')
    cmd.color('gray','pred_R')
    cmd.color('orange','pred_L')
    cmd.multisave(complex_fn, 'all', format='pdb')

# First define what those regions are.
# superimpose receptor chains and calculate rmsd for ligand
def calculate_rmsd(self, pdb_id, gt_pdb_fn, pred_pdb_fn, complex_fn, verbose=False):
        ''' Calculate rmsd between gt and pred
              superimpose pred onto gt (with receptor only)
              assume chain_names[-1] is the target chain (e.g. peptide or idr)
        '''
        rmsds = []
        aa_utils.load_and_select \
            (self.interface_dist, gt_pdb_fn, pred_pdb_fn,
             self.chain_ids[pdb_id], backbone=self.backbone,
             remove_hydrogen=self.remove_hydrogen)
        aa_utils.superimpose_receptors(complex_fn)
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