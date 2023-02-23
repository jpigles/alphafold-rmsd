from pymol import cmd
import pandas as pd
import csv

pdb_df = pd.read_csv('./data/proteins_pdb_best.tsv', sep='\t').astype('object')

def replace_commas(region):
    cmd_region = region
    if ',' in cmd_region:
        cmd_region = cmd_region.replace(',', '+')
        
    return cmd_region


def load_and_select(gt_fn, pred_fn, region1, region2):
    # Load and select native and pred pdbs
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

    # save two objects after superimposing region1
    cmd.color('gray','native')
    cmd.color('red','pred')
    for obj in ['native','pred']:
        cmd.color('yellow', f'{obj}_interface_R')
        cmd.color('blue',f'{obj}_interface_L')

    rmsds = [round(rmsd,3) for rmsd in rmsds]
    if verbose: print(rmsds)
    return rmsds

rmsd_info = []
for i in range(len(pdb_df)):
    # Define pdb, filenames, region1, region2
    pdb = pdb_df.loc[i, 'PDB ID']
    uniprot = pdb_df.loc[i, 'Uniprot_ID']
    region_1 = replace_commas(pdb_df.loc[i, 'region_1'])
    region_2 = replace_commas(pdb_df.loc[i, 'region_2'])
    percent_reg1 = pdb_df.loc[i, 'Percent residues in region_1']
    percent_reg2 = pdb_df.loc[i, 'Percent residues in region_2']
    gt_fn = f'./data/input/RCSB/pdbs_trim/{pdb}.pdb'
    pred_fn = f'./data/output/RCSB_af_full/af_trim/{pdb}.fasta/ranked_0.pdb'
    complex_fn = f'./data/output/RCSB_af_full/complex/{pdb}.pdb'

    rmsds = calculate_rmsd(gt_fn, pred_fn, complex_fn, region_1, region_2)

    rmsd_dic = {'UniProt': uniprot,
                'PDB': pdb,
                'complex_rmsd': rmsds[0],
                'region1_rmsd': rmsds[2],
                'region2_rmsd': rmsds[1],
                'Percent residues in region_1': percent_reg1,
                'Percent residues in region_2': percent_reg2}
    rmsd_info.append(rmsd_dic)

with open('./data/rmsds.tsv', 'w') as file:
    fields = ['Uniprot', 'PDB', 'complex_rmsd', 'region1_rmsd', 'region2_rmsd', 'Percent residues in region_1', 'Percent residues in region_2']
    writer = csv.DictWriter(file, fieldnames=fields, delimiter='\t')
    
    writer.writeheader()
    for item in rmsd_info:
        writer.writerow(item)