from pymol import cmd
import pandas as pd
import csv

pdb_df = pd.read_csv('./data/proteins_pdb_best.tsv', sep='\t').astype('object')

def replace_commas(region):
    '''If a region has several sections, do two things
        First, create full_region which just replaces the comma with + (from ###-###,###-### to ###-###+###-###)
        Second, create anchor_region, which is just the first section (from 111-111,222-222 to 111-111) which will be used to superimpose
        If no commas, full_region and anchor_region are the same'''
    full_region = region.strip()
    anchor_region = region.strip()
    if ',' in region:
        full_region = region.replace(',', '+')
        anchor_region = region.split(',')[0]
        
    return full_region, anchor_region


def load_and_select(gt_fn, pred_fn, region_1, region_2, region1_anchor, region2_anchor):
    # Load and select native and pred pdbs
    cmd.delete('all')
    cmd.load(gt_fn, 'native')
    cmd.load(pred_fn, 'pred')

    for obj in ['native','pred']:
        # select region1 and region2
        cmd.select(f'{obj}_1', f'{obj} and resi {region_1}')
        cmd.select(f'{obj}_2', f'{obj} and resi {region_2}')
        cmd.select(f'{obj}1_anchor', f'{obj} and resi {region1_anchor}')
        cmd.select(f'{obj}2_anchor', f'{obj} and resi {region2_anchor}')
        cmd.color('red', f'{obj}_1')
        cmd.color('green',f'{obj}_2')

def superimpose_region(region_num):
    # superimpose given region and calculate rmsd
    try:
        super = cmd.super(f'native{region_num}_anchor',f'pred{region_num}_anchor')
    except pymol.CmdException:
        print(f'Region {region_num} missing')
        return False

    cmd.color('purple','native_2')
    cmd.color('yellow','native_1')
    cmd.color('blue','pred_2')
    cmd.color('orange','pred_1')
    

def calculate_rmsd(gt_pdb_fn, pred_pdb_fn, complex_fn, region_1, region_2, region1_anchor, region2_anchor, verbose=False):
    '''Calculate rmsd between gt and pred regions and whole proteins
        Region1 is autoinhibitory region, region2 is domain
    '''
    
    # Rmsds are in order of whole protein, region2, region1
    rmsds = []
    load_and_select \
        (gt_pdb_fn, pred_pdb_fn,
        region_1, region_2, 
        region1_anchor, region2_anchor)

    # Superimpose region2 (domains) and calculate rmsd for whole protein and only region2. Save complex based on region2 alignment.
    superimpose_region(2)
    cmd.multisave(complex_fn, 'all', format='pdb')
    # Calculate rmsd of whole protein
    rmsd = cmd.rms_cur('native','pred')
    rmsds.append(rmsd)
    # Calculate rmsd of region2
    rmsd = cmd.rms_cur('native_2', 'pred_2')
    rmsds.append(rmsd)

    # Superimpose region1 (autoinhibitory regions) and calculate rmsd
    region1_sup = superimpose_region(1)
    if region1_sup == False:
        rmsd = 0
    else:
        rmsd = cmd.rms_cur('native_1', 'pred_1')
    rmsds.append(rmsd)

    # save two objects after superimposing region1
    cmd.color('gray','native')
    cmd.color('red','pred')
    for obj in ['native','pred']:
        cmd.color('yellow', f'{obj}_1')
        cmd.color('blue',f'{obj}_1')

    rmsds = [round(rmsd,3) for rmsd in rmsds]
    if verbose: print(rmsds)
    return rmsds

rmsd_info = []
for i in range(len(pdb_df)):
    # Define pdb, filenames, region1, region2
    pdb = pdb_df.loc[i, 'PDB ID']
    uniprot = pdb_df.loc[i, 'Uniprot_ID']
    region_1, region1_anchor = replace_commas(pdb_df.loc[i, 'region_1'])
    region_2, region2_anchor = replace_commas(pdb_df.loc[i, 'region_2'])
    percent_reg1 = pdb_df.loc[i, 'Percent residues in region_1']
    percent_reg2 = pdb_df.loc[i, 'Percent residues in region_2']
    gt_fn = f'./data/input/RCSB/pdbs_trim/{pdb}.pdb'
    pred_fn = f'./data/output/RCSB_af_full/af_trim/{pdb}.fasta/ranked_0.pdb'
    complex_fn = f'./data/output/RCSB_af_full/complex/{pdb}.pdb'
    
    print(f'Trying {pdb}...')
    rmsds = calculate_rmsd(gt_fn, pred_fn, complex_fn, region_1, region_2, region1_anchor, region2_anchor)
    
    rmsd_dic = {'UniProt': uniprot,
                'PDB': pdb,
                'complex_rmsd': rmsds[0],
                'region1_rmsd': rmsds[2],
                'region2_rmsd': rmsds[1],
                'Percent residues in region_1': percent_reg1,
                'Percent residues in region_2': percent_reg2}
    print('Success! Writing rmsds')
    rmsd_info.append(rmsd_dic)

with open('./data/rmsds.tsv', 'w') as file:
    fields = ['UniProt', 'PDB', 'complex_rmsd', 'region1_rmsd', 'region2_rmsd', 'Percent residues in region_1', 'Percent residues in region_2']
    writer = csv.DictWriter(file, fieldnames=fields, delimiter='\t')
    
    writer.writeheader()
    for item in rmsd_info:
        writer.writerow(item)