from pymol import cmd
import pandas as pd
import csv

pdb_df = pd.read_csv('./data/proteins_pdb_best.tsv', sep='\t').astype('object')

def create_region_dict(region, region_num):
    '''Create a dictionary containing an ID for every region in the domain
    For instance, if domain 1 has 123-222,333-444, then make dict {1.0: 123-222+333-444, 1.1: 123-222, 1.2: 333-444}.
    #.0 always contains the full number of regions.'''
    full_region = region.strip()
    # Each dict will have at least one entry.
    region_dict = {f'{region_num}.0': full_region}
    if ',' in region:
        full_region = region.replace(',', '+')
        # Substitute full_region in the dict to one with + in place of ,
        region_dict[f'{region_num}.0'] = full_region
        regions_list = region.split(',')
        for i in range(len(regions_list)):
            # Make subregions 1-indexed to preserve #.0 as full.
            subregion = i + 1
            region_dict[f'{region_num}.{subregion}'] = regions_list[i]

    return region_dict


def load_and_select(gt_fn, pred_fn, region_1, region_2):
    # Load and select native and pred pdbs
    cmd.delete('all')
    cmd.load(gt_fn, 'native')
    cmd.load(pred_fn, 'pred')

    for obj in ['native','pred']:
        # select region1 and region2
        for key in region_1.keys:
            # example: native_1.1, native and resi 111-222
            resi_range = region_1[key]
            cmd.select(f'{obj}_{key}', f'{obj} and resi {resi_range}')
        for key in region_2.keys:
            resi_range = region_2[key]
            cmd.select(f'{obj}_{key}', f'{obj} and resi {resi_range}')

def align_and_calculate(align_reg_key, comp_region_key):
    # superimpose aligned region and calculate two rmsds: One for aligned region and one for complementary region (for example, rmsd for "aligned" region 1.1 
    # and complementary region 2.0)
    rmsds = []
    try:
        align = cmd.align(f'native_{align_reg_key}', f'pred_{align_reg_key}')
        rmsd = cmd.cur_rms(f'native_{align_reg_key}', f'pred_{align_reg_key}')
        rmsds.append(rmsd)
        rmsd = cmd.cur_rms(f'native_{comp_region_key}', f'pred_{comp_region_key}')
        rmsds.append(rmsd)
        return rmsds

    except pymol.CmdException:
        print(f'Region {align_reg_key} missing')
        rmsds = [0, 0]
        return rmsds

def calculate_rmsd(gt_pdb_fn, pred_pdb_fn, complex_fn, region_1, region_2):
    '''Calculate rmsd between gt and pred regions and whole proteins
        Region1 is autoinhibitory region, region2 is domain
    '''
    
    # Rmsds is a dict of variable size with format {'complex_rmsd': ####, '1.0_aligned': ###, '1.0_comp': ###, ...}. Size based on number subregions.
    rmsds = {}
    load_and_select \
        (gt_pdb_fn, pred_pdb_fn,
        region_1, region_2)

    # Superimpose entirety of proteins
    align = cmd.align('native', 'pred')
    cmd.multisave(complex_fn, 'all', format='pdb')
    # Calculate rmsd of whole protein
    rmsd = cmd.rms_cur('native','pred')
    rmsds['complex_rmsd'] = rmsd

    # Superimpose each region and calculate rmsds
    for key in region_1.keys:
        if len(region_1.keys) > 1 and '.0' in key:
            continue
        two_rmsds = align_and_calculate(key, '2.0')
        rmsds[key + '_aligned'] = two_rmsds[0]
        rmsds[key + '_comp'] = two_rmsds[1]

    for key in region_2.keys:
        if len(region_2.keys) > 1 and '.0' in key:
            continue
        two_rmsds = align_and_calculate(key, '1.0')
        rmsds[key + '_aligned'] = two_rmsds[0]
        rmsds[key + '_comp'] = two_rmsds[1]

    rmsds = [round(rmsd,3) for rmsd in rmsds.values]
    return rmsds

rmsd_info = []
for i in range(len(pdb_df)):
    # Define pdb, filenames, region1, region2
    pdb = pdb_df.loc[i, 'PDB ID']
    uniprot = pdb_df.loc[i, 'Uniprot_ID']
    region_1_dict = create_region_dict(pdb_df.loc[i, 'region_1'], 1)
    region_2_dict = create_region_dict(pdb_df.loc[i, 'region_2'], 2)
    percent_reg1 = pdb_df.loc[i, 'Percent residues in region_1']
    percent_reg2 = pdb_df.loc[i, 'Percent residues in region_2']
    gt_fn = f'./data/input/RCSB/pdbs_trim/{pdb}.pdb'
    pred_fn = f'./data/output/RCSB_af_full/af_trim/{pdb}.fasta/ranked_0.pdb'
    complex_fn = f'./data/output/RCSB_af_full/complex/{pdb}.pdb'
    
    print(f'Trying {pdb}...')
    rmsds = calculate_rmsd(gt_fn, pred_fn, complex_fn, region_1_dict, region_2_dict)
    
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