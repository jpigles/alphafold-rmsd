from Bio.PDB.PDBList import PDBList
import pandas as pd
import numpy as np
from os.path import join
from pymol import cmd
import utils
import shutil
import os


###########
# Rule pdb_ids functions

def get_pdb_ids(df):
    '''Retrieves PDB IDs for each protein in the dataframe in the form of ID.chain (e.g. 1A2K.A)'''

    # Create a column to store the PDB IDs for each protein
    df['pdb'] = ''

    for i in range(len(df)):
        # Define UniProt ID and URL
        uniprot_id = df.loc[i, 'uniprot']
        url = 'https://search.rcsb.org/rcsbsearch/v2/query'
    
        pdb_ids = utils.query_rcsb(uniprot_id, url)

        # If received NaN from query, then drop row. Else, prune chains.
        if type(pdb_ids) == float:
            df = df.drop(index=[i])
        else:
            pdb_ids_pruned = utils.prune_extra_chains(pdb_ids)
            df.loc[i, 'pdb'] = pdb_ids_pruned.strip()

    df.reset_index(drop=True, inplace=True)
    return df

def download_pdb_files(df, path):
    '''Downloads the PDB files for each protein in the dataframe and saves them in a directory with the Uniprot ID.'''

    for i in range(len(df)):
        uniprot = df.loc[i, 'uniprot']
        uniprot_path = path + uniprot + '/'
        
        # Try to make a new directory with the gene name. If such a directory
        # already exists then continue
        try:
            os.mkdir(uniprot_path)
        except:
            print('Directory already exists.')

        pdb_ids_chains = df.loc[i, 'pdb']
        
        # Remove chains from the PDB IDs
        pdb_ids_no_chains = utils.remove_chains(pdb_ids_chains)

        # A PDB list object that allows to download PDB files
        pdbl = PDBList(verbose=False)

        print('Downloading structures for %s' % uniprot)

        # Retrieve the PDB file from the PDB and save to the directory with the gene name
        pdbl.download_pdb_files(pdb_ids_no_chains, pdir=uniprot_path, file_format='mmCif')

def correct_offset(df, path):

    offsets = []

    for i in range(len(df)):
        # Designate values for retrieval
        uniprot = df.loc[i, 'uniprot']
        pdb_id = df.loc[i, 'pdb']
        chain = df.loc[i, 'chain']
  
        # Designate file locations. Note that we will be overwriting the CIF files
        cif_path = path + uniprot + '/' + pdb_id + '.cif'

        # Get the offset
        offset = utils.get_offset(cif_path, pdb_id)
        offsets.append(offset)

        # Fix the offset
        fixed_pdb = utils.fix_offset(pdb_id, cif_path, chain, offset)
        print(fixed_pdb)

    # Add column with offset values
    df.insert(len(df.columns), 'label_offset', offsets)

    return df


###########
# Rule interface_analysis functions

def find_domain_completeness(df, path):

    # Create a new empty dataframe to store the results
    df_domain = pd.DataFrame(columns = ['gene_name', 'uniprot', 'protein_length', 'region_1', 'region_2', 'region_1_len', 
                                 'region_2_len', 'pdb', 'pdb_length', 'resolution',
                                 'model', 'chain', 'label_offset', 'pdb residues in region_1', 'pdb residues in region_2', 
                                 'percent_region_1', 'percent_region_2'])
    
    # Convert the domain region strings to ranges or lists of ranges
    df = utils.region_search_range(df)

    for i in range(len(df)):
    
        # Define values for retrieval
        region_1_res = df.loc[i, 'region_1 search']
        region_2_res = df.loc[i, 'region_2 search']
        pdb = df.loc[i, 'pdb']
        uniprot = df.loc[i, 'uniprot']
        path_uniprot = path + uniprot + '/'
        chain = df.loc[i, 'chain']

        # Get structure and dictionary objects
        structure, mmcif_dict = utils.get_structure_dict(pdb, path_uniprot)

        if mmcif_dict['_exptl.method'][0] == 'X-RAY DIFFRACTION':
            resolution = float(mmcif_dict["_refine.ls_d_res_high"][0])
        elif mmcif_dict['_exptl.method'][0] == 'SOLUTION NMR':
            resolution = np.nan

        # Count number of residues in each region
        count_res_reg_1, count_res_reg_2, count_res, model_id = utils.count_domain_residues(region_1_res, region_2_res, structure, chain)

        # Calculate the percentage of residues in each region
        percent_reg_1, percent_reg_2 = utils.calculate_domain_completeness(region_1_res, region_2_res, count_res_reg_1, count_res_reg_2)

        df_domain_part_1 = pd.DataFrame({'gene_name': df.loc[i, 'gene_name'],
                                'uniprot': df.loc[i, 'uniprot'],
                                'protein_length': df.loc[i, 'protein_length'],
                                'region_1': df.loc[i, 'region_1'],
                                'region_2': df.loc[i, 'region_2'],
                                'region_1_len': len(region_1_res),
                                'region_2_len': len(region_2_res),
                                'pdb': pdb, 
                                'pdb_length': count_res, 
                                'resolution': resolution, 
                                'model': model_id, 
                                'chain': chain,
                                'label_offset': df.loc[i, 'label_offset'],
                                'pdb residues in region_1': count_res_reg_1,
                                'pdb residues in region_2': count_res_reg_2,
                                'percent_region_1': percent_reg_1,
                                'percent_region_2': percent_reg_2}, index=[0])
        
        df_domain = df_domain.concat(df_domain_part_1, ignore_index=True)

    return df_domain

def save_domain_quality_files(df, path1, path2, path3, path4, path5):
    '''
    Save several copies of the dataframe with different filters based on the percentage of residues in the inhibitory and active domains.
    '''
    df.to_csv(path1, sep='\t', index=False)

    df_prot_both_80 = df.loc[(df['Percent residues in region_1'] > 80.0) & (df['Percent residues in region_2'] > 80.0)].reset_index(drop = True)
    df_prot_both_80.to_csv(path2, sep='\t', index=False)

    df_prot_1_80 = df.loc[(df['Percent residues in region_1'] > 80.0)].reset_index(drop = True)
    df_prot_1_80.to_csv(path3, sep='\t', index=False)

    df_prot_2_80 = df.loc[(df['Percent residues in region_2'] > 80.0)].reset_index(drop = True)
    df_prot_2_80.to_csv(path4, sep='\t', index=False)

    df_prot_both_60 = df.loc[(df['Percent residues in region_1'] > 60.0) & (df['Percent residues in region_2'] > 60.0)].reset_index(drop = True)
    df_prot_both_60.to_csv(path5, sep='\t', index=False)

    return list(df_prot_both_80, df_prot_1_80, df_prot_2_80, df_prot_both_60)

def copy_best_files(df, inpath, outpath):
    for i in range(len(df)):
        uniprot = df.loc[i, 'uniprot']
        pdb = df.loc[i, 'pdb']
        all_pdbs_path = join(inpath, uniprot, pdb + '.cif')
        best_pdbs_path = join(outpath, pdb + '.cif')
        shutil.copyfile(all_pdbs_path, best_pdbs_path)

    return f'Successfully copied files into {outpath}'

def get_interfaces(df, path):
    '''
    Finds the residues involved in the interface between the two domains.
    '''

    # Define columns for new stats
    df['PDB Mutations'] = ''
    df['Interacting residue pairs'] = ''
    df['Interface Residues'] = ''
    df['Number Interface Residues'] = ''

    df = utils.region_search_range(df)

    for i in range(len(df)):

        # Define values for retrieval 
        region_1_res = df.loc[i, 'region_1 search']
        region_2_res = df.loc[i, 'region_2 search']
        pdb = df.loc[i, 'pdb']
        chain = df.loc[i, 'chain']
        model = df.loc[i, 'model']

        # Get structure and dictionary objects
        structure, mmcif_dict = utils.get_structure_dict(pdb, path)

        # Get mutations
        df.loc[i, 'PDB Mutations'] = mmcif_dict['_entity.pdbx_mutation'][0]

        # Get residues in domains for Neighborsearch
        atoms_ns = utils.get_domain_residues(region_1_res, region_2_res, structure, model, chain)

        # Get interacting residues
        interacting_pairs, interface_res, len_interface_res = utils.domain_neighborsearch(df, region_1_res, region_2_res, atoms_ns)

        df.loc[i, 'Interacting residue pairs'] = interacting_pairs
        df.loc[i, 'Interface Residues'] = interface_res
        df.loc[i, 'Number Interface Residues'] = len_interface_res

    return df

def largest_interface(df):
    ''' 
    Go through all the proteins in df_prot and determine the pdb file with the greatest 
    number of interface residues between the region_1 and the region_2
    '''

    # Convert the int values to numpy.int64 (this is required for the .idxmax method to be applied)
    df.loc[:, 'Number Interface Residues'] = pd.to_numeric(df['Number Interface Residues'])

    # Make a new column to flag the rows to keep
    df['Keep'] = ''

    # Get all the Uniprot_IDs
    proteins = set(df['Uniprot_ID'])

    # Iterate through all the proteins
    for protein in proteins:
        
        print('Determining the interface residues for', protein)
        
        # Grab all the instances of each protein in df_prot
        df_temp = df.loc[df['Uniprot_ID'] == protein]
        
        # Grab the intances of each protein in df_prot with no mutations
        df_temp_no_mut = df_temp.loc[df_temp['PDB Mutations'] == '?']
        
        # Grab all the intances of each protein with mutations
        df_temp_mut = df_temp.loc[df_temp['PDB Mutations'] != '?']
        
        # From the PDB entries with no mutations, select the one with the largest 
        # number of residues at the interface
        if len(df_temp_no_mut) > 0 :
            # Get the index of the max value in the Number Interface Residues column
            max_index_no_mut = df_temp_no_mut['Number Interface Residues'].idxmax(skipna = True)
            
            # Set the row with the max index to True in the Keep column of df_prot
            df.loc[max_index_no_mut, 'Keep'] = True
        
        # If all the PDB entries have mutations, select the one with the largest
        # number of residues at the interface
        elif len(df_temp_mut) > 0:
            # Get the index of the max value in the Number Interface Residues column
            max_index_mut = df_temp_mut['Number Interface Residues'].idxmax(skipna = True)
            
            # Set the row with the max index to True in the Keep column of df_prot
            df.loc[max_index_mut, 'Keep'] = True

    # Make a new df with the proteins with the highest number of interface residues
    df_prot_keep = df.loc[df['Keep'] == True]

    df_prot_result = df.copy()
    df_prot_keep_result = df_prot_keep.copy()

    df_prot_result.loc[:,'Interface Residues'] = df_prot_result['Interface Residues'].apply(utils.to_string)    

    df_prot_keep_result.loc[:, 'Interface Residues'] = df_prot_keep_result['Interface Residues'].apply(utils.to_string)

    df_prot_keep_result.loc[:, 'Interacting residue pairs'] = df_prot_keep_result['Interacting residue pairs'].apply(utils.to_string)

    df_prot_keep_result = df_prot_keep_result.dropna(subset = ['Uniprot_ID']).reset_index(drop = True)
    df_prot_keep_result = df_prot_keep_result.drop(['region_1 search', 'region_2 search', 'Keep'], axis = 'columns')

    return df_prot_keep_result

def calculate_rmsd(gt_pdb_fn, pred_pdb_fn, complex_fn, region_1, region_2):
    '''Calculate rmsd between gt and pred regions and whole proteins
        Region1 is autoinhibitory region, region2 is domain
    '''
    
    # Rmsds is a dict of variable size with format {'complex_rmsd': ####, '1.0_aligned': ###, '1.0_comp': ###, ...}. Size based on number subregions.
    rmsds = {}
    utils.load_and_select \
        (gt_pdb_fn, pred_pdb_fn,
        region_1, region_2)

    # Superimpose entirety of proteins
    align = cmd.align('native', 'pred')
    cmd.multisave(complex_fn, 'all', format='pdb')
    # Calculate rmsd of whole protein
    rmsd = cmd.rms_cur('native','pred')
    rmsds['complex_rmsd'] = round(rmsd, 3)

    # Superimpose each region and calculate rmsds
    for key in region_1:
        if len(region_1) > 1 and '.0' in key:
            continue
        two_rmsds = utils.align_and_calculate(key, '2.0')
        rmsds[key + '_aligned'] = two_rmsds[0]
        rmsds[key + '_comp'] = two_rmsds[1]

    for key in region_2:
        if len(region_2) > 1 and '.0' in key:
            continue
        two_rmsds = utils.align_and_calculate(key, '1.0')
        rmsds[key + '_aligned'] = two_rmsds[0]
        rmsds[key + '_comp'] = two_rmsds[1]

    return rmsds

def get_rmsds(df):
    '''
    Calculate rmsds for each protein in df, aligning first on the autoinhibitory region (region 1) and then on the active region (region 2). Regions with
    multiple subregions are aligned and calculated separately, and then the average is taken.
    '''
    rmsd_info = []
    for i in range(len(df)):
        # Define pdb, filenames, region1, region2
        pdb = df.loc[i, 'PDB ID']
        uniprot = df.loc[i, 'Uniprot_ID']
        region_1_dict = utils.create_region_dict(df.loc[i, 'region_1'], 1)
        region_2_dict = utils.create_region_dict(df.loc[i, 'region_2'], 2)
        percent_reg1 = df.loc[i, 'Percent residues in region_1']
        percent_reg2 = df.loc[i, 'Percent residues in region_2']
        gt_fn = f'./data/input/RCSB/pdbs_trim/{pdb}.pdb'
        pred_fn = f'./data/output/RCSB_af_full/af_trim/{pdb}.fasta/ranked_0.pdb'
        complex_fn = f'./data/output/RCSB_af_full/complex/{pdb}.pdb'

        print(f'Trying {pdb}...')
        rmsds = calculate_rmsd(gt_fn, pred_fn, complex_fn, region_1_dict, region_2_dict)

        # Define default values for columns to retain number of columns per row
        rmsd_dic = {'UniProt': uniprot,
                    'PDB': pdb,
                    'complex_rmsd': 0,
                    '1.0_aligned': 0,
                    '1.0_comp': 0,
                    '1.1_aligned': 0,
                    '1.1_comp': 0,
                    '1.2_aligned': 0,
                    '1.2_comp': 0,
                    '2.0_aligned': 0,
                    '2.0_comp': 0,
                    '2.1_aligned': 0,
                    '2.1_comp': 0,
                    '2.2_aligned': 0,
                    '2.2_comp': 0,
                    '2.3_aligned': 0,
                    '2.3_comp': 0,
                    'Percent residues in region_1': percent_reg1,
                    'Percent residues in region_2': percent_reg2}

        for key in rmsds:
            if key in rmsd_dic:
                rmsd_dic[key] = rmsds[key]

        print('Success! Writing rmsds')
        rmsd_info.append(rmsd_dic)

    final_rmsds = utils.get_region_averages(rmsd_info)
    return final_rmsds