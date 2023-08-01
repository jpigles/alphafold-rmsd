from pdbecif.mmcif_io import CifFileReader, CifFileWriter
from Bio.PDB.PDBList import PDBList
import pandas as pd
import numpy as np
from os.path import join
from pymol import cmd
import pymol
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
        
        # Try to make a new directory with the gene name. If such a directory
        # already exists then continue
        utils.uniprot_dirs(path, uniprot)

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
        offset = utils.get_offset(cif_path, pdb_id, uniprot)
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
        fn = f'{pdb}.cif'

        print('Analyzing %s' % pdb)

        # Get structure and dictionary objects
        structure, mmcif_dict = utils.get_structure_dict(pdb, fn, path_uniprot)

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
        
        df_domain = pd.concat([df_domain, df_domain_part_1], axis=0, ignore_index=True)

    return df_domain

def save_domain_quality_files(df, path1, path2, path3, path4, path5):
    '''
    Save several copies of the dataframe with different filters based on the percentage of residues in the inhibitory and active domains.
    '''
    df_list = []
    df.to_csv(path1, sep='\t', index=False)

    df_prot_both_80 = df.loc[(df['percent_region_1'] > 80.0) & (df['percent_region_2'] > 80.0)].reset_index(drop = True)
    df_prot_both_80.to_csv(path2, sep='\t', index=False)
    df_list.append(df_prot_both_80)

    df_prot_1_80 = df.loc[(df['percent_region_1'] > 80.0)].reset_index(drop = True)
    df_prot_1_80.to_csv(path3, sep='\t', index=False)
    df_list.append(df_prot_1_80)

    df_prot_2_80 = df.loc[(df['percent_region_2'] > 80.0)].reset_index(drop = True)
    df_prot_2_80.to_csv(path4, sep='\t', index=False)
    df_list.append(df_prot_2_80)

    df_prot_both_60 = df.loc[(df['percent_region_1'] > 60.0) & (df['percent_region_2'] > 60.0)].reset_index(drop = True)
    df_prot_both_60.to_csv(path5, sep='\t', index=False)
    df_list.append(df_prot_both_60)

    return df_list

def copy_best_files(df, inpath, outpath):
    # Make sure the output directory exists
    utils.make_dirs([outpath])

    for i in range(len(df)):
        uniprot = df.loc[i, 'uniprot']
        pdb = df.loc[i, 'pdb']

        # Make sure the output directory exists
        utils.uniprot_dirs(outpath, uniprot)

        # Copy the files
        source_pdbs_path = join(inpath, uniprot, pdb + '.cif')
        best_pdbs_path = join(outpath, uniprot, pdb + '.cif')
        shutil.copyfile(source_pdbs_path, best_pdbs_path)

    return f'Successfully copied files into {outpath}'

def get_interfaces(df, path):
    '''
    Finds the residues involved in the interface between the two domains.
    '''

    # Define columns for new stats
    df['pdb_mutations'] = ''
    df['interacting_residue_pairs'] = ''
    df['interface_residues'] = ''
    df['number_interface_residues'] = ''

    df = utils.region_search_range(df)

    for i in range(len(df)):

        # Define values for retrieval 
        region_1_res = df.loc[i, 'region_1 search']
        region_2_res = df.loc[i, 'region_2 search']
        pdb = df.loc[i, 'pdb']
        uniprot = df.loc[i, 'uniprot']
        chain = df.loc[i, 'chain']
        model = df.loc[i, 'model']
        fn = f'{pdb}_{uniprot}.cif'

        # Get structure and dictionary objects
        structure, mmcif_dict = utils.get_structure_dict(pdb, fn, path)

        # Get mutations
        df.loc[i, 'pdb_mutations'] = mmcif_dict['_entity.pdbx_mutation'][0]

        # Get residues in domains for Neighborsearch
        atoms_ns = utils.get_domain_residues(region_1_res, region_2_res, structure, model, chain)

        # Get interacting residues
        interacting_pairs, interface_res, len_interface_res = utils.domain_neighborsearch(region_1_res, region_2_res, atoms_ns)

        df.loc[i, 'interacting_residue_pairs'] = interacting_pairs
        df.loc[i, 'interface_residues'] = interface_res
        df.loc[i, 'number_interface_residues'] = len_interface_res

    return df

def get_af_interfaces(df, path):
    '''
    Finds the residues involved in the interface between the two domains, 
    specifically for the AlphaFold files.
    '''

    # Define columns for new stats
    df['interacting_residue_pairs'] = ''
    df['interface_residues'] = ''
    df['number_interface_residues'] = ''

    df = utils.region_search_range(df)

    for i in range(len(df)):

        # Define values for retrieval 
        region_1_res = df.loc[i, 'region_1 search']
        region_2_res = df.loc[i, 'region_2 search']
        uniprot = df.loc[i, 'uniprot']
        chain = df.loc[i, 'chain']
        model = df.loc[i, 'model']
        fn = f'F-{uniprot}-F1-model_v3.cif'

        # Get structure and dictionary objects
        structure, mmcif_dict = utils.get_structure_dict(uniprot, fn, path)

        # Get residues in domains for Neighborsearch
        atoms_ns = utils.get_domain_residues(region_1_res, region_2_res, structure, model, chain)

        # Get interacting residues
        interacting_pairs, interface_res, len_interface_res = utils.domain_neighborsearch(region_1_res, region_2_res, atoms_ns)

        df.loc[i, 'interacting_residue_pairs'] = interacting_pairs
        df.loc[i, 'interface_residues'] = interface_res
        df.loc[i, 'number_interface_residues'] = len_interface_res

    return df

def largest_interface(df):
    ''' 
    Go through all the proteins in df_prot and determine the pdb file with the greatest 
    number of interface residues between the region_1 and the region_2
    '''

    # Convert the int values to numpy.int64 (this is required for the .idxmax method to be applied)
    df.loc[:, 'number_interface_residues'] = pd.to_numeric(df['number_interface_residues'])

    # Make a new column to flag the rows to keep
    df['Keep'] = ''

    # Get all the Uniprot_IDs
    proteins = set(df['uniprot'])

    # Iterate through all the proteins
    for protein in proteins:
        
        print('Determining the interface residues for', protein)
        
        # Grab all the instances of each protein in df_prot
        df_temp = df.loc[df['uniprot'] == protein]
        
        # Grab the intances of each protein in df_prot with no mutations
        df_temp_no_mut = df_temp.loc[df_temp['pdb_mutations'] == '?']
        
        # Grab all the intances of each protein with mutations
        df_temp_mut = df_temp.loc[df_temp['pdb_mutations'] != '?']
        
        # From the PDB entries with no mutations, select the one with the largest 
        # number of residues at the interface
        if len(df_temp_no_mut) > 0 :
            # Get the index of the max value in the Number Interface Residues column
            max_index_no_mut = df_temp_no_mut['number_interface_residues'].idxmax(skipna = True)
            
            # Set the row with the max index to True in the Keep column of df_prot
            df.loc[max_index_no_mut, 'Keep'] = True
        
        # If all the PDB entries have mutations, select the one with the largest
        # number of residues at the interface
        elif len(df_temp_mut) > 0:
            # Get the index of the max value in the Number Interface Residues column
            max_index_mut = df_temp_mut['number_interface_residues'].idxmax(skipna = True)
            
            # Set the row with the max index to True in the Keep column of df_prot
            df.loc[max_index_mut, 'Keep'] = True

    # Make a new df with the proteins with the highest number of interface residues
    df_prot_keep = df.loc[df['Keep'] == True]

    df_prot_result = df.copy()
    df_prot_keep_result = df_prot_keep.copy()

    df_prot_result.loc[:,'interface_residues'] = df_prot_result['interface_residues'].apply(utils.to_string)    

    df_prot_keep_result.loc[:, 'interface_residues'] = df_prot_keep_result['interface_residues'].apply(utils.to_string)

    df_prot_keep_result.loc[:, 'interacting_residue_pairs'] = df_prot_keep_result['interacting_residue_pairs'].apply(utils.to_string)

    df_prot_keep_result = df_prot_keep_result.dropna(subset = ['uniprot']).reset_index(drop = True)
    df_prot_keep_result = df_prot_keep_result.drop(['region_1 search', 'region_2 search', 'Keep'], axis = 'columns')

    return df_prot_keep_result

def trim_cifs(df, gt_path_in, gt_path_out, pred_path_in, pred_path_out):

    trim_values = []
    for i in range(len(df)):
        
        # Define parameters for selecting files
        uniprot = df.loc[i, 'uniprot']
        pdb = df.loc[i, 'pdb']
        chain = df.loc[i, 'chain']
        fn = f'{uniprot}/{pdb}.cif'
        gt_fn = gt_path_in + fn
        gt_fn_out = gt_path_out + fn
        pred_fn = pred_path_in + f'F-{uniprot}-F1-model_v3.cif'
        pred_fn_out = pred_path_out + fn

        # Make sure the uniprot directory exists
        utils.uniprot_dirs([gt_path_out, pred_path_out], uniprot)

        # Check if trimmed files already exist to save time
        if os.path.isfile(gt_fn_out) and os.path.isfile(pred_fn_out):
            print(f'{pdb} already trimmed')

            cfr = CifFileReader()

            # Get gt dataframe
            gt_obj = cfr.read(gt_fn, output='cif_dictionary', ignore=['_struct_conn'])
            gt_all_chains = pd.DataFrame.from_dict(gt_obj[pdb.upper()]['_atom_site'])
            gt = gt_all_chains[gt_all_chains['label_asym_id'] == chain].reset_index(drop=True)

            # Get gt_trim dataframe
            gt_trim_obj = cfr.read(gt_fn_out, output='cif_dictionary')
            gt_trim = pd.DataFrame.from_dict(gt_trim_obj[pdb.upper()]['_atom_site'])

            # Get pred dataframe
            pred_obj = cfr.read(pred_fn, output='cif_dictionary')
            pred = pd.DataFrame.from_dict(pred_obj[f'AF-{uniprot}-F1']['_atom_site'])

            # Get pred_trim dataframe
            pred_trim_obj = cfr.read(pred_fn_out, output='cif_dictionary')
            pred_trim = pd.DataFrame.from_dict(pred_trim_obj[f'AF-{uniprot}-F1']['_atom_site'])

        else:
            
            print(f'Trying {pdb}...')

            # Initiate reader object
            cfr = CifFileReader()

            # Create dataframe with gt atoms in desired chain
            gt_obj = cfr.read(gt_fn, output='cif_dictionary', ignore=['_struct_conn'])
            gt_all_chains = pd.DataFrame.from_dict(gt_obj[pdb.upper()]['_atom_site'])
            gt = gt_all_chains[gt_all_chains['label_asym_id'] == chain].reset_index(drop=True)

            # Create dataframe with pred atoms (pred file only contains our desired chain)
            pred_obj = cfr.read(pred_fn, output='cif_dictionary')
            pred = pd.DataFrame.from_dict(pred_obj[f'AF-{uniprot}-F1']['_atom_site'])

            print('Length of gt: ' + str(len(gt)) + ', Length of pred:' + str(len(pred)))

            # Find common atoms between files
            print(f'Comparing files for {pdb}...')
            atoms_pred, extra_atoms_gt = utils.compare_atoms(gt, pred)

            # Trim the files
            gt_trim, pred_trim = utils.drop_unshared_atoms(gt, pred, atoms_pred, extra_atoms_gt)

            print('Length of gt_trim: ' + str(len(gt_trim)) + ', Length of pred_trim: ' + str(len(pred_trim)))
            
            # Convert back to mmCIF-like dictionary
            gt_dict = gt_trim.to_dict(orient='list')
            gt_obj[pdb.upper()]['_atom_site'] = gt_dict

            pred_dict = pred_trim.to_dict(orient='list')
            pred_obj[f'AF-{uniprot}-F1']['_atom_site'] = pred_dict

            # Check whether the trimmed files are the same length
            assertion = utils.assert_equal_size(gt_trim, pred_trim)

            if assertion == True:
                print('Trimmed files are the same length')
            else:
                break
            
            if len(gt_trim) == 0:
                print(f'No common atoms found for {pdb}. Removing from dataframe...')
            else:
                print(f'Success! Creating trimmed files for {pdb}...')
                # Write trimmed files
                CifFileWriter(gt_fn_out).write(gt_obj)
                CifFileWriter(pred_fn_out).write(pred_obj)


        # Compile some information on the trimmed files
        trim_values_dict = utils.trim_stats(pdb, gt, gt_trim, pred, pred_trim)
        trim_values.append(trim_values_dict)
    
    # Add trim values to dataframe
    df_trim = pd.DataFrame(trim_values)
    df = df.merge(df_trim, on = 'pdb')

    # Drop any files that have no common atoms.
    df = df[df['gt_len'] != 0].reset_index(drop=True)


    return trim_values, df

def calculate_rmsd(gt_pdb_fn, pred_pdb_fn, complex_fn, region_1, region_2):
    '''Calculate rmsd between gt and pred regions and whole proteins
        Region1 is autoinhibitory region, region2 is domain
    '''
    
    # Rmsds is a dict of variable size with format {'complex_rmsd': ####, '1.0_aligned': ###, '1.0_comp': ###, ...}. Size based on number subregions.
    rmsds = {}
    utils.load_and_select \
        (gt_pdb_fn, pred_pdb_fn,
        region_1, region_2)

    try:
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
    except pymol.CmdException:
        print('Pymol error')

    return rmsds

def get_rmsds(df, gt_path, pred_path, complex_path):
    '''
    Calculate rmsds for each protein in df, aligning first on the autoinhibitory region (region 1) and then on the active region (region 2). Regions with
    multiple subregions are aligned and calculated separately, and then the average is taken.
    '''
    rmsd_info = []
    for i in range(len(df)):
        # Define pdb, filenames, region1, region2
        pdb = df.loc[i, 'pdb']
        uniprot = df.loc[i, 'uniprot']
        fn = f'{uniprot}/{pdb}.cif'
        region_1_dict = utils.create_region_dict(df.loc[i, 'region_1'], 1)
        region_2_dict = utils.create_region_dict(df.loc[i, 'region_2'], 2)
        percent_reg1 = df.loc[i, 'percent_region_1']
        percent_reg2 = df.loc[i, 'percent_region_2']
        gt_fn = gt_path + fn
        pred_fn = pred_path + fn
        complex_fn = complex_path + f'{pdb}_{uniprot}.pdb'

        if not os.path.isfile(gt_fn) and not os.path.isfile(pred_fn):
            print('No files found for ' + pdb + '. Skipping...')

            rmsd_dic = {'uniprot': uniprot,
                    'pdb': pdb,
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
                    'percent_region_1': percent_reg1,
                    'percent_region_2': percent_reg2}
            
            rmsd_info.append(rmsd_dic)

        else:
            print(f'Trying {pdb}...')
            rmsds = calculate_rmsd(gt_fn, pred_fn, complex_fn, region_1_dict, region_2_dict)

            # Define default values for columns to retain number of columns per row
            rmsd_dic = {'uniprot': uniprot,
                        'pdb': pdb,
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
                        'percent_region_1': percent_reg1,
                        'percent_region_2': percent_reg2}

            for key in rmsds:
                if key in rmsd_dic:
                    rmsd_dic[key] = rmsds[key]

            print('Success! Writing rmsds')
            rmsd_info.append(rmsd_dic)

    final_rmsds = utils.get_region_averages(rmsd_info)
    return final_rmsds

def two_state_proteins(df):
    # Find proteins with both open and closed conformations
    print('filtering for proteins with both open and closed conformations...')

    # Create categories for faster computation
    df_cat = df.astype({'uniprot': 'category', 'conformation': 'category'})

    # Find number of counts of open and closed for each UniProt ID.
    conf_cat = df_cat.groupby(['uniprot', 'conformation'], observed=True).size().reset_index(name='counts')

    # Find UniProt IDs with both open and closed conformations
    uniprot_cat = conf_cat.groupby('uniprot').size().reset_index(name='counts')

    # Create list of UniProt IDs with both open and closed conformations
    two_conf = uniprot_cat[uniprot_cat['counts'] == 2]['uniprot'].tolist()

    return two_conf

def calculate_disorder(df):
    # Calculate percent disorder of region 1 for each protein

    print('Calculating disorder...')

    df = utils.region_search_range(df)

    for i in range(len(df)):
        uniprot = df.loc[i, 'uniprot']
        region_1_res = df.loc[i, 'region_1 search']
        disorder_residues = []

        with open(f'./data/disorder_stats/sp{uniprot}.fasta.espritz', 'r') as f:
            residues = f.readlines()

            for residue in residues:
                if residue[0] == 'D':
                    disorder_residues.append(residues.index(residue) + 1)
        
        common_residues = utils.common_member(region_1_res, disorder_residues)
        percent_disorder = (len(common_residues) / len(region_1_res)) * 100

        df.loc[i, 'percent_disorder_1'] = round(percent_disorder, 3)
    
    df.drop(columns=['region_1 search'], inplace=True)
    return df

def mean_plddt(df, path):
    # Calculate mean plDDT for our region of interest.

    print('Calculating mean plDDT...')

    # Turn region ranges into list of residues
    df = utils.region_search_range(df).reset_index(drop=True)

    for i in range(len(df)):
        uniprot = df.loc[i, 'uniprot']
        fn = f'F-{uniprot}-F1-model_v3.cif'
        region_1_res = df.loc[i, 'region_1 search']
        region_2_res = df.loc[i, 'region_2 search']

        # Read in AF cif file
        cfr = CifFileReader()
        cif_obj = cfr.read(path + fn, output='cif_wrapper')
        cif_data = list(cif_obj.values())[0]

        # Create lists of residue plDDT values.
        region_1_plddt = []
        region_2_plddt = []

        for n in range(len(cif_data._ma_qa_metric_local.label_seq_id)):
            if int(cif_data._ma_qa_metric_local.label_seq_id[n]) in region_1_res:
                region_1_plddt.append(float(cif_data._ma_qa_metric_local.metric_value[n]))
            elif int(cif_data._ma_qa_metric_local.label_seq_id[n]) in region_2_res:
                region_2_plddt.append(float(cif_data._ma_qa_metric_local.metric_value[n]))

        # Calculate mean plDDT for each region
        region_1_mean = np.mean(region_1_plddt)
        region_2_mean = np.mean(region_2_plddt)

        df.loc[i, 'region_1_mean_plddt'] = round(region_1_mean, 3)
        df.loc[i, 'region_2_mean_plddt'] = round(region_2_mean, 3)

    df.drop(columns=['region_1 search', 'region_2 search'], inplace=True)
    return df

def mean_paes(df, path):
    # Calculate the average pae for region 1 to region 1, region 2 to region 2, and region 1 to region 2

    print('Calculating mean pae...')

    for i in range(len(df)):
        uniprot = df.loc[i, 'uniprot']
        fn = f'F-{uniprot}-F1-predicted_aligned_error_v3.json'
        region_1 = df.loc[i, 'region_1']
        region_2 = df.loc[i, 'region_2']

        # Region bounds are in the format [start, end] for each region. Regions with multiple sections look like [[start, end], [start, end], ...]
        reg1_bounds = utils.region_bounds(region_1)
        reg2_bounds = utils.region_bounds(region_2)

        reg1_array = np.array(reg1_bounds)
        reg2_array = np.array(reg2_bounds)

        # Read in json file
        prot_array = utils.pae_from_json(path, fn)

        '''
        We want means of reg1 compared against reg1, reg1 compared against reg2, and reg2 compared against reg2.
        '''

        mean11 = utils.calculate_pae_mean(prot_array, reg1_array, reg1_array)
        mean12 = utils.calculate_pae_mean(prot_array, reg1_array, reg2_array)
        mean22 = utils.calculate_pae_mean(prot_array, reg2_array, reg2_array)
        

        df.loc[i, 'mean_pae_1_1'] = round(mean11, 3)
        df.loc[i, 'mean_pae_1_2'] = round(mean12, 3)
        df.loc[i, 'mean_pae_2_2'] = round(mean22, 3)
    
    return df