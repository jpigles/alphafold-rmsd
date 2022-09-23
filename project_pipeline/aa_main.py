
import os
import project_pipeline.aa_utils as aa_utils
import pickle
import argparse
import configargparse

from pymol import cmd
from project_pipeline.aa_parser import parse_args
from os.path import join, exists

def process_input(args):
    ''' Combine multi chains with poly-g linker
          and merge multiple fasta files into one.
        Also store residue id where each single chain starts
          and chain name of each protein
    '''
    #if args.remove_x:


    chain_start_resid_ids = aa_utils.generate_fasta_from_pdb\
        (args.input_pdb_dir, args.linker_fasta_dir, args.pdb_ids, args.n_g,
         remove_x=args.remove_x)

    gt_chain_bd_resid_ids = aa_utils.read_bd_resid_id_all \
        (args.input_pdb_dir, args.pdb_ids)

    with open(args.chain_start_resid_ids_fn, 'wb') as fp: #...data/input/idr_84/poly_g_6/chain_start_ids.pkl
        pickle.dump(chain_start_resid_ids, fp)

    with open(args.gt_chain_bd_resid_ids_fn, 'wb') as fp: #...data_input/idr_84/poly_g_6/gt_chain_bd_ids.pkl
        pickle.dump(gt_chain_bd_resid_ids, fp)

    # get chain names and store locally
    chain_names = aa_utils.get_chain_identifiers_all \
        (args.pdb_ids, args.input_pdb_dir) #Another dictionary of chain IDs associated with PDB files.
    with open(args.chain_names_fn, 'wb') as fp:
        pickle.dump(chain_names, fp)

    '''
    We end up with a dictionary of starting IDs (the number of the sequence that the chain starts at) for each chain in each PDB file, 
    border residues for each chain in each PDB file, and the identity of each chain (A, B, ...) in each PDB file, in binary files. Then we also get the
    fasta files that contain the "full" sequence patched together with poly-G linkers.
    '''


def process_output(args):
    ''' remove linker from predicted pdb and reorder
          residue number as per gt pdb'
        calculate rmsd and export to csv file
    '''
    #Is gt ground-truth? 

    with open(args.chain_names_fn, 'rb') as fp:
        chain_names = pickle.load(fp)

    with open(args.chain_start_resid_ids_fn, 'rb') as fp:
        chain_start_ids = pickle.load(fp)

    with open(args.gt_chain_bd_resid_ids_fn, 'rb') as fp:
        gt_chain_bd_ids = pickle.load(fp)

    rmsds = [args.rmsd_names]
    for pdb_id in args.pdb_ids:
        print() #?
        print(pdb_id, chain_names[pdb_id], gt_chain_bd_ids[pdb_id])

        dir = join(args.output_dir, pdb_id + '.fasta') #data/output/idr_84_af_full/poly_g_6/{PDB_ID}.fasta
        if not exists(dir):
            print(f'directory {dir} doesn\'t exist')
            continue

        pred_fn = join(dir, args.pred_fn_str) #? output/idr_84_af_full/poly_g_6/{PDB_ID}.fasta.ranked_0.pdb? This is the predicted filename, so in my case it will be my Alphafold file
        gt_pdb_fn = join(args.input_pdb_dir, pdb_id+ '.pdb') #.cif
        pred_removed_linker_fn = join(dir, args.removed_linker_fn_str) #... {PDB_ID}.fasta.ranked_0_removed_linker.pdb?
        aa_utils.remove_linker(pred_fn, pred_removed_linker_fn, args.n_g,
                            chain_names[pdb_id], chain_start_ids[pdb_id],
                            gt_chain_bd_ids[pdb_id])

        # prune unknown aa from sequences
        if args.prune_unknown:
            assert(False)

        complex_fn = join(args.output_dir, pdb_id + '.fasta', args.complex_fn_str)
        cur_rmsds = aa_utils.calculate_rmsd\
            (gt_pdb_fn, pred_removed_linker_fn, complex_fn,
             chain_names[pdb_id], backbone=args.backbone,
             remove_hydrogen=args.remove_hydrogen)
        cur_rmsds.insert(0, pdb_id)
        rmsds.append(cur_rmsds)

    aa_utils.write_to_csv(rmsds, args.rmsd_fn)
    cmd.quit()

def assert_fasta(args):
    ''' check if aa sequence from processed prediction pdb match
          sequence from gt pdb
    '''
    for id in args.pdb_ids:
        pdb_fn = join(args.input_pdb_dir, id + '.pdb')
        pred_fn = join(args.output_dir, id + '.fasta', args.pred_fn)
        seq_1 = aa_utils.read_residue_from_pdb(pdb_fn)
        seq_2 = aa_utils.read_residue_from_pdb(pred_fn)
        #print(seq_1)
        #print(seq_2)
        assert(seq_1 == seq_2)

if __name__ == "__main__":

    parser = configargparse.ArgumentParser() #Creates the ArgumentParser object using a package that has more support for config files.
    config = parse_args(parser) #This is the parse_args() method, so it will inspect the command line and convert the arguments into the proper types. Also note that this uses the function from parser.py and not from the argparse package.
    # This config now contains additional command line args, hardcoded arg values, and a defined path.
    args = argparse.Namespace(**config) #** is unpacking config. Normally it would be a dictionary, but this stores everything in Namespace in a different way. So when we run $python main.py --config config/test.ini, it knows how to interpret those commands. I think.

    for operation in args.operations:
        if operation == 'process_input':
            process_input(args)
        elif operation == 'process_output':
            process_output(args)
        elif operation == 'assert_fasta':
            assert_fasta(args)
