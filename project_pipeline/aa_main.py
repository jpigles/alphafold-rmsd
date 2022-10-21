
import project_pipeline.aa_module as aa_module
import argparse
import configargparse

from parser import parse_args


def process_input(args):
    ''' Combine multi chains with poly-g linker
          and merge multiple fasta files into one.
        Also store residue id where each single chain starts
          and chain name of each protein
    '''
    pipeline = aa_module.pipeline('init_input_procs', args)
    pipeline.process_input()


def process_input_from_fasta(args):
    ''' Combine multi chains with poly-g linker
          and merge multiple fasta files into one.
        Also store residue id where each single chain starts
          and chain name of each protein
    '''
    pipeline = aa_module.pipeline('init_input_procs_fasta', args)
    pipeline.process_input_from_fasta()

    '''
    We end up with a dictionary of starting IDs (the number of the sequence that the chain starts at) for each chain in each PDB file, 
    border residues for each chain in each PDB file, and the identity of each chain (A, B, ...) in each PDB file, in binary files. Then we also get the
    fasta files that contain the "full" sequence patched together with poly-G linkers.
    '''


def process_output(args):
    ''' Remove linker from predicted pdb and reorder
          residue number as per gt pdb
        Calculate rmsd and export to csv file
    '''
    pipeline = aa_module.pipeline('init_output_procs', args)
    pipeline.process_output()


if __name__ == "__main__":

    parser = configargparse.ArgumentParser() #Creates the ArgumentParser object using a package that has more support for config files.
    config = parse_args(parser) #This is the parse_args() method, so it will inspect the command line and convert the arguments into the proper types. Also note that this uses the function from parser.py and not from the argparse package.
    # This config now contains additional command line args, hardcoded arg values, and a defined path.
    args = argparse.Namespace(**config) #** is unpacking config. Normally it would be a dictionary, but this stores everything in Namespace in a different way. So when we run $python main.py --config config/test.ini, it knows how to interpret those commands. I think.

    for operation in args.operations:
        if operation == 'process_input':
            if args.from_fasta:
                process_input_from_fasta(args)
            else:
                process_input(args)
        elif operation == 'process_output':
            process_output(args)
        elif operation == 'assert_fasta':
            assert_fasta(args)