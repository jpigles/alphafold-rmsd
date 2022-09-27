import pymol

filenames = ['../data/input/alphafold_files_cif/F-A0A0R4I961-F1-model_v3.cif']

for filename in filenames:
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(filename)
        pymol.cmd.save(filename.replace('cif', 'pdb'))