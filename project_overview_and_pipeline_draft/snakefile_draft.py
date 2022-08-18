
#Establishing the end directory
rule all:
    input:
        # Here I place my final analysis. What is my final analysis?

rule pdb_files:
    input:
        #List of desired PDB files.
    output:
        #Directory containing all of the PDB files.
    script:
        #scripts/get_pdb_files.py or something similar
        'Note: this must be updated because RCSB and/or Uniprot have changed their API.'
    
rule best_pdb_files:
    input:
        #Directory containing all of the PDB files, as well as file with information on proteins.
    output:
        #In Jorge's script it's a .tsv file, so I guess stick with that.
    script:
        #scripts/best_pdb.py. Should be fine.

rule interface_residues:
    input:
        #Takes in the .tsv file generated from the last rule.
    output:
        #Creates another .tsv file containing the interface residues (see Jorge's sop)
    script:
        #scripts/get_interface_all.py. Should be fine.

rule compare_alphafold:
    #This is where it gets interesting. What I should do is, I should get interface proteins from both the PDB and the AlphaFold files.
    #Or I could compare the AlphaFold and PDB files directly.
        #If I do this, then I could generate a list of residues for each protein where AlphaFold and PDB
        # do not match, then compare that to the interface residues and see if there's any overlap.
    #What would be the output of such analysis? A measure of how well AlphaFold matches autoinhibited structures versus its own measure of correctness?
    #And then we could follow del Alamo et al (March 2022) and look at alternative conformational states and see how AlphaFold does there.
    #This could potentially lead to the revelation that sampling multiple conformations allows AlphaFold to capture all possible states of autoinhibited proteins.
