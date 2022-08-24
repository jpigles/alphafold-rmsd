import pandas as pd


df_samples = pd.read_csv('data/proteins_pdb.csv', sep = ',')
#change sample names to use Uniprot accession codes, not gene names.
sample_names = list(df_samples['Uniprot_ID'])
#Establishing the end directory
rule all:
    input:
        # Here I place my final analysis. What is my final analysis?

rule pdb_files:
    input:
        'data/proteins_pdb.csv'    
    output:
        directory(expand('data/structures/{protein}', protein = sample_names))
    script:
        'get_pdb_files_aa.py'
        # 'Note: this must be updated because RCSB and/or Uniprot have changed their API.'
    
rule best_pdb_files:
    input:
        'data/proteins_pdb.csv'
        expand('data/structures/{protein}', protein = sample_names) 
    output:
        'data/proteins_pdb_best.tsv'
    script:
        'scripts/best_pdb_aa.py'

rule interface_residues:
    input:
        #Takes in the .tsv file generated from the last rule.
    output:
        #Creates another .tsv file containing the interface residues (see Jorge's sop)
    script:
        #scripts/get_interface_all.py. Should be fine.

rule compare_structures:
    """
    2 ways to compare
    1st: Directly compare the coordinates between the AF and PDB files. This gives us a percentage of residues that are within a certain distance of each other.
    Then we look at how many of the interface residues are within that cutoff distance. Or we look at the average distance between residues for each protein. Graph that against Alphafold confidence score.
    2nd: A follow-up analysis to determine the interface residues in the AF and PDB files separately (using Jorge's script) then compare. Acts as a secondary verification.
    
    So: what is the cutoff distance?
        Consult "Advances and pitfalls of protein structural alignment", Hitomi & Holm, https://doi.org/10.1016/j.sbi.2009.04.003. Possibilities: least-squares superimposition
        (rigid or flexible) or distance difference matrix.
            If we use LSS, what residues do we center on? Active sites? Do I have to write a program to align them, or can I use existing software? https://doi.org/10.1093/nar/gkaa443 (FATCAT), https://doi.org/10.1093/protein/11.9.739 (CE)
            PDB has a page with pairwise structure alignment resources: https://www.rcsb.org/docs/tools/pairwise-structure-alignment#:~:text=Structure%20alignment%20is%20a%20valuable,by%20standard%20sequence%20alignment%20techniques. 
            I don't understand what distance difference matrices are."""
