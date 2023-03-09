# autoinhibition_and_alphafold2
This is my first project in my Masters. This aims to evaluate the performance of AlphaFold2 on autoinhibitory interfaces in proteins exhibiting autoinhibitory functionality.

Please note that pdb file 7w39 was manually removed due to file size limitations when converting from mmcif to pdb. In addition, 3vd8, 4d7r, 7apj, and 2vgq were removed for being chimeric proteins that did not contain one of the inhibitory regions for the desired protein. 3nhn and 1qcf had an error in author-assigned sequence number where sequence number 115 was skipped but the corresponding residue was not, increasing the auth seq numbers of all residues after 114 by 1. This error was manually fixed.
For every pdb file with NMR as the method of collection, the first model was chosen for rmsd calculations.

There are quite a few things I may need to go back and fix, and could potentially add for reference.
    - 3nhn is missing residues 425-514, which are defined in its 2nd range in proteins_pdb_best.tsv (alongside 76-136), but these missing residues were not taken into account by the quality control test. Uncertain why. I may have to go modify that code. In addition, that range is perhaps not what should have been included. In SRC, another kinase with SH2 and SH3 domains, we were more concerned with the C-terminus, which in HCK is 515-526.