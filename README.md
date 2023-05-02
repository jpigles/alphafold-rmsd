# autoinhibition_and_alphafold2
This is my first project in my Masters. This aims to evaluate the performance of AlphaFold2 on autoinhibitory interfaces in proteins exhibiting autoinhibitory functionality.

Please note that pdb file 7w39 was manually removed due to file size limitations when converting from mmcif to pdb. In addition, 3vd8, 4d7r, 7apj, and 2vgq were removed for being chimeric proteins that did not contain one of the inhibitory regions for the desired protein. 3nhn and 1qcf had an error in author-assigned sequence number where sequence number 115 was skipped but the corresponding residue was not, increasing the auth seq numbers of all residues after 114 by 1. This error was manually fixed.
For every pdb file with NMR as the method of collection, the first model was chosen for rmsd calculations.

There are quite a few things I may need to go back and fix, and could potentially add for reference.
    - 3nhn is missing residues 425-514, which are defined in its 2nd range in proteins_pdb_best.tsv (alongside 76-136), but these missing residues were not taken into account by the quality control test. Uncertain why. I may have to go modify that code. In addition, that range is perhaps not what should have been included. In SRC, another kinase with SH2 and SH3 domains, we were more concerned with the C-terminus, which in HCK is 515-526.
    - I could include experimental structures with one region present and one region missing. It would still give a (partial) estimation of AlphaFold2's efficacy.
    - Could add 2BDW (C. elegans CamKII). It's got the full structure of the protein in autoinhibited form. We have 2wel, which is human.
    - LRR domain may not actually be involved in autoinhibition in the NLRP3 inflammasome (6NPY). 6NPY is the inactive conformation, uncertain if it's autoinhibited. 
    - 1wao (P53041) was missed and it contains a nearly-whole protein.
    - P63086 (rat MAPK1) is entirely missing from my selection of proteins, even though it shares 99.4% sequence similarity with human MAPK1 (P28482).
        - As an aside to that, it's hard to define whether a structure for P28482 is inhibited or autoinhibited.
    - 3JWN (A0A0R4I961) appears to contain the autoinhibited form of FimH, but it is missing from my list. The only file I have for FimH is 6GU0, but there are several more.
    - 6LOJ has the active form of IpaH9.8-LRR (compared to autoinhibited full-length 6LOL).
    - Should look into 5E7J (active form of O00571).
    - 2R09 and 2R0D are missing for O08967.
    - 4CKI is missing for P07949.
    - 4CKJ is in the "open", active position, while the AF prediction is in the "closed", autoinhibited position (10.1016/j.molcel.2014.01.015)
    - 1PKG (active) is missing from P10721.
    - The autoinhibitory region of Heat shock cognate 71 kDa protein (P11142) acts on its nucleolar targeting signal and not on its active site (very interesting). It's difficult to say whether the protein is in the autoinhibited or active forms based on the given structures, although we may assume that it is autoinhibited because nucleolar localization due to heat shock is initiated by phosphorylation, possibly of Thr265. In our models, we only have a phosphate in the active site at the phosphate binding pocket.
    - 2Y1M and 2Y1N are missing for P22681, even though they're the full protein :/. 
    - I am entirely missing P63086 (MAPK1 for Rattus norvegicus) and P63085 (MAPK1 for Mus musculus). I really think I should have these...
    - CRKII (P46108) and CrkL (P46109) are weird. The binding site of the SH2 domain in CrkL is occluded compared to CRKII, but phosphorylation leads to further inhibition (same as CRKII). Meanwhile, the binding site of the SH3 domain is entirely open in CrkL and can interact with ligands at any time. In this sense, it is hard to characterize CrkL as "autoinhibited" versus "active".
    - 2AYN is missing for P54578.
    - 2J0L is missing for Q00944
    - P49137 is missing 1NY3 and 1NXK (active forms) and 1KWP (autoinhibited form).


THIS IS IMPORTANT. I realize now that the reason I'm missing so many files is due to the prune_ids step, which simply prunes any files where our protein of interest is present in more than one chain. But that problem would be very simple for us to fix, so I think I'm going to have to go back through and devise a way to include those files while restricting the chains used to just one. It would also give me an opportunity to re-organize my code, because much of it is a mess right now. 