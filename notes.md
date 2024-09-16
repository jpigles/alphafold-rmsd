# autoinhibition_and_alphafold2
This is my first project in my Masters. This aims to evaluate the performance of AlphaFold2 on autoinhibitory interfaces in proteins exhibiting autoinhibitory functionality.

Please note that pdb file 7w39 was manually removed due to file size limitations when converting from mmcif to pdb. In addition, 3vd8, 4d7r, 7apj, and 2vgq were removed for being chimeric proteins that did not contain one of the inhibitory regions for the desired protein. 3nhn and 1qcf had an error in author-assigned sequence number where sequence number 115 was skipped but the corresponding residue was not, increasing the auth seq numbers of all residues after 114 by 1. This error was manually fixed.
For every pdb file with NMR as the method of collection, the first model was chosen for rmsd calculations.

Bash command to download all necessary alphafold files
```
while read p; do gcloud storage cp "$p" ~/file/path/; done <af_file_locations_single_domains.csv 
```

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
    - 4CKI is missing for P07949. Also, the inhibitory section appears to be wrong. Instead of 900-909, it's more like 730-740 (see 24560924 and the open GRL vs closed GRL). Although, I'm thoroughly confused, because I'm seeing conflicting information about the same conformation. 4CKJ is in the "open", active position with its GRL, but 2IVT is supposedly in the active conformation while still having the GRL be in the closed conformation.(10.1016/j.molcel.2014.01.015). I will probably go with the "closed" conformer (and the AF model) being auto-inhibited and the "open" conformer being active, but 2IVT can be considered "active" (see 2pvf for active comparison and 2psq for apo/autoinhibited). In fact, I think the 4CKJ authors got it wrong; both the open and closed GRL are active, while the AF model shows the autoinhibited form (which is taken by 2psq)
    - 1PKG (active) is missing from P10721.
    - The autoinhibitory region of Heat shock cognate 71 kDa protein (P11142) acts on its nucleolar targeting signal and not on its active site (very interesting). It's difficult to say whether the protein is in the autoinhibited or active forms based on the given structures, although we may assume that it is autoinhibited because nucleolar localization due to heat shock is initiated by phosphorylation, possibly of Thr265. In our models, we only have a phosphate in the active site at the phosphate binding pocket.
    - 2Y1M and 2Y1N are missing for P22681, even though they're the full protein :/. 
    - I am entirely missing P63086 (MAPK1 for Rattus norvegicus) and P63085 (MAPK1 for Mus musculus). I really think I should have these...
    - CRKII (P46108) and CrkL (P46109) are weird. The binding site of the SH2 domain in CrkL is occluded compared to CRKII, but phosphorylation leads to further inhibition (same as CRKII). Meanwhile, the binding site of the SH3 domain is entirely open in CrkL and can interact with ligands at any time. In this sense, it is hard to characterize CrkL as "autoinhibited" versus "active".
    - 2AYN is missing for P54578.
    - 2J0L is missing for Q00944
    - P49137 is missing 1NY3 and 1NXK (active forms) and 1KWP (autoinhibited form).
    - 4D7Q is chimeric.
    - Can also pull 5UPD for Q9BZ95.
    - Can add Q96L73 with autoinhibited structure 3OOI (PMID 21196496). 
    - Can add Q96028 with autoinhibited structure 5LSU (PMID 27571355).
    - For Q9UM73, 5IUI is in the "dfg-in" conformation (appears to be the autoinhibited conformation) and 5IUG is in the "dfg-out" conformation. 5IUG is missing density for region 1, but it's not too bad: ~ 70%.

------
Notes for correcting mmcif files.
- 2ptk has the incorrect db_align_begin. It's 80 when it should be 81.
- 2rgn is chimeric. Exclude.
- 3ig3, db_align_beg is wrong, auth_seq_align_beg is correct. Change to that.
- 6yr5 only has 14 amino acids. Removed from proteins_by_pdb.
- 2d9f has the second isoform as per Uniprot, so it is one amino acid off from the Alphafold file (AF G183 -> 2d9f G183,S184). db_align_beg is 84, seq_align_beg is 1.
- 1a0n seq_align_beg = 1, db_align_beg = 80.
- 6amv contains a non-canonical sequence, P00519-2, but this difference only extends until residue 45. Afterwards they're the same. seq_align_beg = 46, db_align_beg = 27.
- 6amw same as 6amv.
- 1hct, db_align_beg = 144.
- In the case of 7epu, it's using the sequence information for the wrong chain. ~~I'll have to fix my code to account for that. This may fix some other ones.~~ This has been done
- 2mf9 has sequence Q14318-2.
- 7tfc had an issue where the first chain returned in the query was 7tfc.AA, which meant that the chain was selected as 'A' due to slicing [:4]. Resolved by using .split() to select pdb id and chain.
- For 2m0v, I discovered that the calculations for percent completeness in each region was counting the region in every entity in the NMR files, thus greatly over-representing the completeness (it was 2000 in some cases). I changed the script to only select model 0 in each file, but I'm not sure if this will cause problems later.
- 4m8t, the author sequence is correct, the database is off by one. 
- 4mao, same as 4m8t.
- 6nr0, author sequence is correct.
- 6zei _very annoyingly_ has two different proteins under chain A, so even if I use the chain to find the index of the seq_begin, I will _still_ get the wrong seq_begin. So I switched to using the UniProt ID to get the index. Hopefully this is more robust. And actually, besides that fact, we only have 16 amino acids from our protein of interest, so it's pretty useless anyway. It appears to be a chimeric protein. Unfortunately, it will always get past the percent counting step, but what can you do.
- 3bhh is completely correct, and the full structure is there, so I'm not sure what's going on with the trim_values Looks like it got corrected twice for some reason. 
- 1na6, author_seq is correct.
- 2ptk, author_seq is correct.
- 2xp2, auth_seq is correct. db_align_beg is _way_ off.
- 2dx1, both db_align_beg and auth_seq are wrong. Difference between db and seq_align_beg is 129.
- 4ped, auth_seq is correct.
- 2yfx, auth_seq is correct.
- 1p14, db_align_beg should start at 1004.
- 2ojj, db_align_beg should start at 42.
- 2oji, same as 2ojj.
- 2ojg, same as 2ojj.
- 7us1, auth_seq is correct, but I also had to fix pdbx_db_accession to match the Uniprot id (S4X0T1 to Q9JK66)
- 1s9j, auth_seq is correct.
- 2v7o, wrong Uniprot id in pdbx_db_accession (apparently an older entry that was merged into the current entry) (Q8N4I3 to Q13555)
- 2vgq is chimeric, remove. (Could I automate removal of chimeric entries?)
- 2wel, db_align_beg should start at 13 for Q13557. 
- 2c0i, db_align_beg should be 81.
- 2c0o, same as 2c0i.
- 5tj4 is chimeric, remove.
- Not sure what went wrong with 2w4o. 
- 1qkr, auth_seq is correct. 
- 4r7h, auth_seq is correct. 
- 7m8e, incorrect Uniprot ID. C3TR27 should be P60240.
- 5tib is chimeric, remove
- 5tj2 is chimeric, remove.
- 4d7r is chimeric, remove.
- 2src, db_align_beg should start on 86.
- 1ksw, db_align_beg should start on 86.
- 2h8h, db_align_beg should start on 2.
- 4wxx, auth_seq is correct.
- 1jpa, db_align_beg should start on 587. 
- 1fmk, db_align_beg should start on 86.
- 6nif is a fusion protein, remove.
- 6nbs, seq_align_beg should start on 8.
- 1opl, db_align_beg should start on -18.
- 2fo0 has the same issue as 6AMV. On the first row, change db_align_beg to 43 and the Uniprot ID from P00519-2 to P00519.
- 5d7q, auth_seq is correct.
- 1tr2, auth_seq is correct.
- 1st6, auth_seq is correct.
- 1gnv, problem is that auth_seq correctly skips 9 residues, but seq_align_beg does not. Auth_seq goes from 74 to 84, seq_align_beg goes from 74 to 75. Not certain if there's a real fix to that. db_align_beg is also incorrect. I think I'll go with the larger fragment. Fixed db_align_beg on first row to 117.
- 3h2u, auth_seq is correct. Also fixed B4DTM7 to P18206
- 6nji, like 1gnv, skips a few residues, but the second fragment is very small. First row db_align_beg starts at 380.
- 6njj fixed as per 6nji.
- 6njh fixed as per 6nji.
- 6pst, last row incorrect Uniprot ID fixed (to P00579)
- 1g83, auth_seq is correct.
- 2ayo, db_align_beg starts at 91.
- 2ayn, same as 2ayo.
- 6t58 is a fusion protein, best to eliminate it.
- 2xkx is just... one row per amino acid? What? Eliminate it.
- 3fi7, same issue as 1gnv. Second fragment is larger, keep it. Second row db_align_beg starts at 64.
- 6fek, auth_seq (which is correct) skips from 825 to 841 but seq_align_beg goes 126 to 127. Uncertain how or if I should fix this.
- 3krj, fusion protein, best to eliminate.
- 4zyn, seq_align_beg skips from 78 to 100, but auth_seq skips from 73 to 140. Change db_align_beg to 46.
- 7syf is a fusion protein, remove
- 7am4, auth_seq skips from 75 to 85 but seq_align_beg stays consistent (75 to 76). Change db_align_beg to 10.
- 7am8, same as 7am4.
- 7am3, same as 7am4.
- 5ox2, same as 7am4.
- 7am5, same as 7am4.
- 7am7, same as 7am4.
- 7am6, same as 7am4.
- 1dui, same as 7am4.
- 1sue, same as 7am4.
- 1gns, similar issue. seq_align_beg goes 71 to 72 while auth_seq goes 74 to 84. Change db_align_beg to 13.
- 1sua, same as 7am4.
- 7apj is chimeric, remove.
- 3t6p, seq_align_beg goes from 99 to 113, auth_seq goes from 362 to 386 (263 difference to 273 difference). Change db_align_beg to 275.
- 6lvs, seq_align_beg goes 144 to 156, auth_seq goes 219 to 240. I wrote a short script to fix all of the weird skip files listed above. db_align_beg changed to 101.
- 7dwb has 7 added amino acids, 173 to 180. According to paper, it's an inserted Strep-tag. 

For Espritz, I used the Disprot prediction type and the Best Sw decision threshold.


Uniprot access numbesr to add:
Q922S4 (PDE2A_MOUSE) (for comparison to O76074)
