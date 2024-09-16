This document contains the guide to collecting, creating, and analyzing the interfaces of AlphaFold2 structural predictions.

## Setting up your virtual environment
First, you'll have to (install Anaconda)[https://conda.io/projects/conda/en/latest/user-guide/install/index.html] if you don't already have it.
Navigate to the folder containing this repository.

```
$cd ~/some/directory/path/autoinhibition_protein_data
```

Create a new environment using the environment.yml file.

```
$ conda env create -f rmsd_snek.yml
```

Activate the environment

```
$ conda activate myenv
```

## Downloading AlphaFold2 structures
Once you have the list of proteins of interest, navigate to project_pipeline and run scripts/download_alphafold.py, making sure to pass it the file.

```
$ cd ./project_pipeline
$ python scripts/download_alphafold.py -f data/proteins_of_interest.csv
```

It will download the appropriate AlphaFold2 structures into the data folder and create a file with the proteins it did not find an AlphaFold2 structure for.

## Creating AlphaFold2 structures
For each protein that requires a novel AlphaFold2 structure, you will first need to generate a multiple sequence alignment (MSA) for it in the format of an .a3m file, which you will do using the (ColabFold notebook)[https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=G4yBrceuFbf3].

First, navigate to [uniprot.org] and search using the UniProt accession ID of your protein. Copy the sequence and paste it into the ```query_sequence``` field in the ColabFold notebook (I recommend also changing the ```jobname``` to the UniProt ID). Run the notebook. Click the folder icon to the left. A folder should appear with the jobname. As the MSA is one of the first steps that the notebook runs, the .a3m file should appear relatively quickly. Stop the notebook once the a3m file appears, download it, and place it into a designated MSA folder. Repeat this process for all your proteins.

You will now be using these a3m files to create the structures.

### Copying a repository