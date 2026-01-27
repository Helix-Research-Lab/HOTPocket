## The Human Omnibus of Targetable Pockets
#### Kristy Carpenter & Russ Altman

![HOTPocket logo](hotpocket.png)

This repository accompanies the manuscript "The Human Omnibus of Targetable Pockets." HOTPocket = "The Human Omnibus of Targetable Pockets", aka the dataset. *hotpocketNN* = the method.

### Citation
> K. A. Carpenter and R. B. Altman, “The Human Omnibus of Targetable Pockets,” J. Cheminformatics, vol. 17, no. 1, p. 180, Dec. 2025, doi: 10.1186/s13321-025-01125-x.


In addition to citing our manuscript, please cite **all seven** constituent methods directly:
> B. Huang and M. Schroeder, “LIGSITEcsc: predicting ligand binding sites using the Connolly surface and degree of conservation,” BMC Struct. Biol., vol. 6, no. 19, 2006, doi: 10.1186/1472-6807-6-19.
>
> R. Krivák and D. Hoksza, “P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure,” J. Cheminformatics, vol. 10, no. 1, p. 39, Dec. 2018, doi: 10.1186/s13321-018-0285-8.
>
> V. Le Guilloux, P. Schmidtke, and P. Tuffery, “Fpocket: An open source platform for ligand pocket detection,” BMC Bioinformatics, vol. 10, no. 1, p. 168, Dec. 2009, doi: 10.1186/1471-2105-10-168.
>
> A. Meller et al., “Predicting locations of cryptic pockets from single protein structures using the PocketMiner graph neural network,” Nat. Commun., vol. 14, no. 1, p. 1177, Mar. 2023, doi: 10.1038/s41467-023-36699-3.
>
> P. A. Ravindranath and M. F. Sanner, “AutoSite: an automated approach for pseudo-ligands prediction—from ligand-binding sites identification to predicting key ligand atoms,” Bioinformatics, vol. 32, no. 20, pp. 3142–3149, Oct. 2016, doi: 10.1093/bioinformatics/btw367.
>
> W. Tian, C. Chen, X. Lei, J. Zhao, and J. Liang, “CASTp 3.0: computed atlas of surface topography of proteins,” Nucleic Acids Res., vol. 46, no. W1, pp. W363–W367, Jul. 2018, doi: 10.1093/nar/gky473.
>
> S. Wang, J. Xie, J. Pei, and L. Lai, “CavityPlus 2022 Update: An Integrated Platform for Comprehensive Protein Cavity Detection and Property Analyses with User-friendly Tools and Cavity Databases,” J. Mol. Biol., vol. 435, no. 14, p. 168141, Jul. 2023, doi: 10.1016/j.jmb.2023.168141.


### Installation
1. Clone this repository: `git clone https://github.com/Helix-Research-Lab/HOTPocket`
2. Download the data from [Zenodo](https://zenodo.org/records/17653326) and move all desired files into this directory (`HOTPocket`)
3. Install dependencies:
	1. Set up `hotpocket` environment (required): `conda env create --name hotpocket --file=environment.yml`
	2. Install dependencies for ESM2 (required to run *hotpocketNN* on new structures): follow instructions [here](https://github.com/facebookresearch/esm)
4. Run setup script: `./setup.sh`

### Usage
#### Browsing HOTPocket
We provide the 2.4 million predicted pockets across the whole human proteome as described in the accompanying manuscript in the `proteome_hotpocket_embs_NN_score.csv` file available on [Zenodo](https://zenodo.org/records/16891050). An alternate version computing using Feature Set C (see manuscript for feature set definitions) is available as `proteome_hotpocket_comb_NN_score.csv`, also on Zenodo. For more information on how to use these dataframes, see the [data dictionary](data_dictionary.md). 

:construction: Coming soon: Website for easier browsing! :construction:

#### Running *hotpocketNN* on a novel structure(s)
If your structure of interest does not have precomputed *hotpocketNN* pockets, you can compute them yourself following these steps:
1. Get your structure file(s), sequence(s), and input dataframe

    You must have PDB structure files of any proteins you would like to run *hotpocketNN* on. Generate a fasta file containing the exact sequence represented in the PDB file using the `pdbs_to_fasta` function in `utils.py`.

    For the next step, you will need a dataframe of information about your structure. This should take a similar form as `data/proteome_structure_files_canonical.csv`; there needs to be `uniprot id`, `structure type`, `structure id`, and `structure file` columns. See the [data dictionary](data_dictionary.md) for descriptions of these columns. Note that the `structure file` column in the version of `proteome_structure_files_canonical.csv` from Zenodo has `nan` for all rows -- to use this dataframe, please fill in the path to your local copies of the structure files. 

2. Run constituent pocket-finding methods

    `assembly/make_db.py` contains all functionality needed to run the 7 constituent pocket-finding methods on a novel structure. For AutoSite, CASTp, CavitySpace, Fpocket, LIGSITEcs, and P2Rank, you must use the `hotpocket` environment. For PocketMiner, you must be using the `pocketminer` environment. You may either run the script from the command line with the corresponding flag (e.g. for AutoSite, `python make_db.py --autosite`, for PocketMiner, `python make_db.py --pocketminer`). You can also run the functions individually from the console with your desired parameters (e.g. if you do not need to build the directory in `../data`, you can call the function with `makedir=False`; we have also made some of these functions more easily parallelizable).

    You may also obtain BioLiP2 annotations for your structures of interest by running `python make_db.py --knowndb` or calling the `make_known_df()` function, passing along the location of your (perhaps more recent than our) `BioLiP.txt` file.

4. Get ESM2 embedding

    Use the `esm-extract` functionality as described in the [ESM2 respository](https://github.com/facebookresearch/esm) to generate ESM2 embeddings. You must have the `esmfold` environment activated. Note that you may need to increase the `--truncation_seq_length` parameter if your sequence is longer than the default.

    Run: `esm-extract esm2_t36_3B_UR50D [FASTA] [OUTPUT DIRECTORY] --include per_tok --truncation_seq_length [MAXIMUM SEQUENCE LENGTH]`

5. Run *hotpocketNN*

    *hotpocketNN* will take the candidate pockets predicted by the constituent pocket-finding methods and the ESM2 embeddings and generate a score, with which you can filter and/or rank the output pockets.

    Navigate into the `assembly` directory and run `python run_hotpocketnn.py --in_strucid [INPUT DATAFRAME FROM STEP 1] --embdir [EMBEDDINGS DIRECTORY FROM STEP 3] --pdbdir [DIRECTORY CONTAINING STRUCTURES FROM STEP 1] --outdir [OUTPUT DIRECTORY] --errfile [FILENAME TO WRITE ERRORS]`

6. Visualize

    Optionally, visualize the output pockets using the functions contained within `evaluation/get_pymol_pockets.py`.

Note that this whole process can be done with as few as one structure or as many as a whole other proteome. In the large-scale case, parallelization is imperative to get results in this lifetime. If you need help parallelizing these steps, please reach out and we can provide additional scripts.

#### Evaluating pocket predictions
The scripts in the `evaluation` directory contain all functions needed to repeat the analyses conducted in the accompanying manuscript. Specifically:
- `dcccriterion.py` is for calculating the DCCcriterion as shown in Table 4.
- `get_pymol_pockets.py` is for generating Pymol selection commands that can be used to create visualizations similar to Figure 2, Figure 3, Figure 5, Figure 9, and Figure 10.
- `ligand_distance_histograms.py` is for generating histograms of the distance between each pocket center of mass and the nearest atom of a relevant ligand as shown in Figure 6.
- `naive_roc.py` is for generating the ROC curve as shown in Figure 4. 
- `roc_prc.py` is for generating the ROC and PRC curves as shown in Figure 7 and Figure 8,

#### Tuning the *hotpocketNN* model
If you would like to explore the *hotpocketNN* model beyond loading the state dictionary and running it in `assembly/run_hotpocketnn.py`, to try out additional hyperparameters, or to replicate the model tests and hyperparameter tuning described in the manuscript, all required functions are contained in the scripts within the `model` directory. `cnn_sweep.py` and `nn_sweep.py` contain everything needed to run the CNN and NN models, respectively, and are set up for tracking with [W&B](https://wandb.ai). `logistic_regression.py` contains everything needed to run all variants of the logistic regression model.

### Questions?
Please raise an Issue in this repository if anything seems broken or if there is additional functionality you would like to see.
