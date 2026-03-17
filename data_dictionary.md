### Data dictionary
This document contains descriptions of columns found in the output files of *hotpocketNN* -- namely, `proteome_hotpocket_embs_NN_score.csv` and `proteome_hotpocket_comb_NN_score.csv`. The columns in `proteome_structure_files_canonical.csv` and all `pred_*_pockets.csv` dataframes that share these names also share these descriptions.

- *`uniprot id`* - the UniProt ID for the protein in question (can have multiple Uniprot IDs per structure) (string, *e.g.* `Q6N022`)
- *`structure type`* - whether this structure is an experimentally-determined structure from the Protein Data Bank or a computationally-predicted structure from the AlphaFold2 Protein Structure Database; must either be `PDB` for the former case or `AF2` for the latter case (string, `AF2` or `PDB`)
- *`structure id`* - the identifier for the structure in question, for PDB structures this is the PDB ID and for AF2 structures this is the UniProt ID (string, *e.g.* `1ubq` or `Q6N022`). This does not need to match the PDB filename.
- *`structure file`* - the full path to the PDB structure file (string, *e.g.* `/path/to/structure.pdb`). The filename can be different from the `structure id`.
- *`pocket id`* - the identifier for the predicted pocket on this structure (string, *e.g.* `29` or `pocket29`)
- *`pocket res`* - a string representation of all residues in this predicted pocket; residues are delimited by spaces and take the form [CHAIN]-[ONE LETTER AMINO ACID CODE][AMINO ACID INDEX] (string, *e.g.* `A-M101 A-I103 A-H104 A-A111`)
- *`pocket length`* - the number of residues in this predicted pocket (int, *e.g.* `4`)
- *`contributing method name`* - the pocket-finding method that originally predicted this pocket, will be one of the 7 constituent pocket finding methods (string, one of `autosite`, `castp`, `cavity`, `fpocket`, `ligsite`, `pocketminer`, `prank`)
- *`score`* - if the contributing method gives a score, then this column will provide that score, otherwise it will be `NaN` (the contributing methods for which we track scores are `fpocket`, `pocketminer`, and `prank`) (float, *e.g.* `0.23` or `NaN`)
- *`embs NN score`* - the score from the *hotpocketNN* model using Feature Set B (ESM2 embeddings only) (float ranging between `0` and `1`)
- *`comb NN score`* - the score from the *hotpocketNN* model using Feature Set C (ESM2 embeddings combined with constituent method predictions) (float ranging between `0` and `1`)
- *`avg confidence`* - if this pocket is from an AF2-predicted structure, then this column contains the average pLDDT score over all residues in the pocket; if this pocket is from a PDB structure, then this column always contains the value `100.0` (float ranging between `0` and `100`)