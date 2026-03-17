import os
import sys
import re
import time
import argparse
import subprocess
import glob
import pickle
import pandas as pd
import numpy as np
from localpdb import PDB
from tqdm import tqdm
from Bio.PDB import PDBParser
from io import StringIO

sys.path.insert(0, "..")
from utils import extract_resids, get_pairwise_distances, get_neighbors, AA_3TO1
from clean_n_uniprot import process_chain_df


# loads the uniprot canonical fasta file and extracts all uniprot accession ids
#
# params:
# uniprot_fasta_path (str) - path to uniprot canonical fasta file
#
# returns
# accs (list of str) - list of uniprot accessions within the fasta file
def get_canonical_accs(uniprot_fasta_path):
    accs = []
    for line in open(uniprot_fasta_path, "r").read().split("\n"):
        if len(line) < 4:
            continue
        if line[0] == ">" and line[3] == "|":
            acc = line.split("|")[1]
            accs.append(acc)
    return accs


# if the structure filepath exists, adds structure information to growing data frame
#
# params:
# prot (str) - uniprot id
# filepath (str) - path to structure file
# structure_type (str) - PDB or AF2
# structure_id (str) - PDB or AF2 id
# dfdict (dict) - current dataframe in dictionary form
#
# returns:
# dfdict (dict) - updated dataframe in dictionary form
def add_structure(prot, filepath, structure_type, structure_id, dfdict):
    try:
        if os.path.isfile(filepath):
            dfdict["uniprot id"].append(prot)
            dfdict["structure type"].append(structure_type)
            dfdict["structure id"].append(structure_id)
            dfdict["structure file"].append(filepath)
        else:
            print("WARNING: Missing %s structure for %s" % (structure_type, prot))
    except Exception as e:
        print("An unexpected exception occurred: " + str(e))
    return dfdict


# downloads all castp data to folder
#
# params:
# df (pandas DataFrame) - dataframe of proteome structure information
#
# returns:
# None
def castp_download(df):
    df = df[df["structure type"] == "PDB"]
    if not os.path.isdir("../data"):
        os.mkdir("../data")

    if not os.path.isdir("../data/castp"):
        os.mkdir("../data/castp")

    for pdb in df["structure id"].unique():
        if os.path.isdir(os.path.join("../data/castp", pdb)):
            continue
        url = "sts.bioe.uic.edu/castp/data/pdb/%s/%s/processed/%s.zip" % (
            pdb[1:3],
            pdb,
            pdb,
        )
        wget_parts = ["wget", url, "-q"]
        result = subprocess.run(wget_parts)
        time.sleep(0.2)

        if result.returncode != 0:
            continue
        os.mkdir(os.path.join("../data/castp", pdb))
        unzip_parts = [
            "unzip",
            "%s.zip" % pdb,
            "-d",
            os.path.join("../data/castp", pdb),
        ]
        subprocess.run(unzip_parts)
        subprocess.run(["rm", "%s.zip" % pdb])

# batch downloads all PDB and AF2 CavitySpace pockets
def setup_cavityspace():
    os.system("python ../data/cavityspace/cavityspace_downloader.py -i ../data/cavityspace/AF-all_cavities_index.pkl -k cavity -o ../data/cavityspace/AF-all-cavities")
    os.system("python ../data/cavityspace/cavityspace_downloader.py -i ../data/cavityspace/hrefPDB-all_cavities_index.pkl -k cavity -o ../data/cavityspace/PDB-all-cavities")


# removes any hetatm records from the pdb file for safe pocketfinding
#
# params:
# pdbpath (str) - path to the pdb file
#
# returns:
# None
def strip_hetatm(pdbpath):
    pattern = re.compile(r"HETATM")
    lines = open(pdbpath, "r").read().split("\n")
    newpdb = open(pdbpath, "w")
    for l in lines:
        if not pattern.search(l):
            newpdb.write(l)
            newpdb.write("\n")
    newpdb.close()


# moves structure to specified directory and prepares it for pocket finding
#
# params:
# fpath (str) - original structure filepath
# dirname (str) - folder to copy structure to
#
# returns:
# new_fpath (str) - moved structure filepath
def move_process_struct(fpath, dirname):
    # copy structure file to folder
    cp_parts = ["cp", fpath, dirname]
    subprocess.run(cp_parts)

    # localpdb stores pdbs as .ent.gz, need to unzip to run
    fname = os.path.split(fpath)[1]
    new_fpath = os.path.join(dirname, fname)
    if fname[-3:] == ".gz":
        unzip_parts = ["gunzip", new_fpath]
        subprocess.run(unzip_parts)
        fname = fname[:-3]
        new_fpath = os.path.join(dirname, fname)

    # fpocket only runs on .pdb, not .ent, even tho it's the same file type
    if fname[-4:] == ".ent":
        tmp_fpath = new_fpath[:-4] + ".pdb"
        rename_parts = ["mv", new_fpath, tmp_fpath]
        subprocess.run(rename_parts)
        new_fpath = tmp_fpath
        fname = fname[:-4] + ".pdb"

    # make sure pocketfinding isn't affected by ligands in pdb files!
    strip_hetatm(new_fpath)

    return new_fpath


# for fpocket, ligsite, autosite, and pocketminer -- copies each structure file into a folder in the
# data directory named with the method name, runs the method on the structures
# (this is the step before assembling the standardized dataframe)
#
# params:
# df (pandas DataFrame) - dataframe with all proteome structure files
# method (str) - method to run, either fpocket, ligsite, autosite, or pocketminer
# start (int) - if running program in multiple parallel batches, the row index of df at which to start this batch
# n (int) - if running program in multiple parallel batches, the number of rows in df to include in this batch; if this is set to any negative number, will process the whole df
# autosite_script_relative_path (str) - where the run_autosite.sh script is located relative to the current working directory, only necessary for running autosite
#
# returns:
# None
def run_pocket_program(df, method, start=0, n=-1, autosite_script_relative_path="."):
    if not os.path.isdir("../data"):
        os.mkdir("../data")

    assert method in ["fpocket", "ligsite", "pocketminer", "autosite"]
    dirname = os.path.join("../data", method)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    if method == "pocketminer":
        fpaths = []
        outnames = []
        protdirs = []

    if n > -1:
        stop = min(start + n, len(df))
        df = df[start:stop]

    for fpath in tqdm(df["structure file"]):
        fname = os.path.split(fpath)[1]
        # Use os.path.splitext to properly handle filenames with multiple dots
        if "%s_out" % os.path.splitext(fname)[0] in os.listdir(dirname):
            continue

        tmp_fpath = move_process_struct(fpath, dirname)
        fname = os.path.split(tmp_fpath)[1]

        # run method on structure file
        if method == "fpocket":
            method_parts = ["fpocket", "-f", tmp_fpath]
            subprocess.run(method_parts)

            # remove copied structure file to save space -- we don't need it anymore
            rm_parts = ["rm", tmp_fpath]
            subprocess.run(rm_parts)
        elif method == "ligsite":
            # ligsite doesn't automatically create a folder for outputs, so we
            # will make one for organizational purposes
            protdir = os.path.join(dirname, "%s_out" % os.path.splitext(fname)[0])
            if os.path.isdir(
                protdir
            ):  # check again, sometimes already done with parallel processes
                continue
            os.mkdir(protdir)
            mv_parts = ["mv", tmp_fpath, protdir]
            subprocess.run(mv_parts)
            wd = os.getcwd()
            os.chdir(protdir)
            method_parts = ["lcs", "-i", fname]
            subprocess.run(method_parts)
            os.chdir(wd)
        elif method == "pocketminer":
            fname_base = os.path.splitext(fname)[0]
            protdir = os.path.join(dirname, "%s_out" % fname_base)
            if not os.path.isdir(protdir):
                os.mkdir(protdir)
            outname = fname_base
            fpaths.append(tmp_fpath)
            outnames.append(outname)
            protdirs.append(protdir)
        elif method == "autosite":
            fname_base = os.path.splitext(fname)[0]
            wd = os.getcwd()
            os.chdir(dirname)
            method_parts = [
                "bash",
                os.path.join(autosite_script_relative_path, "run_autosite.sh"),
                fname,
                "%s_out" % fname_base,
            ]
            subprocess.run(method_parts)

            # remove copied structure file to save space -- we don't need it anymore
            rm_parts = [
                "rm",
                fname,
                "%s_out_AutoSiteSummary.log" % fname_base,
                "%s_out.pdbqt" % fname_base,
                "%s_out_truncated.pdb" % fname_base,
            ]
            subprocess.run(rm_parts)
            os.chdir(wd)

    if method == "pocketminer":
        sys.path.append("../data/gvp/src")
        import pocketminer

        pocketminer.main(fpaths, outnames, protdirs)
        for idx, fpath in enumerate(fpaths):
            mv_parts = ["mv", fpath, protdirs[idx]]
            subprocess.run(mv_parts)

# for either PDB or AF2 structures in df, creates the input file and runs the p2rank program
# note that you may need to load java for this -- on stanford sherlock, the command is `ml java/21`
#
# params:
# df (pandas DataFrame) - dataframe with all proteome structure files
# struc_type (str) - either "PDB" or "AF2"
#
# returns:
# None
def run_prank(df, struc_type):
    if not struc_type in df["structure type"].tolist():
        return

    subdf = df[df["structure type"] == struc_type].copy()
    ds_file_name = "../data/prank/%s_all.ds" % struc_type.lower()
    outdir = "../data/prank/%s_outs" % struc_type.lower()
    subdf[["structure file"]].to_csv(ds_file_name, index=False, header=False)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    af2_flag = ""
    if struc_type == "AF2":
        af2_flag = " -c alphafold"
    
    os.system("prank predict -o %s -threads 8%s %s" % (outdir, af2_flag, ds_file_name))


# creates folder for prank, creates input file, and runs the p2rank program
#
# params:
# df (pandas DataFrame) - dataframe with all proteome structure files
#
# returns:
# None
def make_prank_dir(df):
    if not os.path.isdir("../data"):
        os.mkdir("../data")
    if not os.path.isdir("../data/prank"):
        os.mkdir("../data/prank")

    run_prank(df, "AF2")
    run_prank(df, "PDB")


# takes pdb file with ONLY pocket residues and converts to list of res ids
# every res in pdb file will be returned in res id list
#
# params:
# row (Pandas DataFrame row) - row of dataframe with pocket info
# pocket_lines (dict) - OPTIONAL, dict containing pdb lines for each pocket, forgoes need to have individual pocket pdb files
#
# returns:
# resids (str) - list of pocket res ids delimited by spaces
def get_pocket_res(row, parser, pocket_lines=None):
    if pocket_lines == None:
        structure = parser.get_structure("struc", row["pocket path"])
    else:
        lines = pocket_lines[row["structure id"]][row["pocket id"]]
        pdb_txt = "\n".join(lines)
        pdb_io = StringIO(pdb_txt)
        structure = parser.get_structure("struc", pdb_io)

    resid_list, _ = extract_resids("", structure=structure)
    resids = " ".join(resid_list)
    return resids


# given a central coordinate for a pocket, returns the resids for all residues within r angstroms away from that coordinate
#
# params:
# pdb_path (str) - path to pdb file from which to get pocket lines
# cen (str) - pdb line with coordinates of pocket center
# r (int) - radius of pocket, in angstroms
#
# returns:
# resids (str) - standard form of all resids within radius r of cen
def get_lines_within_radius(pdb_path, cen, r=5):
    parser = PDBParser(QUIET=True)
    pdb_io = StringIO(cen)
    structure = parser.get_structure("center", pdb_io)
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    cen_coords = atom.get_coord()

    resids_list = extract_resids(pdb_path, cen=cen_coords, r=r)[0]
    resids = " ".join(resids_list)
    return resids


# pocketminer clustering as described in pocketminer paper
# finds 10 closest residues to central residue, averages scores, reports as pocket if over threshold
#
# params:
# pocket_scores (np array) - scores for all individual residues
# pairwise_distances (np array) - all pairwise c-alpha distances
# thresh (float) - threshold for hotspot, must be between 0 and 1
# n_neighbors (int) - number of neighboring residues to include in hotspot calculation
#
# returns:
# hotspots (list) - list of hotspot residue indices
# hotspot_scores (list) - list of scores corresponding to each element of the hotspots list
def get_hotspot_scores(pocket_scores, pairwise_distances, thresh, n_neighbors=10):
    hotspots = []
    hotspot_scores = []
    for idx in range(len(pocket_scores)):
        neighbor_idx = get_neighbors(idx, pairwise_distances, n_neighbors)
        neighborhood_scores = [pocket_scores[n] for n in neighbor_idx]
        neighborhood_avg = np.mean(neighborhood_scores)
        if neighborhood_avg > thresh:
            if (
                not set(neighbor_idx) in hotspots
                and not neighborhood_avg in hotspot_scores
            ):
                hotspots.append(set(neighbor_idx))
                hotspot_scores.append(neighborhood_avg)
    return hotspots, hotspot_scores


# obtains fpocket-generated druggabilty scores for fpocket-predicted pockets
#
# params:
# info_fname (str) - path to fpocket info output for structure (*_info.txt)
#
# returns:
# scores (dict) - dictionary with keys as the pocket ids, and values as the druggability scores of each pocket id
def get_fpocket_scores(info_fname):
    scores = dict()
    f = open(info_fname, "r")
    current_pocket = -1
    for line in f.read().split("\n"):
        if len(line) == 0:
            continue
        spl = line.split("\t")
        if len(spl) == 1:
            spl = line.split()
            if spl[0] == "Pocket" and spl[2] == ":":
                current_pocket = int(spl[1])
        elif len(spl) == 3 and spl[0] == "":
            if spl[1] == "Druggability Score : ":
                if current_pocket == -1:
                    print("ERROR: not inside pocket but found score!")
                    print(line)
                elif current_pocket in scores.keys():
                    print(
                        "ERROR: multiple druggability scores detected for single pocket"
                    )
                    print("Pocket %d" % current_pocket)
                else:
                    scores[current_pocket] = float(spl[2])
        else:
            print("unexpected line:")
            print(line)
    f.close()
    return scores


# to be used when processing prank outputs
# prank uses chain and residue number in outputs
# this function will standardize to chain, res name, and res number
#
# params:
# row (pandas row) - row with pocket information (used in apply)
# res_col (str) - name of column with original pocket residue information
# all_resids (list) - list of all resids, extracted from pdb file
#
# returns:
# resid (str) - standardized format of pocket res list
def prank_to_resid(row, res_col, all_resids):
    original_pocket_resids = row[res_col].strip().split(" ")
    new_pocket_resids = []
    for res in all_resids:
        chain = res.split("-")[0]
        resid = "-".join(res.split("-")[1:])
        resname = resid[0]
        resnum = int(resid[1:])
        prank_rep = "%s_%d" % (chain, resnum)
        if prank_rep in original_pocket_resids:
            new_pocket_resids.append(res)
    return " ".join(new_pocket_resids)


# cleans empty pockets and pockets assigned to wrong uniprot id
#
# params:
# df (pandas DataFrame) - dataframe of predicted pockets to clean
# memo_filepath (str) - path to memo pickle
#
# returns:
# df (pandas DataFrame) - cleaned dataframe of predicted pockets
def clean_df(df, memo_filepath="../data/pdb_chain_map_memo.p"):
    if os.path.isfile(memo_filepath):
        memo = pickle.load(open(memo_filepath, "rb"))
    else:
        memo = {}
    df = df.dropna(subset=["pocket res"])
    df, memo = process_chain_df(df, memo=memo)
    pickle.dump(memo, open(memo_filepath, "wb"))
    return df


# creates df with paths to pdb and af2 structure files
#
# params:
# af2db_path (str) - path to af2 db
# localpdb_path (str) - path to localpdb db
# uniprot_path (str) - path to uniprot directory with .fasta and .idmapping files
#
# returns
# None
def make_files_df(af2db_path, localpdb_path, uniprot_path="../data/uniprot"):
    if not os.path.isdir("../data"):
        os.mkdir("../data")

    # load in all of uniprot human proteome
    uniprot_fasta_path = os.path.join(uniprot_path, "UP000005640_9606.fasta")
    uniprot_id_path = os.path.join(uniprot_path, "UP000005640_9606.idmapping")
    uniprot_canonical = get_canonical_accs(uniprot_fasta_path)
    uniprot_ids = pd.read_csv(
        uniprot_id_path, sep="\t", names=["uniprot id", "id type", "external id"]
    )

    # set up localpdb
    lpdb = PDB(db_path=localpdb_path)

    # setting up new df
    dfdict = {
        "uniprot id": [],
        "structure type": [],
        "structure id": [],
        "structure file": [],
    }

    # for each protein
    for prot in tqdm(uniprot_canonical):
        # find the AF2 structure -- this is the uniprot id in the leftmost column
        af2_fname = "AF-%s-F1-model_v1.pdb" % prot
        af2_filepath = os.path.join(af2db_path, af2_fname)
        dfdict = add_structure(prot, af2_filepath, "AF2", prot, dfdict)

        # find any pdb ids -- second column will say PDB
        pdb_ids = uniprot_ids.loc[uniprot_ids["uniprot id"] == prot]
        pdb_ids = pdb_ids.loc[pdb_ids["id type"] == "PDB"]
        pdb_ids = pdb_ids["external id"]

        if len(pdb_ids) == 0:
            continue

        for pdb in pdb_ids:
            try:
                pdb_filepath = lpdb.entries.loc[pdb.lower()]["pdb_fn"]
            except:
                print("WARNING: no PDB structure found for %s" % pdb)
            else:
                dfdict = add_structure(prot, pdb_filepath, "PDB", pdb.lower(), dfdict)
                (pocket_scores[idx])
    df = pd.DataFrame(dfdict)
    df.to_csv("../data/proteome_structure_files_canonical.csv", sep=",")


# creates df of known binding pockets from biolip
#
# params:
# biolip_path (str) - path to biolip dataframe
#
# returns:
# None
def make_known_df(biolip_path="../data/BioLiP.txt"):
    if not os.path.isdir("../data"):
        os.mkdir("../data")

    # load biolip df
    cols = [
        "pdb",
        "chain",
        "resolution",
        "binding site number",
        "ligand id",
        "ligand chain",
        "ligand serial number",
        "binding site residues, pdb numbering",
        "binding site residues, reindexed numbering",
        "cat res, pdb numbering",
        "cat res, reindexed numbering",
        "EC number",
        "GO terms",
        "binding affinity, literature",
        "binding affinity, binding moad",
        "binding affinity, pdbbind",
        "binding affinity, bindingdb",
        "uniprot",
        "pubmed",
        "ligand res number",
        "sequence",
    ]
    biolip_df = pd.read_csv(biolip_path, sep="\t", names=cols)

    # condense df into just the columns we want
    df = biolip_df[
        [
            "uniprot",
            "pdb",
            "chain",
            "binding site residues, pdb numbering",
            "binding site residues, reindexed numbering",
            "ligand id",
        ]
    ]
    df.to_csv("../data/known_pockets.csv", sep=",", index=False)


# creates dataframe with autosite predictions
#
# params:
# df (pandas DataFrame) - dataframe of proteome structure information
# makedir (bool) - flag for whether or not to construct autosite data folder
# start (int) - if creating dataframe in multiple parallel batches, the row index of df at which to start this batch
# n (int) - if creating dataframe in multiple parallel batches, the number of rows in df to include in this batch; if this is set to any negative number, will process the whole df
# autosite_script_relative_path (str) - where the run_autosite.sh script is located relative to the current working directory
#
# returns:
# None
def make_autosite_df(
    df, makedir=True, start=0, n=-1, autosite_script_relative_path="."
):
    if not os.path.isdir("../data"):
        os.mkdir("../data")

    if makedir:
        run_pocket_program(
            df, "autosite", autosite_script_relative_path=autosite_script_relative_path
        )

    # lists for autosite df
    uniprot_id = []
    structure_type = []
    structure_id = []
    pocket_id = []
    res_lists = []

    if n > -1:
        stop = min(start + n, len(df))
        df = df[start:stop]

    uniprot_ids = df["uniprot id"].tolist()
    structure_types = df["structure type"].tolist()
    structure_ids = df["structure id"].tolist()
    for row_idx, fpath in enumerate(df["structure file"]):
        fname = os.path.split(fpath)[1]
        fname_base = os.path.splitext(fname)[0]
        out_dir = os.path.join("../data", "autosite", "%s_out" % fname_base)
        if not os.path.isdir(out_dir):
            continue
        pocket_npys = [f for f in os.listdir(out_dir) if f[-4:] == ".npy"]
        pocket_npys = [f for f in pocket_npys if "_fp_" in f]
        structure = None
        for pocket_npy in pocket_npys:
            current_res_list = set()
            pid = int(pocket_npy.split("_")[-1][:-4])
            fps = np.load(os.path.join(out_dir, pocket_npy))
            for coord in fps:
                resids, structure = extract_resids(
                    fpath, cen=coord, r=5, structure=structure
                )
                current_res_list.update(resids)
            uniprot_id.append(uniprot_ids[row_idx])
            structure_type.append(structure_types[row_idx])
            structure_id.append(structure_ids[row_idx])
            pocket_id.append(pid)
            res_lists.append(" ".join(list(current_res_list)))

    autosite_df = pd.DataFrame(
        {
            "uniprot id": uniprot_id,
            "structure type": structure_type,
            "structure id": structure_id,
            "pocket id": pocket_id,
            "pocket res": res_lists,
        }
    )

    autosite_df = clean_df(autosite_df)

    if n == -1:
        autosite_df.to_csv("../data/pred_autosite_pockets.csv", sep=",", index=False)
    else:
        autosite_df.to_csv(
            "../data/pred_autosite_pockets_%d_%d.csv" % (start, stop),
            sep=",",
            index=False,
        )


# creates dataframe with castp predictions
#
# params:
# df (pandas DataFrame) - dataframe of proteome structure information
# makedir (bool) - flag for whether or not to construct castp data folder
#
# returns:
# None
def make_castp_df(df, makedir=True):
    if not os.path.isdir("../data"):
        os.mkdir("../data")

    if makedir:
        castp_download(df)

    # lists for castp df
    uniprot_id = []
    structure_type = []
    structure_id = []
    pocket_id = []

    pocket_lines = dict()
    uniprot_ids = df["uniprot id"].tolist()
    structure_types = df["structure type"].tolist()
    for row_idx, struc_id in enumerate(df["structure id"]):
        pockets_path = os.path.join(
            "../data", "castp", struc_id, "%s.poc" % struc_id
        )
        info_path = "%sInfo" % pockets_path
        if not os.path.isfile(pockets_path) or not os.path.isfile(info_path):
            continue
        info_df = pd.read_csv(info_path, sep="\t")
        pockets_pdb = open(pockets_path, "r").read().split("\n")
        pocket_lines[struc_id] = dict()
        for p in info_df["ID"]:
            pattern = re.compile(r"\s\s*%s\s*POC" % p)
            pocket_lines[struc_id][p] = [
                line for line in pockets_pdb if pattern.search(line)
            ]
            uniprot_id.append(uniprot_ids[row_idx])
            structure_type.append(structure_types[row_idx])
            structure_id.append(struc_id)
            pocket_id.append(p)

    castp_df = pd.DataFrame(
        {
            "uniprot id": uniprot_id,
            "structure type": structure_type,
            "structure id": structure_id,
            "pocket id": pocket_id,
        }
    )

    # get pocket residues out of castp
    if len(castp_df) > 0:
        parser = PDBParser(QUIET=True)
        castp_df["pocket res"] = castp_df.apply(
            get_pocket_res, args=(parser, pocket_lines), axis=1
        )

    castp_df = clean_df(castp_df)

    castp_df.to_csv("../data/pred_castp_pockets.csv", sep=",", index=False)


# creates dataframe with cavityspace predictions
#
# params:
# df (pandas DataFrame) - dataframe of proteome structure information
# makedir (bool) - flag for whether or not to construct castp data folder
# cavity_download_parent_dir (str) - path to data folder that contains cavityspace download
#
# returns:
# None
def make_cavityspace_df(df, makedir=False, cavity_download_parent_dir=".."):
    if not os.path.isdir("../data"):
        os.mkdir("../data")

    if makedir:
        setup_cavityspace()
        
    # lists for cavityspace df
    uniprot_id = []
    structure_type = []
    structure_id = []
    pocket_id = []
    pocket_path = []

    dir_contents = dict()
    cavityspace_dir = os.path.join(cavity_download_parent_dir, "data", "cavityspace")
    dirname = dict()
    dirname["AF2"] = os.path.join(cavityspace_dir, "AF-all-cavities")
    dirname["PDB"] = os.path.join(cavityspace_dir, "PDB-all-cavities")
    dir_contents["AF2"] = os.listdir(dirname["AF2"])
    dir_contents["PDB"] = os.listdir(dirname["PDB"])

    uniprot_ids = df["uniprot id"].tolist()
    structure_types = df["structure type"].tolist()
    for row_idx, struc_id in enumerate(df["structure id"]):
        struc_type = structure_types[row_idx]
        if struc_type == "AF2":
            filename = "AF-%s-F1-model_v1" % struc_id
        elif struc_type == "PDB":
            filename = "%s_%s_" % (uniprot_ids[row_idx], struc_id.upper())
        else:
            print("ERROR: unexpected structure type")
            return

        pocket_results = [f for f in dir_contents[struc_type] if filename in f]
        for pocket_tar in pocket_results:
            if pocket_tar[-7:] == ".tar.gz" and "_cavity_result_" in pocket_tar:
                unzip_parts = [
                    "tar",
                    "-xf",
                    os.path.join(dirname[struc_type], pocket_tar),
                    "-C",
                    dirname[struc_type],
                ]
                subprocess.run(unzip_parts)
            else:
                continue
            p = pocket_tar[:-7].split("_")[-1]
            assert p.isnumeric()
            pocket_pdb_name = os.path.join(
                dirname[struc_type],
                "%s_%s.pdb" % ("_".join(pocket_tar.split("_")[:-2]), p),
            )
            if not os.path.isfile(pocket_pdb_name):
                continue
            uniprot_id.append(uniprot_ids[row_idx])
            structure_type.append(struc_type)
            structure_id.append(struc_id)
            pocket_id.append(p)
            pocket_path.append(pocket_pdb_name)
        if len(pocket_results) > 0:
            surface_files = glob.glob(
                os.path.join(dirname[struc_type], "%s*_surface_*.pdb" % filename)
            )
            box_files = glob.glob(
                os.path.join(dirname[struc_type], "%s*_box_*.txt" % filename)
            )
            vacant_files = glob.glob(
                os.path.join(dirname[struc_type], "%s*_vacant_*.pdb" % filename)
            )
            subprocess.run(["rm"] + surface_files + box_files + vacant_files)

    cavityspace_df = pd.DataFrame(
        {
            "uniprot id": uniprot_id,
            "structure type": structure_type,
            "structure id": structure_id,
            "pocket id": pocket_id,
            "pocket path": pocket_path,
        }
    )
    if len(cavityspace_df) == 0:
        cavityspace_df["pocket res"] = None
    else:
        # get pocket residues out of castp
        cavityspace_df["pocket res"] = cavityspace_df.apply(
            get_pocket_res, args=(parser,), axis=1
        )
    cavityspace_df = cavityspace_df.drop(columns=["pocket path"])

    cavityspace_df = clean_df(cavityspace_df)

    cavityspace_df.to_csv("../data/pred_cavityspace_pockets.csv", sep=",", index=False)


# creates dataframe with fpocket predictions
#
# params:
# df (pandas DataFrame) - dataframe of proteome structure information
# makedir (bool) - flag for whether or not to construct fpocket data folder
# start (int) - if creating dataframe in multiple parallel batches, the row index of df at which to start this batch
# n (int) - if creating dataframe in multiple parallel batches, the number of rows in df to include in this batch; if this is set to any negative number, will process the whole df
#
# returns:
# None
def make_fpocket_df(df, makedir=True, start=0, n=-1):
    if not os.path.isdir("../data"):
        os.mkdir("../data")

    if makedir:
        run_pocket_program(df, "fpocket", start=start, n=n)

    if n > -1:
        stop = min(start + n, len(df))
        df = df[start:stop]

    # lists for fpocket df
    uniprot_id = []
    structure_type = []
    structure_id = []
    pocket_id = []
    pocket_path = []
    pocket_score = []

    uniprot_ids = df["uniprot id"].tolist()
    structure_types = df["structure type"].tolist()
    structure_ids = df["structure id"].tolist()
    for row_idx, fpath in enumerate(tqdm(df["structure file"])):
        fname = os.path.split(fpath)[1]
        # Use os.path.splitext to properly handle filenames with multiple dots
        fname_base = os.path.splitext(fname)[0]
        prot_dir = os.path.join("../data/fpocket", "%s_out" % fname_base)
        pockets_dir = os.path.join(prot_dir, "pockets")
        if not os.path.isdir(pockets_dir):
            continue
        pockets = os.listdir(pockets_dir)
        pockets = [f for f in pockets if f[-4:] == ".pdb"]

        # Use filename (not structure_id) since fpocket outputs files based on input filename
        info_fname = os.path.join(prot_dir, "%s_info.txt" % fname_base)
        scores_dict = get_fpocket_scores(info_fname)
        for p in pockets:
            pid = p.split("_")[0]
            uniprot_id.append(uniprot_ids[row_idx])
            structure_type.append(structure_types[row_idx])
            structure_id.append(structure_ids[row_idx])
            pocket_id.append(pid)
            pocket_path.append(os.path.join(pockets_dir, p))
            pocket_score.append(scores_dict[int(pid[6:])])

    fpocket_df = pd.DataFrame(
        {
            "uniprot id": uniprot_id,
            "structure type": structure_type,
            "structure id": structure_id,
            "pocket id": pocket_id,
            "pocket path": pocket_path,
            "score": pocket_score,
        }
    )

    # get pocket residues out of fpocket
    parser = PDBParser(QUIET=True)
    fpocket_df["pocket res"] = fpocket_df.apply(get_pocket_res, args=(parser,), axis=1)
    fpocket_df = fpocket_df.drop(columns=["pocket path"])

    fpocket_df = clean_df(fpocket_df)

    if n == -1:
        fpocket_df.to_csv(
            "../data/pred_fpocket_pockets_scored.csv", sep=",", index=False
        )
    else:
        fpocket_df.to_csv(
            "../data/pred_fpocket_pockets_scored_%d_%d.csv" % (start, stop),
            sep=",",
            index=False,
        )


# creates dataframe with ligsite predictions
#
# params:
# df (pandas DataFrame) - dataframe of proteome structure information
# makedir (bool) - flag for whether or not to construct ligsite data folder
#
# returns:
# None
def make_ligsite_df(df, makedir=True):
    ligsite_radius = 8

    if not os.path.isdir("../data"):
        os.mkdir("../data")

    if makedir:
        run_pocket_program(df, "ligsite")

    # lists for ligsite df
    uniprot_id = []
    structure_type = []
    structure_id = []
    pocket_id = []
    pocket_res = []

    uniprot_ids = df["uniprot id"].tolist()
    structure_types = df["structure type"].tolist()
    structure_ids = df["structure id"].tolist()
    for row_idx, fpath in enumerate(df["structure file"]):
        fname = os.path.split(fpath)[1]
        fname_base = os.path.splitext(fname)[0]
        out_dir = os.path.join("../data", "ligsite", "%s_out" % fname_base)
        pockets_path = os.path.join(out_dir, "pocket.pdb")
        pdb_path = os.path.join(out_dir, "%s.pdb" % fname_base)
        if (
            not os.path.isdir(out_dir)
            or not os.path.isfile(pockets_path)
            or os.stat(pockets_path).st_size == 0
        ):
            continue
        pocket_centers = open(pockets_path, "r").read().split("\n")
        for idx, cen in enumerate(pocket_centers):
            if len(cen) < 54:
                continue
            resids = get_lines_within_radius(pdb_path, cen, r=ligsite_radius)
            uniprot_id.append(uniprot_ids[row_idx])
            structure_type.append(structure_types[row_idx])
            structure_id.append(structure_ids[row_idx])
            pocket_id.append(idx)
            pocket_res.append(resids)

    ligsite_df = pd.DataFrame(
        {
            "uniprot id": uniprot_id,
            "structure type": structure_type,
            "structure id": structure_id,
            "pocket id": pocket_id,
            "pocket res": pocket_res,
        }
    )

    ligsite_df = clean_df(ligsite_df)

    ligsite_df.to_csv("../data/pred_ligsite_pockets.csv", sep=",", index=False)


# creates dataframe with pocketminer predictions
#
# params:
# df (pandas DataFrame) - dataframe of proteome structure information
# makedir (bool) - flag for whether or not to construct pocketminer data folder
# start (int) - if creating dataframe in multiple parallel batches, the row index of df at which to start this batch
# n (int) - if creating dataframe in multiple parallel batches, the number of rows in df to include in this batch; if this is set to any negative number, will process the whole df
#
# returns:
# None
def make_pocketminer_df(df, makedir=True, start=0, n=-1):
    if not os.path.isdir("../data"):
        os.mkdir("../data")

    if makedir:
        run_pocket_program(df, "pocketminer", start=start, n=n)

    # lists for pocket df
    uniprot_id = []
    structure_type = []
    structure_id = []
    pocket_id = []
    pocket_res = []
    scores = []

    if n > -1:
        stop = min(start + n, len(df))
        df = df[start:stop]

    uniprot_ids = df["uniprot id"].tolist()
    structure_types = df["structure type"].tolist()
    structure_ids = df["structure id"].tolist()
    for row_idx, fpath in enumerate(tqdm(df["structure file"])):
        fname = os.path.split(fpath)[1]
        fname_base = os.path.splitext(fname)[0]
        pockets_dir = os.path.join(
            "../data/pocketminer", "%s_out" % fname_base
        )
        if not os.path.isdir(pockets_dir):
            print("dir not found")
            continue
        pocket_scores_file = os.path.join(
            pockets_dir, "%s-preds.npy" % fname_base
        )
        if not os.path.isfile(pocket_scores_file):
            print("scores not found")
            continue
        pocket_scores = np.load(pocket_scores_file)
        pocket_scores = pocket_scores.flatten()
        seq_resids, structure = extract_resids(fpath)
        if len(pocket_scores) != len(seq_resids):
            print("%d != %d" % (len(pocket_scores), len(seq_resids)))
            continue
        pairwise_distances = get_pairwise_distances(structure, seq_resids)
        hotspots, hotspot_scores = get_hotspot_scores(
            pocket_scores, pairwise_distances, thresh=0.7
        )
        for idx, hotspot in enumerate(hotspots):
            resids_list = [seq_resids[i] for i in hotspot]
            resids = " ".join(resids_list)

            uniprot_id.append(uniprot_ids[row_idx])
            structure_type.append(structure_types[row_idx])
            structure_id.append(structure_ids[row_idx])
            pocket_id.append(idx)
            pocket_res.append(resids)
            scores.append(hotspot_scores[idx])

    pocketminer_df = pd.DataFrame(
        {
            "uniprot id": uniprot_id,
            "structure type": structure_type,
            "structure id": structure_id,
            "pocket id": pocket_id,
            "pocket res": pocket_res,
            "score": scores,
        }
    )

    pocketminer_df = clean_df(pocketminer_df)

    if n == -1:
        pocketminer_df.to_csv(
            "../data/pred_pocketminer_pockets_scored.csv", sep=",", index=False
        )
    else:
        pocketminer_df.to_csv(
            "../data/pred_pocketminer_pockets_scored_%d_%d.csv" % (start, stop),
            sep=",",
            index=False,
        )


# creates dataframe with prank predictions
#
# params:
# df (pandas DataFrame) - dataframe of proteome structure information
# makedir (bool) - flag for whether or not to construct prank data folder
# start (int) - if creating dataframe in multiple parallel batches, the row index of df at which to start this batch
# n (int) - if creating dataframe in multiple parallel batches, the number of rows in df to include in this batch; if this is set to any negative number, will process the whole df
#
# returns:
# None
def make_prank_df(df, makedir=False, start=0, n=-1):
    if not os.path.isdir("../data"):
        os.mkdir("../data")

    if makedir:
        make_prank_dir(df)

    uniprot_id = []
    structure_type = []
    structure_id = []
    pocket_id = []

    prank_df = pd.DataFrame(
        {
            "uniprot id": [],
            "structure type": [],
            "structure id": [],
            "pocket id": [],
            "pocket res": [],
            "score": [],
        }
    )

    if n > -1:
        stop = min(start + n, len(df))
        df = df[start:stop]

    uniprot_ids = df["uniprot id"].tolist()
    structure_types = df["structure type"].tolist()
    structure_ids = df["structure id"].tolist()
    for row_idx, fpath in enumerate(df["structure file"]):
        fname = os.path.split(fpath)[1]
        subdir = "%s_outs" % structure_types[row_idx].lower()
        out_fname = os.path.join(
            "../data", "prank", subdir, "%s_predictions.csv" % fname
        )
        if not os.path.isfile(out_fname):
            continue
        pockets_df = pd.read_csv(out_fname)
        if len(pockets_df) == 0:
            continue
        name_col = pockets_df.columns[0]
        score_col = pockets_df.columns[2]
        res_col = pockets_df.columns[-2]
        assert name_col.strip() == "name"
        assert score_col.strip() == "score"
        assert res_col.strip() == "residue_ids"
        pockets_df["uniprot id"] = uniprot_ids[row_idx]
        pockets_df["structure type"] = structure_types[row_idx]
        pockets_df["structure id"] = structure_ids[row_idx]
        pockets_df["pocket id"] = pockets_df.apply(
            lambda row: row[name_col].strip(), axis=1
        )
        pockets_df["score"] = pockets_df[score_col]
        all_resids, _ = extract_resids(fpath)
        pockets_df["pocket res"] = pockets_df.apply(
            prank_to_resid, args=(res_col, all_resids), axis=1
        )
        pockets_df = pockets_df[
            [
                "uniprot id",
                "structure type",
                "structure id",
                "pocket id",
                "pocket res",
                "score",
            ]
        ]
        prank_df = pd.concat([prank_df, pockets_df])

    prank_df = clean_df(prank_df)

    if n == -1:
        prank_df.to_csv("../data/pred_prank_pockets_scored.csv", sep=",", index=False)
    else:
        prank_df.to_csv(
            "../data/pred_prank_pockets_scored_%d_%d.csv" % (start, stop),
            sep=",",
            index=False,
        )


def main(args):
    if args.proteome_df_path:
        df = pd.read_csv(args.proteome_df_path, sep=",")
    if args.filedb:
        make_files_df(args.af2db_path, args.localpdb_path)
    if args.knowndb:
        make_known_df()
    if args.castp:
        make_castp_df(df, makedir=True)
    if args.fpocket:
        make_fpocket_df(df, makedir=True, start=args.start, n=args.n)
    if args.pocketminer:
        make_pocketminer_df(df, makedir=True, start=args.start, n=args.n)
    if args.ligsite:
        make_ligsite_df(df, makedir=True)
    if args.cavityspace:
        make_cavityspace_df(df, makedir=True)
    if args.autosite:
        make_autosite_df(df, makedir=True, start=args.start, n=args.n)
    if args.prank:
        make_prank_df(df, makedir=True, start=args.start, n=args.n)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--proteome_df_path",
        type=str,
        help="location of dataframe with proteome structure files",
    )
    parser.add_argument(
        "--af2db_path",
        type=str,
        help="location of AlphaFold2 human structures directory, only needed for creating proteome structure files dataframe",
    )
    parser.add_argument(
        "--localpdb_path",
        type=str,
        help="location of localpdb directory, only needed for creating proteome structure files dataframe",
    )
    parser.add_argument(
        "--filedb",
        action="store_true",
        help="flag to create proteome_structure_files_canonical.csv",
    )
    parser.add_argument(
        "--knowndb", action="store_true", help="flag to create known_pockets.csv"
    )
    parser.add_argument(
        "--castp", action="store_true", help="flag for assembling castp info"
    )
    parser.add_argument(
        "--cavityspace",
        action="store_true",
        help="flag for assembling cavityspace info",
    )
    parser.add_argument(
        "--fpocket", action="store_true", help="flag for assembling fpocket info"
    )
    parser.add_argument(
        "--pocketminer",
        action="store_true",
        help="flag for assembling pocketminer info",
    )
    parser.add_argument(
        "--ligsite", action="store_true", help="flag for assembling ligsite info"
    )
    parser.add_argument(
        "--autosite", action="store_true", help="flag for assembling autosite info"
    )
    parser.add_argument(
        "--prank", action="store_true", help="flag for assembling prank info"
    )
    parser.add_argument(
        "--start",
        type=int,
        default=0,
        help="for running dataframe assembly in parallel, starting row",
    )
    parser.add_argument(
        "-n",
        type=int,
        default=-1,
        help="for running dataframe assembly in parallel, number of rows in batch",
    )
    args = parser.parse_args()
    main(args)
