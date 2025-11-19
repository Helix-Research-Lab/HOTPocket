#!/bin/bash

mkdir data
cd data
git clone https://github.com/Mickdub/gvp.git
cd gvp
git checkout pocket_pred
conda env create -f pocketminer.yml
cp ../../pocketminer.py src
cd ../..

if [ -f "all_pocket_preds.tar.gz" ]; then
    tar -xvf all_pocket_preds.tar.gz
    mv all_pocket_preds/* data
    rm -rf all_pocket_preds
fi
if [ -f "hotpocket_csvs.tar.gz" ]; then
    tar -xvf hotpocket_csvs.tar.gz
    mv hotpocket_csvs/* data
    rm -rf hotpocket_csvs
fi
if [ -f "model_checkpoints.tar.gz" ]; then
    tar -xvf model_checkpoints.tar.gz
    mv model_checkpoints data
fi
if [ -f "model_inputs.tar.gz" ]; then
    tar -xvf model_inputs.tar.gz
    mv model_inputs data
fi
if [ -f "proteome_structure_files_canonical.csv.tar.gz" ]; then
    tar -xvf proteome_structure_files_canonical.csv.tar.gz
    mv proteome_structure_files_canonical.csv data
fi
if [ -f "uniprot.tar.gz" ]; then
    tar -xvf uniprot.tar.gz
    mv uniprot data
fi
if [ -f "BioLiP.txt" ]; then
    mv BioLiP.txt data
fi

echo "setup complete!"
