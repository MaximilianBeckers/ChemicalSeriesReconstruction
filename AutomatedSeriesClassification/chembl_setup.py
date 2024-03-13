#
#  Copyright (c) 2021, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: 
#
#     * Redistributions of source code must retain the above copyright 
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials provided 
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
#       nor the names of its contributors may be used to endorse or promote 
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Created by Maximilian Beckers, December 2021

import gzip
import os
import pickle
import pandas as pd
import numpy as np

import click
import wget
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdSubstructLibrary
from tqdm import tqdm

RDLogger.DisableLog("rdApp.warning")

#-------------------------------------------------------------        
def make_rdkit_substr_lib(path):
        
    print("Preparing RDKit substructure library ...");
    chembl_version = 27;
    sdf_path = os.path.join(path, f'chembl_{chembl_version}.sdf.gz');
        
    # Now make the substructure library
    sss_path = os.path.join(path, f'chembl{chembl_version}_sssdata.pkl')
    if not os.path.exists(sss_path):
        click.echo(f'RDKit Version: {rdBase.rdkitVersion}')
        data = []

        with gzip.GzipFile(sdf_path) as gz:
            suppl = Chem.ForwardSDMolSupplier(gz)
            for mol in tqdm(suppl, desc=f'Processing ChEBML {chembl_version}', unit_scale=True):
                                
                if mol is None or mol.GetNumAtoms() > 50:
                    continue
                    
                mol = Chem.AddHs(mol);
                fp = Chem.PatternFingerprint(mol)
                smi = Chem.MolToSmiles(mol)
                data.append((smi, fp))

        click.echo(f'Outputting to {sss_path}')
        with open(sss_path, 'wb') as file:
            mols = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
            fps = rdSubstructLibrary.PatternHolder()
            for smi,fp in data:
                mols.AddSmiles(smi)
                fps.AddFingerprint(fp)
            library = rdSubstructLibrary.SubstructLibrary(mols,fps)
            pickle.dump(library, file, protocol=pickle.HIGHEST_PROTOCOL)
    
    
#-------------------------------------------------------------        
def make_chembl_smiles(path):
    
    chembl_version = 27;
    print("Preparing smiles input file for Arthor ...");
        
    sdf_path = os.path.join(path, f'chembl_{chembl_version}.sdf.gz');
        
    # Now make the substructure library
    sss_path = os.path.join(path, f'chembl_{chembl_version}.smi')
    print(f'RDKit Version: {rdBase.rdkitVersion}');
    data = []

    with gzip.GzipFile(sdf_path) as gz:
        suppl = Chem.ForwardSDMolSupplier(gz)
        for mol in tqdm(suppl, desc=f'Processing ChEBML {chembl_version}', unit_scale=True):
                                
            if mol is None or mol.GetNumAtoms() > 50:
                continue
                    
            mol = Chem.AddHs(mol);
            smi = Chem.MolToSmiles(mol)
            data.append(smi)

    print(f'Outputting to {sss_path}')
    df_out = pd.DataFrame();
    df_out["canonical_smiles"] = data;
    df_out["chembl_id"] = np.arange(len(data));
    
    df_out.to_csv(sss_path, sep=" ", index=False);
    
    
#--------------------------------------
def download_chembl(path):

    chembl_version = 27;
    print("Downloading ChEMBL data ...");
    
    chembl_url = (
        f'https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/'
        f'chembl_{chembl_version}/chembl_{chembl_version}.sdf.gz'
    )

    sdf_path = os.path.join(path, f'chembl_{chembl_version}.sdf.gz')
    if not os.path.exists(sdf_path):
        wget.download(chembl_url, out=path)
    
    return sdf_path;