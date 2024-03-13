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
# Created by Maximilian Beckers, December 2021, initial Code by Franziska Kruger et al. J. Chem. Inf. Model. 2020, 60, 6, 2888–2902

import numpy as np
import pandas as pd
import os
from AutomatedSeriesClassification import cluster_utils
from rdkit import Chem
from rdkit.Chem import AllChem

#-----------------------------------------------
def read_data(filename, seperator, smiles_column):
    
    df = pd.read_csv(filename, sep=seperator);
    
    num_compounds = len(df[smiles_column].to_list());
    #num_compounds= 10000;
    
    print("");
    print("****************************");
    print("******* Reading data *******");
    print("****************************");

    print("Number of compounds: " + repr(num_compounds));

    fp_list = [];
    compound_indices = [];
    smiles_list = [];
    #project_id = [];
    mol_list = [];

    #clean the data and get fingerprints from the compounds
    for tmp_compound in range(num_compounds):

        try:
            if ("*" in df[smiles_column].iloc[tmp_compound]) | ("." in df[smiles_column].iloc[tmp_compound]):
                continue;      
            
            tmp_mol = Chem.MolFromSmiles(df[smiles_column].iloc[tmp_compound]);
            
            if tmp_mol.GetNumHeavyAtoms() < 2:
                continue;
            
            #tmp_mol = Chem.AddHs(tmp_mol);
            tmp_fp = AllChem.GetMorganFingerprintAsBitVect(tmp_mol, 2);
            fp_list.append(tmp_fp);
            compound_indices.append(tmp_compound);
            mol_list.append(tmp_mol);
            smiles_list.append(df['Structure'][tmp_compound]);
                #project_id.append(tmp_pid);
        except:
            alu = 0; 

        num_printout = 10;
        if (tmp_compound % int((num_compounds/num_printout))) == 0:
            progress = 100*tmp_compound/float(num_compounds);
            print('{:.2f}% finished ...'.format(progress));

    num_compounds = len(compound_indices);
    print("Number of compounds after cleaning: " + repr(num_compounds));

    df_filtered = df.iloc[compound_indices, :];

    # set up project database for arthor substructure matching        
    proj_db = cluster_utils.make_project_db(smiles_list);

    df = 0;
    
    return mol_list, proj_db, df_filtered, smiles_list;
