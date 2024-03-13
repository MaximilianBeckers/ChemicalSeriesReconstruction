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

import unittest
import sys
import os

sys.path.append('..')
from ChemicalSeriesReconstruction import ChemicalSeriesReconstruction
import pandas as pd

#--------------------------------------------------------------------------
class TestChemicalSeriesReconstruction(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        curr_path = os.path.dirname(os.path.abspath(__file__));

        df= pd.read_csv(os.path.join(curr_path, "test_smiles.csv"));
        smiles_list = list(df["Structure"].to_numpy()[:1000]);

        min_cluster_size = 10;
        flimit = 0.001;
        scaffolds = None;
        dates = [];
        size_sliding_window = None;
        jaccard_similarity_threshold = None;
        
        cls.given_series_data = pd.read_csv(os.path.join(curr_path, "series_data.csv"));
        cls.given_mcs_data = pd.read_csv(os.path.join(curr_path, "mcs_data.csv"));
        cls.series = ChemicalSeriesReconstruction(smiles_list, min_cluster_size, flimit, scaffolds, dates, size_sliding_window, jaccard_similarity_threshold);
        
    def test_classes(self):
        self.assertListEqual(self.series.series_data["Class"].to_list(), self.given_series_data["Class"].to_list());
        
        
    def test_mcs(self):
        self.assertListEqual(self.series.mcs_data["Scaffold ID"].to_list(), self.given_mcs_data["Scaffold ID"].to_list());
    
#-----------------------------------------------------------------------------
class TestChemicalSeriesReconstructionJaccard(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        
        curr_path = os.path.dirname(os.path.abspath(__file__));
          
        df= pd.read_csv(os.path.join(curr_path, "test_smiles.csv"));
        smiles_list = list(df["Structure"].to_numpy()[:1000]);

        min_cluster_size = 10;
        flimit = 0.001;
        scaffolds = None;
        dates = [];
        size_sliding_window = None;
        jaccard_similarity_threshold = 0.5;
        
        cls.given_series_data = pd.read_csv(os.path.join(curr_path, "series_data_jaccard.csv"));
        cls.given_mcs_data = pd.read_csv(os.path.join(curr_path, "mcs_data_jaccard.csv"));
        cls.series = ChemicalSeriesReconstruction(smiles_list, min_cluster_size, flimit, scaffolds, dates, size_sliding_window, jaccard_similarity_threshold);
        
    def test_classes(self):
        self.assertListEqual(self.series.series_data["Class"].to_list(), self.given_series_data["Class"].to_list());
        
        
    def test_mcs(self):
        self.assertListEqual(self.series.mcs_data["Scaffold ID"].to_list(), self.given_mcs_data["Scaffold ID"].to_list());
        
#-------------------------------------------------------------------------------
class TestChemicalSeriesReconstructionJaccard(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        
        curr_path = os.path.dirname(os.path.abspath(__file__));
        
        df= pd.read_csv(os.path.join(curr_path, "test_smiles.csv"));
        smiles_list = list(df["Structure"].to_numpy()[:1000]);

        min_cluster_size = 10;
        flimit = 0.001;
        scaffolds = None;
        size_sliding_window = 365;
        jaccard_similarity_threshold = 0.5;
        
        #generate timestamps
        import datetime
        from random import randrange

        dates = [];
        current = datetime.datetime(2013, 9, 20,13, 00);
        for i in range(len(smiles_list)):
            current = current + datetime.timedelta(minutes=100);
            dates.append(current)
        
        
        cls.given_series_data = pd.read_csv(os.path.join(curr_path, "series_data_time.csv"));
        cls.given_mcs_data = pd.read_csv(os.path.join(curr_path, "mcs_data_time.csv"));
        cls.series = ChemicalSeriesReconstruction(smiles_list, min_cluster_size, flimit, scaffolds, dates, size_sliding_window, jaccard_similarity_threshold);
        
    def test_classes(self):
        self.assertListEqual(self.series.series_data["Class"].to_list(), self.given_series_data["Class"].to_list());
        
    def test_active_classes(self):
        self.assertListEqual(self.series.series_data["Class (Active)"].to_list(), self.given_series_data["Class (Active)"].to_list());
            
    def test_dates(self):
        self.assertListEqual(self.series.series_data["Date"].to_list(), self.given_series_data["Date"].to_list());
        
    def test_mcs(self):
        self.assertListEqual(self.series.mcs_data["Scaffold ID"].to_list(), self.given_mcs_data["Scaffold ID"].to_list());