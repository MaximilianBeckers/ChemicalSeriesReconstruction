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

import numba
import numpy as np
import math
from numba import cuda, prange, float32

NUM_FP = 2000;
NB_CLEANING = 256;
NB = 1024;
    
#------------------------------------------------------------
@cuda.jit
def get_dist_matrix_gpu(fp_array, dist_matrix):
    
    x, y = cuda.grid(2);
    
    bit_sum = cuda.local.array(shape=NB, dtype=float32);
    
    if (x < fp_array.shape[0]) and (y < fp_array.shape[0]):
    
        for i in range(bit_sum.shape[0]):
            bit_sum[i] = fp_array[x, i] + fp_array[y, i];
        
    
    if (x < fp_array.shape[0]) and (y < fp_array.shape[0]):
        bitwise_and = 0.0;
        for i in bit_sum:
            if i == 2:
                bitwise_and = bitwise_and + 1.0;

        bitwise_or = 0.0;
        for i in bit_sum:
            if i > 0.0:
                bitwise_or = bitwise_or + 1.0;

        dist_matrix[x,y] = 1.0 - (bitwise_and / bitwise_or);
    
#------------------------------------------------------------
@cuda.jit
def get_dists_to_sample_gpu(fp_sample, fp_array, dist_matrix):
    
    x = cuda.grid(1);
    
    bit_sum = cuda.local.array(shape=NB, dtype=float32);
    
    if x < fp_array.shape[0]:
    
        for i in range(bit_sum.shape[0]):
            bit_sum[i] = fp_sample[i] + fp_array[x, i];
    
        bitwise_and = 0.0;
        for i in bit_sum:
            if i == 2:
                bitwise_and = bitwise_and + 1.0;

        bitwise_or = 0.0;
        for i in bit_sum:
            if i > 0.0:
                bitwise_or = bitwise_or + 1.0;

        dist_matrix[x] = 1.0 - (bitwise_and / bitwise_or);    

#--------------------------------------------------------------
@cuda.jit
def find_outliers_GPU(fp_array, outlier_compds, cutoff, lim_num_neighbors):
        
    x = cuda.grid(1);
        
    dists = cuda.local.array(shape=NUM_FP, dtype=float32);
    bit_sum = cuda.local.array(shape=NB_CLEANING, dtype=float32);
        
    if x < fp_array.shape[0]:

        num_neighbors = 0;
        fp_1 = fp_array[x];
        for tmp_fp_2 in range(fp_array.shape[0]):
            
            for i in range(bit_sum.shape[0]):
                bit_sum[i] = fp_1[i] + fp_array[tmp_fp_2, i];

            bitwise_and = 0.0;
            for i in bit_sum:
                if i == 2:
                    bitwise_and = bitwise_and + 1.0;

            bitwise_or = 0.0;
            for i in bit_sum:
                if i > 0.0:
                    bitwise_or = bitwise_or + 1.0;

            tmp_dist = 1.0 - (bitwise_and / bitwise_or);
            
            if tmp_dist < (1-cutoff):
                num_neighbors = num_neighbors+1;
            
            if num_neighbors >= lim_num_neighbors:
                outlier_compds[x] = False; 
                break;