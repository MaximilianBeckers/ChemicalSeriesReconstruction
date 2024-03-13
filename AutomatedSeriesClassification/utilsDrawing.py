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

from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Chem
import numpy as np
import matplotlib.pyplot as plt
import re

#---------------------------------------------------------------------------------
def moltosvg(mol, molSize=(450,250), kekulize=True):
    
    mc = Chem.Mol(mol.ToBinary());
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
            
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    # the MolDraw2D code is not very good at the moment at dealing with atom queries,
    # this is a workaround until that's fixed.
    # The rendering is still not going to be perfect because query bonds are not properly indicated
    opts = drawer.drawOptions()
    for atom in mc.GetAtoms():
        if atom.HasQuery() and atom.DescribeQuery().find('AtomAtomicNum')!=0:
            opts.atomLabels[atom.GetIdx()]=atom.GetSmarts()
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    
    return svg.replace('svg:','')


#-------------------------------------------------------------------------------
def SvgsToGrid(svgs, labels, svgsPerRow=5, molSize=(250,150), fontSize=12):
    
    matcher = re.compile(r'^(<.*>\n)(<rect .*</rect>\n)(.*)</svg>',re.DOTALL) 
    hdr='' 
    ftr='</svg>' 
    rect='' 
    nRows = len(svgs)//svgsPerRow; 
    if len(svgs)%svgsPerRow : nRows+=1;
    blocks = ['']*(nRows*svgsPerRow);
    labelSizeDist = fontSize*5;
    fullSize=(svgsPerRow*(molSize[0]+molSize[0]/10.0),nRows*(molSize[1]+labelSizeDist));
    print(fullSize);
    
    print(labels[0]);
    print(labels[1]);
    
    count=0
    for svg,name in zip(svgs,labels):
        h,r,b = matcher.match(svg).groups()
        if not hdr: 
            hdr = h.replace("width='"+str(molSize[0])+"px'","width='%dpx'"%fullSize[0])
            hdr = hdr.replace("height='"+str(molSize[1])+"px'","height='%dpx'"%fullSize[1])
        if not rect: 
            rect = r
        legend = '<text font-family="sans-serif" font-size="'+str(fontSize)+'px" text-anchor="middle" fill="black">\n'
        legend += '<tspan x="'+str(molSize[0]/2.)+'" y="'+str(molSize[1]+fontSize*2)+'">'+name.split('|')[0]+'</tspan>\n'
        if len(name.split('|')) > 1:
            legend += '<tspan x="'+str(molSize[0]/2.)+'" y="'+str(molSize[1]+fontSize*3.5)+'">'+name.split('|')[1]+'</tspan>\n'
        legend += '</text>\n'
        blocks[count] = b + legend
        count+=1

    for i,elem in enumerate(blocks): 
        row = i//svgsPerRow 
        col = i%svgsPerRow 
        elem = rect+elem 
        blocks[i] = '<g transform="translate(%d,%d)" >%s</g>'%(col*(molSize[0]+molSize[0]/10.0),row*(molSize[1]+labelSizeDist),elem) 
    res = hdr + '\n'.join(blocks)+ftr 
    return res 


#-------------------------------------------------------------------------------
def barplot_vertical(fig, ax, dict_input,ylabel,legend,xticks,legend_x_pos):
    N=len(dict_input)
    N2=len(list(dict_input.values())[0])
    ind=np.arange(N)
    width=0.35
    if N2<=20:
        clist=list(np.arange(0,20,2))+list(np.arange(1,21,2))
    else:
        clist=list(np.arange(0,N2//2*2+2,2))+list(np.arange(1,N2//2*2+1,2))
    colors=plt.cm.tab20(clist[0:len(list(dict_input.values())[0])])
    #p=np.zeros(N)
    bottom=np.zeros(N)
    for i in range(len(list(dict_input.values())[0])):
        ax.bar(ind,[v[i] for v in dict_input.values()], width, bottom=bottom, color=colors[i])
        bottom = bottom + np.array([v[i] for v in dict_input.values()])
    ax.set_ylabel(ylabel,fontsize=22)
    ax.set_xticks(ind)
    ax.set_xticklabels(xticks, rotation='vertical')
    if legend!=[]:
        ax.legend(legend, fontsize=20, ncol=4,bbox_to_anchor=(legend_x_pos, -0.3))
    plt.setp(ax.get_xticklabels(), fontsize=18)
    plt.setp(ax.get_yticklabels(), fontsize=18)
    return fig,ax