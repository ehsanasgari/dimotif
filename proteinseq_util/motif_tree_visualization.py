#! /usr/bin/python

# -*- coding: utf-8 -*-
__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu"
__project__ = "DIMOTIF 2018"
__website__ = "llp.berkeley.edu/dimotif"


import sys
sys.path.append('../')
import ete3
import random
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace, TextFace, RectFace, random_color, ProfileFace
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from PyQt4 import QtGui
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
from proteinseq_util.motif_properties import MotifProperties

def get_color_gradient(self):
    cNorm  = colors.Normalize(vmin=0, vmax=1)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plt.get_cmap('coolwarm'))
    color_scale = []
    for scale in np.linspace(0, 1, 201):
            b=[int(x) for x in scalarMap.to_rgba(scale)[:3]]
            hex_color = '#%02x%02x%02x' % (b[0],b[1],b[2])
            [r,g,b,a] = scalarMap.to_rgba(scale, bytes=True)
            color_scale.append( QtGui.QColor( r, g, b, a ) )
    return color_scale

class VisualizeTreeOfMotifs(object):
    '''


    '''
    def __init__(self, nwk_format, motifs):
        self.nwk=nwk_format
        ete3.ProfileFace.get_color_gradient=get_color_gradient
        self.MotifProp = MotifProperties()
        motif2struct, self.propVec=self.MotifProp.get_motifs_pss_biophys(motifs)
        self.motif2struct=dict()
        convert={'e':'Mostly participates in beta ladder', 'g':'Mostly participates in helix', 'h':'Mostly participates in helix', 'n':'Mostly participates in loop, hydrogen bonded turn, irregular structure', 's':'Mostly participates in loop, hydrogen bonded turn, irregular structure', 't':'Mostly participates in loop, hydrogen bonded turn, irregular structure'}
        for x,y in motif2struct.items():
            self.motif2struct[x]=convert[y]

        self.struct2color_dic={'Mostly participates in beta ladder':'#ffeb54', 'Mostly participates in helix':'#fc5901', 'Mostly participates in loop, hydrogen bonded turn, irregular structure':'#7ce9ff'}


    def create_circle(self, filename, title, ignore_branch_length=True):

        plt.clf()
        axis_font = {'size':'29'}
        plt.rc('xtick', labelsize=0.1)
        plt.rc('ytick', labelsize=0.1)
        plt.rc({'font.size':0.2})
        plt.rc('text', usetex=True)

        # legend creation
        if self.motif2struct and self.struct2color_dic:
            leg=[]
            for cls,color in self.struct2color_dic.items():
                leg.append(mpatches.Patch(color=color, label=str(cls)))

        t = Tree(self.nwk)
        node_to_keep=[]
        for l in t.iter_leaves():
            if len(l.name)>1:
                node_to_keep.append(l.name)
        t.prune(tuple(node_to_keep))
        # iterate over tree leaves only
        for l in t.iter_leaves():
            ns = NodeStyle()
            if self.motif2struct and self.struct2color_dic:
                ns["bgcolor"] =self.struct2color_dic[self.motif2struct[l.name]] if l.name in self.motif2struct and self.motif2struct[l.name] in self.struct2color_dic else 'white'
            # Gray dashed branch lines
            #ns["hz_line_type"] = 1
            #ns["hz_line_color"] = "#cccccc"
            #

            l.img_style = ns
            if self.propVec:
                if l.name in self.propVec:
                    l.add_features(profile = self.propVec[l.name])
                    l.add_features(deviation = [0 for x in range(len(self.propVec[l.name]))])
                    l.add_face(ProfileFace(max_v=1, min_v=-1, center_v=0, width=200, height=40, style='heatmap', colorscheme=2), column=3, position='aligned')
            l.name="  "+l.name+"  "

        # Create an empty TreeStyle
        ts = TreeStyle()

        # Set our custom layout function
        ts.layout_fn = VisualizeTreeOfMotifs.layout

        # Draw a tree
        ts.mode = "c"

        # We will add node names manually
        ts.show_leaf_name = False
        # Show branch data
        ts.show_branch_length = False
        ts.show_branch_support = False
        ts.force_topology=ignore_branch_length
        ts.title.add_face(TextFace(title, fsize=20, ftype='Times'), column=15)

        # legend creation
        if self.motif2struct and self.struct2color_dic:
            for k , (cls, col) in enumerate(self.struct2color_dic.items()):
                x=RectFace(15,15, 'black', col)
                #x.opacity=0.5
                ts.legend.add_face(x, column=8)
                ts.legend.add_face(TextFace(' '+str(cls)+'   ', fsize=30,ftype='Times'), column=9)
        ts.legend.add_face(TextFace('--- Properties vector order ↓', fsize=30,ftype='Times'), column=10)
        for y in ['Mean molecular weight of amino acids','Mean flexibility of amino acids','Mean DIWV instability index of sequence','Mean surface accessibility of amino acids','Mean KD hydrophobicity','Mean hydrophilicity of amino acids']:
            x=RectFace(5,5, 'black', 'black')
            #x.opacity=0.5
            ts.legend.add_face(x, column=11)
            ts.legend.add_face(TextFace(' '+y+'   ', fsize=25,ftype='Times'), column=12)

        t.render(filename+'.pdf',tree_style=ts,dpi=5000)


    @staticmethod
    def layout(node):
        if node.is_leaf():
            # Add node name to laef nodes
            N = AttrFace("name", fsize=14, fgcolor="black")
            faces.add_face_to_node(N, node, 0)
        if "weight" in node.features:
            # Creates a sphere face whose size is proportional to node's
            # feature "weight"
            C = CircleFace(radius=node.weight, color="RoyalBlue", style="sphere")
            # Let's make the sphere transparent
            C.opacity = 0.3
            # And place as a float face over the tree
            faces.add_face_to_node(C, node, 0, position="float")
