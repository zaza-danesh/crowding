#!/usr/bin/env python3

import argparse, os, glob, math
import matplotlib.pyplot as plt
import numpy as np

# base_dir = "/home/zahedeh/CROWDING/"
# big_dir  = "/home/spoel/wd/ub"
protein_list = [ "2kim", "2k57", "ub" ]

COL = ["b", "r", "g", "c"]
C=0
for protein in protein_list:
    C = C+1
    c=-1
    for nprotein in 2, 8:
        if (protein == "2kim" and nprotein == 64) or (protein == "2k57" and nprotein == 64 ) :
            continue
        c=c+1        
        for rep in 1, 2, 3:
            r,y = [],[]
            with open("../RESULTS/%s/rdf_%s%s_%s.xvg" % (protein, protein, nprotein, rep), "r") as infile:
                for line in infile:
                    if line.find("@") < 0 and line.find("#") < 0:
                        data = line.split()
                        r.append(float(data[0]))
                        y.append(float(data[1]))
            # if rep == 1:
            #     plt.plot(r, y, color=COL[c], label=nprotein)
            # else:
            #     plt.plot(r, y, color=COL[c])
                plt.plot(r, y, color=COL[c], label= (nprotein if rep == 1 else None))
            plt.legend()
        plt.title(protein)    
        plt.xlabel('distance ($nm$)')
        plt.ylabel('radial distribution function')
        plt.xlim([1.5, 4.5])    
        # plt.ylim([0, 1.5])    
    plt.savefig('../rdf_%s.pdf' % (protein) )
    plt.close()
