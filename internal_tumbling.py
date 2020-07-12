#!/usr/bin/env python3
import argparse, os, subprocess, glob, math
import numpy as np
import matplotlib.pyplot as plt
fig, ax = plt.subplots(3)
[h, w] = fig.get_size_inches()
fig.set_figheight(h*1.5)
fig.set_figwidth(w*1.5)
col = ['blue', 'red', 'green', 'orange']
i=0
fs=14
for protein in [ '2k57', 'ub',  '2kim']:
    i=i+1
    plt.subplot(3,1,i)    
    c=-1
    for n in [1, 2, 4, 8]:
        c=c+1        
        file = "../RESULTS/%s/%s%s_tau_e.xvg" % (protein, protein, str(n))
        v =[]
        with open(file, 'r') as infile:
            for line in infile:
                line = line.strip('\n')
                line = line.replace('[', '').replace(']', '').replace(',','')
                for e in line.split(' '):
                    if e !='':
                        v.append(float(e)/1000.)
        l = len(v[0: int(len(v)/2)])
        plt.subplot(3,1,i)
        plt.plot(np.arange(1, l), np.array(v[0:l-1]),  color=col[c], label=str(n))
        plt.fill_between(np.arange(1, l), np.array(v[0:l-1])-np.array(v[l:-1]) , np.array(v[0:l-1])+np.array(v[l:-1]) , color = col[c], alpha=0.5)
        if i == 1:
            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0., fontsize=fs)            
            plt.text(56, 90, 'A', fontsize=fs)
            plt.yticks(fontsize=fs)
            plt.xticks(fontsize=fs)
        if i == 2:
            plt.text(78, 250, 'B', fontsize=fs)
            plt.ylabel('$\\tau_e (ns)$', fontsize=fs+4)
            plt.yticks(fontsize=fs)
            plt.xticks(fontsize=fs)
        if i == 3:
            plt.text(104, 145, 'C', fontsize=fs)
            plt.yticks(fontsize=fs)
            plt.xticks(fontsize=fs)
            plt.xlabel('Residue', fontsize=fs)

plt.savefig('../fig_tau_e.pdf', bbox_inches="tight")


#             plt.xticks(fontsize=LG_FONTSIZE)
#             plt.yticks(fontsize=LG_FONTSIZE)
#             if i == 0:
#                 plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0., fontsize=LG_FONTSIZE)            
#                 plt.text(58, .8, CPT[i], fontsize=CPT_FONTSIZE)
#             if i == 1:
#                 plt.ylabel('$S^2$', fontsize=LBL_FONTSIZE)
#                 plt.text(80, .8, CPT[i], fontsize=CPT_FONTSIZE)
#             if i == 2:
#                 plt.xlabel('Residue', fontsize=LBL_FONTSIZE)
#                 plt.text(107, .8, CPT[i], fontsize=CPT_FONTSIZE)
#             plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
#         plt.savefig('../fig_S2.pdf', bbox_inches='tight')



