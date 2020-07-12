#!/usr/bin/env python3
import argparse, os, subprocess, glob, math
import numpy as np
import matplotlib.pyplot as plt
fig, ax = plt.subplots(4)
[h, w] = fig.get_size_inches()
fig.set_figheight(h*1.5)
fig.set_figwidth(w*1.5)
col = ['blue', 'red', 'green', 'orange']
i=0
fs=14
for protein in [ '2k57', 'ub',  '2kim', '2kim_test']:
    i=i+1
    plt.subplot(4,1,i)    
    c=-1
    for n in [1, 2, 4, 8]:
        c=c+1        
        if protein != '2kim_test':
            xvg = "../RESULTS/%s/%s%s_rmsf.*.*chain*xvg" % (protein, protein, str(n))
        else:
            xvg = "../RESULTS/%s/%s%s_rmsf.chain*xvg" % ('2kim', '2kim', str(n))            
        xvg_files = glob.glob(xvg)
        V =[]
        for file in xvg_files:
            v = []
            with open(file, 'r') as infile:
                for line in infile:
                    if line.find('#') < 0 and line.find('@') < 0:
                        v.append(np.power(float(line.split()[1]), 2))
            V.append(v)
        plt.plot(np.arange(1, len(v)+1), np.mean(V, axis=0), color = col[c], label=str(n))
        plt.fill_between(np.arange(1, len(v)+1), np.mean(V, axis=0) + np.std(V, axis=0)/math.sqrt(len(V)), np.mean(V, axis=0) - np.std(V, axis=0)/math.sqrt(len(V)), color = col[c], alpha=0.3)
        plt.xlim([0, len(v)])
    if protein == '2k57': 
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0., fontsize=fs)
        plt.ylim([0, 0.1])
        plt.yticks([ 0.05, 0.1], fontsize=fs)
        plt.text(len(v)+1.5, 0.1, 'A', fontsize=fs+4)
        plt.xticks(fontsize=fs)
    elif protein == 'ub':
        plt.ylim([0, 0.05])
        plt.yticks([ 0.025, 0.05], fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.text(len(v)+3, 0.05, 'B', fontsize=fs+4)
    else:
        
        if i == 3:
            plt.text(len(v)+4, 0.8, 'C', fontsize=fs+4)
            plt.ylim([0, 0.8])
            plt.yticks([ 0.4, 0.8], fontsize=fs)
            plt.xticks(fontsize=fs)
            plt.text(-15, 1.05, 'MSF ($nm^2$)', fontsize=fs+4, rotation=90)
        if i == 4:
            plt.ylim([0, 0.2])
            plt.yticks([ 0.1, 0.2], fontsize=fs)
            plt.xlabel('Residue', fontsize=fs)
            plt.text(len(v)+4, 0.2, 'D', fontsize=fs+4)
            plt.xticks(fontsize=fs)
plt.savefig('../fig_rmsf.pdf', bbox_inches='tight')



