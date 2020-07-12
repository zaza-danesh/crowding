#!/usr/bin/env python3
import argparse, os, subprocess
import matplotlib.pyplot as plt
import math
import numpy as np

# property_list = ["$\\eta$", "$D_P$", "$\\tau_M$" ]
box_length = [7. , 8. , 9.]
LG_FONTSIZE=13
CPT_FONTSIZE=20
markers = ['o', 'v','*', 'd', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p' , 'h', 'H', '+', 'x', 'D',  '|', '_', 'P', 'X', ',', '.']
colors = ('b', 'g', 'c', 'm', 'y', 'k', 'r')
phi=0.00
scale_factor =(1.0-0.5*phi**2) / ((1-phi)**3 ) 
phi = np.arange(0.0, 0.3, 0.01)
fs=16
vol =[]
with open('vol_fraction.dat', 'r') as infile:
    for line in infile:
        v = [float(i) for i in line.split()]
        vol.append(v)

for prop in ["viscosity.dat", "diff_coeff.dat","tau.dat", "water_diff_coeff.dat"]:
    i=-1
    for protein in ["2k57", "ub", "2kim"]:
        i = i+1
        avg_vec, std_vec = [], []
        with open(prop, 'r') as infile:
            for line in infile:
                if protein in line:
                    if line.split()[1] == "1":
                        avg_norm = float(line.split()[2])
                        std_norm = float(line.split()[3])
                        if prop == "viscosity.dat":
                            avg_norm = 1 #0.855*0.65 #water viscosity
                            std_norm = 0
                    avg = float(line.split()[2])
                    std = float(line.split()[3])
                    avg_vec.append(avg/avg_norm)
                    std_vec.append((avg/avg_norm)*math.sqrt( (std/avg)**2 + (std_norm/avg_norm)**2 ) )
        # if prop == "tau.dat":
        #     print(protein, avg_vec, std_vec)
        plt.errorbar(vol[i][:], avg_vec[0:4], std_vec[0:4], label=protein, marker=markers[i], color=colors[i], markersize=10, linestyle='--')

    if prop == "viscosity.dat":
#        y_label = '$\\eta/\\eta_{water}$'
        y_label = '$\\eta$ mPa s'
        x_label= ' '
        plt.xlabel(x_label, fontsize=fs)        
        plt.ylabel(y_label, fontsize=fs)
        plt.xticks([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3], fontsize=fs-2)
        plt.yticks([0, 2, 4, 6, 8, 10, 12, 14], fontsize=fs-2)
        plt.xlim([0.00, 0.3])
        plt.text(-0.05, 15,'A', fontsize=fs+4)
        plt.savefig('../fig_protein_visc.pdf', bbox_inches='tight')
        plt.close()
        # plt.show()
    elif prop == "diff_coeff.dat":
        plt.plot(phi, scale_factor*(((1-phi)**3 ) / (1.0-0.5*phi**2)) , color='gray', label='HS')
        y_label = '$D_{P}/D_{P,_0}$'
        x_label= ' '
        plt.legend(loc=1, fontsize=LG_FONTSIZE)
        plt.xlabel(x_label, fontsize=fs)        
        plt.ylabel(y_label, fontsize=fs)
        plt.xticks([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3], fontsize=fs-2)
        plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6], fontsize=fs-2)
        plt.xlim([0.00, 0.3])
        plt.text(-0.05, 1.65,'B', fontsize=fs+4)
        plt.savefig('../fig_protein_diff.pdf', bbox_inches='tight')
        plt.close()
        # plt.show()
    elif prop == "tau.dat":
        plt.ylabel('$\\tau_M / \\tau_{M,0} $', fontsize=fs)
        x_label= 'Protein Vol. Fraction'        
        plt.xlabel(x_label, fontsize=fs)
        plt.xticks([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3], fontsize=fs-2)
        plt.yticks([0, 5, 10, 15, 20, 25], fontsize=fs-2)
        # plt.ylim([0, 25])
        plt.xlim([0.00, 0.3])
        plt.text(-0.05, 26.5, 'C', fontsize=fs+4)
        plt.savefig('../fig_protein_tumbling.pdf', bbox_inches='tight')
        plt.close()
        # plt.show()
    elif prop == "water_diff_coeff.dat":
        plt.plot(phi, scale_factor*(((1-phi)**3 ) / (1.0-0.5*phi**2)) , color='gray', label='HS')
        y_label = '$D_{W}/D_{W,_0}$'
        x_label= 'Protein Vol. Fraction'        
        plt.xlabel(x_label, fontsize=fs)        
        plt.ylabel(y_label, fontsize=fs)
        plt.xticks([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3], fontsize=fs-2)
        plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], fontsize=fs-2)
        # plt.ylim([0, 1.05])
        plt.xlim([0.00, 0.3])
        plt.text(-0.05, 1.02, 'D', fontsize=fs+4)
        plt.savefig('../fig_water_diff.pdf', bbox_inches='tight')
        plt.close()
    # plt.show()



        



# for i in range(0, len(file_list)):        
#     for j in range(0, len(protein_list)):
#         avg_vec, std_vec = [], []
#         # print(prop, protein_list[j])
#         with open(prop, 'r') as infile:
#             for line in infile:
#                 if protein_list[j] in line:
#                     if line.split()[1] == "1":
#                         avg_norm = float(line.split()[2])
#                         std_norm = float(line.split()[3])
#                         if prop == "viscosity.dat":
#                             avg_norm = 0.855 #water viscosity
#                             std_norm = 0
#                     # avg.append(float(line.split()[2]))
#                     avg = float(line.split()[2])
#                     std = float(line.split()[3])
#                     avg_vec.append(avg/avg_norm)
#                     std_vec.append((avg/avg_norm)*math.sqrt( (std/avg)**2 + (std_norm/avg_norm)**2 ) )
#         if prop == "tau.dat":
#             print(protein_list[j], avg_vec, std_vec)
        
#         plt.errorbar([0.04, 0.08, 0.14, 0.25], avg_vec[0:4], std_vec[0:4], label=protein_list[j], marker=markers[j], color=colors[j], markersize=10, linestyle='--')
