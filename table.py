#!/usr/bin/env python3
import argparse, os, subprocess

property_list = ["$\\eta$", "$D_P$", "$\\tau_M$" ]
protein_list = ["2k57", "ub", "2kim"]
file_list = ["viscosity.dat", "diff_coeff.dat","tau.dat"]
box_length = [7. , 8. , 9.]
exp= [5.1, 4.4 ,8.05]
#correction_factor = (5.5*10**(-17))  #spatial unit in cm and the correction needed to be divided by viscosity in cm

with open('../tbl_diff.tex', 'w') as table:
    table.write("\\begin{tabular}{llcccccc}"+"\n")
    table.write("\\hline"+"\n")
    table.write(" &  Protein & 1 copy & 2 copies & 4 copies & 8 copies & 64 copies & {\\em in vitro}  \\\\"+"\n")
    table.write("\\hline"+"\n")
    for i in range(0, len(property_list)):        
        table.write("\\hline"+"\n")
        for j in range(0, len(protein_list)):
            propvec = []
            with open(file_list[i], 'r') as infile:
                for line in infile:
                    if protein_list[j] in line:
                        propvec.append(' &')
                        propvec.append(line.split()[2])
                        propvec.append('$\\pm$ ')
                        propvec.append(line.split()[3])
            if len(propvec) < 20:
                propvec.append('&')                
                propvec.append('-')
            if j == (len(protein_list) - 1 ) /2 :
                if file_list[i] == "tau.dat":
                    table.write(property_list[i]+" &  " + protein_list[j]+ ' '.join(map(str, propvec)) + " &  " +str(exp[j])+" \\\\"+"\n")
                else:
                    table.write(property_list[i]+" &  " + protein_list[j]+ ' '.join(map(str, propvec)) + " &  " +  " -\\\\"+"\n")
            else:
                if file_list[i] == "tau.dat":                    
                    table.write(" &  " + protein_list[j]+ ' '.join(map(str, propvec)) +" &  " + str(exp[j])+" \\\\"+"\n")
                else:
                    table.write(" &  " + protein_list[j]+ ' '.join(map(str, propvec)) + " &  " +" -\\\\"+"\n")
    table.write("\\hline"+"\n")
    table.write("\\end{tabular}"+"\n")

