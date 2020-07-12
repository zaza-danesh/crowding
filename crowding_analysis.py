#!/usr/bin/env python3
import argparse, os, subprocess
import numpy as np
import re
# import matplotlib.pyplot 
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tkinter
import scipy.optimize
from scipy.integrate import simps
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit


def func1(x, a, b):
    return a*a + (1.-a*a)*np.exp(-x/b)    


mu_0 = 1.2566370614*1e-6 # N/A^2
h = 0.39903132 #kJ/mol ps
wN = 60.8
wH = 600.13
CSA=-157
gammaN = -27.126
gammaH = 267.522
r_NH=0.101
c=wN*wN*CSA*CSA/3
r_NH6 = r_NH*r_NH*r_NH*r_NH*r_NH*r_NH
d=mu_0*mu_0*h*h*gammaH*gammaH*gammaN*gammaN/(4*np.pi*np.pi*r_NH6)

exp= [5.1, 4.4 ,8.05]

LBL_FONTSIZE=17
TCK_FONTSIZE=14
LG_FONTSIZE=13
CPT_FONTSIZE=20
linestyles = ['_', '-', '--', ':']
markers = ['o', 'v','*', 'd', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p' , 'h', 'H', '+', 'x', 'D',  '|', '_', 'P', 'X', ',', '.']
colors = ('b', 'g', 'c', 'm', 'y', 'k', 'r')
CPT=['A', 'B', 'C', 'D']
for m in Line2D.markers:
    try:
        if len(m) == 1 and m != ' ':
            markers.append(m)
    except TypeError:
        pass
styles = markers




def split_list(a_list):
    half = int(len(a_list)/2)    
    return a_list[:half], a_list[half:]

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-traj_list", "--traj_list", help="",   type=str,  default=None)
    parser.add_argument("-internal_order_parameter", "--internal_order_parameter", help="",   type=str,  default=None)
    parser.add_argument("-internal_tumbling", "--internal_tumbling", help="",   type=str,  default=None)
    parser.add_argument("-rmsf_plot", "--rmsf_plot", help="",   type=str,  default=None)
    parser.add_argument("-hbond_pp_plot", "--hbond_pp_plot", help="",   type=str,  default=None)
    

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    diff_corr = [(5.5*10**(-17))/7., (5.5*10**(-17))/8., (5.5*10**(-17))/9.]  #spatial unit in cm and the correction needed to be divided by viscosity in cm

    path='../RESULTS/'
    args  = parseArguments()
    traj_directory, protein_list = [], []
    phi=0.04
    scale_factor =(1.0-0.5*phi**2) / ((1-phi)**3 ) 
    phi = np.arange(0.04, 0.26, 0.01)

    with open('traj_list.txt') as fn:
        for line in fn:
            line_split = line.split("/")
            traj_directory.append(line.strip())
            protein_list.append(line_split[0])


    if args.rmsf_plot:
        if args.rmsf_plot == "rmsf":
            unit = ""            
            y_protein=[]
            plt.figure(figsize=(5,4))
            for i in range(0, len(traj_directory)+1):
                y_copy=[]
                c=-1
                for j in 1, 2, 4, 8:
                    c=c+1
                    if i < 3:
                        md_name =protein_list[i]+str(j)
                    elif i == 3:
                        md_name =protein_list[i-1]+str(j)
                    y_rep=[]
                    for rep in 1, 2, 3:
                        y_chain = []
                        for chain in range(1, j+1):
                            y=[]
                            if i < 3:
                                with open(path+traj_directory[i]+md_name+'_'+args.rmsf_plot+'.'+str(rep)+'.chain'+str(chain)+'.xvg') as fn:
                                    for line in fn:
                                        if line.find("@") < 0 and line.find("#") < 0:
                                            data = line.split()
                                            y.append(float(data[1]))
                            elif i == 3:
                                with open(path+traj_directory[i-1]+md_name+'_'+args.rmsf_plot+'.chain'+str(chain)+'.xvg') as fn:
                                    for line in fn:
                                        if line.find("@") < 0 and line.find("#") < 0:
                                            data = line.split()
                                            y.append(float(data[1]))
                            y_copy.append(np.power(y,2))
                    if i == 0:
                        plt.subplot(411)
                    elif i == 1:
                        plt.subplot(412)
                    elif i == 2:
                        plt.subplot(413)
                    else:
                        plt.subplot(414)
                    plt.plot(np.mean(y_copy, axis=0), color=colors[c], label=str(j))
                if i == 0:
                    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0., fontsize=LG_FONTSIZE)            
                    plt.text(60, .35, CPT[i], fontsize=CPT_FONTSIZE)
                if i == 1:
                    plt.ylabel('MSF($nm^2$)', fontsize=LBL_FONTSIZE)
                    plt.text(84, .35, CPT[i], fontsize=CPT_FONTSIZE)
                if i == 2:
                    plt.text(112, .35, CPT[i], fontsize=CPT_FONTSIZE)
                if i == 3:
                    plt.xlabel('Residue', fontsize=LBL_FONTSIZE)
                    plt.text(112, .35, CPT[i], fontsize=CPT_FONTSIZE)
                plt.ylim([0, 0.45])
                plt.xticks(fontsize=TCK_FONTSIZE)
                plt.yticks([0, 0.2, 0.4], fontsize=TCK_FONTSIZE)
                plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
            

            plt.savefig('../fig_rmsf.pdf', bbox_inches='tight')

    if args.internal_order_parameter:
        plt.figure(figsize=(5,4))        
        for i in range(0, len(traj_directory)):
            c=-1
            for j in 1, 2, 4, 8:
                c=c+1
                params=[]
                md_name =protein_list[i]+str(j)
                with open(path+traj_directory[i]+md_name+'_order_parameter.xvg') as fn:
                    for line in fn:
                        line = line.strip('\n')
                        line = line.replace('[', '').replace(']', '').replace(',','')
                        for e in line.split(' '):
                            if e != '':
                                params.append(float(e))
                l=len(params[0:int(len(params)/2)])
                if i == 0:
                    plt.subplot(311)
                elif i == 1:
                    plt.subplot(312)
                else:
                    plt.subplot(313)
                plt.errorbar(list(range(1,l)), np.array(params[0:l-1])**2, np.array(params[l:-1])**2, color=colors[c], label=str(j))
            # plt.legend(fontsize=LG_FONTSIZE)
            plt.xticks(fontsize=LG_FONTSIZE)
            plt.yticks(fontsize=LG_FONTSIZE)
            if i == 0:
                plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0., fontsize=LG_FONTSIZE)            
                plt.text(58, .8, CPT[i], fontsize=CPT_FONTSIZE)
            if i == 1:
                plt.ylabel('$S^2$', fontsize=LBL_FONTSIZE)
                plt.text(80, .8, CPT[i], fontsize=CPT_FONTSIZE)
            if i == 2:
                plt.xlabel('Residue', fontsize=LBL_FONTSIZE)
                plt.text(107, .8, CPT[i], fontsize=CPT_FONTSIZE)
            plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
        plt.savefig('../fig_S2.pdf', bbox_inches='tight')

    if args.internal_tumbling:
        plt.figure(figsize=(5,4))
        for i in range(0, len(traj_directory)):
            c=-1
            for j in 1, 2, 4, 8:
                c=c+1
                params=[]
                md_name =protein_list[i]+str(j)
                with open(path+traj_directory[i]+md_name+'_tau_e.xvg') as fn:
                    for line in fn:
                        line = line.strip('\n')
                        line = line.replace('[', '').replace(']', '').replace(',','')
                        for e in line.split(' '):
                            if e != '':
                                params.append(float(e)/1000)
                l=len(params[0:int(len(params)/2)])
                if i == 0:
                    plt.subplot(311)
                elif i == 1:
                    plt.subplot(312)
                else:
                    plt.subplot(313)
                plt.plot(list(range(1,l)), np.array(params[0:l-1]), color=colors[c], label=str(j))
            plt.xticks(fontsize=LG_FONTSIZE)
            plt.yticks(fontsize=LG_FONTSIZE)
            if i == 0:
                plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0., fontsize=LG_FONTSIZE)            
                plt.text(58, 50, CPT[i], fontsize=CPT_FONTSIZE)
            if i == 1:
                plt.ylabel('$\\tau_e (ns)$', fontsize=LBL_FONTSIZE)
                plt.text(80, 130, CPT[i], fontsize=CPT_FONTSIZE)
            if i == 2:
                plt.xlabel('Residue', fontsize=LBL_FONTSIZE)
                plt.text(107, 80, CPT[i], fontsize=CPT_FONTSIZE)
            plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
        plt.savefig('../fig_tau_e.pdf', bbox_inches='tight')




    # if args.order_parameter:
    #     tau_M_protein = []
    #     for i in range(0, len(traj_directory)):
    #         S2_copy = []
    #         count=-1
    #         for j in 1, 2, 4, 8:
    #             count=count+1
    #             md_name =protein_list[i]+str(j)
    #             S2_rep = []
    #             for rep in 1, 2, 3:
    #                 S2_chain = []
    #                 for chain in range(1, j+1):
    #                     params, S2 = [],[]
    #                     # with open('data.dat') as fn:
    #                     with open(path+traj_directory[i]+md_name+'_'+'rotacf_noaver.'+str(rep)+'.'+'chain'+str(chain)+'.init') as fn:
    #                         c=0
    #                         for line in fn:
    #                             line = line.strip('\n')
    #                             line = line.replace('[', '').replace(']', '').replace(',','')
    #                             for e in line.split(' '):
    #                                 if e != '':
    #                                     params.append(float(e))
    #                     p_best = params
    #                     for ii in range(0, int(len(p_best)/2)):
    #                         S2.append(p_best[2*ii])
    #                     S2_rep.append(S2)
    #                     plt.plot(S2)
    #             # plt.legend()
    #             plt.title(md_name)
    #             plt.ylim([0,1])
    #             plt.ylabel('$S^2$')
    #             plt.xlabel('Non-Water Mass Fraction')
    #             plt.show()



    if args.hbond_pp_plot:
        if args.hbond_pp_plot == "hb_num_pp":            
            y_label = "$\\Delta$ protein protein HB"
            fig_name='../fig_protein_hb_num_pp.pdf'        
        y_protein=[]
        for i in range(0, len(traj_directory)):
            y_copy=[]
            for j in 1, 2, 4, 8:
                md_name =protein_list[i]+str(j)
                y=[]
                for rep in 1, 2, 3:
                    with open(path+traj_directory[i]+md_name+'_'+args.hbond_pp_plot+'.'+str(rep)+'.xvg') as fn:
                        for line in fn:
                            if line.find("@") < 0 and line.find("#") < 0:
                                data = line.split()
                                y.append(float(data[1]))

                y_copy.append(np.mean(y)/j)
            plt.plot([4/100, 8/100, 14/100, 25/100], [y_copy[0]-y_copy[0], y_copy[1] - y_copy[0], y_copy[2]-y_copy[0], y_copy[3]-y_copy[0] ], label=protein_list[i], marker=markers[i], color=colors[i], markersize=10, linestyle='--')
        plt.legend(loc=2, fontsize=LG_FONTSIZE)
        plt.xlabel('Non-Water Mass Fraction', fontsize=LBL_FONTSIZE)        
        plt.ylabel(y_label, fontsize=LBL_FONTSIZE)
        plt.xticks(fontsize=TCK_FONTSIZE)
        plt.yticks(fontsize=TCK_FONTSIZE)
        plt.text(-0.02, 2, 'A', fontsize=CPT_FONTSIZE)
        # plt.show()
        plt.savefig(fig_name, bbox_inches='tight')
        plt.close()
                
 