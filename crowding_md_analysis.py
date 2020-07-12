#!/usr/bin/env python3

import argparse, os, glob, math
import numpy as np

base_dir = "/home/zahedeh/CROWDING/"
big_dir  = "/home/spoel/wd/ub"
protein_list = [ "2kim", "2k57", "ub" ]

corr_fac_diff = 6.44*10**(-10)

corr_fac_tau = 4.28*10**(-4)

def tpr_name(prot, nprotein):
    if nprotein == 64:
        if prot == "ub":
            return("%s/%d/ub64_1.tpr" % (big_dir, 1))
        else:
            return None
    elif prot == "ub":
        return ("%s/ff_benchmark/ffam-ws/%s%d_prod.1.tpr" % ( base_dir, prot, nprotein ) )
    else:
        return ("%s/ffam-ws/%s/%s%d_prod.1.tpr" % ( base_dir, prot, prot, nprotein ) )

def edr_name(prot, nprotein, replica):
    if nprotein == 64:
        if prot == "ub":
            return("%s/%d/ub64_%d.edr" % (big_dir, replica, replica))
        else:
            return None
    elif prot == "ub":
        return ("%s/ff_benchmark/ffam-ws/%s%d_prod.%d.edr" % ( base_dir, prot, nprotein, replica ) )
    else:
        return ("%s/ffam-ws/%s/%s%d_prod.%d.edr" % ( base_dir, prot, prot, nprotein, replica ) )

def xtc_name(prot, nprotein, replica):
    if nprotein == 64:
        if prot == "ub":
            return ("%s/%d/traj_comp.xtc" % ( big_dir, replica ) )
        else:
            return None
    elif prot == "ub":
        return ("%s/ff_benchmark/ffam-ws/%s%d_prod.%d.xtc" % ( base_dir, prot, nprotein, replica ) )
    else:
        return ("%s/ffam-ws/%s/%s%d_prod.%d.xtc" % ( base_dir, prot, prot, nprotein, replica ) )

def ndx_name(prot, nprotein):
    return ("prot.%d.ndx" % nprotein)
        
def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-gmx", "--gromacs_command",  type=str,  default="gmx")
    parser.add_argument("-make_ndx", "--make_ndx",  action="store_true")
    parser.add_argument("-hbond_intra", "--hbond_intra",  action="store_true")
    parser.add_argument("-hbond_pp", "--hbond_pp",  action="store_true")
    parser.add_argument("-hbond_pw", "--hbond_pw",  action="store_true")
    parser.add_argument("-pairdist", "--pairdist",  action="store_true")
    parser.add_argument("-hb_analyze", "--hb_analyze", action="store_true")
    parser.add_argument("-msd", "--msd", help="",  action="store_true")
    parser.add_argument("-rdf", "--rdf", help="",  action="store_true")
    parser.add_argument("-visco", "--visco", action="store_true")
    parser.add_argument("-diff_coeff", "--diff_coeff", action="store_true")
    parser.add_argument("-extract_visco", "--extract_visco", action="store_true")
    parser.add_argument("-tau", "--tau", action="store_true")
    parser.add_argument("-water_diff_coeff", "--water_diff_coeff", action="store_true")

    args = parser.parse_args()
    return args

def aver_std(hbsum):
    hbs = 0
    hbsum2 = 0
    for h in hbsum:
        hbs += h
        hbsum2 += h*h
    aver = hbs/len(hbsum)
    std = math.sqrt(hbsum2/len(hbsum)-aver*aver)
    return [ aver, std]

def write_histo(histofile, histo, nprotein):
    scale = 0
    for x in range(1,nprotein):
        scale += nprotein-x
        with open(histofile, "w") as hfn:
            for n in sorted(histo):
                hfn.write("%d  %f\n" % (n, float(histo[n])/(3.0*scale)))

def hbpp_analyze(prot, nprotein, tab):
    hbaversum = []
    histo = {}
    for hpp in glob.glob(("hbnum_pp_%d*xvg" % nprotein)):
        hbsum = []
        with open(hpp, "r") as infile:
            for line in infile:
                if not (line.startswith("#") or line.startswith("@")):
                    words = line.split()
                    if float(words[0]) >= Banalyze:
                        count = int(line.split()[1])
                        if count in histo:
                            histo[count] += 1
                        else:
                            histo[count] = 1
                        hbsum.append(count)
        if len(hbsum) > 0:
            [ aver, std ] = aver_std(hbsum)
            hbaversum.append(aver)
    if len(hbaversum) > 0:
        [ aver, std ] = aver_std(hbaversum)
        tab.write("%s  %d  %.1f  %.1f\n" % ( prot, nprotein, aver, std+0.05 ))
        histofile = ("%s_%d_hbpp_histo.xvg" % ( prot, nprotein ))
        write_histo(histofile, histo, nprotein)
    else:
        tab.write("No data for %s %s\n" % ( prot, nprotein ))

def hb_analyze(template, prot, nprotein, tab):
    hbaversum = []
    for hpp in glob.glob(("%s_%d*xvg" % ( template, nprotein))):
        hbsum = []
        with open(hpp, "r") as infile:
            for line in infile:
                if not (line.startswith("#") or line.startswith("@")):
                    words = line.split()
                    if float(words[0]) >= Banalyze:
                        hbsum.append(float(line.split()[1]))
        if len(hbsum) > 0:
            [ aver, std ] = aver_std(hbsum)
            hbaversum.append(aver)
    if len(hbaversum) > 0:
        [ aver, std ] = aver_std(hbaversum)
        tab.write("%s  %d  %.1f  %.0f\n" % ( prot, nprotein, aver, std+0.5 ))
    else:
        tab.write("No data for %s %s\n" % ( prot, nprotein ))


def extract_diff(prot, nprotein, tab, visco, box_size):
    diffaversum = []
    # np.random.randint(2,1)
    # for replica in [1, 2, 3]:
    #     diffsum = []
    #     for nblock in range(0,10):            
    #         with open("msd_%s%s_%s_%s.xvg" % (prot, nprotein, replica, nblock), "r") as infile:
    #             for line in infile:
    #                 if 'Protein]' in line:
    #                     diff = float(line.split()[4]) + corr_fac_diff/(visco*box_size)
    #                     diffsum.append(diff)
    #     if len(diffsum) > 0:            
    #         [aver, std] = aver_std(diffsum)
    #         diffaversum.append(aver)
    
    for i in range(0, 10):
        diffsum = []
        for boot in range(1, 31):
            replica = np.random.randint(3,size=1)[0]+1
            nblock = np.random.randint(10,size=1)[0]
            with open("msd_%s%s_%s_%s.xvg" % (prot, nprotein, replica, nblock), "r") as infile:
                for line in infile:
                    if 'Protein]' in line:
                        diff = float(line.split()[4]) + corr_fac_diff/(visco*box_size)
                        diffsum.append(diff)
        if len(diffsum) > 0:            
            [aver, std] = aver_std(diffsum)
            diffaversum.append(aver)

    print(len(diffaversum))
    if len(diffaversum) > 0:
        [ aver, std ] = aver_std(diffaversum)
        tab.write("%s  %d  %.3f  %.3f\n" % ( prot, nprotein, aver, std ))
        # print( prot, nprotein, aver, std )
    else:
        tab.write("No data for %s %s\n" % ( prot, nprotein ))

def extract_water_diff(prot, nprotein, tab, visco, box_size):
    diffaversum = []
    for replica in [1, 2, 3]:
        diffsum = []
        # for nblock in range(0,10):            
        with open("%s%s_msd_Water.%s.xvg" % (prot, nprotein, replica), "r") as infile:
            for line in infile:
                if 'Water]' in line:
                    diff = float(line.split()[4]) + corr_fac_diff/(visco*box_size)
                    diffsum.append(diff)
        if len(diffsum) > 0:
            [aver, std] = aver_std(diffsum)
            diffaversum.append(aver)
    if len(diffaversum) > 0:
        [ aver, std ] = aver_std(diffaversum)
        

        tab.write("%s  %d  %.3f  %.3f\n" % ( prot, nprotein, aver, std ))
        print( prot, nprotein, aver, std )
    else:
        tab.write("No data for %s %s\n" % ( prot, nprotein ))

def extract_tau(prot, nprotein, tab, visco, box_size):
    tauaversum = []
    for replica in [1, 2, 3]:
        tauchain = []
        for chain in range(1, nprotein+1):
            params = []
            with open('../../CROWDING/ffam-ws/%s/%s%d_rotacf_noaver.%d.chain%d.best' % (prot, prot, nprotein, replica, chain)) as infile:
                c=0
                for line in infile:
                    line = line.strip('\n')
                    line = line.replace('[', '').replace(']', '')
                    for e in line.split(' '):
                        if e != '':
                            params.append(float(e))
            p_best = params
            temp = 1./(p_best[-1])
            # print(prot, nprotein, p_best[-1], temp)
            temp = temp + corr_fac_tau/(visco*box_size*box_size*box_size)
            temp = 1./(temp)
            tauchain.append(temp)

        if len(tauchain) > 0:
            [aver, std] = aver_std(tauchain)
            # temp = 1./aver
            # temp = temp + 
            # aver = 1./temp            
            tauaversum.append(aver)

    if len(tauaversum) > 0:
        [ aver, std ] = aver_std(tauaversum)
        tab.write("%s  %d  %.1f  %.1f\n" % ( prot, nprotein, aver/1000, std/1000 ))
    else:
        tab.write("No data for %s %s\n" % ( prot, nprotein ))


def extract_visco(prot, nprotein, tab):
    viscoaver = []
    for replica in [1, 2, 3]:
        with open("evis_%s%s_%s.xvg" % (prot, nprotein, replica), "r") as infile:
            for line in infile:
                if line.find("@") < 0 and line.find("#") < 0:
                    data = line.split()
                    v = float(data[4])
        viscoaver.append(v*1000)        
    if len(viscoaver) > 0:
        [aver, std] = aver_std(viscoaver)
        tab.write("%s  %d  %.3f  %.3f\n" % ( prot, nprotein, aver, std ))        
        # print( prot, nprotein, aver, std )
        return aver    
    else:
        tab.write("No data for %s %s\n" % ( prot, nprotein ))

def do_visco():
    for protein in protein_list:
        os.chdir(protein)
        for j in 1, 2, 4, 8, 64:
            mydir = str(j)
            os.mkdir(mydir)
            os.chdir(mydir)
            myjob = "vis_" + protein + str(j) + ".job"
            print(myjob)
            job = open(myjob, "w")
            job.write("#!/bin/sh\n")
            job.write("#SBATCH -t 71:00:00\n")
            job.write("#SBATCH -n 1\n")
            mylogfile = myjob[:-4] + ".out"
            job.write("#SBATCH -o %s\n" % mylogfile)
            for rep in 1, 2, 3:
                md_name = protein + str(j) + "_" + str(rep)
                edr = edr_name(protein, j, rep)
                if not edr:
                    continue
                base_name = base_dir+md_name

                vis  = "vis_" + md_name + ".xvg"
                exvg = "ener_" + md_name + ".xvg"
                if not os.path.isfile(vis):
                    job.write(("%s energy -f %s -vis %s -o %s\n\n" %
                               ( gromacs_command, 
                                 edr, vis, exvg)))
                    job.write("mv evisco.xvg %s\n" % ( "e" + vis))
                    job.write("mv eviscoi.xvg %s\n" % ( "ei" + vis))
            job.close()
            os.system(("sbatch %s" % myjob))
            os.chdir("..")
        os.chdir("..")
        
def do_msd():
    for protein in protein_list:
        for j in 1, 2, 4, 8, 64:
            os.chdir(protein)
            for rep in 1, 2, 3:
                md_name = protein + str(j) + "_" + str(rep)
                myjob = "msd_" + md_name + ".job"
                xtc = xtc_name(protein, j, rep)
                tpr = tpr_name(protein, j)
                if not xtc or not tpr:
                    continue
                ndx = ndx_name(protein, j)
                job = open(myjob, "w")
                job.write("#!/bin/sh\n")
                job.write("#SBATCH -t 48:00:00\n")
                job.write("#SBATCH -n 1\n")
                base_name = base_dir+md_name
                nblocks = 10
                simlength = 1000
                
                msd  = "msd_" + md_name + "_full.xvg"
                dmol = "diffmol_" + md_name + "_full.xvg"
                if not os.path.isfile(msd) or not os.path.isfile(dmol):
                        job.write(("echo 1 | %s msd -tu ns -f %s -s %s -o %s -mol %s -trestart %d\n\n" %
                               ( gromacs_command, 
                                 xtc, tpr, msd, dmol,
                                 simlength+1)))
                for t in range(nblocks):
                    msd  = "msd_" + md_name + "_" + str(t) + ".xvg"
                    dmol = "diffmol_" + md_name + "_" + str(t) + ".xvg"
                    if not os.path.isfile(msd) or not os.path.isfile(dmol):
                        job.write(("echo 1 | %s msd -tu ns -f %s -s %s -o %s -b %s -e %s -mol %s -trestart %d\n\n" %
                               ( gromacs_command, 
                                 xtc, tpr, msd,
                                 t*simlength/nblocks, 
                                 (t+1)*simlength/nblocks,
                                 dmol,
                                 simlength+1)))
                job.close()
                os.system("sbatch " + myjob)
            os.chdir("..")

def do_rdf():
    for protein in protein_list:
        for j in 1, 2, 4, 8, 64:
            os.chdir(protein)
            for rep in 1, 2, 3:
                md_name = protein + str(j) + "_" + str(rep)
                myjob = "rdf_" + md_name + ".job"
                xtc = xtc_name(protein, j, rep)
                tpr = tpr_name(protein, j)
                if not xtc or not tpr:
                    continue
                ndx = ndx_name(protein, j)
                job = open(myjob, "w")
                job.write("#!/bin/sh\n")
                job.write("#SBATCH -t 48:00:00\n")
                job.write("#SBATCH -n 1\n")
                base_name = base_dir+md_name
                # nblocks = 10
                selrpos = "mol_com"
                seltype = "whole_mol_com"
                sel = "Protein"
                ref = "Protein"
                dt = 0.01
                begin = 500
                # simlength = 1000
                
                rdf  = "rdf_" + md_name + ".xvg"
                # dmol = "diffmol_" + md_name + "_full.xvg"
                if not os.path.isfile(rdf):
                    job.write(("%s rdf -tu ns -f %s -s %s -o %s --selrpos %s --seltype %s -sel %s -ref %s -b %s -dt %s \n\n" %
                           ( gromacs_command, 
                             xtc, tpr, rdf, 
                              selrpos, seltype, sel, ref, begin, dt)))                
                job.close()
                os.system("sbatch " + myjob)
            os.chdir("..")


if __name__ == '__main__':
    args  = parseArguments()
    gromacs_command = args.gromacs_command
    B  = 0
    E  = 1000
    dt = 0.1
    Banalyze = 500
    L = [9, 7, 8 ] 

    if args.msd:
        do_msd()
        exit(1)

    if args.rdf:
        do_rdf()
        exit(1)

    tabpp = None
    tabpw = None
    tabhb = None
    tabdiff = None
    tabWdiff = None
    tabvisco = None
    tabtau = None
    if args.hb_analyze:
        tabpp = open("hbnum_pp.dat", "w")
        tabpw = open("hbnum_pw.dat", "w")
        tabhb = open("hbnum_intra.dat", "w")
    if args.diff_coeff:
        tabdiff = open("diff_coeff.dat", "w")
        tabvisco = open("viscosity.dat", "w")        
    if args.tau:
        tabtau = open("tau.dat", "w")
        tabvisco = open("viscosity.dat", "w")        

    if args.water_diff_coeff:
        tabWdiff = open("water_diff_coeff.dat", "w")
        tabvisco = open("viscosity.dat", "w")
        
    if args.extract_visco:
    	tabvisco = open("viscosity.dat", "w")
    if args.visco:
        do_visco()
    # for prot in protein_list:
    for index, prot in enumerate(protein_list):
        os.makedirs(prot, exist_ok=True)
        os.chdir(prot)
        # for nprotein in [ 1, 2, 4, 8]:
        for nprotein in [ 1, 2, 4, 8, 64 ]:
            tpr = tpr_name(prot, nprotein)
            if not tpr:
                continue
            pndx = ndx_name(prot, nprotein)
            if args.make_ndx:
                cmd  = ("%s make_ndx -f %s -o %s < ../split_chain.txt" % ( gromacs_command, tpr, pndx ) )
                os.system(cmd)              
            elif args.hbond_intra:
                for replica in [ 1, 2, 3 ]:
                    jobfn = ("hb_intra_%s.%d.%d.sh" % (prot, nprotein, replica ))
                    job = open(jobfn, "w")
                    job.write("#!/bin/sh\n")
                    job.write("#SBATCH -c 4\n")
                    job.write("#SBATCH -t 24:00:00\n")
                    for ch1 in range(nprotein):
                        chain1 = ("Protein_chain_%d" % ch1)
                        if nprotein == 1:
                            chain1 = 0
                        hbnum = ("hbnum_intra_%d_%d.xvg" % ( nprotein, replica ))
                        cmd = ("echo %s %s | %s hbond -tu ns -dt %s -b %s -e %s -s %s -n %s -f %s -num %s" % ( chain1, chain1, gromacs_command, dt, B, E, tpr, pndx, xtc, hbnum ))
                        job.write("%s\n" % cmd)
                    job.close()
                    os.system("sbatch " + jobfn)
            elif args.pairdist:
                for replica in [ 1, 2, 3 ]:
                    xtc = xtc_name(prot, nprotein, replica)
                    jobfn = ("pairdist_%s.%d.%d.sh" % (prot, nprotein, replica ))
                    job = open(jobfn, "w")
                    job.write("#!/bin/sh\n")
                    job.write("#SBATCH -c 4\n")
                    job.write("#SBATCH -t 24:00:00\n")                    
                    for ch1 in range(nprotein):
                        ch1=ch1+1
                        chain1 = ("chain%d" % ch1)
                        for ch2 in range(ch1+1, nprotein+1):
                            chain2 = ("chain%d" % ch2)
                            output = ("dist_%d_%d_%d_%d.xvg" % ( nprotein, replica, ch1, ch2))
                            cmd = ("%s pairdist -tu ns -dt %s -b %s -e %s -s %s -n %s -f %s -ref %s -sel %s -selrpos whole_mol_com -o %s " % ( gromacs_command, dt, B, E, tpr, pndx, xtc, chain1, chain2, output ))
                            job.write("%s\n" % cmd)
                    job.close()
                    os.system("sbatch " + jobfn)
            elif args.hbond_pp:
                for replica in [ 1, 2, 3 ]:
                    xtc = xtc_name(prot, nprotein, replica)
                    jobfn = ("hbpp_%s.%d.%d.sh" % (prot, nprotein, replica ))
                    job = open(jobfn, "w")
                    job.write("#!/bin/sh\n")
                    job.write("#SBATCH -c 4\n")
                    job.write("#SBATCH -t 24:00:00\n")
                    for ch1 in range(nprotein):
                        ch1=ch1+1
                        chain1 = ("chain%d" % ch1)
                        for ch2 in range(ch1+1, nprotein+1):
                            chain2 = ("chain%d" % ch2)
                            hbnum = ("hbnum_pp_%d_%d_%d_%d.xvg" % ( nprotein, replica, ch1, ch2))
                            cmd = ("echo %s %s | %s hbond -tu ns -dt %s -b %s -e %s -s %s -n %s -f %s -num %s" % ( chain1, chain2, gromacs_command, dt, B, E, tpr, pndx, xtc, hbnum ))
                            job.write("%s\n" % cmd)
                    job.close()
                    os.system("sbatch " + jobfn)
            elif args.hbond_pw:
                for replica in [ 1, 2, 3 ]:
                    xtc = xtc_name(prot, nprotein, replica)
                    jobfn = ("hbpw_%s.%d.%d.sh" % (prot, nprotein, replica ))
                    job = open(jobfn, "w")
                    job.write("#!/bin/sh\n")
                    job.write("#SBATCH -c 4\n")
                    job.write("#SBATCH -t 24:00:00\n")
                    for ch1 in range(nprotein):
                        chain1 = ("chain%d" % ch1)
                        if nprotein == 1:
                            chain1 = 0
                        chain2 = nprotein
                        hbnum = ("hbnum_pw_%d_%d_%d.xvg" % ( nprotein, replica, ch1))
                        hbac = ("hbac_pw_%d_%d_%d.xvg" % ( nprotein, replica, ch1))
                        cmd = ("echo %s %s | %s hbond -dt %f -tu ns -b %s -e %s -s %s -n %s -f %s -num %s" % ( chain1, chain2, gromacs_command, dt, B, E, tpr, pndx, xtc, hbnum ))
                        job.write("%s\n" % cmd)
                    job.close()
                    os.system("sbatch " + jobfn)
            elif args.hb_analyze:
                hbpp_analyze(prot, nprotein, tabpp)
                hb_analyze("hbnum_pw", prot, nprotein, tabpw)
                hb_analyze("hbnum_intra", prot, nprotein, tabhb)

            elif args.diff_coeff:
                box_size = L[index]
                if nprotein == 64:
                    box_size = 15.8
                    # print(box_size)
                mydir = str(nprotein)
                os.chdir(mydir)
                visco = extract_visco(prot, nprotein, tabvisco)
                os.chdir("..")
                extract_diff(prot, nprotein, tabdiff, visco, box_size)
            elif args.water_diff_coeff:
                box_size = L[index]
                if nprotein == 64:
                    box_size = 15.8                    
                    # print(box_size)
                mydir = str(nprotein)
                os.chdir(mydir)
                visco = extract_visco(prot, nprotein, tabvisco)
                os.chdir("..")
                extract_water_diff(prot, nprotein, tabWdiff, visco, box_size)

            elif args.tau:
                box_size = L[index]
                if nprotein == 64:
                    box_size = 15.8
                    # print(box_size)
                mydir = str(nprotein)
                os.chdir(mydir)
                visco = extract_visco(prot, nprotein, tabvisco)
                os.chdir("..")
                extract_tau(prot, nprotein, tabtau, visco, box_size)                                  
            elif args.extract_visco:
                mydir = str(nprotein)
                os.chdir(mydir)
                extract_visco(prot, nprotein, tabvisco)
                os.chdir("..")
        os.chdir("..")
    for tab in [ tabpp, tabpw, tabhb, tabdiff, tabvisco, tabtau, tabWdiff]:
        if tab:
            tab.close()

        
