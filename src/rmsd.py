#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap, align, rms
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import re
import argparse

# input_check
def input_check(*, deffnm, psf, xtc, selfile, csv, pdf):
    args = locals()
    deffnm = args.pop("deffnm")
    psf = args.pop("psf")
    xtc = args.pop("xtc")
    selfile = args.pop("selfile")
    
    ### file names
    if deffnm is None:
        if psf is None:
            print("Error: Specify a PSF file name or use default file name.")
            sys.exit(1)
        if xtc is None:
            print("Error: Specify a XTC file name or use default file name.")
            sys.exit(1)
        if args["csv"] is None:
            args["csv"] = xtc[:-4]
        if args["pdf"] is None:
            args["pdf"] = xtc[:-4] + "_rmsd.pdf"
    else:
        if psf is None:
            psf = deffnm + ".psf"
        if xtc is None:
            xtc = deffnm + ".xtc"
        if args["csv"] is None:
            args["csv"] = deffnm
        if args["pdf"] is None:
            args["pdf"] = deffnm + ".pdf"
    
    if args["csv"][-4:] == ".csv":
        args["csv"] = args["csv"][:-4]
    
    ### selections
    seles = []
    s = open(selfile, "r")
    for l in s:
        seles.append(l)
    s.close()
    args["seles"] = seles
    
    ### options
    args["u"] = mda.Universe(psf, xtc)
    u_sels = []
    for sel in seles:
        u_sels.append(args["u"].select_atoms(sel))
    
    ### check all inputs
    print("input check")
    print("psf: " + psf)
    print("xtc: " + xtc)
    print("csv: " + args["csv"], "_00.csv")
    print("pdf: " + args["pdf"])
    print("selections:")
    for i in range(len(seles)):
        print(seles[i] + ": " + str(u_sels[i].n_atoms) + " atoms selected.")
    
    return args

# analysis
def rmsd_analysis(*, u, seles):
    ### analysis
    rmsd_df = {}
    for sel in seles:
        align.AlignTraj(u, u, select=sel, in_memory=True).run()
        rmsd_df[sel]=diffusionmap.DistanceMatrix(u, select=sel).run()
    
    return rmsd_df

# output
def output(*, csv, pdf, rmsd_df):
    ### csv
    for i in range(len(rmsd_df)):
        c_i = csv + "_" + str(i).zfill(2) + ".csv"
        key = list(rmsd_df.keys())[i]
        c = open(c_i, "w")
        for j in range(len(rmsd_df[key].results.dist_matrix)):
            for k in range(len(rmsd_df[key].results.dist_matrix[j])):
                if k == 0:
                    c.write(str(rmsd_df[key].results.dist_matrix[j][k]))
                else:
                    c.write("," + str(rmsd_df[key].results.dist_matrix[j][k]))
            c.write("\n")
        c.close()
    
    ### pdf
    p = PdfPages(pdf)
    for i in rmsd_df.keys():
        plt.figure()
        plt.imshow(rmsd_df[i].results.dist_matrix, cmap="viridis")
        plt.xlabel("frame")
        plt.ylabel("frame")
        plt.title(i)
        plt.colorbar(label=r'RMSD ($\AA$)')
        p.savefig()
    p.close()

# parser
def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("--deffnm", type=str, default=None, help="Default file name for input and output files which are not specified. Input a file name without extention. If not specified, PDB file name would be the default file name. default: None")
    parser.add_argument("--psf", type=str, required=True, help="Input PSF file. Specify a file name or use default file name. required")
    parser.add_argument("--xtc", type=str, required=True, help="Input XTC file. Specify a file name or use default file name. required")
    parser.add_argument("--selfile", type=str, required=True, help="Input TXT file which contains the selection lines. required.")
    parser.add_argument("--csv", type=str, default=None, help="Output csv file. Specify a file name. default: None")
    parser.add_argument("--pdf", type=str, default=None, help="Output pdf file. Specify a file name. default: None")
    args = parser.parse_args()
    
    return args

# main
def main():
    input_args = vars(parse_args())
    args = input_check(**input_args)
    output_args = {"csv": args.pop("csv"), "pdf": args.pop("pdf")}
    output_args["rmsd_df"] = rmsd_analysis(**args)
    output(**output_args)
    
    return 0

if __name__ == "__main__":
    main()