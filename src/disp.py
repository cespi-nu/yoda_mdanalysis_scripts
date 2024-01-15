#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import re
import argparse

# input_check
def input_check(*, deffnm, psf, xtc, selfile, edr, start, stop, step, csv, pdf):
    args = locals()
    deffnm = args.pop("deffnm")
    csv = args.pop("csv")
    psf = args.pop("psf")
    xtc = args.pop("xtc")
    selfile = args.pop("selfile")
    edr = args.pop("edr")
    start = args.pop("start")
    stop = args.pop("stop")
    step = args.pop("step")
    
    ### file names
    if deffnm is None:
        if psf is None:
            print("Error: Specify a PSF file name or use default file name.")
            sys.exit(1)
        if xtc is None:
            print("Error: Specify a XTC file name or use default file name.")
            sys.exit(1)
        if csv is None:
            csv = xtc[:-4]
        if args["pdf"] is None:
            args["pdf"] = xtc[:-4] + "_rmsd.pdf"
    else:
        if psf is None:
            psf = deffnm + ".psf"
        if xtc is None:
            xtc = deffnm + ".xtc"
        if csv is None:
            csv = deffnm
        if args["pdf"] is None:
            args["pdf"] = deffnm + ".pdf"
    
    if csv[-4:] == ".csv":
        csv = csv[:-4]
    
    args["pos"] = csv + "_pos"
    args["disp"] = csv + "_disp.csv"
    
    ### selections
    seles = []
    s = open(selfile, "r")
    for l in s:
        seles.append(l)
    s.close()
    args["seles"]= seles
    
    ### options
    args["u"] = mda.Universe(psf, xtc)
    u_sels = {}
    for sel in seles:
        u_sels[sel] = args["u"].select_atoms(sel)
    
    if start < 0 or start > len(args["u"].trajectory):
        print("Error: Start frame is out of range.")
        sys.exit(1)
    if stop != -1:
        if stop < start or stop > len(args["u"].trajectory):
            print("Error: Stop frame is out of range.")
            sys.exit(1)
    else:
        stop = len(args["u"].trajectory)
    args["frames"] = list(range(start, stop, step))
    
    ### box size
    edr_data = mda.auxiliary.EDR.EDRReader(edr)
    full_frames = len(edr_data.get_data("Box-X")["Time"])
    tbin = full_frames // len(args["u"].trajectory)
    size_df = edr_data.get_data("Box-X")["Box-X"]
    args["size"] = size_df[start*tbin:stop*tbin:step*tbin]
    
    ### check all inputs
    print("input check")
    print("psf: " + psf)
    print("xtc: " + xtc)
    print("pos: " + args["pos"] + "_00.csv")
    print("disp: " + args["disp"])
    print("pdf: " + args["pdf"])
    print("selections:")
    for i in range(len(seles)):
        print(seles[i] + ": " + str(u_sels[list(u_sels.keys())[i]].n_atoms) + " atoms selected.")
    
    return args

# analysis
def disp_analysis(*, u, seles, frames, size):
    ### positions
    pos_df = {}
    for sel in seles:
        pos_df[sel] = []
    for f in frames:
        u.trajectory[f]
        for sel in seles:
            pos_df[sel].append(u.select_atoms(sel).center_of_geometry())
    
    ### displacement
    disp_df = {}
    for sel in pos_df.keys():
        s_disp = 0
        s_pos = pos_df[sel][0]
        for f in range(len(pos_df[sel])):
            pos = pos_df[sel][f]
            d = np.linalg.norm(pos - s_pos)
            if d > size[f] / 2:
                pos2 = pos
                if abs(pos2[0] - s_pos[0]) > size[f] / 2:
                    if pos2[0] > s_pos[0]:
                        pos2[0] -= size[f]
                    else:
                        pos2[0] += size[f]
                if abs(pos2[1] - s_pos[1]) > size[f] / 2:
                    if pos2[1] > s_pos[1]:
                        pos2[1] -= size[f]
                    else:
                        pos2[1] += size[f]
                d = np.linalg.norm(pos2 - s_pos)
            s_disp += d
            s_pos = pos
        disp_df[sel] = s_disp
    
    return (pos_df, disp_df)

# output
def output(*, pos, disp, pdf, pos_df, disp_df, frames, size):
    ### pos
    for i in range(len(disp_df)):
        pos_i = str(pos) + "_" + str(i).zfill(2) + ".csv"
        p = open(pos_i, "w")
        p.write("frame,x,y,z")
        for j in range(len(pos_df[list(pos_df.keys())[i]])):
            p.write("\n" + str(frames[j]) + "," + str(pos_df[list(pos_df.keys())[i]][j][0]) + "," + str(pos_df[list(pos_df.keys())[i]][j][1]) + "," + str(pos_df[list(pos_df.keys())[i]][j][2]))
        p.close()
    
    ### disp
    d = open(disp, "w")
    for sel in disp_df.keys():
        d.write(sel + "," + str(disp_df[sel]) + "\n")
    d.close()
    
    ### pdf
    p = PdfPages(pdf)
    for i in range(len(pos_df)):
        pos_x = []
        pos_y = []
        for j in range(len(pos_df[list(pos_df.keys())[i]])):
            pos_x.append(pos_df[list(pos_df.keys())[i]][j][0])
            pos_y.append(pos_df[list(pos_df.keys())[i]][j][1])
        
        ### plot in whole system
        plt.figure()
        plt.scatter(pos_x, pos_y, s=15, c=frames, cmap="cool")
        plt.colorbar()
        plt.title(list(pos_df.keys())[i])
        plt.xlabel("x")
        plt.xlim(-20, max(size)+20)
        plt.ylabel("y")
        plt.ylim(-10, max(size)+20)
        plt.plot(pos_x, pos_y, "-", alpha=0.3)
        p.savefig()
        
        ### zoom
        plt.figure()
        plt.scatter(pos_x, pos_y, s=15, c=frames, cmap="cool")
        plt.colorbar()
        plt.title(list(pos_df.keys())[i])
        plt.xlabel("x")
        plt.ylabel("y")
        plt.plot(pos_x, pos_y, "-", alpha=0.3)
        p.savefig()
    p.close()

# parser
def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("--deffnm", type=str, default=None, help="Default file name for input and output files which are not specified. Input a file name without extention. If not specified, PDB file name would be the default file name. default: None")
    parser.add_argument("--psf", type=str, required=True, help="Input PSF file. Specify a file name or use default file name. required")
    parser.add_argument("--xtc", type=str, required=True, help="Input XTC file. Specify a file name or use default file name. required")
    parser.add_argument("--selfile", type=str, required=True, help="Input TXT file containing the selection lines. required.")
    parser.add_argument("--edr", type=str, required=True, help="Input EDR file containing box size information for calculation of the displacement. required.")
    parser.add_argument("--start", type=int, default=0, help="Start frame. default: 0")
    parser.add_argument("--stop", type=int, default=-1, help="Stop frame. default: -1")
    parser.add_argument("--step", type=int, default=1, help="Step frame. default: 1")
    parser.add_argument("--csv", type=str, default=None, help="Output csv files containig position / displacement of selected residues. Specify a file name. default: None")
    parser.add_argument("--pdf", type=str, default=None, help="Output pdf file. Specify a file name. default: None")
    args = parser.parse_args()
    
    return args

# main
def main():
    input_args = vars(parse_args())
    args = input_check(**input_args)
    output_args = {"pos": args.pop("pos"), "disp": args.pop("disp"), "pdf": args.pop("pdf"), "frames": args["frames"], "size": args["size"]}
    df = disp_analysis(**args)
    output_args["pos_df"] = df[0]
    output_args["disp_df"] = df[1]
    output(**output_args)
    
    return 0

if __name__ == "__main__":
    main()