#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os
import re
import argparse

# utils
def dist(p1, p2):
    return np.linalg.norm(p1 - p2)

def pl2ps(pl):
    ps = "protein and segid " + str(pl[0]) + " and resid " + str(pl[1])
    if len(pl) == 3:
        ps += " and name " + str(pl[2])
    return ps

# input check
def input_check(*, deffnm, psf, xtc, lipids, lcsv, start, stop, step, csv, pdf, r):
    ### output files
    spos = xtc.find("*")
    if deffnm is None:
        if csv is None:
            if spos != -1:
                csv = xtc[:spos] + "_nint.csv"
            else:
                csv = xtc[:-4] + "_nint.csv"
        if pdf is None:
            if spos != -1:
                pdf = xtc[:spos] + "_nint.pdf"
            else:
                pdf = xtc[:-4] + "_nint.pdf"
    else:
        if csv is None:
            csv = deffnm + ".csv"
        if pdf is None:
            pdf = deffnm + ".pdf"
    
    ### input files
    if psf is None:
        print("Error: Specify a PSF file name.")
        sys.exit(1)
    if xtc is None:
        print("Error: Specify a XTC file name or a series file name.")
        sys.exit(1)
    
    xtcs = []
    if spos != -1:
        xtc = xtc[:spos] + "01" + xtc[spos + 1:]
        i = 1
        while os.path.isfile(xtc):
            xtcs.append(xtc)
            i = i + 1
            xtc = xtc[:spos] + str(i).zfill(2) + xtc[spos + 2:]
    else:
        xtcs.append(xtc)
    
    ### make universe
    us = []
    for x in xtcs:
        us.append(mda.Universe(psf, x))
    
    ### frames
    frames = []
    frms = []
    for u in us:
        if start < 0 or start > len(u.trajectory):
            print("Error: Start frame is out of range.")
            sys.exit(1)
        else:
            first = start
        if stop != -1:
            if stop < start or stop > len(u.trajectory):
                print("Error: Stop frame is out of range.")
                sys.exit(1)
            else:
                last = stop
        else:
            last = len(u.trajectory)
        if len(frames) > 0:
            if stop == -1:
                first = 1
        frms.append((first, last, step))
        frames.append(range(first, last, step))
    
    ### lipids
    l_ls = lipids
    if lcsv is not None:
        lf = open(lcsv, "r")
        for line in lf:
            l = line.split(",")[0]
            if re.match(r"([0-9]+)", l) is None:
                continue
            if int(l) not in l_ls:
                l_ls.append(int(l))
        lf.close()
    
    ### input check
    print("Input files")
    print(f"Input PSF file: {psf}")
    print(f"Input XTC file: {xtcs}\n")
    print("Output files")
    print(f"Output CSV file: {csv}")
    print(f"Output PDF file: {pdf}\n")
    print("Options")
    for i in range(len(us)):
        print(f"Frames(universe {i}): {frms[i][0]}-{frms[i][1]} (step: {frms[i][2]}) {len(frames[i])} frames")
    print(f"lipids: {l_ls}")
    print(f"border interaction radius: {r}\n")
    
    ### output args
    args = {"us": us, "frames": frames, "l_ls": l_ls, "r": r, "csv": csv, "pdf": pdf}
    
    return args

# ananlysis
def analysis(*, us, frames, l_ls, r):
    ### list all protein atoms to analyse
    nint_ls = {}
    for l in l_ls:
        nint_ls[l] = []
    for i in range(len(us)):
        for f in frames[i]:
            us[i].trajectory[f]
            for l in l_ls:
                upn = us[i].select_atoms(f"protein and name N* and around {r} (resname POPC and resid {l} and name P)", updating=True)
                for a in upn.atoms:
                    if (a.segid, a.resid, a.name) not in nint_ls[l]:
                        nint_ls[l].append((a.segid, a.resid, a.name))
        print(f"universe {i+1} is analysed!")
    
    print("all protein atoms close to lipids' P are listed!")
    
    ### making data frame
    nint_df = {}
    for l in l_ls:
        nint_df[l] = {}
        for a in nint_ls[l]:
            nint_df[l][a] = {}
    
    print("data frame is generated!")
    
    ### analysis
    frame = 0
    for i in range(len(us)):
        for f in frames[i]:
            us[i].trajectory[f]
            frame += 1
            for l in l_ls:
                ulp = us[i].select_atoms(f"resname POPC and resid {l} and name P")
                for pl in nint_ls[l]:
                    upn = us[i].select_atoms(pl2ps((pl)))
                    nint_df[l][pl][frame] = dist(ulp.atoms[0].position, upn.atoms[0].position)
        print(f"universe {i+1} is analysed!")
    
    print("analysis is done!")
    
    return nint_df

### output
def output(*, nint_df, csv, pdf):
    ### csv
    c = open(csv, "w")
    c.write("resid,segid,name,frame,distance\n")
    for l in nint_df.keys():
        for s,r,n in nint_df[l].keys():
            for f,d in nint_df[l][(s,r,n)].items():
                c.write(f"{l},{s},{r},{n},{f},{d}\n")
    c.close()
    
    ### pdf
    p = PdfPages(pdf)
    for l in nint_df.keys():
        if nint_df[l] == {}:
            continue
        fig,ax = plt.subplots(tight_layout=True, figsize=(50,20))
        tag = []
        data = []
        for s,r,n in nint_df[l].keys():
            ##### tag, frames, data
            tag.append(f"{s}_{r}_{n}")
            frames = nint_df[l][(s,r,n)].keys()
            data.append(list(nint_df[l][(s,r,n)].values()))
        fig.colorbar(ax.matshow(data), ax=ax)
        ax.set_title(f"lipid{l}-protein distance")
        ax.matshow(data, cmap="viridis")
        ax.set_xlabel("frame")
        ax.set_xticks(range(0, len(frames), len(frames) // 10))
        ax.set_xticklabels(range(0, len(frames), len(frames) // 10))
        ax.set_ylabel("protein atom")
        ax.set_yticks(range(len(tag)))
        ax.set_yticklabels(tag)
        ax.set_aspect(5)
        p.savefig()
    p.close()

# parser
def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("--deffnm", type=str, default=None, help="Default file name for output files which are not specified. Input a file name without extention. If not specified, XTC file name would be the default file name. default: None")
    parser.add_argument("--psf", type=str, default=None, help="Input PSF file. Specify a file name. default: None")
    parser.add_argument("--xtc", type=str, default=None, help="Input XTC file. Specify a file name or a series name. If file name contains '\*', this program works for all serial files. default: None")
    parser.add_argument("--lipids", type=int, nargs="*", default=[], help="The residue number of lipids to analyse. Use this option or --lcsv option to specify lipids to analyse. default: None")
    parser.add_argument("--lcsv", type=str, required=True, help="Input CSV file containing the residue number of lipids (the summary of lipids_inside.py output is compatible). If --lipids option is specified together, both of their content would be used. The first line would be discarded. required.")
    parser.add_argument("--start", type=int, default=0, help="Start frame. default: 0")
    parser.add_argument("--stop", type=int, default=-1, help="Stop frame. default: -1")
    parser.add_argument("--step", type=int, default=1, help="Step frame. default: 1")
    parser.add_argument("--csv", type=str, default=None, help="Output csv file. Specify a file name or use default file name (--deffnm). default: None")
    parser.add_argument("--pdf", type=str, default=None, help="Output pdf file. Specify a file name or use default file name (--deffnm). default: None")
    parser.add_argument("-r", "--r", type=float, default=5, help="The border interaction radius for analysis. default: 5.0")
    args = parser.parse_args()
    
    return args

# main
def main():
    input_args = vars(parse_args())
    args = input_check(**input_args)
    out_args = {"csv": args.pop("csv"), "pdf": args.pop("pdf")}
    out_args["nint_df"] = analysis(**args)
    output(**out_args)
    
    return 0

if __name__ == "__main__":
    main()