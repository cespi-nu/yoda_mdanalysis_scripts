#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import argparse

# utils
def f2uf(us, frames):
    uf = {}
    frame = 0
    for i in range(len(us)):
        for f in frames[i]:
            uf[frame] = (i,f)
            frame += 1
    
    return uf

def resid2sel(resis):
    
    return f"protein and resid {resi[0]}:{resi[1]}"

# imput check
def input_check(*, deffnm, psf, xtc, vdwr, csv, pdf, start, end, step, s1, s2, s3, s4, r, cutoff):
    ### files
    if deffnm is None:
        if psf is None:
            print("Error: Specify a PSF file name or use default file name.")
            sys.exit(1)
        if xtc is None:
            print("Error: Specify a XTC file name or use default file name.")
            sys.exit(1)
        if csv is None:
            csv = f"{xtc[:-4]}_fenhole.csv"
        if pdf is None:
            pdf = f"{xtc[:-4]}_fenhole.pdf"
    else:
        if psf is None:
            psf = f"{deffnm}.psf"
        if xtc is None:
            xtc = f"{deffnm}.xtc"
        if csv is None:
            csv = f"{deffnm}.csv"
        if pdf is None:
            pdf = f"{deffnm}.pdf"
    
    spos = xtc.find("*")
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
    
    if vdwr is None:
        if os.path.isfile("/amorphous/yoda/resources/data/vdwr.csv"):
            vdwr = "/amorphous/yoda/resources/data/vdwr.csv"
        else:
            print("Error: Specify a CSV file containing van der waals radius of atoms.")
            sys.exit(1)
    vdwr_dict = {}
    v = open(vdwr, "r")
    for line in v:
        l = line.split(",")
        if l[0] == "atom":
            continue
        vdwr_dict[l[0]] = float(l[1])
    v.close()
    
    print("Input files")
    print(f"PSF file: {psf}")
    print(f"XTC file: {xtcs}")
    print(f"VDWR file: {vdwr}")
    print("Output files")
    print(f"CSV file: {csv}")
    print(f"PDF file: {pdf}")
    
    ### make universe
    us = []
    for x in xtcs:
        us.append(mda.Universe(psf, x))
    
    ### options
    ##### frames
    frames = []
    frms = []
    for u in us:
        if start < 0 or start > len(u.trajectory):
            print("Error: Start frame is out of range.")
            sys.exit(1)
        else:
            first = start
        if end != -1:
            if end < start or end > len(u.trajectory):
                print("Error: End frame is out of range.")
                sys.exit(1)
            else:
                last = end
        else:
            last = len(u.trajectory)
        if len(frames) > 0:
            if end == -1:
                first = 1
        frms.append((first, last, step))
        frames.append(range(first, last, step))
    ##### frame <-> universe, frames
    uf = {}
    frame = 0
    for i in range(len(us)):
        for f in frames[i]:
            uf[frame] = (i,f)
            frame += 1
    
    ##### subunit pairs
    adj = {}
    ##### protein segment argument
    seg_arg = {}
    for s in segids:
        useg = us[0].select_atoms(f"segid {s} and name CA")
        us[0].trajectory[0]
        seg_pos = useg.center_of_geometry()
        seg_arg[s] = np.arctan2(seg_pos[1] - u_cen[1], seg_pos[0] - u_cen[0]) * 180 / np.pi
        if seg_arg[s] <0:
            seg_arg[s] += 360
    seg_sorted = sorted(seg_arg.items(), key=lambda x:x[1])
    seg_order = [s[0] for s in seg_sorted]
    seg_adj = {}
    for i in range(len(seg_order)):
        seg_adj[seg_order[i]] = seg_order[i-1]
    
    ### helices
    helices = {1: s1, 2: s2, 3: s3, 4: s4}
    
    
    print("Options")
    for i in range(len(us)):
        print(f"Frames(universe {i}): {frms[i][0]}-{frms[i][1]} (steps: {frms[i][2]}) {len(frames[i])} frames")
    print(f"Subunits order: {su_order}")
    print("Helices")
    for i in range(4):
        print(f"Helix {i+1}: {helices[i+1][0]}:{helices[i+1][1]}")
    print(f"radius: {r}")
    print(f"cutoff: {cutoff}")
    
    
    ### output args
    args = {"us": us, "frames": frames, "csv": csv, "pdf": pdf, "uf": uf, "adj": seg_adj, "helices": helices, "r": r, "cutoff": cutoff, "vdwr": vdwr_dict}
    
    return args

# analysis
def analysis(*, u, frames, adj, helices, r, cutoff, vdwr):
    pos = {}
    frame = 0
    for s in helices.keys():
        pos[s] = {}
        for r in range(helices[s][0], helices[s][1] + 1):
            pos[s][r] = {}
    for i in range(len(u)):
        for f in frames[i]:
            
    
    return 

# output
### position
def output_pos(pos_df, frames, position):
    c = open(position, "w")
    c.write("frame,subunit,residue,x,y,z\n")
    
    for f in frames:
        for sp in pos_df[f].keys():
            for rp in pos_df[f][sp].keys():
                c.write(str(f) + "," + str(sp) + "," + str(rp) + "," + str(pos_df[f][sp][rp][0]) + "," + str(pos_df[f][sp][rp][1]) + "," + str(pos_df[f][sp][rp][2]) + "\n")
    c.close()
### fenestration distance
def output_dist(fen_df, frames, distance):
    c = open(distance, "w")
    c.write("subunit pair,residue pair")
    for f in frames:
        c.write("," + str(f))
    c.write("\n")
    
    out_df = {}
    for sp in fen_df.keys():
        out_df[sp] = {}
        for rp in fen_df[sp].keys():
            out_df[sp][rp] = []
            for d in fen_df[sp][rp]:
                out_df[sp][rp].append(d)
    
    for sp in out_df.keys():
        for rp in out_df[sp].keys():
            c.write(sp + "," + rp)
            for f in range(len(frames)):
                c.write("," + str(out_df[sp][rp][f]))
            c.write("\n")
    c.close()
### summary
def output_summary(fen_df, summary):
    c = open(summary, "w")
    c.write("subunit pair,residue pair,mean,std,min,max\n")
    
    for sp in fen_df.keys():
        for rp in fen_df[sp].keys():
            c.write(sp + "," + rp)
            c.write("," + str(np.mean(fen_df[sp][rp])))
            c.write("," + str(np.std(fen_df[sp][rp])))
            c.write("," + str(np.min(fen_df[sp][rp])))
            c.write("," + str(np.max(fen_df[sp][rp])))
            c.write("\n")
    c.close()
### pdf
def output_pdf(fen_df, frames, pdf):
    p = PdfPages(pdf)
    for sp in fen_df.keys():
        mat = []
        for rp in fen_df[sp].keys():
            mat.append(fen_df[sp][rp])
        fig, ax = plt.subplots()
        fig.colorbar(ax.matshow(mat), ax=ax)
        ax.set_title(sp)
        ax.matshow(mat)
        ax.set_xlabel("Frame")
        ax.set_xticks(range(0,len(frames),int(len(frames)/10)))
        ax.set_xticklabels(range(0,len(frames),int(len(frames)/10)))
        ax.set_ylabel("residue pair")
        ax.set_yticks(range(len(fen_df[sp].keys())))
        ax.set_yticklabels(fen_df[sp].keys())
        ax.set_aspect(10)
        p.savefig(fig)
    p.close()

# parse arguments
def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("--psf", type=str, default=None, help="Input PSF file.")
    parser.add_argument("--xtc", type=str, default=None, help="Input XTC file.")
    parser.add_argument("--vdwr", type=str, default=None, help="Input CSV file containing van der waals radius of atoms.")
    parser.add_argument("--csv", type=str, default=None, help="Output CSV file.")
    parser.add_argument("--pdf", type=str, default=None, help="Output PDF file.")
    parser.add_argument("--start", type=int, default=0, help="Start frame.")
    parser.add_argument("--end", type=int, default=-1, help="End frame.")
    parser.add_argument("--step", type=int, default=1, help="Step frame.")
    parser.add_argument("--s1", type=int, nargs='2', help="Input the start and end residue number of the S1 helix with a blank.")
    parser.add_argument("--s2", type=int, nargs='2', help="Input the start and end residue number of the S2 helix with a blank.")
    parser.add_argument("--s3", type=int, nargs='2', help="Input the start and end residue number of the S3 helix with a blank.")
    parser.add_argument("--s4", type=int, nargs='2', help="Input the start and end residue number of the S4 helix with a blank.")
    parser.add_argument("--r", type=int, default=2.4, help="Input the border radii of spheres for analysis. default: 2.4 [A] (the phosphate ion radius)")
    parser.add_argument("--space", type=float, default=0.5, help="The space between grids. default: 0.5 [A]")
    parser.add_argument("--cutoff", type=float, default=10.0, help="Input the cutoff distance to stop analysis. default: 10.0 [A]")
    parser.add_argument("--deffnm", type=str, default=None, help="Prefix for output CSV and PDF.")
    
    args = parser.parse_args()
    return args

def main():
    ### parse arguments
    input_args = vars(parse_args())
    ### input check
    args = check_input(**input_args)
    ### analysis
    pos_df = pos(args["u"], args["frames"], args["su_pairs"], args["resi_pairs"])
    fen_df = fen(pos_df, args["su_pairs"], args["resi_pairs"])
    ### output
    output_pos(pos_df, args["frames"], args["position"])
    output_dist(fen_df, args["frames"], args["distance"])
    output_summary(fen_df, args["summary"])
    output_pdf(fen_df, args["frames"], args["pdf"])
    
    return 0

if __name__ == "__main__":
    main()