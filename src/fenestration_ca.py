#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import argparse

# imput check
def check_input(*, psf, xtc, position, distance, summary, pdf, start, end, step, su_order, resi_A, resi_B, deffnm):
    # args list
    args = locals()
    deffnm = args.pop("deffnm")
    psf = args.pop("psf")
    xtc = args.pop("xtc")
    # check input and output files
    if psf is None:
        print("Error: --psf is required.")
        sys.exit(1)
    if xtc is None:
        print("Error: --xtc is required.")
        sys.exit(1)
    if deffnm is None:
        if args["position"] is None:
            args["position"] = xtc[:-4] + "_fen-ca_distance.csv"
        if args["distance"] is None:
            args["distance"] = xtc[:-4] + "_fen-ca_position.csv"
        if args["summary"] is None:
            args["summary"] = xtc[:-4] + "_fen-ca_summary.csv"
        if args["pdf"] is None:
            args["pdf"] = xtc[:-4] + "_fen-ca.pdf"
    else:
        if args["position"] is None:
            args["position"] = deffnm + "_position.csv"
        if args["distance"] is None:
            args["distance"] = deffnm + "_distance.csv"
        if args["summary"] is None:
            args["summary"] = deffnm + "_summary.csv"
        if args["pdf"] is None:
            args["pdf"] = deffnm + ".pdf"
    
    # options
    resi_A = args.pop("resi_A")
    resi_B = args.pop("resi_B")
    start = args.pop("start")
    end = args.pop("end")
    step = args.pop("step")
    # universe
    args["u"] = mda.Universe(psf, xtc)
    # frames
    if start < 0 or start > len(args["u"].trajectory):
        print("Error: Start frame is out of range.")
        sys.exit(1)
    if end != -1:
        if end < start or end > len(args["u"].trajectory):
            print("Error: end frame is out of range.")
            sys.exit(1)
    else:
        end = len(args["u"].trajectory)
    args["frames"] = list(range(start, end, step))
    # residue pairs
    args["resi_pairs"] = list(zip(resi_A, resi_B))
    rp_pri = "resi_pairs:("
    for i in args["resi_pairs"]:
        rp_pri += str(i) + ","
    rp_pri = rp_pri[:-1] + ")"
    print("resi_pairs:" + rp_pri)
    # subunit pairs
    args["su_pairs"] = []
    for i in range(len(su_order)):
        if len(su_order[i]) == 1:
            su_order = "PRO" + su_order[i]
        args["su_pairs"].append((su_order[i-1], su_order[i]))
    
    # check all inputs
    files = {"psf":psf, "xtc":xtc, "posision":args["position"], "distance":args["position"], "summary":args["summary"], "pdf":args["pdf"]}
    options = {"frames":args["frames"], "subunit pairs":args["su_pairs"], "residue pairs":args["resi_pairs"]}
    print("Input and output files:\n")
    for k, v in files.items():
        print(k + ":" + v + "\n")
    print("Options:\n")
    for k, v in options.items():
        print(str(k) + ":" + str(v) + "\n")
    
    return args

# analysis
### position
def pos(u, frames, su_pairs, resi_pairs):
    pos_df = {}
    for f in frames:
        f_dict = {}
        u.trajectory[f]
        for sp1,sp2 in su_pairs:
            sp2_dict = {}
            for r1,r2 in resi_pairs:
                u1 = u.select_atoms("name CA and resid " + str(r1) + " and segid " + sp2)
                sp2_dict[r1] = u1.positions[0]
                u2 = u.select_atoms("name CA and resid " + str(r2) + " and segid " + sp2)
                sp2_dict[r2] = u2.positions[0]
            f_dict[sp2] = sp2_dict
        pos_df[f] = f_dict
    
    return pos_df
### fenestration distance
def fen(pos_df, su_pairs, resi_pairs):
    def dist(p1, p2):
        return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)
    
    fen_df = {}
    for s1,s2 in su_pairs:
        sp = s1 + "-" + s2
        fen_df[sp] = {}
        for r1,r2 in resi_pairs:
            rp = str(r1) + "-" + str(r2)
            fen_df[sp][rp] = []
    for f in pos_df.keys():
        for s1,s2 in su_pairs:
            sp = s1 + "-" + s2
            for r1,r2 in resi_pairs:
                rp = str(r1) + "-" + str(r2)
                d = dist(pos_df[f][s1][r1], pos_df[f][s2][r2])
                fen_df[sp][rp].append(d)
    
    return fen_df

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
    parser.add_argument("--position", type=str, default=None, help="Output CSV file of residue positions.")
    parser.add_argument("--distance", type=str, default=None, help="Output CSV file of residue distances.")
    parser.add_argument("--summary", type=str, default=None, help="Output CSV file of summary.")
    parser.add_argument("--pdf", type=str, default=None, help="Output PDF file.")
    parser.add_argument("--start", type=int, default=0, help="Start frame.")
    parser.add_argument("--end", type=int, default=-1, help="End frame.")
    parser.add_argument("--step", type=int, default=1, help="Step frame.")
    parser.add_argument("--su_order", type=str, nargs='*', default=None, help="Input the subunit name in order with a blank.")
    parser.add_argument("--resi_A", type=str, nargs="*", default=None, help="Residue pairs of residue numbers to measure distances.")
    parser.add_argument("--resi_B", type=str, nargs="*", default=None, help="Residue pairs of residue numbers to measure distances.")
    parser.add_argument("--deffnm", type=str, default="fen_ca", help="Prefix for CSV and PDF.")
    
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