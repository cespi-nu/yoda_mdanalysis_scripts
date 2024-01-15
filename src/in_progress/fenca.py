#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import argparse

# utils
def srp2sel(sp, rp):
    return (f"segid {sp[0]} and resid {rp[0]} and name CA", f"segid {sp[1]} and resid {rp[1]} and name CA")

# input check
def input_check(*, deffnm, psf, xtc, data, summary, pdf, start, stop, step, rps, resi_A, resi_B):
    ### files
    if deffnm is None:
        if psf is None:
            print("Error: Specify a PSF file name or use default file name.")
            sys.exit(1)
        if xtc is None:
            print("Error: Specify a XTC file name or use default file name.")
            sys.exit(1)
        if data is None:
            data = f"{xtc[:-4]}_lin.csv"
        if summary is None:
            summary = f"{xtc[:-4]}_summary.csv"
        if pdf is None:
            pdf = f"{xtc[:-4]}.pdf"
    else:
        if psf is None:
            psf = f"{deffnm}.psf"
        if xtc is None:
            xtc = f"{deffnm}.xtc"
        if data is None:
            data = f"{deffnm}_lin.csv"
        if summary is None:
            summary = f"{deffnm}_summary.csv"
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
    
    print("Input files")
    print(f"PSF file: {psf}")
    print(f"XTC file: {xtcs}")
    print("Output files")
    print(f"Data CSV file: {data}")
    print(f"Summary CSV file: {summary}")
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
    ##### frame <-> universe, frames
    uf = {}
    frame = 0
    for i in range(len(us)):
        for f in frames[i]:
            uf[frame] = (i,f)
            frame += 1
    
    ##### subunit pairs
    adj = {}
    if su_order is None:
        su = []
        u_protein1 = us[0].select_atoms("protein and name CA")
        for s in u_protein1.segids:
            if s not in su:
                su.append(s)
        su_order = su_orderrted(su)
    else:
        su_order = []
        for s in su_order:
            if len(s) == 1:
                su_order.append("PRO" + s)
            else:
                su_order.append(s)
    for i in range(len(su_order)-1):
        adj[su_order[i]] = su_order[i+1]
    adj[su_order[-1]] = su_order[0]
    
    rpairs = []
    if rps is None:
        if resi_A is None or resi_B is None:
            print("Error: Specify a CSV file of residue pairs or input residue pairs in --resi_A and --resi_B.")
            sys.exit(1)
        else:
            if len(resi_A) != len(resi_B):
                print("Error: Specify the same number of residue pairs.")
                sys.exit(1)
            else:
                for i in range(len(resi_A)):
                    rpairs.append((resi_A[i], resi_B[i]))
    else:
        r = open(rps, "r")
        for l in r:
            if re.match(r"([0-9]+),([0-9]+)", l) is False:
                continue
            rpairs.append((int(l.split(",")[0]), int(l.split(",")[1])))
    
    print("Options")
    for i in range(len(us)):
        print(f"Frames(universe {i}): {frms[i][0]}-{frms[i][1]} (steps: {frms[i][2]}) {len(frames[i])} frames")
    print(f"Subunit pairs: {su_order}")
    print(f"Residue pairs: (A-side, B-side) = {rpairs}")
    
    
    ### output args
    args = {"us": us, "frames": frames, "summary": summary, "data": data, "pdf": pdf, "uf": uf, "adj": adj, "rpairs": rpairs}
    
    return args

# analysis
def analysis(*, us, frames, adj, rpairs):
    rpdist_df = {}
    frame = 0
    for i in range(len(us)):
        for f in frames[i]:
            for sp in adj.items():
                for rp in rpairs:
                    sel1, sel2 = srp2sel()
                    u1 = us[i].select_atoms(sel1)
                    u2 = us[i].select_atoms(sel2)
                    rpdist_df[f"{sp}-{adj[sp]}"][f"{rp[0]}-{rp[1]}"][frame] = np.linalg.norm(u1.positions[0] - u2.positions[0])
            frame += 1
    
    return rpdist_df


# output
def output(*, rpdist_df, frames, data, summary, pdf, uf):
    ### data
    d = open(data, "w")
    d.write("subunit pair,residue pair,universe,frame,distance\n")
    for sp in rpdist_df.keys():
        for rp in rpdist_df[sp].keys():
            for f in rpdist_df[sp][rp].keys():
                d.write(f"{sp},{rp},{uf[f][0]},{uf[f][1]},{rpdist_df[sp][rp][f]}\n")
    d.close()
    
    ### summary
    s = open(summary, "w")
    s.write("subunit pair,residue pair,mean,std,min,min frame,max,max frame\n")
    summary_df = {}
    for sp in rpdist_df.keys():
        summary_df[sp] = {}
        for rp in rpdist_df[sp].keys():
            summary_df[sp][rp] = {"mean": np.mean(rpdist_df[sp][rp]), "std": np.std(rpdist_df[sp][rp])}
            memo = {min: (100,-1), max: (0,-1)}
            for f in rpdist_df[sp].keys():
                if rpdist_df[sp][rp][f] < memo[min][0]:
                    memo[min] = (rpdist_df[sp][rp][f], f)
                if rpdist_df[sp][rp][f] > memo[max][0]:
                    memo[max] = (rpdist_df[sp][rp][f], f)
            summary_df[sp][rp]["min"] = memo[min][0]
            summary_df[sp][rp]["max"] = memo[max][0]
            s.write(f"{sp},{rp},{np.mean(rpdist_df[sp][rp])},{np.std(rpdist_df[sp][rp])},{memo[min][0]},{uf[memo[min][1]]},{memo[max][0]},{uf[memo[max][1]]}\n")
    s.close()
    
    ### pdf
    p = PdfPages(pdf)
    ##### summary
    fig, ax = plt.subplots(ncols=5, nrows=3, tight_layout=True)
    for i in range(len(summary_df.keys())):
        sp = summary_df.keys()[i]
        fig.suptitle(sp)
        rplabels = reversed(summary_df[sp].keys())
        rpdata = []
        rpstd = []
        for rp in rplabels:
            rpdata.append(summary_df[sp][rp]["mean"])
            rpstd.append(summary_df[sp][rp]["std"])
        rpbar = reversed(rpdata)
        rperr = reversed(rpstd)
        ax[i//5, i%5].barh(rplabels, rpbar, xerr=rperr)
    ##### data
    spddata = {}
    for sp in summary_df.keys():
        fig, ax = plt.subplots()
        spddata[sp] = []
        tag = []
        for rp in summary_df[sp].keys():
            spddata[sp].append(rpdist_df[sp][rp].values())
            tag.append(rp)
        ax.matshow(spddata[sp], cmap="viridis")
        ax.set_xlabel("frame")
        ax.set_xticks(range(0, len(rpdist_df[sp][rp].keys()), len(rpdist_df[sp][rp].keys()) // 10))
        ax.set_xticklabels(range(0, len(rpdist_df[sp][rp].keys()), len(rpdist_df[sp][rp].keys()) // 10))
        ax.set_ylabel("protein atom")
        ax.set_yticks(range(len(tag)))
        ax.set_yticklabels(tag)
        ax.set_aspect(5)
    p.savefig()
    
    return 0


# parse arguments
def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("--psf", type=str, default=None, help="Input PSF file.")
    parser.add_argument("--xtc", type=str, default=None, help="Input XTC file. Specify a file name or a series name. If file name contains '\*', this program works for all serial files. default: None")
    parser.add_argument("--data", type=str, default=None, help="Output CSV file of residue-pair distances.")
    parser.add_argument("--summary", type=str, default=None, help="Output CSV file of summary.")
    parser.add_argument("--pdf", type=str, default=None, help="Output PDF file.")
    parser.add_argument("--start", type=int, default=0, help="Start frame. default: 0")
    parser.add_argument("--end", type=int, default=-1, help="End frame. default: -1")
    parser.add_argument("--step", type=int, default=1, help="Step frame. default: 1")
    parser.add_argument("--rps", type=str, default=None, help="Input CSV file of residue pairs. The first column would be the residue numbers of subunit A, and the second column would be the residue numbers of subunit B of subunit pair AB.")
    parser.add_argument("--resi_A", type=int, nargs="*", default=None, help="Residue pairs of residue numbers to measure distances (A-side of subunit pair AB).")
    parser.add_argument("--resi_B", type=int, nargs="*", default=None, help="Residue pairs of residue numbers to measure distances (B-side of subunit pair AB).")
    parser.add_argument("--deffnm", type=str, default=None, help="Prefix for CSV and PDF.")
    
    args = parser.parse_args()
    return args

def main():
    ### parse arguments
    input_args = vars(parse_args())
    ### input check
    args = check_input(**input_args)
    ### analysis
    analysis_args = {"us": args["us"], "frames": args["frames"], "adj": args["adj"], "rpairs": args["rpairs"]}
    rpdist_df = analysis(**analysis_args)
    ### output
    output_args = {"rpdist_df": rpdist_df, "frames": args["frames"], "data": args["data"], "summary": args["summary"], "pdf": args["pdf"], "uf": args["uf"]}
    output(**output_args)
    
    return 0

if __name__ == "__main__":
    main()