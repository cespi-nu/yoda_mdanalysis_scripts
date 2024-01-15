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
    
    # check input and output files
    if args["psf"] is None:
        print("Error: --psf is required.")
        sys.exit(1)
    if args["xtc"] is None:
        print("Error: --xtc is required.")
        sys.exit(1)
    if deffnm is None:
        if args["position"] is None:
            args["position"] = args["xtc"][:-4] + "_fen-ca_distance.csv"
        if args["distance"] is None:
            args["distance"] = args["xtc"] + "_fen-ca_position.csv"
        if args["summary"] is None:
            args["summary"] = args["xtc"][:-4] + "_fen-ca_summary.csv"
        if args["pdf"] is None:
            args["pdf"] = args["xtc"][:-4] + "_fen-ca.pdf"
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
    
    files = {"psf":args["psf"], "xtc":args["xtc"], "posision":args["position"], "distance":args["position"], "summary":args["summary"], "pdf":args["pdf"]}
    options = {"start":args["start"], "end":args["end"], "step":args["step"], "su_order":args["su_order"]}
    print("Input and output files:\n")
    for k, v in files.items():
        print(k + ":" + v + "\n")
    print("Options:\n")
    for k, v in options.items():
        print(str(k) + ":" + str(v) + "\n")
    
    args["resi_pairs"] = list(zip(resi_A, resi_B))
    rp_pri = "resi_pairs:("
    for i in args["resi_pairs"]:
        rp_pri += str(i) + ","
    rp_pri = rp_pri[:-1] + ")"
    print("resi_pairs:" + rp_pri)
    
    return args

# analysis and output
def fen_hole(*, psf, xtc, position, distance, summary, pdf, start, end, step, su_order, resi_pairs):
    ### make universe
    u = mda.Universe(psf, xtc)
    
    ### make frame list, subunits pair list and sorted residues list
    if start < 0 and start > len(u.trajectory):
        print("Error: start frame is out of range.")
        sys.exit(1)
    if end == -1:
        end = len(u.trajectory)
    if end <= start and end > len(u.trajectory):
        print("Error: end frame is out of range.")
    fls = range(start, end, step)
    
    su_pairs = []
    for i in range(len(su_order)):
        su_pairs.append((su_order[i-1], su_order[i]))
    
    ### make data frame and write all positions of residues to pos csv file
    pos_df = {}
    p = open(position, "w")
    p.write("frame,resid,chain,x,y,z\n")
    for ts in fls:
        pos_dict_ts = {}
        u.trajectory[ts]
        for r1,r2 in resi_pairs:
            u1 = u.select_atoms("name CA and resid " + r1)
            u2 = u.select_atoms("name CA and resid " + r2)
            for n in range(len(u1.segids)):
                pos_dict_ts[(u1.resids[n], u1.segids[n])] = (u1.positions[n][0], u1.positions[n][1], u1.positions[n][2])
                p.write(str(ts) + "," + r1 + "," + str(u1.segids[n]) + "," + str(u1.positions[n][0]) + "," + str(u1.positions[n][1]) + "," + str(u1.positions[n][2]) + "\n")
            for n in range(len(u2)):
                pos_dict_ts[(u2.resids[n], u2.segids[n])] = (u2.positions[n][0], u2.positions[n][1], u2.positions[n][2])
                p.write(str(ts) + ", " + r2 + ", " + str(u2.segids[n]) + ", " + str(u2.positions[n][0]) + ", " + str(u2.positions[n][1]) + ", " + str(u2.positions[n][2]) + "\n")
        pos_df[ts] = pos_dict_ts
    p.close()
    
    ### distance analysis
    def measure_dist(pos1, pos2):
        return np.sqrt((pos1[0]-pos2[0])**2 + (pos1[1]-pos2[1])**2 + (pos1[2]-pos2[2])**2)
    
    fen_df = {}
    for ts in fls:
        fen_ts = []
        for s1,s2 in su_pairs:
            fen_sp = []
            for r1,r2 in resi_pairs:
                fen_sp.append(measure_dist(pos_df[ts][(int(r1), s2)], pos_df[ts][(int(r2), s1)]))
            fen_ts.append(fen_sp)
        fen_df[ts] = fen_ts
    
    ### write distances to dist csv file
    spairs = []
    for s1, s2 in su_pairs:
        spairs.append(str(s1) + "-" + str(s2))
    rpairs = []
    for r1, r2 in resi_pairs:
        rpairs.append(str(r1) + "-" + str(r2))
    
    d = open(distance, "w")
    d.write("frame, subunit pair, residue pair, distance\n")
    for i in range(len(fls)):
        j = 0
        for s in spairs:
            k = 0
            for r1,r2 in resi_pairs:
                d.write(str(fls[i]) + "," + str(s) + "," + str(r1) + "-" + str(r2) + "," + str(fen_df[i][j][k]) + "\n")
                k += 1
            j += 1
    d.close()
    
    ### output summary csv file
    sp_df = []
    for sp in range(len(spairs)):
        sp_df.append([])
        for rp in range(len(resi_pairs)):
            sp_df[sp].append([])
    for f in range(len(fls)):
        for sp in range(len(spairs)):
            for rp in range(len(resi_pairs)):
                sp_df[sp][rp].append(fen_df[fls[f]][sp][rp])
    
    sp_data = []
    for sp in range(len(spairs)):
        sp_data.append([])
    for sp in range(len(spairs)):
        for rp in range(len(resi_pairs)):
            sp_data[sp].append((np.mean(sp_df[sp][rp]), np.std(sp_df[sp][rp]), np.min(sp_df[sp][rp]), np.max(sp_df[sp][rp])))
    
    s = open(summary, "w")
    s.write("subunit pairs, residue pairs, mean, std, min, max\n")
    for sp in range(len(spairs)):
        for rp in range(len(resi_pairs)):
            s.write(str(spairs[sp]) + "," + str(rpairs[rp]) + "," + str(sp_data[sp][rp][0]) + "," + str(sp_data[sp][rp][1]) + "," + str(sp_data[sp][rp][2]) + "," + str(sp_data[sp][rp][3]) + "\n")
    s.close()
    
    ### plot and write pdf file
    p = PdfPages(pdf)
    rp_labels = []
    for rp1,rp2 in resi_pairs:
        rp_labels.append(str(rp1) + "-" + str(rp2))
    
    ##### summary
    fig, axs = plt.subplots(nrows=3, ncols=5, tight_layout=True, figsize=(20, 14))
    fig.suptitle("summary")
    for ax in axs.flat:
        ax.set(xlabel="residue pair", ylabel="distance [A]")
    for sp in range(len(spairs)):
        
        sp_mean, sp_std, sp_min, sp_max = zip(*sp_data[sp])
        axs[sp//5, sp%5].bar(range(len(resi_pairs)), sp_mean, yerr=sp_std, capsize=10)
        axs[sp//5, sp%5].scatter(range(len(resi_pairs)), sp_min, color="black")
        axs[sp//5, sp%5].scatter(range(len(resi_pairs)), sp_max, color="black")
        axs[sp//5, sp%5].set_xticks(range(len(resi_pairs)))
        axs[sp//5, sp%5].set_xticklabels(rp_labels, rotation=45)
        axs[sp//5, sp%5].set_ylim(0, max(sp_max)//5*5+5)
        axs[sp//5, sp%5].set_title(spairs[sp])
    p.savefig(fig)
    ##### raw data
    for i in range(len(fls)):
        fig, axs = plt.subplots(nrows=3, ncols=5, tight_layout=True, figsize=(20, 14))
        for ax in axs.flat:
            ax.set(xlabel="residue pair", ylabel="distance [A]")
        for sp in range(len(spairs)):
            axs[sp//5, sp%5].plot(fen_df[i][sp], range(len(resi_pairs)))
            axs[sp//5, sp%5].set_title(spairs[sp])
            axs[sp//5, sp%5].set_yticks(range(len(resi_pairs)))
            axs[sp//5, sp%5].set_yticklabels(rp_labels)
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
    # analysis and output
    fen_hole(**args)

if __name__ == "__main__":
    main()