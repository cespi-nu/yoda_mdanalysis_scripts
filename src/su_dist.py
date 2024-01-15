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
def f2uf(us, frames):
    uf = {}
    frame = 0
    for i in range(len(us)):
        for f in frames[i]:
            uf[frame] = (i,f)
            frame += 1
    
    return uf



# input check
def input_check(*, deffnm, psf, xtc, start, stop, step, csv, pdf):
    ### output files
    spos = xtc.find("*")
    if deffnm is None:
        if csv is None:
            if spos != -1:
                csv = xtc[:spos] + "_sudist.csv"
            else:
                csv = xtc[:-4] + "_sudist.csv"
        if pdf is None:
            if spos != -1:
                pdf = xtc[:spos] + "_sudist.pdf"
            else:
                pdf = xtc[:-4] + "_sudist.pdf"
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
    
    ### su_order
    adj = {}
    su = []
    us[0].trajectory[0]
    u_protein1 = us[0].select_atoms("protein and name CA")
    prot_cen = u_protein1.center_of_geometry()
    for s in u_protein1.segids:
        if s not in su:
            su.append(s)
    sargs = {}
    for s in su:
        u_seg = us[0].select_atoms(f"segid {s} and name CA")
        seg_cen = u_seg.center_of_geometry()
        seg_vec = np.array(seg_cen) - np.array(prot_cen)
        seg_arg = np.arctan2(seg_vec[1], seg_vec[0]) * 180 / np.pi
        if seg_arg < 0:
            seg_arg += 360
        sargs[s] = seg_arg
    su_ordered = sorted(sargs.items(), key=lambda x:x[1])
    su_order = []
    for s,a in su_ordered:
        su_order.append(s)
    for i in range(len(su_order)-1):
        adj[su_order[i]] = su_order[i+1]
    adj[su_order[-1]] = su_order[0]
    
    ### uf
    uf = f2uf(us, frames)
    
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
    print(f"Subunit order: {su_order}\n")
    
    ### output args
    args = {"us": us, "frames": frames, "uf": uf, "adj": adj, "csv": csv, "pdf": pdf}
    
    return args

# ananlysis
def analysis(us, frames, adj):
    frame = 0
    su_sels = {}
    su_pos = {}
    for s in adj.keys():
        su_pos[s] = {}
    for i in range(len(us)):
        su_sels[i] = {}
        for s in adj.keys():
            su_sels[i][s] = us[i].select_atoms(f"segid {s}")
        for f in frames[i]:
            us[i].trajectory[f]
            for s in adj.keys():
                su_pos[s][frame] = su_sels[i][s].center_of_geometry()
            frame += 1
    
    su_dist = {}
    labels = []
    for s in adj.keys():
        label = s[3]+adj[s][3]
        su_dist[label] = []
        labels.append(label)
        for f in su_pos[s].keys():
            pos1 = np.array(su_pos[s][f][:2])
            pos2 = np.array(su_pos[adj[s]][f][:2])
            su_dist[s[3]+adj[s][3]].append(np.linalg.norm(pos1-pos2))
    
    return (su_pos, su_dist, labels)

# output
def output(*, csv, pdf, su_pos, su_dist, adj, labels, uf):
    frame = len(uf)
    
    c = open(csv, "w")
    c.write("universe,frame,subunits,distance\n")
    for f in range(frame):
        for sp in su_dist.keys():
            c.write(f"{uf[f][0]},{uf[f][1]},{sp},{su_dist[sp][f]}\n")
    c.close()
    
    p = PdfPages(pdf)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cens = {}
    for s in adj.keys():
        cens[s] = [[],[]]
    for s in adj.keys():
        for f in range(frame):
            cens[s][0].append(su_pos[s][f][0])
            cens[s][1].append(su_pos[s][f][1])
    for s in adj.keys():
        a = ax.scatter(cens[s][0], cens[s][1], label=s, c=range(frame), s=5, cmap="jet")
        if s == "PROA":
            fig.colorbar(a)
    for j, label in enumerate(adj.keys()):
        ax.annotate(label, (su_pos[label][frame-1][0]+2, su_pos[label][frame-1][1]+2))
    p.savefig(fig)
    
    ### distance
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for sp in labels:
        ax.plot(range(frame), su_dist[sp], label=sp)
    ax.legend(labels)
    p.savefig(fig)
    
    p.close()
    
    return 0

# parser
def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("--deffnm", type=str, default=None, help="Default file name for output files which are not specified. Input a file name without extention. If not specified, XTC file name would be the default file name. default: None")
    parser.add_argument("--psf", type=str, default=None, help="Input PSF file. Specify a file name. default: None")
    parser.add_argument("--xtc", type=str, default=None, help="Input XTC file. Specify a file name or a series name. If file name contains '\*', this program works for all serial files. default: None")
    parser.add_argument("--start", type=int, default=0, help="Start frame. default: 0")
    parser.add_argument("--stop", type=int, default=-1, help="Stop frame. default: -1")
    parser.add_argument("--step", type=int, default=1, help="Step frame. default: 1")
    parser.add_argument("--csv", type=str, default=None, help="Output csv file. Specify a file name or use default file name (--deffnm). default: None")
    parser.add_argument("--pdf", type=str, default=None, help="Output pdf file. Specify a file name or use default file name (--deffnm). default: None")
    args = parser.parse_args()
    
    return args

# main
def main():
    input_args = vars(parse_args())
    args = input_check(**input_args)
    analysis_args = {"us": args["us"], "frames": args["frames"], "adj": args["adj"]}
    su_pos, su_dist, labels = analysis(**analysis_args)
    output_args = {"csv": args["csv"], "pdf": args["pdf"], "su_pos": su_pos, "su_dist": su_dist, "adj": args["adj"], "labels": labels, "uf": args["uf"]}
    output(**output_args)
    
    return 0

if __name__ == "__main__":
    main()