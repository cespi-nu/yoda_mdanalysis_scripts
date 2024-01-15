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

# input check
def input_check(*, deffnm, psf, xtc, start, stop, step, r1, r2, atom, zcen, upmax, csv):
    ### files
    if deffnm is None:
        if psf is None:
            print("Error: Specify a PSF file name or use default file name.")
            sys.exit(1)
        if xtc is None:
            print("Error: Specify a XTC file name or use default file name.")
            sys.exit(1)
        if csv is None:
            csv = f"{xtc[:-4]}_lin"
    else:
        if psf is None:
            psf = f"{deffnm}.psf"
        if xtc is None:
            xtc = f"{deffnm}.xtc"
        if csv is None:
            csv = deffnm
    
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
    
    if csv[-4:] == ".csv":
        csv = csv[:-4]
    
    pymol = f"{csv}.py"
    summary = f"{csv}_summary.csv"
    data = f"{csv}_data.csv"
    
    print("Input files")
    print(f"PSF file: {psf}")
    print(f"XTC file: {xtcs}")
    print("Output files")
    print(f"Summary csv file: {summary}")
    print(f"Data csv file: {data}")
    print(f"Pymol script file: {pymol}")
    
    ### make universe
    us = []
    for x in xtcs:
        us.append(mda.Universe(psf, x))
    
    ### options
    if r1 <= 0:
        print("Error: Specify the cutoff radius [--r1].")
        sys.exit(1)
    if r2 is None:
        r2 = r1
    elif r2 <= 0 or r2 > r1:
        print("Error: Specify the border radius [--r2].")
        sys.exit(1)
    
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
    
    uf = {None: (None, None)}
    frame = 0
    for i in range(len(us)):
        for f in frames[i]:
            uf[frame] = (i,f)
            frame += 1
    
    if atom is None:
        atom = ["C13", "C14", "C15", "C218", "C316"]
    else:
        if atom[0] == "all":
            atom = u.select_atoms("resname POPC and resid 1").atoms.names
    
    if zcen is None:
        if upmax != -1:
            l_upper = us[0].select_atoms(f"resname POPC and resid 1:{upmax}")
            l_lower = us[0].select_atoms(f"resname POPC and resid {upmax+1}:100000")
            zcen = (l_upper.center_of_geometry()[2] + l_lower.center_of_geometry()[2]) / 2
            nl = upmax + l_lower.n_residues
        else:
            protein = us[0].select_atoms("protein")
            zcen = protein.center_of_geometry()[2]
            l_upper = us[0].select_atoms(f"resname POPC and prop z >= {zcen} and name P")
            lipid = us[0].select_atoms("resname POPC")
            nl = lipid.n_residues
            upmax = l_upper.n_residues
    
    print("Options")
    print(f"Lipids to analyse\n  upper leaflet: {upmax}\n  lower leaflet: {nl - upmax}")
    print(f"Atom to analyse: {atom}")
    print(f"Cutoff radius: {r1}")
    print(f"Border radius: {r2}")
    for i in range(len(us)):
        print(f"Frames(universe {i}): {frms[i][0]}-{frms[i][1]} (steps: {frms[i][2]}) {len(frames[i])} frames")
    print(f"Border z-coordinate: {zcen}")
    
    
    ### output args
    args = {"us": us, "frames": frames, "r1": r1, "r2": r2, "atom": atom, "summary": summary, "data": data, "pymol": pymol, "nl": nl, "upmax": upmax, "uf": uf, "psf": psf, "xtcs": xtcs,}
    
    return args

# analysis
def analysis(*, us, frames, r1, r2, atom, nl, upmax):
    ### segids
    segids = []
    u_pca = us[0].select_atoms("protein and name CA")
    for pca in u_pca.segids:
        if pca not in segids:
            segids.append(pca)
    
    ### making data frame
    lin_df = {}
    for l in range(nl):
        lin_df[l+1] = {"arg": {}, "leaflet": None, "fen": {}, "first": None}
        for a in atom:
            lin_df[l+1][a] = {}
    for s in segids:
        lin_df[s] = {"arg": {}}
    
    ### data collection
    frame = 0
    for i in range(len(us)):
        u_protein = us[i].select_atoms("protein")
        u_lipid = us[i].select_atoms("resname POPC")
        for f in frames[i]:
            us[i].trajectory[f]
            u_cen = u_protein.center_of_geometry()
            ### lipid atom radius
            lins = []
            for a in atom:
                ula = u_lipid.select_atoms(f"name {a}")
                for la in ula.atoms:
                    rad2 = (la.position[0] - u_cen[0]) ** 2 + (la.position[1] - u_cen[1]) ** 2
                    if rad2 <= r1 ** 2:
                        l = la.resid
                        lins.append(l)
                        rad = np.sqrt(rad2)
                        lin_df[l][a][frame] = rad
                        if lin_df[l]["first"] is None:
                            if rad <= r2:
                                lin_df[l]["first"] = frame
            ### lipid center argument
            for l in lins:
                ul = u_lipid.select_atoms(f"resid {l}")
                l_cen = ul.center_of_geometry()
                lin_df[l]["arg"][frame] = np.arctan2(l_cen[1] - u_cen[1], l_cen[0] - u_cen[0]) * 180 / np.pi
                if lin_df[l]["arg"][frame] <0:
                    lin_df[l]["arg"][frame] += 360
                if l <= upmax:
                    lin_df[l]["leaflet"] = "upper"
                else:
                    lin_df[l]["leaflet"] = "lower"
            ### frame update
            frame += 1
    
    ### fenestration
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
    ##### determine fenestration
    for f in range(frame):
        for l in range(nl):
            if lin_df[l+1]["leaflet"] is None:
                continue
            else:
                for s,a in seg_sorted:
                    if f in lin_df[l+1]["arg"].keys():
                        if lin_df[l+1]["arg"][f] <= a:
                            lin_df[l+1]["fen"][f] = seg_adj[s][3] + s[3]
                            break
                        else:
                            lin_df[l+1]["fen"][f] = seg_adj[list(seg_adj.keys())[0]][3] + list(seg_adj.keys())[0][3]
    
    
    return lin_df

# output
def output(lin_df, summary, data, pymol, uf, psf, xtcs):
    ### csv summary
    minr_df = {}
    for l in lin_df.keys():
        if type(l) is int:
            minr = [None, -1, 1000]
            for a in lin_df[l].keys():
                if a == "arg" or a == "leaflet" or a == "fen" or a == "first":
                    continue
                for f,r in lin_df[l][a].items():
                    if r < minr[2]:
                        minr = [a, f, r]
            minr_df[l] = minr
    s = open(summary, "w")
    s.write("resid,fenestration,leaflet,first universe,first frame,universe,frame,name,min r\n")
    for l in minr_df.keys():
        if minr_df[l][1] != -1:
            s.write(f"{l},{lin_df[l]['fen'][minr_df[l][1]]},{lin_df[l]['leaflet']},{uf[lin_df[l]['first']][0]+1},{uf[lin_df[l]['first']][1]},{uf[minr_df[l][1]][0]+1},{uf[minr_df[l][1]][1]},{minr_df[l][0]},{minr_df[l][2]}\n")
    s.close()
    
    ### csv data
    d = open(data, "w")
    d.write("resid,atom,universe,frame,radius,fenestration,leaflet\n")
    for l in lin_df.keys():
        if type(l) is int:
            if lin_df[l]["leaflet"] is None:
                continue
            else:
                for a in lin_df[l].keys():
                    if a == "arg" or a == "leaflet" or a == "fen" or a == "first":
                        continue
                    if lin_df[l][a] == {}:
                        continue
                    else:
                        for f,r in lin_df[l][a].items():
                            d.write(f"{l},{a},{uf[f][0]+1},{uf[f][1]},{r},{lin_df[l]['fen'][f]},{lin_df[l]['leaflet']}\n")
    d.close()
    
    ### pymol
    pm = open(pymol, "w")
    pm.write("from pymol import cmd\nimport os\n\n")
    fst = None
    for xtc in xtcs:
        if "/" in xtc:
            name = str(xtc.split("/")[-1])
        else:
            name = str(xtc)
        name = name.split("_")[0]
        if fst is None:
            fst = name
        pm.write(f"cmd.load('{psf}', '{name}')\n")
        pm.write(f"cmd.load_traj('{xtc}', '{name}')\n")
    pm.write("cmd.select('protein', 'polymer.protein')\n")
    resi = []
    for l in minr_df.keys():
        if minr_df[l][1] != -1:
            if l not in resi:
                resi.append(l)
    resi_sorted = sorted(resi)
    for r in resi_sorted:
        pm.write(f"cmd.select('l{r}', selection='resn POPC and resi {r}')\n")
    pm.write("cmd.hide('everything')\n")
    pm.write("cmd.show('cartoon', 'polymer.protein')\n")
    pm.write("cmd.show('sticks', 'l*')\n")
    pm.write("cmd.dss()\n")
    pm.write("cmd.color('white', 'resn POPC and name C*')\n")
    pm.write("cmd.spectrum('segi', selection='polymer.protein and name C*')\n")
    pm.write("cmd.disable('*')\n")
    pm.write(f"cmd.enable('{fst}')\n")
    pm.close()
    
    return 0

# parser
def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("--deffnm", type=str, default=None, help="Default file name for input and output files which are not specified. Input a file name without extention. If not specified, XTC file name would be the default file name. default: None")
    parser.add_argument("--psf", type=str, default=None, help="Input PSF file. Specify a file name. default: None")
    parser.add_argument("--xtc", type=str, default=None, help="Input XTC file. Specify a file name or a series name. If file name contains '\*', this program works for all serial files. default: None")
    parser.add_argument("--start", type=int, default=0, help="Start frame. default: 0")
    parser.add_argument("--stop", type=int, default=-1, help="Stop frame. default: -1")
    parser.add_argument("--step", type=int, default=1, help="Step frame. default: 1")
    parser.add_argument("--r1", type=float, required=True, help="The cutoff radius for analysis. required.")
    parser.add_argument("--r2", type=float, default=None, help="The border radius for pore entrance. If not specified, the r1 would be used.")
    parser.add_argument("--atom", type=str, nargs="*", default=None, help="Atom names to analyse. Input 'all' to analyse all atoms. default: C13, C14, C15, C218, C316")
    parser.add_argument("--zcen", type=float, default=None, help="Input border z-coordinate of lipids' leaflets. default: None")
    parser.add_argument("--upmax", type=int, default=-1, help="The number of lipids in upper leaflet. default: -1")
    parser.add_argument("--csv", type=str, default=None, help="Output csv file. Specify a file name or use default file name (--deffnm). default: None")
    args = parser.parse_args()
    
    return args

# main
def main():
    input_args = vars(parse_args())
    args = input_check(**input_args)
    analysis_args = {"us": args["us"], "frames": args["frames"], "r1": args["r1"], "r2": args["r2"], "atom": args["atom"], "nl": args["nl"], "upmax": args["upmax"]}
    output_args = {"lin_df": analysis(**analysis_args), "summary": args["summary"], "data": args["data"], "pymol": args["pymol"], "uf": args["uf"], "psf": args["psf"], "xtcs": args["xtcs"]}
    output(**output_args)
    
    return 0

if __name__ == "__main__":
    main()