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

# util
def lu_decomp(mat):
    n = len(mat)
    mat_l = mat.copy()
    mat_u = mat.copy()
    i = 1
    while i < n:
        l = mat_l[i:n,i-1]
        ut = mat_l[i-1,i:n] / mat_l[i-1,i-1]
        mat_la = mat_l[i:n,i:n]
        mat_lut = np.dot(l,ut)
        mat_l[i:n,i:n] = mat_la - mat_lut
        mat_l[i-1,i:n] = 0
        mat_u[i-1,i-1] = 1
        mat_u[i:n,i-1] = 0
        mat_u[i-1,i:n] = ut
        i += 1
    mat_u[n-1,0:n-1] = 0
    mat_u[n-1,n-1] = 1
    
    return (mat_l, mat_u)

def rgrs_plane(points):
    ### data load
    points_data = {}
    points_data["n"] = len(points)
    points_data["x"] = sum([p[0] for p in points])
    points_data["y"] = sum([p[1] for p in points])
    points_data["z"] = sum([p[2] for p in points])
    points_data["x2"] = sum([p[0] ** 2 for p in points])
    points_data["y2"] = sum([p[1] ** 2 for p in points])
    points_data["xy"] = sum([p[0] * p[1] for p in points])
    points_data["xz"] = sum([p[0] * p[2] for p in points])
    points_data["yz"] = sum([p[1] * p[2] for p in points])
    
    ### matrix
    points_mat_a = np.matrix([[points_data["n"], points_data["x"], points_data["y"]],
                              [points_data["x"], points_data["x2"], points_data["xy"]],
                              [points_data["y"], points_data["xy"], points_data["y2"]]])
    points_mat_b = np.matrix([[points_data["z"]],
                              [points_data["xz"]],
                              [points_data["yz"]]])
    ##### LU decomposition
    points_mat_l, points_mat_u = lu_decomp(points_mat_a)
    ##### solve y
    points_mat_y = np.matrix([[0],[0],[0]], dtype=float)
    for i in range(len(points_mat_l)):
        points_mat_y[i,0] = (points_mat_b[i,0] - sum(points_mat_l[i,0:i]*points_mat_y[0:i,0])) / points_mat_l[i,i]
    ##### solve x
    points_mat_x = np.matrix([[0],[0],[0]], dtype=float)
    for i in range(len(points_mat_u)-1,-1,-1):
        points_mat_x[i,0] = points_mat_y[i,0] - sum(points_mat_u[i,i+1:]*points_mat_x[i+1:,0])
    
    return points_mat_x

def diff_plane(point, plane):
    return point[2] - (plane[0,0] + plane[1,0] * point[0] + plane[2,0] * point[1])

def resids_selection(resids):
    selection = "protein and ("
    for r in resids:
        selection += f"resid {r} or "
    return selection[:-4] + ")"

# input check
def input_check(*, deffnm, psf, xtc, start, stop, step, scatom, selection, csv, pdf):
    ### files
    if deffnm is None:
        if psf is None:
            print("Error: Specify a PSF file name or use default file name.")
            sys.exit(1)
        if xtc is None:
            print("Error: Specify a XTC file name or use default file name.")
            sys.exit(1)
        if csv is None:
            csv = f"{xtc[:-4]}.csv"
        if pdf is None:
            pdf = f"{xtc[:-4]}.pdf"
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
    
    ### make universe
    us = []
    for x in xtcs:
        us.append(mda.Universe(psf, x))
    
    ### options
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
    
    ### scatom
    sca = {}
    if scatom is None:
        scatom = "sc_atom.csv"
    if os.path.isfile(scatom):
        continue
    else:
        s = open(scatom, "w")
        s.write("amino acid,atoms\n")
        s.write("ALA,CB\n")
        s.write("CYS,SG\n")
        s.write("ASP,OD1,OD2\n")
        s.write("GLU,OE1,OE2\n")
        s.write("PHE,CZ\n")
        s.write("GLY,CA\n")
        s.write("HIS,CE1,NE2\n")
        s.write("HSD,CE1,NE2\n")
        s.write("HSE,CE1,NE2\n")
        s.write("ILE,CD\n")
        s.write("LYS,NZ\n")
        s.write("LEU,CD1,CD2\n")
        s.write("MET,CE\n")
        s.write("ASN,OE1,NE2\n")
        s.write("PRO,CG\n")
        s.write("GLN,OE1,NE2\n")
        s.write("ARG,NH1,NH2\n")
        s.write("SER,OG\n")
        s.write("THR,OG1,CG2\n")
        s.write("VAL,CG1,CG2\n")
        s.write("TRP,NE1,CZ2,CH2,CZ3\n")
        s.write("TYR,OH\n")
        s.close()
    s = open(scatom, "r")
    for line in s:
        if re.match("amino acid", line):
            continue
        else:
            l = line.split(",")
            sca[l[0]] = l[1:]
    s.close()
    
    ### residue selection
    sel = []
    s = open(selection, "r")
    for line in s:
        sel.append(line.split(","))
    s.close()
    
    ### check inputs
    print("Input files")
    print(f"PSF file: {psf}")
    print(f"XTC file: {xtcs}")
    print(f"SC atom file: {scatom}")
    print("Output files")
    print(f"CSV file: {csv}")
    print(f"PDF file: {pdf}")
    print("Options")
    for i in range(len(us)):
        print(f"Frames(universe {i}): {frms[i][0]}-{frms[i][1]} (step: {frms[i][2]}) {len(frames[i])} frames")
    print(f"Residues A: {sel[0]}")
    print(f"residue B: {sel[1]}")
    
    
    ### output args
    args = {"us": us, "frames": frames, "scatom": sca, "selection": sel, "csv": csv, "pdf": pdf}
    
    return args

# analysis
def fen_sc(*, us, frames, scatom, selection):
    ### sc_atom
    sc_atom = {}
    for rn, an in scatom.items():
        sc_atom[rn] = an
    
    ### 
    fen_df = {}
    frame = 0
    for i in range(len(us)):
        fen_df[i] = {}
        as_a = us[i].select_atoms(resids_selection(selection[0]))
        as_b = us[i].select_atoms(resids_selection(selection[1]))
        for f in frames[i]:
            us[i].trajectory[f]
            ### subunit A
            a_rgrs = rgrs_plane(as_a.positions)
            a_x, a_y, a_z = as_a.positions
            a_xmin = min(a_x)
            a_xmax = max(a_x)
            a_ymin = min(a_y)
            a_ymax = max(a_y)
            a_zmin = min(a_z)
            a_zmax = max(a_z)
            ### subunit B
            b_rgrs = rgrs_plane(as_b.positions)
            b_x, b_y, b_z = as_b.positions
            b_xmin = min(b_x)
            b_xmax = max(b_x)
            b_ymin = min(b_y)
            b_ymax = max(b_y)
            b_zmin = min(b_z)
            b_zmax = max(b_z)
            ### dataframe
            fen_df[i][f] = {"a_rg": a_rgrs, "a_stats": [a_xmin, a_xmax, a_ymin, a_ymax, a_zmin, a_zmax], "b_rg": b_rgrs, "b_stats": [b_xmin, b_xmax, b_ymin, b_ymax, b_zmin, b_zmax]}
            
    
    return 

# output
def output(df, frames, selection, csv, pdf):
    ### csv
    c = open(csv, "w")
    c.write("universe,frame,leaflet,aprx plane,average diff,std diff\n")
    for i in range(len(df.keys())):
        for f in df[i].keys():
            for l in df[i][f].keys():
                c.write(f"{i},{f},{l},z = {df[i][f][l]['plane_constants'][0,0]} + {df[i][f][l]['plane_constants'][1,0]}x + {df[i][f][l]['plane_constants'][2,0]},{df[i][f][l]['avrg_diff']},{df[i][f][l]['std_diff']}\n")
    
    
    ### pdf
    ##### make dataframe for plot
    plot3d_data = {}
    for i in df.keys():
        for f in df[i].keys():
            plot3d_data[f] = {}
            for l in df[i][f].keys():
                plot3d_data[f][l] = {"x": [], "y": [], "z": [], "plane": df[i][f][l]["plane_constants"], "avrg_diff": df[i][f][l]["avrg_diff"], "std_diff": df[i][f][l]["std_diff"]}
                for p in df[i][f][l]["positions"]:
                    plot3d_data[f][l]["x"].append(p[0])
                    plot3d_data[f][l]["y"].append(p[1])
                    plot3d_data[f][l]["z"].append(p[2])
    ##### plot
    p = PdfPages(pdf)
    lcolor = {"upper": "blue", "lower": "red"}
    lmarker = {"upper": "o", "lower": "x"}
    sumloc = {"upper": 0.05, "lower": -0.05}
    for f in plot3d_data.keys():
        fig = plt.figure()
        ax = fig.add_subplot(121, projection="3d")
        for l in plot3d_data[f].keys():
            fig.suptitle(f"Frame {f}")
            ax.scatter(plot3d_data[f][l]["x"], plot3d_data[f][l]["y"], plot3d_data[f][l]["z"], c=lcolor[l], marker=lmarker[l], label=f"{l} leaflet")
            x = np.linspace(min(plot3d_data[f][l]["x"]), max(plot3d_data[f][l]["x"]), 100)
            y = np.linspace(min(plot3d_data[f][l]["y"]), max(plot3d_data[f][l]["y"]), 100)
            x, y = np.meshgrid(x, y)
            z = plot3d_data[f][l]["plane"][0,0] + plot3d_data[f][l]["plane"][1,0] * x + plot3d_data[f][l]["plane"][2,0] * y
            ax.plot_surface(x, y, z, alpha=0.5, color=lcolor[l])
            ax.text2D(0.12, sumloc[l], f"{l} leaflet\nregression plane:\nz = {plot3d_data[f][l]['plane'][0,0]:.2f} + {plot3d_data[f][l]['plane'][1,0]:.3f}x + {plot3d_data[f][l]['plane'][2,0]:.3f}y\ndiff avrg: {plot3d_data[f][l]['avrg_diff']:.3f}\ndiff std: {plot3d_data[f][l]['std_diff']:.3f}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.legend(loc="upper right")
        p.savefig()
        plt.close()
    p.close()
    
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
    parser.add_argument("--selection", type=str, required=True, help="Input CSV file containing the resids of protein for analysis. The first/second row is for the selection of A/B subunit in fenestration between subunit A and B.")
    parser.add_argument("--scatom", type=str, default=None, help="Input CSV file containing the name of atoms for in analysis. If not specified, a file with default value will be generated and used. default: None")
    parser.add_argument("--csv", type=str, default=None, help="Output CSV file. Specify a file name or use default file name (--deffnm). default: None")
    parser.add_argument("--pdf", type=str, default=None, help="Output PDF file. Specify a file name or use default file name (--deffnm). default: None")
    args = parser.parse_args()
    
    return args

# main
def main():
    input_args = vars(parse_args())
    args = input_check(**input_args)
    analysis_args = {"us": args["us"], "frames": args["frames"], "selection": args["selection"]}
    output_args = {"df": lipids_surface(**analysis_args), "frames": args["frames"], "selection": args["selection"], "csv": args["csv"], "pdf": args["pdf"]}
    output(**output_args)
    
    return 0

if __name__ == "__main__":
    main()