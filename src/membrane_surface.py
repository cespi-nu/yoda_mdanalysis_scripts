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
def f2uf(us, frames):
    uf = {}
    frame = 0
    for i in range(len(us)):
        for f in frames[i]:
            uf[frame] = (i,f)
            frame += 1
    
    return uf

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

# input check
def input_check(*, deffnm, psf, xtc, start, stop, step, upmax, csv, pdf):
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
    
    selection = {}
    u_lipid_p = us[0].select_atoms("resname POPC and name P")
    if upmax == -1:
        ucen = u_lipid_p.center_of_geometry()
        u_lipid_upper = us[0].select_atoms(f"resname POPC and name P and prop z >= {ucen[2]}")
        selection["upper"] = (u_lipid_upper.n_residues)
    else:
        selection["upper"] = (upmax)
    selection["lower"] = (u_lipid_p.n_residues)
    
    ### uf
    uf = f2uf(us, frames)
    
    ### check inputs
    print("Input files")
    print(f"PSF file: {psf}")
    print(f"XTC file: {xtcs}")
    print("Output files")
    print(f"CSV file: {csv}")
    print(f"PDF file: {pdf}")
    print("Options")
    for i in range(len(us)):
        print(f"Frames(universe {i}): {frms[i][0]}-{frms[i][1]} (step: {frms[i][2]}) {len(frames[i])} frames")
    print(f"Lipids in upper leaflet: {selection['upper']}")
    print(f"Lipids in lower leaflet: {selection['lower'] - selection['upper']}")
    
    
    ### output args
    args = {"us": us, "frames": frames, "selection": selection, "csv": csv, "pdf": pdf, "uf": uf}
    
    return args

# analysis
def lipids_surface(*, us, frames, selection):
    ### data collection
    df = {}
    frame = 0
    for i in range(len(us)):
        df[i] = {}
        u_lipid_upper = us[i].select_atoms(f"resname POPC and name P and resid 1:{selection['upper']}")
        u_lipid_lower = us[i].select_atoms(f"resname POPC and name P and resid {selection['upper']+1}:{selection['lower']}")
        usels = {"upper": u_lipid_upper, "lower": u_lipid_lower}
        for f in range(len(frames[i])):
            df[i][frame] = {"upper": {"plane_constants": None, "avrg_diff": None, "std_diff": None, "positions": []}, "lower": {"plane_constants": None, "avrg_diff": None, "std_diff": None, "positions": []}}
            us[i].trajectory[frames[i][f]]
            for l in selection.keys():
                df[i][frame][l]["positions"] = usels[l].positions
                df[i][frame][l]["plane_constants"] = (rgrs_plane(usels[l].positions))
                ##### data analysis
                diff = []
                for r in range(usels[l].n_residues):
                    d = diff_plane(df[i][frame][l]["positions"][r], df[i][frame][l]["plane_constants"])
                    diff.append(d)
                df[i][frame][l]["avrg_diff"] = (np.mean(diff))
                df[i][frame][l]["std_diff"] = (np.std(diff))
            
            frame += 1
    
    return df

# output
def output(df, uf, selection, csv, pdf):
    ### csv
    c = open(csv, "w")
    c.write("universe,universe,frame,leaflet,reg a,reg b,reg c,average diff,std diff\n")
    for i in range(len(df.keys())):
        for f in df[i].keys():
            for l in df[i][f].keys():
                c.write(f"{i},{uf[f][0]},{uf[f][1]},{l},{df[i][f][l]['plane_constants'][0,0]},{df[i][f][l]['plane_constants'][1,0]},{df[i][f][l]['plane_constants'][2,0]},{df[i][f][l]['avrg_diff']},{df[i][f][l]['std_diff']}\n")
    
    
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
            ax.text2D(0.15, sumloc[l], f"{l} leaflet\nregression plane:\nz = {plot3d_data[f][l]['plane'][0,0]:.2f} + {plot3d_data[f][l]['plane'][1,0]:.3f}x + {plot3d_data[f][l]['plane'][2,0]:.3f}y\ndiff avrg: {plot3d_data[f][l]['avrg_diff']:.3f}\ndiff std: {plot3d_data[f][l]['std_diff']:.3f}")
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
    parser.add_argument("--upmax", type=int, default=-1, help="The max residue number of upper leaflet. If not specified, This program will automatically classify by z-coordinate, could be inaccurate. default: -1")
    parser.add_argument("--csv", type=str, default=None, help="Output csv file. Specify a file name or use default file name (--deffnm). default: None")
    parser.add_argument("--pdf", type=str, default=None, help="Output pdf file. Specify a file name or use default file name (--deffnm). default: None")
    args = parser.parse_args()
    
    return args

# main
def main():
    input_args = vars(parse_args())
    args = input_check(**input_args)
    analysis_args = {"us": args["us"], "frames": args["frames"], "selection": args["selection"]}
    output_args = {"df": lipids_surface(**analysis_args), "uf": args["uf"], "selection": args["selection"], "csv": args["csv"], "pdf": args["pdf"]}
    output(**output_args)
    
    return 0

if __name__ == "__main__":
    main()