#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os
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
def input_check(*, deffnm, psf, xtc, start, stop, step, top, bot, csv, pdf):
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
    
    top_sel = f"protein and (resid {top.pop(0)}"
    for r in top:
        top_sel += f" or resid {r}"
    top_sel += ")"
    bot_sel = f"protein and (resid {bot.pop(0)}"
    for r in bot:
        bot_sel += f" or resid {r}"
    bot_sel += ")"
    
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
    
    ### output args
    args = {"us": us, "frames": frames, "top": top_sel, "bot": bot_sel, "csv": csv, "pdf": pdf, "uf": uf}
    
    return args

# analysis
def analysis(*, us, frames, top, bot):
    tilt_df = {}
    frame = 0
    
    for i in range(len(us)):
        utop = us[i].select_atoms(top)
        ubot = us[i].select_atoms(bot)
        for f in frames[i]:
            us[i].trajectory[f]
            top_cen = utop.center_of_geometry()
            bot_cen = ubot.center_of_geometry()
            axis_vect = top_cen - bot_cen
            axis_tilt = np.arccos(axis_vect[2] / np.linalg.norm(axis_vect)) * 180 / np.pi
            tilt_direction = np.arctan2(axis_vect[1], axis_vect[0]) * 180 / np.pi
            if axis_vect[1] < 0:
                tilt_direction = - tilt_direction
            tilt_df[frame] = {"vect": axis_vect, "tilt": axis_tilt, "direction": tilt_direction}
            frame += 1
    
    return tilt_df

# output
def output(*, tilt_df, csv, pdf, uf):
    ### csv
    c = open(csv, "w")
    c.write("universe,frame,axis_vect_x,axis_vect_y,axis_vect_z,axis_tilt,tilt_direction\n")
    for f in tilt_df.keys():
        c.write(f"{uf[f][0]},{uf[f][1]},{tilt_df[f]['vect'][0]},{tilt_df[f]['vect'][1]},{tilt_df[f]['vect'][2]},{tilt_df[f]['tilt']},{tilt_df[f]['direction']}\n")
    c.close()
    
    ### pdfabout:blank#blocked
    tilt = []
    for f in tilt_df.keys():
        tilt.append(tilt_df[f]["tilt"])
    direction = []
    for f in tilt_df.keys():
        direction.append(tilt_df[f]["direction"])
    p = PdfPages(pdf)
    fig, ax = plt.subplots(ncols=1, nrows=2, tight_layout=True)
    ax[0].plot(list(tilt_df.keys()), tilt)
    ax[0].set_xlabel("Frame")
    ax[0].set_ylabel("Tilt angle (degree)")
    ax[0].set_ylim([0, 180])
    ax[0].set_yticks([0, 30, 60, 90, 120, 150, 180])
    ax[0].set_title("Pore axis tilt angle")
    ax[1].plot(list(tilt_df.keys()), direction)
    ax[1].set_xlabel("Frame")
    ax[1].set_ylabel("Tilt direction (degree)")
    ax[1].set_ylim([-180, 180])
    ax[1].set_yticks([-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180])
    ax[1].set_title("Pore axis tilt direction")
    p.savefig()
    plt.close()
    p.close()
    
    return 0

# parser
def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("--deffnm", type=str, default=None, help="Default file name for input and output files which are not specified. Input a file name without extention. If not specified, PDB file name would be the default file name. default: None")
    parser.add_argument("--psf", type=str, default=None, help="Input PSF file. Specify a file name or use default file name. default: None")
    parser.add_argument("--xtc", type=str, default=None, help="Input XTC file. Specify a file name or use default file name. default: None")
    parser.add_argument("--start", type=int, default=0, help="Start frame. default: 0")
    parser.add_argument("--stop", type=int, default=-1, help="Stop frame. default: -1")
    parser.add_argument("--step", type=int, default=1, help="Step frame. default: 1")
    parser.add_argument("--top", type=str, nargs='+', help="Input the residues for top-selection with a blank. required.")
    parser.add_argument("--bot", type=str, nargs='+', help="Input the residues for bottom-selection with a blank. required.")
    parser.add_argument("--csv", type=str, default=None, help="Output csv file. Specify a file name or set this to True to generate csv output. default: None")
    parser.add_argument("--pdf", type=str, default=None, help="Output pdf file. Specify a file name or set this to True to generate pdf output. default: None")
    args = parser.parse_args()
    
    return args

# main
def main():
    input_args = vars(parse_args())
    args = input_check(**input_args)
    analysis_args = {"us": args["us"], "frames": args["frames"], "top": args["top"], "bot": args["bot"]}
    output_args = {"tilt_df": analysis(**analysis_args), "csv": args["csv"], "pdf": args["pdf"], "uf": args["uf"]}
    output(**output_args)
    
    return 0

if __name__ == "__main__":
    main()