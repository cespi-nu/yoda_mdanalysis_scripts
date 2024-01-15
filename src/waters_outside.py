#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse


def waters_outside(tpr, trj, opng=False, ocsv=False, time_unit='ps', title=False, start_step=0, end_step=False, stride=1):
    # Load the trajectory
    u = mda.Universe(tpr, trj)
    
    # prepare list for trajectory iteration
    waters_outside = []
    times = []
    frame = []
    
    # trajectory iteration
    # times
    for ts in u.trajectory:
        times.append(u.trajectory.time)
        frame.append(u.trajectory.frame)
    if end_step > len(times):
        end_step = -1
    times = np.array(times)
    frame = np.array(frame)
    if end_step:
        times_sliced = times[start_step:end_step:stride]
        frame_sliced = frame[start_step:end_step:stride]
    else:
        times_sliced = times[start_step::stride]
        frame_sliced = frame[start_step::stride]
    if time_unit == "ns":
        times_sliced = times_sliced / 1000
    
    for i in frame_sliced:
        # xy center of the pore
        u.trajectory[i]
        acc_x = 0
        acc_y = 0
        ca_44 = u.select_atoms("resid 44 and name CA")
        for x in ca_44.positions[:,0]:
            acc_x = acc_x + x
        vcen_x = acc_x/7
        for y in ca_44.positions[:,1]:
            acc_y = acc_y + y
        vcen_y = acc_y/7

        # z center of the membrane
        acc_z = 0
        lipid_p = u.select_atoms("resname POPC and name P")
        for z in lipid_p.positions[:,2]:
            acc_z = acc_z + z
        vcen_z = acc_z/lipid_p.n_atoms 

        # lipid-liquid interface
        u_lipid_p = u.select_atoms(f"resname POPC and name P and prop z > {vcen_z}")
        l_lipid_p = u.select_atoms(f"resname POPC and name P and prop z < {vcen_z}")
        acc_u = 0
        for z in u_lipid_p.positions[:,2]:
            acc_u = acc_u + z
        vup_z = acc_u/u_lipid_p.n_atoms - 5
        acc_l = 0
        for z in l_lipid_p.positions[:,2]:
            acc_l = acc_l + z
        vlow_z = acc_l/l_lipid_p.n_atoms + 5

        # prepare for water molecule selection
        waters_at_membrane_hight = u.select_atoms(f"resname TIP3 and name OH2 and segid TIP3 and prop z > {vlow_z-5} and prop z < {vup_z+5}")
        waters_outside_frame = []

        # water molecule selection
        for w in waters_at_membrane_hight:
            x = w.position[0]
            y = w.position[1]
            if (x-vcen_x)**2 + (y-vcen_y)**2 > 900:
                waters_outside_frame.append([w.resid, w.segment])
        waters_outside.append(waters_outside_frame)
    
    # waters number
    waters_number = []
    for i in waters_outside:
        waters_number.append(len(i))
    
    # set output file base name
    outfile_base = os.path.splitext(trj)[0]
    
    # plot
    plt.plot(times_sliced, waters_number)
    figname = outfile_base
    if title:
        figname = title
    plt.title(figname)
    pngout = outfile_base + "_waters_outside.png"
    if opng:
        pngout = opng
    plt.savefig(pngout)
    
    # output csv
    csvout = outfile_base + "_waters_outside.csv"
    if ocsv:
        csvout = ocsv
    with open(csvout, "w") as f:
        i = 0
        for frame in waters_outside:
            f.write(f"Frame {frame_sliced[i]}\n")
            for w in waters_outside[i]:
                f.write(f"{w[0]},{str(w[1])[9:-1]}\n")
            i = i + 1

# parse arguments
def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("--trajfile", type=str, required=True, help="Input trajectory file. required.")
    parser.add_argument("--tprfile", type=str, default="step5_input.psf", help="Input topology file. default: step5_input.psf")
    parser.add_argument("--outpng", type=str, default=False, help="Output .png file.")
    parser.add_argument("--outcsv", type=str, default=False, help="Output .csv file.")
    parser.add_argument("--title", type=str, default=False, help="figure title.")
    parser.add_argument("--time_unit", type=str, choices=["ps", "ns"], default="ps", help="Time unit of output plot. ps or ns.")
    parser.add_argument("--start_step", type=int, default=0, help="Start frame number.")
    parser.add_argument("--end_step", type=int, default=False, help="End frame number.")
    parser.add_argument("--stride", type=int, default=1, help="Stride of frame number.")
    args = parser.parse_args()
    return args

# main function
def main():
    args = parse_args()
    waters_outside(
        tpr=args.tprfile,
        trj=args.trajfile,
        opng=args.outpng,
        ocsv=args.outcsv,
        time_unit=args.time_unit,
        title=args.title,
        start_step=args.start_step,
        end_step=args.end_step,
        stride=args.stride
        )

if __name__ == "__main__":
    main()
