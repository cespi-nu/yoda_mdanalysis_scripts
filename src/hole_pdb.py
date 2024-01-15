#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
from MDAnalysis.analysis import hole2
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os
import argparse

# check input parameters
def input_check(*, deffnm=None, dot_den=15, n=1, csv=False, pdf=False, vmd=False, pdb=None, out=None, sph=None, vdwr=None, end_r=22, sample=0.2, cp=None, cv=[0,0,1], ignore_res=None, random=None, keep=False):
    ### make dictionary of arguments
    args = locals()
    deffnm = args.pop("deffnm")
    
    ### check if vmd output is required
    verify_dict = {"Y":True, "y":True, "N":False, "n":False, "":True, "yes":True, "no":False, "YES":True, "NO":False}
    if args["vmd"]:
        if args["keep"] is False:
            print("Warning: VMD output requires keep_files option to be True.\n")
            verify = input("Continue by setting keep_files to True? (Y/n): ")
            if verify_dict[verify]:
                print("keep_files is set to True.\n")
                args["keep"] = True
            else:
                sys.exit("Exiting...")
    
    ### set input / output file names and options
    ##### input / output file names
    if not deffnm:
        if args["pdb"] is None:
            sys.exit("Error: No PDB file is specified.")
        else:
            if args["csv"] == "True":
                args["csv"] = args["pdb"][:-4] + ".csv"
            if args["pdf"] == "True":
                args["pdf"] = args["pdb"][:-4] + ".pdf"
            if args["out"] is None or args["out"] == "True": 
                args["out"] = args["pdb"][:-4] + ".out"
            if args["sph"] is None or args["sph"] == "True":
                args["sph"] = args["pdb"][:-4] + ".sph"
            if args["keep"]:
                if args["vmd"] == "True":
                    args["vmd"] = args["pdb"][:-4] + ".vmd"
    else:
        if args["pdb"] is None:
            args["pdb"] = deffnm + ".pdb"
        if args["csv"] == "True":
            args["csv"] = deffnm + ".csv"
        if args["pdf"] == "True":
            args["pdf"] = deffnm + ".pdf"
        if args["out"] is None or args["out"] == "True": 
            args["out"] = deffnm + ".out"
        if args["sph"] is None or args["sph"] == "True":
            args["sph"] = deffnm + ".sph"
        if args["keep"]:
            if args["vmd"] == "True":
                args["vmd"] = deffnm + ".vmd"
    ##### options
    default_ignore_res = ["SOL", "WAT", "TIP", "HOH", "K  ", "NA ", "CL "]
    if args["ignore_res"] is None:
        args["ignore_res"] = default_ignore_res
    
    ### verify input / output files names and options
    ##### input / output file names
    print("Input PDB file: " + str(args["pdb"]) + "\n")
    print("Hole output csv: " + str(args["csv"]) + "\n")
    print("Hole output pdf: " + str(args["pdf"]) + "\n")
    print("keep_files option: " + str(args["keep"]) + "\n")
    if args["keep"]:
        print("Output files are going to be kept.\nOutput OUT file: " + str(args["out"]) + "\nOutput SPH file: " + str(args["sph"]) + "\n") # actual argument name for hole is "logfile" and "sphpdb"
    else:
        print("Output log file (*.out) and SPH-style PDB file (*.sph) are not going to be kept.\n")
    if args["vmd"]:
        print("The VMD script file for VMD visualization will be generated.\nOutput VMD script file: " + str(args["vmd"]) + "\n")
    else:
        print("The VMD script file for VMD visualization will not be generated.\n")
    print("This program will run hole analysis " + str(args["n"]) + " time(s).\n")
    ##### options
    options = {}
    options["vdwradii"] = args["vdwr"] # actual argument name for hole is "radius"
    options["end_radius"] = args["end_r"] # actual argument name for hole is "endrad"
    options["cpoint"] = args["cp"]
    options["cvect"] = args["cv"]
    options["sample"] = args["sample"]
    options["ignore_residues"] = args["ignore_res"]
    options["random"] = args["random"]
    print("Options for hole analysis;\n")
    for val in options.keys():
        print(val + ": " + str(options[val]) + "\n")
    if args["vmd"]:
        print("Options for vmd visualization;\n")
        print("dot_density: " + str(args["dot_den"]) + "\n")
    
    return args

# hole analysis on a pdb
def hole_analysis(*, n, dot_den, csv, pdf, vmd, pdb, out, sph, vdwr, end_r, sample, cp, cv, ignore_res, random, keep):
    ### setting arguments for hole analysis and output
    args = locals()
    n = args.pop("n")
    dot_den = args.pop("dot_den")
    csv = args.pop("csv")
    pdf = args.pop("pdf")
    vmd = args.pop("vmd")
    pdb = args.pop("pdb")
    hole_args = ["outfile", "sphpdb_file", "vdwradii_file", "end_radius", "sample", "cpoint", "cvect", "ignore_residues", "random_seed", "keep_files"]
    args_vals = args.values()
    hole_options = dict(zip(hole_args, args_vals))
    
    ### run hole analysis
    holes = []
    if n == 1:
        h = hole2.hole(pdb, **hole_options)
        holes.append(h)
    else:
        out_prefix = args["out"][:-4]
        sph_prefix = args["sph"][:-4]
        for i in range(n):
            hole_options["outfile"] = str(out_prefix) + "_" + str(i+1) + ".out"
            hole_options["sphpdb_file"] = str(sph_prefix) + "_" + str(i+1) + ".sph"
            h = hole2.hole(pdb, **hole_options)
            holes.append(h)
    
    ### results for output
    output_args = {"holes": holes, "sph": args["sph"], "csv": csv, "pdf": pdf, "vmd": vmd, "n": n, "dot_den": dot_den}
    return output_args

def hole_output(*, holes, sph, csv, pdf, vmd, n, dot_den):
    ### making a data frame of hole results
    results=[]
    for i in range(n):
        pore_axis = holes[i][0].rxn_coord
        pore_radius = holes[i][0].radius
        result_n = (i, pore_axis, pore_radius)
        results.append(result_n)
    
    ### csv output
    if csv:
        c = open(csv, "w")
        for i in range(n):
            c.write("run," + str(i+1) + "\n")
            c.write("pore axis,pore_radius\n")
            for j in range(len(results[i][1])):
                c.write(str(results[i][1][j]) + "," + str(results[i][2][j]) + "\n")
        c.close()
    
    ### pdf output
    if pdf:
        if n == 1:
            p = PdfPages(pdf)
            plt.plot(results[0][2], results[0][1])
            plt.xlabel("pore radius (A)")
            plt.ylabel("pore axis (A)")
            p.savefig()
        else:
            p = PdfPages(pdf)
            for i in range(n):
                if i % 4 == 0:
                    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(9, 16))
                    for ax in axs.flat:
                        ax.set(xlabel="pore radius (A)", ylabel="pore axis (A)")
                    axs[0, 0].plot(results[i][2], results[i][1])
                if i % 4 == 1:
                    axs[0, 1].plot(results[i][2], results[i][1])
                if i % 4 == 2:
                    axs[1, 0].plot(results[i][2], results[i][1])
                if i % 4 == 3:
                    axs[1, 1].plot(results[i][2], results[i][1])
                    p.savefig()
            if n % 4 != 0:
                p.savefig()
        p.close()
    
    ### vmd config file output
    if vmd:
        if n == 1:
            hole2.create_vmd_surface(filename=vmd, sphpdb=sph, dot_density=dot_den)
        else:
            sph_prefix = sph[:-4]
            vmd_prefix = vmd[:-4]
            for i in range(n):
                sph = sph_prefix + "_" + str(i+1) + ".sph"
                vmd = vmd_prefix + "_" + str(i+1) + ".vmd"
                hole2.create_vmd_surface(filename=vmd, sphpdb=sph, dot_density=dot_den)

# parse arguments
def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("--deffnm", type=str, default=None, help="Default file name for input and output files which are not specified. Input a file name without extention. If not specified, PDB file name would be the default file name. default: None")
    parser.add_argument("--dot_den", type=int, default=15, help="Dot density for VMD visualization. Setting to 5-35 is recommended. default: 15")
    parser.add_argument("--n", type=int, default=1, help="Number of times to run hole analysis. default: 1")
    parser.add_argument("--csv", type=str, default=False, help="Output csv file. Specify a file name or set this to True to generate csv output. default: False")
    parser.add_argument("--pdf", type=str, default=False, help="Output pdf file. Specify a file name or set this to True to generate pdf output. default: False")
    parser.add_argument("--vmd", type=str, default=False, help="Output vmd config file. Specify a file name or set this to True to generate vmd config file. default: False")
    parser.add_argument("--pdb", type=str, default=None, help="Input PDB file. Specify a file name or use default file name. default: None")
    parser.add_argument("--out", type=str, default=None, help="Output log (*.out) file. Specify a file name or use default file name. default: None")
    parser.add_argument("--sph", type=str, default=None, help="Output SPH-style PDB (*.sph) file. Specify a file name or use default file name. default: None")
    parser.add_argument("--vdwr", type=str, default=None, help="Input file of van der Waals radii for each atom. If not specified, the 'simple2.rad' file would be used. default: None")
    parser.add_argument("--end_r", type=float, default=22, help="Set the radius [A] to end hole analysis. default: 22")
    parser.add_argument("--sample", type=float, default=0.2, help="Set the sampling distance [A] for the hole analysis. default: 0.2")
    parser.add_argument("--cp", type=float, nargs=3, default=None, help="Input the center point of the pore of [x, y, z] in order with a blank, or 'center_of_geometry'. If not specified, automatically the center point would be automatically searched.  default: None")
    parser.add_argument("--cv", type=float, nargs=3, default=[0,0,1], help="Input the pore axis vector of [x, y, z] in order with a blank, or 'None'. If this is 'None', Hole will automatically search. default: [0,0,1]")
    parser.add_argument("--ignore_res", type=str, nargs="*", default=None, help="Input residue names (< 3 letters only) with a blank to set atoms not included for pore radius calculation. If not specified, water molecules and ions would be ignored. default: ['SOL', 'WAT', 'TIP', 'HOH', 'K', 'NA', 'CL']")
    parser.add_argument("--random", type=int, default=None, help="Input a random seed to set the random number generator. If not specified, the random seed would be automatically set. default: None")
    parser.add_argument("--keep", action='store_true', help="Use this option to keep output files of hole analysis. If you want to generate VMD config file, this option must be used.")
    args = parser.parse_args()
    
    return args

# main
def main():
    ### parse argument and set input parameters
    args=parse_args()
    input_keys = ["deffnm", "dot_den", "n", "csv", "pdf", "vmd", "pdb", "out", "sph", "vdwr", "end_r", "sample", "cp", "cv", "ignore_res", "random", "keep"]
    input_values = [args.deffnm, args.dot_den, args.n, args.csv, args.pdf, args.vmd, args.pdb, args.out, args.sph, args.vdwr, args.end_r, args.sample, args.cp, args.cv, args.ignore_res, args.random, args.keep]
    input_args = dict(zip(input_keys, input_values))
    ### input check
    hole_args = input_check(**input_args)
    ### run hole analysis
    output_args = hole_analysis(**hole_args)
    ### generate output files
    hole_output(**output_args)

if __name__ == "__main__":
    main()