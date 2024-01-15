#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
from MDAnalysis.analysis import hole2
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import argparse

# check input parameters
def input_check(*, deffnm, n=1, psf=None, xtc=None, sele="protein", vdwr=None, cp=None, cv=[0,0,1], sample=0.2, end_r=22, ignore_res=None, out_prefix=None, start=0, end=-1, step=1, random=None, csv=False, pdf=False, vmd=False, dot_den=15):
    ### make dictionary of arguments
    args = locals()
    deffnm = args.pop("deffnm")
    if args["vmd"][-4:] == ".vmd":
        args["vmd"] = args["vmd"][:-4]
    
    ### set input / output file names and options
    ##### input / output file names
    if not deffnm:
        if args["xtc"] is None:
            sys.exit("Error: No XTC file is specified.")
        else:
            if args["psf"] is None:
                args["psf"] = args["xtc"][:-4] + ".psf"
            if args["out_prefix"] is None:
                args["out_prefix"] = args["xtc"][:-4]
                if args["csv"] == "True":
                    args["csv"] = args["xtc"][:-4] + ".csv"
                if args["pdf"] == "True":
                    args["pdf"] = args["xtc"][:-4] + ".pdf"
                if args["vmd"] == "True":
                    args["vmd"] = args["xtc"][:-4]
            else:
                if args["csv"] == "True":
                    args["csv"] = args["out_prefix"] + ".csv"
                if args["pdf"] == "True":
                    args["pdf"] = args["out_prefix"] + ".pdf"
                if args["vmd"] == "True":
                    args["vmd"] = args["out_prefix"]
    else:
        if args["psf"] is None:
            args["psf"] = deffnm + ".psf"
        if args["xtc"] is None:
            args["xtc"] = deffnm + ".xtc"
        if args["out_prefix"] is None:
            args["out_prefix"] = deffnm
        if args["csv"] == "True":
            args["csv"] = deffnm + ".csv"
        if args["pdf"] == "True":
            args["pdf"] = deffnm + ".pdf"
        if args["vmd"] == "True":
            args["vmd"] = deffnm
    
        ##### ignore_residues option
    default_ignore_res = ["SOL", "WAT", "TIP", "HOH", "K  ", "NA ", "CL "]
    if args["ignore_res"] is None:
        args["ignore_res"] = default_ignore_res
    else:
        verify_dict = {"Y":True, "y":True, "N":False, "n":False, "":True, "yes":True, "no":False, "YES":True, "NO":False}
        print("IGNORE_RESIDUES option is set to " + str(args["ignore_res"]) + "\n")
        verify_ignore_res = input("Would you like to add these residues to the default IGNORE_RESIDUES option ['SOL', 'WAT', 'TIP', 'HOH', 'K', 'NA', 'CL']? (Y/n): ")
        if verify_dict[verify_ignore_res]:
            args["ignore_res"] = args["ignore_res"] + default_ignore_res
        else:
            args["ignore_res"] = args["ignore_res"]
    
    ### verify input / output files names and options
    ##### input / output file names
    print("Input psf file: " + str(args["psf"]))
    print("Input xtc file: " + str(args["xtc"]))
    print("Output log files: " + str(args["out_prefix"]) + "_run00_hole000.out")
    print("Output sph files: " + str(args["out_prefix"]) + "_run00_hole000.sph")
    if args["csv"]:
        print("Hole output csv: " + str(args["csv"]))
    else:
        print("CSV output will not be generated.")
    if args["pdf"]:
        print("Hole output pdf: " + str(args["pdf"]))
    else:
        print("PDF output will not be generated.")
    if args["vmd"]:
        print("Output vmd: " + str(args["vmd"]) + "_run00.vmd")
    else:
        print("VMD output will not be generated.")
    print("This program will run hole analysis " + str(args["n"]) + " time(s).")
    ##### options
    options_hole = {}
    options_hole["select"] = args["sele"]
    options_hole["vdwradii_file"] = args["vdwr"]
    options_hole["cpoint"] = args["cp"]
    options_hole["cvect"] = args["cv"]
    options_hole["sample"] = args["sample"]
    options_hole["end_radius"] = args["end_r"]
    options_hole["ignore_residues"] = args["ignore_res"]
    options_run = {}
    options_run["start"] = args["start"]
    options_run["stop"] = args["end"]
    options_run["step"] = args["step"]
    options_run["random"] = args["random"]
    if args["vmd"]:
        dot_density = args["dot_den"]
    print("Options for hole analysis;")
    for val in options_hole.keys():
        print(val + ": " + str(options_hole[val]))
    print("Options for hole run;")
    for val in options_run.keys():
        print(val + ": " + str(options_run[val]))
    print("Options for vmd visualization;")
    print("dot_density: " + str(dot_density))
    
    return args

# hole analysis
def hole_analysis(*, n, universe, sele, vdwr, cp, cv, sample, end_r, ignore_res, out_prefix, start, stop, step, random):
    ### setting arguments for hole analysis and output
    input = locals()
    n = input.pop("n")
    holes = []
    
    hole_opts_keys = ["universe", "select", "vdwradii_file", "cpoint", "cvect", "sample", "end_radius", "ignore_residues", "prefix"]
    hole_opts_vals = [input["universe"], input["sele"], input["vdwr"], input["cp"], input["cv"], input["sample"], input["end_r"], input["ignore_res"], input["out_prefix"]]
    hole_opts = dict(zip(hole_opts_keys, hole_opts_vals))
    run_opts_keys = ["start", "stop", "step", "random_seed"]
    run_opts_vals = [input["start"], input["stop"], input["step"], input["random"]]
    run_opts = dict(zip(run_opts_keys, run_opts_vals))
    
    
    ### run hole analysis
    holes_index = range(n)
    for i in holes_index:
        hole_opts["prefix"] = input["out_prefix"] + "_run{:02d}_".format(i+1)
        h = hole2.HoleAnalysis(**hole_opts)
        h.run(**run_opts)
        holes.append(h)
    return holes

# output
def hole_output(*, holes, csv, pdf, vmd, out_prefix, dot_den):
    
    ### dataframe for hole analysis
    holes_df = []
    for run in range(len(holes)):
        z_all = []
        r_all = []
        z_bin = []
        r_mean = []
        r_std = []
        br_dict = {}
        br_all = []
        for f in range(len(holes[run].results.profiles)):
            zrb = []
            for rxn_coord, radius, cen_line_D in holes[run].results.profiles[f]:
                zrb.append((rxn_coord, radius, rxn_coord//1))
            zls,rls,bls = zip(*zrb)
            z_all.append(zls)
            r_all.append(rls)
            for i in range(len(bls)):
                br_all.append((bls[i],rls[i]))
        for b,r in br_all:
            if b in br_dict:
                br_dict[b].append(r)
            else:
                br_dict[b] = [r]
        br_dict_sorted = dict((b,r) for b, r in sorted(br_dict.items()))
        for b in br_dict_sorted:
            z_bin.append(b)
            r_mean.append(np.mean(br_dict[b]))
            r_std.append(np.std(br_dict[b]))
        z_bin = np.array(z_bin)
        r_mean = np.array(r_mean)
        r_std = np.array(r_std)
        holes_df.append((z_all, r_all, z_bin, r_mean, r_std))
    
    ### csv output
    if csv:
        c = open(csv, "w")
        c.write("run,frame,pore axis,pore_radius\n")
        for r in range(len(holes_df)):
            for f in range(len(holes_df[r][0])):
                for i in range(len(holes_df[r][0][f])):
                    c.write(str(r+1) + "," + str(f+1) + "," + str(holes_df[r][0][f][i]) + "," + str(holes_df[r][1][f][i]) + "\n")
        c.close()
    
    ### pdf output
    if pdf:
        if len(holes_df) == 1:
            p = PdfPages(pdf)
            figs, axs = plt.subplots(nrows=2, ncols=1, tight_layout=True, figsize=(10, 7))
            for ax in axs.flat:
                ax.set_xlabel("pore radius (A)")
                ax.set_ylabel("pore axis (A)")
            for f in range(len(z_all)):
                axs[0].plot(holes_df[0][0][f], holes_df[0][1][f])
            axs[0].set_title("run 1, all frames")
            axs[1].plot(z_bin, r_mean)
            axs[1].set_title("run 1, mean & std")
            axs[1].fill_between(z_bin, r_mean - r_std, r_mean + r_std, alpha=0.25)
            p.savefig(figs)
        else:
            p = PdfPages(pdf)
            for i in range(len(holes_df)):
                z_bin = np.array(holes_df[i][2])
                r_mean = np.array(holes_df[i][3])
                r_std = np.array(holes_df[i][4])
                if i % 5 == 0:
                    fig, axs = plt.subplots(nrows=5, ncols=2, tight_layout=True, figsize=(10, 7))
                    for ax in axs.flat:
                        ax.set(xlabel="pore radius (A)", ylabel="pore axis (A)")
                for f in range(len(holes_df[i][0])):
                    axs[i%5, 0].plot(holes_df[i][0][f], holes_df[i][1][f])
                axs[i%5, 0].set_title("run " + str(i+1) + " all frames")
                axs[i%5, 1].plot(z_bin, r_mean)
                axs[i%5, 1].fill_between(z_bin, r_mean - r_std, r_mean + r_std, alpha=0.25)
                axs[i%5, 1].set_title("run " + str(i+1) + " mean & std")
                if i % 5 == 4:
                    p.savefig(fig)
            if len(holes_df) % 5 != 0:
                p.savefig(fig)
        p.close()
    
    ### vmd config file output
    if vmd:
        vmd_prefix = vmd
        for i in range(len(holes)):
            vmd_file = vmd_prefix + "_run{:02d}.vmd".format(i+1)
            holes[i].create_vmd_surface(filename=vmd_file, dot_density=dot_den)

# parse arguments
def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("--deffnm", type=str, default=None, help="Default file name for input and output files which are not specified. Input a file name without extention. If not specified, PDB file name would be the default file name. default: None")
    parser.add_argument("--n", type=int, default=1, help="Number of times to run hole analysis. default: 1")
    parser.add_argument("--psf", type=str, default=None, help="Input PSF file. Specify a file name or use default file name. default: None")
    parser.add_argument("--xtc", type=str, default=None, help="Input XTC file. Specify a file name or use default file name. default: None")
    parser.add_argument("--sele", type=str, default="protein", help="Atom-selection string for MDAnalysis Universe. If not specified,  default: protein")
    parser.add_argument("--vdwr", type=str, default=None, help="Input file of van der Waals radii for each atom. If not specified, the 'simple2.rad' file would be used. default: None")
    parser.add_argument("--cp", type=float, nargs=3, default=None, help="Input the [x, y, z] of the center point of the pore in order with a blank, or 'center_of_geometry'. If not specified, the center point would be automatically searched. default: None")
    parser.add_argument("--cv", type=float, nargs=3, default=[0,0,1], help="Input the [x, y, z] of the pore axis vector in order with a blank, or 'None'. If this is 'None', Hole will automatically search. default: [0,0,1]")
    parser.add_argument("--sample", type=float, default=0.2, help="Set the sampling distance [A] for the hole analysis. default: 0.2")
    parser.add_argument("--end_r", type=float, default=22, help="Set the radius [A] to end hole analysis. default: 22")
    parser.add_argument("--ignore_res", type=str, nargs="*", default=None, help="Input residue names (< 3 letters only) in order with a blank to set atoms not included for pore radius calculation. If not specified, water molecules and ions would be ignored. default: ['SOL', 'WAT', 'TIP', 'HOH', 'K', 'NA', 'CL']")
    parser.add_argument("--out_prefix", type=str, default=None, help="Input a prefix for output (.out and .sph) files. If not specified, the default file name would be used. default: None")
    parser.add_argument("--start", type=int, default=0, help="Input the starting frame number. default: 0")
    parser.add_argument("--end", type=int, default=-1, help="Input the ending frame number. default: -1")
    parser.add_argument("--step", type=int, default=1, help="Input the step width. default: 1")
    parser.add_argument("--random", type=int, default=None, help="Input a random seed to set the random number generator. If not specified, the random seed would be automatically set. default: None")
    parser.add_argument("--csv", type=str, default=None, help="Output csv file. Specify a file name or set this to True to generate csv output. default: None")
    parser.add_argument("--pdf", type=str, default=None, help="Output pdf file. Specify a file name or set this to True to generate pdf output. default: None")
    parser.add_argument("--vmd", type=str, default=None, help="Output vmd config file. Specify a file name or set this to True to generate vmd config file. default: None")
    parser.add_argument("--dot_den", type=int, default=15, help="Dot density for VMD visualization. Setting to 5-35 is recommended. default: 15")
    args = parser.parse_args()
    
    return args

# main
def main():
    ### parse argument and set input parameters
    input_args=vars(parse_args())
    ### input check
    args=input_check(**input_args)
    ### making universe
    u = mda.Universe(args["psf"], args["xtc"])
    ### run hole analysis
    hole_args_keys = ["n", "universe", "sele", "vdwr", "cp", "cv", "sample", "end_r", "ignore_res", "out_prefix", "start", "stop", "step", "random"]
    hole_args_vals = [args["n"], u, args["sele"], args["vdwr"], args["cp"], args["cv"], args["sample"], args["end_r"], args["ignore_res"], args["out_prefix"], args["start"], args["end"], args["step"], args["random"]]
    hole_args = dict(zip(hole_args_keys, hole_args_vals))
    holes = hole_analysis(**hole_args)
    ### generate output files
    output_args_keys = ["holes", "csv", "pdf", "vmd", "out_prefix", "dot_den"]
    output_args_vals = [holes, args["csv"], args["pdf"], args["vmd"], args["out_prefix"], args["dot_den"]]
    output_args = dict(zip(output_args_keys, output_args_vals))
    hole_output(**output_args)

if __name__ == "__main__":
    main()