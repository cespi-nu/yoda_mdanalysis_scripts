#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("--sel", type=str, nargs='+', required=True, help="Input text to select a residue. required.")
    args = parser.parse_args()
    
    return args

def main():
    args = parseargs()
    print(*args.sel)

if __name__ == "__main__":
    main()