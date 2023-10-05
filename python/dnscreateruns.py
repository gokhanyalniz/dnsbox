#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
from os import makedirs
from shutil import copy

import dns


def main():

    parser = argparse.ArgumentParser(
        description="Create run folders.",
    )
    parser.add_argument("parentdir", type=str, help="path to the parent dir for runs")
    parser.add_argument(
        "Re_i",
        type=float,
    )
    parser.add_argument(
        "Re_f",
        type=float,
    )
    parser.add_argument(
        "Re_delta",
        type=float,
    )

    args = vars(parser.parse_args())
    dnscreateruns(**args)


def dnscreateruns(parentdir, Re_i, Re_f, Re_delta):

    parentdir = Path(parentdir)
    paramsfile = parentdir / "parameters.in"
    icfile = parentdir / "state.000000"
    slurmfile = parentdir / "dns.slurm"

    parameters = dns.readParameters(paramsfile)
    parameters["initiation"]["ic"] = 0
    parameters["initiation"]["i_start"] = 0
    parameters["initiation"]["t_start"] = 0

    Res = np.arange(Re_i, Re_f, Re_delta)
    for iRe in range(len(Res)):
        Re = Res[iRe]
        dirname = f"re{Re:.2f}"
        outdir = parentdir / dirname

        makedirs(outdir)
        copy(icfile, outdir)
        copy(slurmfile, outdir)

        parameters["physics"]["re"] = Re
        dns.writeParameters(parameters, outdir / "parameters.in")

if __name__ == "__main__":
    main()
