#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from subprocess import call
import numpy as np
from os import makedirs
from shutil import copy

import dns


def main():

    parser = argparse.ArgumentParser(
        description="Continue a run from the last state file saved",
    )
    parser.add_argument("rundir", type=str, help="path to the run")
    parser.add_argument("newroot", type=str, help="path to the run")
    parser.add_argument(
        "-i_finish_plus",
        type=int,
        dest="i_finish_plus",
        help="number of time steps to add to i_finish",
    )
    parser.add_argument(
        "--noray",
        action="store_true",
        dest="noray",
    )

    args = vars(parser.parse_args())

    dnsreroot(**args)


def dnscontinue(rundir, newroot, i_finish_plus=None, noray=False):

    rundir = Path(rundir)
    newroot = Path(newroot)
    rundir_out = newroot / rundir.name
    makedirs(rundir_out)

    states = sorted(list(rundir.glob("state.*")))
    if len(states) > 2:
        stateout = states[-2]
    else:
        stateout = states[-1]

    parameters = dns.readParameters(rundir / "parameters.in")
    parameters["initiation"]["ic"] = 0
    parameters["initiation"]["i_start"] = 0
    if i_finish_plus is not None:
        parameters["termination"]["i_finish"] += i_finish_plus

    if noray:
        parameters["physics"]["sigma_r"] = False

    parameters["initiation"]["t_start"] = 0
    dns.writeParameters(parameters, rundir_out / "parameters.in")

    copy(stateout, rundir_out / "state.000000")
    for f in rundir.glob("*.slurm"):
        copy(f, rundir_out / f.name)



if __name__ == "__main__":
    main()
