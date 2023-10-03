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
    parser.add_argument(
        "-script",
        dest="script",
        help="Submission script. If given, submit the job by sbatch script",
    )
    parser.add_argument(
        "--newdir",
        action="store_true",
        dest="newdir",
    )

    args = vars(parser.parse_args())

    dnscontinue(**args)


def dnscontinue(rundir, i_finish_plus=None, noray=False, script=None, newdir=False):

    rundir = Path(rundir)
    if not newdir:
        rundir_out = rundir
    else:
        rundir_out = rundir / "continue"
        makedirs(rundir_out)

    states = sorted(list(rundir.glob("state.*")))
    if len(states) > 2:
        i_final_state = int(states[-2].name[-6:])
        stateout = states[-2]
    else:
        i_final_state = int(states[-1].name[-6:])
        stateout = states[-1]

    parameters = dns.readParameters(rundir / "parameters.in")
    if not newdir:
        parameters["initiation"]["ic"] = i_final_state
        itime_final = i_final_state * parameters["output"]["i_save_fields"]
        parameters["initiation"]["i_start"] = itime_final
    else:
        parameters["initiation"]["ic"] = 0
        parameters["initiation"]["i_start"] = 0
    if i_finish_plus is not None:
        parameters["termination"]["i_finish"] += i_finish_plus

    if noray:
        parameters["physics"]["sigma_r"] = False

    if not newdir:
        stat_file = rundir / "stat.gp"
        if Path.is_file(stat_file):
            stats = np.loadtxt(rundir / "stat.gp")
            times = t_final_state = stats[stats[:, 0] == itime_final]
            if len(times) > 0:
                t_final_state = stats[stats[:, 0] == itime_final][-1][1]
            else:
                t_final_state = itime_final * parameters["time_stepping"]["dt"]
        else:
            t_final_state = itime_final * parameters["time_stepping"]["dt"]

        parameters["initiation"]["t_start"] = t_final_state
    else:
        parameters["initiation"]["t_start"] = 0
    dns.writeParameters(parameters, rundir_out / "parameters.in")

    if not newdir:
        files = list(rundir.glob("*.gp"))
        for file in files:
            with open(file, "r") as f:
                lines = f.readlines()

            with open(file, "w") as f:
                for line in lines:
                    try:
                        i_time = int(re.search(r"\d+", line).group())

                        if i_time > itime_final:
                            break
                        else:
                            f.write(line)

                    except:
                        f.write(line)
    else:
        copy(stateout, rundir_out / stateout.name)
        for f in rundir.glob("*.slurm"):
            copy(f, rundir_out / f.name)

    if not script == None:
        call(["sbatch", script])


if __name__ == "__main__":
    main()
