#!/usr/bin/env python3
import argparse
from pathlib import Path

import dns


def main():

    parser = argparse.ArgumentParser(
        description="Delete states in a directory until the current resumption point.",
    )
    parser.add_argument("rundir", type=str, help="path to the run")

    args = vars(parser.parse_args())

    dnsrmtocontinuation(**args)

    
def dnsrmtocontinuation(rundir):

    rundir = Path(rundir)
    parameters = dns.readParameters(rundir / "parameters.in")
    i_start = parameters["initiation"]["i_start"]
    i_save_fields = parameters["output"]["i_save_fields"]

    istate_start = i_start // i_save_fields
    
    states = sorted(list(rundir.glob("state.*")))
    if len(states) > 0:
        for state in states:
            i_state = int(state.name[-6:]) 
            if i_state < istate_start:
                state.unlink()

if __name__ == "__main__":
    main()
