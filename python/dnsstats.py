#!/usr/bin/env python3
"""Produces plots of statistics of a run.
Figures are saved in a seperate directory within the run directory, and by
default they are shown.

"""

import argparse
from pathlib import Path
from sys import argv

import numpy as np
from matplotlib import pyplot as plt

import dns


def main():
    parser = argparse.ArgumentParser(description="Produce plots of time series.")
    parser.add_argument("runDir", type=str, help="path to the run folder.")
    parser.add_argument(
        "Ni",
        type=int if "--t" not in argv else float,
        help="initial line number in stat.gp, or if --tfilter, initial time.",
    )
    parser.add_argument(
        "Nf",
        type=int if "--t" not in argv else float,
        help="final line number if stat.gp, or if --tfilter, final time.",
    )
    parser.add_argument(
        "--noshow", action="store_true", dest="noshow", help="do not display the plots."
    )
    parser.add_argument(
        "--tfilter",
        action="store_true",
        dest="tfilter",
        help="use initial final time instead of line numbers.",
    )
    parser.add_argument(
        "--tex", action="store_true", dest="tex", help="use LaTeX to render text."
    )
    parser.add_argument(
        "--diet",
        action="store_true",
        dest="diet",
        help="plot only energy and input-dissipation.",
    )
    args = vars(parser.parse_args())

    dnsstats(**args)


def dnsstats(
    runDir, Ni, Nf, tfilter=False, noshow=False, tex=False, diet=False,
):

    dns.setPlotDefaults(tex=tex)

    runDir = Path(runDir)
    figuresDir = dns.createFiguresDir(runDir)

    statsfile = "stat.gp"
    stepsfile = "steps.gp"
    stats = np.loadtxt(runDir / statsfile, ndmin=2)

    if not diet and Path.is_file(runDir / stepsfile):
        steps = np.loadtxt(runDir / stepsfile, ndmin=2)

    mhdfile = "stat_mhd.gp"
    rayfile = "stat_ray.gp"
    fracfile = "stat_frac.gp"

    phasesfile = runDir / "phases.gp"
    if Path.is_file(phasesfile):
        phases = np.loadtxt(phasesfile)
        phase = True
    else:
        phase = False

    # non-dimensionalize based on laminar velocity and its characteristic
    # length scale
    # assumed here qF = 1
    nml = dns.readParameters(runDir / "parameters.in")
    Re = nml["physics"]["Re"]
    try:
        tilt_angle = nml["physics"]["tilt_angle"]
    except:
        tilt_angle = 0

    try:
        Ha = nml["physics"]["Ha"]

        if Path.is_file(runDir / mhdfile):
            mhds = np.loadtxt(runDir / mhdfile, ndmin=2)
            mhd = True
            frac = True
    except:
        Ha = 0

    try:
        sigma_R = nml["physics"]["sigma_R"]

        if Path.is_file(runDir / rayfile):
            rays = np.loadtxt(runDir / rayfile, ndmin=2)
            ray = True
            frac = True
    except:
        sigma_R = 0

    if Path.is_file(runDir / fracfile):
        fracs = np.loadtxt(runDir / fracfile, ndmin=2)
        frac = True
        frac0 = True

    Lx = nml["grid"]["Lx"]
    Lz = nml["grid"]["Lz"]
    nx = nml["grid"]["nx"]
    ny = nml["grid"]["ny"]
    nz = nml["grid"]["nz"]

    if mhd or ray:
        if abs(tilt_angle) > 0:
            title = f"$\\mathrm{{Ha}}={Ha:.1f}$, $\\mathrm{{Re}}={Re:.1f}$, $\\sigma_R={sigma_R:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $\\theta={tilt_angle:.1f}$, $N=({nx},{ny},{nz})$"
        else:
            title = f"$\\mathrm{{Ha}}={Ha:.1f}$, $\\mathrm{{Re}}={Re:.1f}$, $\\sigma_R={sigma_R:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $N=({nx},{ny},{nz})$"
    else:
        if abs(tilt_angle) > 0:
            title = f"$\\mathrm{{Re}}={Re:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $\\theta={tilt_angle:.1f}$, $N=({nx},{ny},{nz})$"
        else:
            title = f"$\\mathrm{{Re}}={Re:.1f}$, $L=({Lx:.1f},{dns.Ly:.1f},{Lz:.1f})$, $N=({nx},{ny},{nz})$"

    # no non-dimensionalization whatsoever, use stat.gp as is
    if tfilter:
        Ni_ = np.transpose(np.nonzero(stats[:, 1] > Ni))[0][0]
        Nf_ = np.transpose(np.nonzero(stats[:, 1] < Nf))[-1][0]
        stats = stats[Ni_:Nf_]

        if mhd:
            Ni_ = np.transpose(np.nonzero(mhds[:, 1] > Ni))[0][0]
            Nf_ = np.transpose(np.nonzero(mhds[:, 1] < Nf))[-1][0]
            mhds = mhds[Ni_:Nf_]

        if ray:
            Ni_ = np.transpose(np.nonzero(rays[:, 1] > Ni))[0][0]
            Nf_ = np.transpose(np.nonzero(rays[:, 1] < Nf))[-1][0]
            rays = rays[Ni_:Nf_]

        if frac0:
            Ni_ = np.transpose(np.nonzero(fracs[:, 1] > Ni))[0][0]
            Nf_ = np.transpose(np.nonzero(fracs[:, 1] < Nf))[-1][0]
            fracs = fracs[Ni_:Nf_]

        if phase:
            Ni_ = np.transpose(np.nonzero(phases[:, 1] > Ni))[0][0]
            Nf_ = np.transpose(np.nonzero(phases[:, 1] < Nf))[-1][0]
            phases = phases[Ni_:Nf_]

        if not diet:
            Ni_ = np.transpose(np.nonzero(steps[:, 1] > Ni))[0][0]
            Nf_ = np.transpose(np.nonzero(steps[:, 1] < Nf))[-1][0]
            steps = steps[Ni_:Nf_]

    amp = np.pi**2
    if mhd:
        amp += 4*Ha**2

    Elam = 1 / 4
    Edotlam = amp / (8 * Re)

    if ray:
        Edotlam += 2 * sigma_R * Elam

    KineticEnergy = stats[:, 2]
    Production = stats[:, 3]
    Dissipation = stats[:, 4]
    normRHS = stats[:, 5]

    timeLabel = "$t$"
    ekinLabel = "$E / E_L$"
    prodLabel = "$I / I_L$"
    dissLabel = "$\\epsilon / \\epsilon_L$"
    normLabel = "$|{{\\bf F}}({{\\bf u}})|$"

    # Turbulent kinetic energy
    figKin, axKin = plt.subplots()
    axKin.set_xlabel(timeLabel)
    axKin.set_ylabel(ekinLabel)
    axKin.plot(stats[:, 1], KineticEnergy / Elam)
    axKin.set_title(title)
    figKin.savefig(figuresDir / "ekin.png")

    # Input
    figin, axin = plt.subplots()
    axin.set_xlabel(timeLabel)
    axin.set_ylabel(prodLabel)
    if not (mhd or ray):
        axin.plot(stats[:, 1], Production / Edotlam)
    else:
        input_forcing = Production - (mhds[:, 2] + rays[:, 2])
        axin.plot(stats[:, 1], input_forcing / Edotlam, label="Body")
    if mhd:
        axin.plot(mhds[:, 1], mhds[:, 2] / Edotlam, label="MHD")
    if ray:
        axin.plot(rays[:, 1], rays[:, 2] / Edotlam, label="Rayleigh")
    if mhd or ray:
        axin.legend()
    axin.set_title(title)
    figin.savefig(figuresDir / "input.png")

    # Dissipation
    figdis, axdis = plt.subplots()
    axdis.set_xlabel(timeLabel)
    axdis.set_ylabel(dissLabel)
    if not (mhd or ray):
        axdis.plot(stats[:, 1], Dissipation / Edotlam)
    else:
        dissip_forcing = Dissipation - (mhds[:, 3] + rays[:, 3])
        axdis.plot(stats[:, 1], dissip_forcing / Edotlam, label="Viscous")
    if mhd:
        axdis.plot(mhds[:, 1], mhds[:, 3] / Edotlam, label="MHD")
    if ray:
        axdis.plot(rays[:, 1], rays[:, 3] / Edotlam, label="Rayleigh")
    if mhd or ray:
        axdis.legend()
    axdis.set_title(title)
    figdis.savefig(figuresDir / "dissip.png")

    if frac:
        figf, axf = plt.subplots()
        axf.set_xlabel(timeLabel)
        axf.set_ylabel("$F$")
        if frac0:
            axf.plot(fracs[:, 1], fracs[:, 2], label="code")

        else:

            try:
                Ry = nml["symmetries"]["Ry"]
            except:
                Ry = 0

            if mhd:
                f_mhd = 2 * (KineticEnergy - mhds[:, 3] / (2*Ha**2/Re))
                if Ry:
                    f_mhd *= ny / (ny - 2)
                axf.plot(stats[:, 1], f_mhd, label="mhd")

            elif ray:
                f_ray = 2 * (KineticEnergy - rays[:, 3] / (2*sigma_R))
                if Ry:
                    f_ray *= ny / (ny - 2)
                axf.plot(stats[:, 1], f_ray, label="ray")

            if mhd or ray:
                axf.legend()
    
        axf.set_title(title)
        figf.savefig(figuresDir / "v2_avg.png")

    # Production-Dissipation plot
    # factor = 0.1
    factor = 0
    figProdDis, axProdDis = plt.subplots()
    axProdDis.set_xlabel(prodLabel)
    axProdDis.set_ylabel(dissLabel)
    axProdDis.plot(Production / Edotlam, Dissipation / Edotlam, color="gray", alpha=0.5)

    minDP = min(0, min(min(Production / Edotlam), min(Dissipation / Edotlam)))
    maxDP = max(max(Production / Edotlam), max(Dissipation / Edotlam))
    unitline = np.linspace(minDP, (1.0 + factor) * maxDP, 100)
    axProdDis.plot(unitline, unitline, "g")
    axProdDis.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    axProdDis.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    # axProdDis.set_xlim(min(0, min(Production / Edotlam)), (1.0 + factor) * max(Production / Edotlam))
    # axProdDis.set_ylim(min(0, min(Dissipation / Edotlam)), (1.0 + factor) * max(Dissipation / Edotlam))
    axProdDis.set_xlim(min(Production / Edotlam), (1.0 + factor) * max(Production / Edotlam))
    axProdDis.set_ylim(min(Dissipation / Edotlam), (1.0 + factor) * max(Dissipation / Edotlam))
    axProdDis.set_title(title)
    figProdDis.savefig(figuresDir / "peps.png")

    if Path.is_file(phasesfile):
        tphases = phases[:, 1]
        dtphases = tphases[1:] - tphases[:-1]
        phases_x = np.unwrap(phases[:, 2])
        phases_z = np.unwrap(phases[:, 3])
        dphases_x = (phases_x[1:] - phases_x[:-1]) / dtphases
        dphases_z = (phases_z[1:] - phases_z[:-1]) / dtphases

        # Assumed Ly=2, then give Lx and Lz in units of Ly

        fig, ax = plt.subplots()
        ax.plot(tphases[:-1], dphases_x * Lx / (4 * np.pi))
        ax.set_xlabel("$t$")
        ax.set_ylabel(f"$\\dot{{\\phi_x}} L_x / (2\\pi)$")
        ax.set_xlim(left=tphases[0], right=tphases[-2])
        ax.set_title(title)
        fig.savefig(figuresDir / f"dphase_x.png", bbox_inches="tight")

        fig, ax = plt.subplots()
        ax.plot(tphases[:-1], dphases_z * Lz / (4 * np.pi))
        ax.set_xlabel("$t$")
        ax.set_ylabel(f"$\\dot{{\\phi_z}} L_z / (2\\pi)$")
        ax.set_xlim(left=tphases[0], right=tphases[-2])
        ax.set_title(title)
        fig.savefig(figuresDir / f"dphase_z.png", bbox_inches="tight")

    if not diet:

        # |RHS|
        figRHS, axRHS = plt.subplots()
        axRHS.set_xlabel(timeLabel)
        axRHS.set_ylabel(normLabel)
        axRHS.plot(stats[:, 1], normRHS)
        axRHS.set_title(title)
        figRHS.savefig(figuresDir / "normrhs.png")

        if Path.is_file(runDir / stepsfile):
            dt = steps[:, 2]
            courant = steps[:, 3]
            err_corr = steps[:, 4]
            ncorr = steps[:, 5]

            # time stepping stats
            figts, (axts1, axts2, axts4, axts3) = plt.subplots(
                nrows=4, ncols=1, sharex=True
            )
            axts3.plot(steps[:, 1], courant)
            axts3.set_xlabel(timeLabel)
            axts3.set_ylabel("Courant number")

            axts4.plot(steps[:, 1], dt)
            axts4.set_xlabel(timeLabel)
            axts4.set_ylabel("$dt$")

            axts2.plot(steps[:, 1][:-1], err_corr[1:])
            axts2.set_ylabel("Corrector rel. err.")

            axts1.plot(steps[:, 1][:-1], ncorr[1:])
            axts1.set_ylabel("Corrector iterations")
            axts1.yaxis.get_major_locator().set_params(integer=True)
            axts1.set_title(title)

            figts.savefig(figuresDir / "steps.png")

    if not noshow:
        plt.show()


if __name__ == "__main__":
    main()
