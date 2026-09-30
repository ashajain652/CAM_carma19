#!/usr/bin/env python3
"""Pre-build CARMA RRTMG optics files for alumina models, with CARMA's own Mie code.

Why: CARMA's default Mie routine (miess, Toon & Ackerman 1981) caps its series at
|m|*x terms (min 135). For alumina, |m| < 1 in LW06/LW07 (next to the reststrahlen
band), so it fails for radii above ~175 um ("miess:: The upper limit for acap is
not enough"). The alumina_d<D>um size models therefore run with
carma_do_optics = .false. and read the files written here.

Mie code: CARMA's bhmie (Bohren & Huffman) with one fix (bhmie_patched.F90): its
downward recurrence for d(n) started too low to converge for large absorbing
particles (up to 0.3% error in Qsca). The patched version matches miess to <= 7e-6
wherever miess works. --mie toon uses CARMA's miess unchanged (for validation).

How: reads NBIN, rmin, vmrat, RHO_ALUMINA, mie_rh and the refidx table from the
model's carma_model_mod.F90, runs carma_optics_driver (CARMA's base/ code built by
build.sh; see carma_optics_driver.F90), and writes one
<model>_CRALUMnn_rrtmg.nc per bin in the same layout as CARMA_CreateOpticsFile
(cam/carma_intr.F90), plus global attributes that record the provenance.

Usage:
  make_alumina_optics.py alumina_d250um [...] [--mie bohren|toon] [--outdir DIR]
  make_alumina_optics.py --all-sizes            # the 10 alumina_d<D>um models

Default outdir: $DIN_LOC_ROOT/atm/cam/physprops/alumina, which is where
build-namelist looks for the alumina_d<D>um files (DIN_LOC_ROOT defaults to
~/MIT/cesm2_port/inputdata; set it, or use --outdir, on other machines).
Needs python3 with numpy and netCDF4.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parent
# This directory is components/cam/tools/carma_optics_offline in a CESM checkout.
MODELS = HERE.parents[1] / "src/physics/carma/models"
DIN = Path(os.environ.get("DIN_LOC_ROOT", Path.home() / "MIT/cesm2_port/inputdata"))
SIZES = [2, 6, 10, 14, 18, 50, 100, 150, 200, 250]
NLW, NSW = 16, 14
SHORTNAME = "CRALUM"
FNUM = r"([0-9.eEdD+-]+)_f"


def parse_model(model):
    src = (MODELS / model / "carma_model_mod.F90").read_text()

    def one(pat):
        m = re.findall(pat, src)
        if len(m) != 1:
            sys.exit(f"{model}: expected one match for {pat!r}, found {len(m)}")
        return m[0]

    p = {
        "nbin": int(one(r"integer, public, parameter\s*::\s*NBIN\s*=\s*(\d+)")),
        "rmin": float(one(r"parameter\s*::\s*rmin\s*=\s*" + FNUM)),
        "vmrat": float(one(r"parameter\s*::\s*vmrat\s*=\s*" + FNUM)),
        "rho": float(one(r"parameter\s*::\s*RHO_ALUMINA\s*=\s*" + FNUM)),
        "rh": float(one(r"mie_rh\(NMIE_RH\)\s*=\s*\(/\s*" + FNUM)),
    }
    table = src[src.index("refidx(:) = (/"):]
    table = table[: table.index("/)") + 2]
    pairs = re.findall(r"\(\s*" + FNUM + r"\s*,\s*" + FNUM + r"\s*\)", table)
    if len(pairs) != NLW + NSW:
        sys.exit(f"{model}: found {len(pairs)} refidx values, expected {NLW + NSW}")
    p["n"] = np.array([float(a) for a, _ in pairs])
    p["k"] = np.array([float(b) for _, b in pairs])
    return p


def run_driver(p, mie):
    exe = HERE / "carma_optics_driver"
    if not exe.exists():
        sys.exit(f"{exe} not found; run build.sh first")
    fmt = lambda a: ", ".join(f"{v:.17g}" for v in a)
    nl = (f"&optics_nl\n nbin = {p['nbin']}\n rmin = {p['rmin']:.17g}\n vmrat = {p['vmrat']:.17g}\n"
          f" rho_elem = {p['rho']:.17g}\n rh = {p['rh']:.17g}\n mie_routine = '{mie}'\n"
          f" n_re = {fmt(p['n'])}\n n_im = {fmt(p['k'])}\n/\n")
    with tempfile.NamedTemporaryFile("w", suffix=".nml", delete=False) as fh:
        fh.write(nl)
        nlpath = fh.name
    res = subprocess.run([str(exe), nlpath], capture_output=True, text=True)
    os.unlink(nlpath)
    if res.returncode != 0 or "MIEFAIL" in res.stdout:
        sys.exit(f"driver failed:\n{res.stdout[-2000:]}\n{res.stderr[-2000:]}")

    nb = p["nbin"]
    out = {"wave": np.zeros(NLW + NSW), "bin": np.zeros((nb, 5)),
           "abs_lw": np.zeros((nb, NLW)), "ext_sw": np.zeros((nb, NSW)),
           "ssa_sw": np.zeros((nb, NSW)), "asm_sw": np.zeros((nb, NSW))}
    for line in res.stdout.splitlines():
        t = line.split()
        if not t:
            continue
        if t[0] == "WAVE":
            out["wave"][int(t[1]) - 1] = float(t[2])
        elif t[0] == "BIN":
            out["bin"][int(t[1]) - 1] = [float(v) for v in t[2:7]]
        elif t[0] == "LW":
            out["abs_lw"][int(t[1]) - 1, int(t[2]) - 1] = float(t[3])
        elif t[0] == "SW":
            ib, iw = int(t[1]) - 1, int(t[2]) - 1
            out["ext_sw"][ib, iw], out["ssa_sw"][ib, iw], out["asm_sw"][ib, iw] = map(float, t[3:6])
    return out


def write_files(model, p, o, mie, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    wave_m = o["wave"] * 1e-2
    paths = []
    for ib in range(p["nbin"]):
        r, rlow, rup, rmass, rho = o["bin"][ib]
        c_name = f"{SHORTNAME}{ib + 1:02d}"
        path = outdir / f"{model}_{c_name}_rrtmg.nc"
        with Dataset(path, "w", format="NETCDF3_CLASSIC") as nc:
            # Same dimensions, variables, attributes and order as CARMA_CreateOpticsFile.
            nc.createDimension("rh_idx", 1)
            nc.createDimension("lw_band", NLW)
            nc.createDimension("sw_band", NSW)
            rhv = nc.createVariable("rh", "f8", ("rh_idx",))
            lwv = nc.createVariable("lw_band", "f8", ("lw_band",))
            swv = nc.createVariable("sw_band", "f8", ("sw_band",))
            rhv.units, lwv.units, swv.units = "fraction", "m", "m"
            rhv.long_name, lwv.long_name, swv.long_name = "relative humidity", "longwave bands", "shortwave bands"
            abs_lw = nc.createVariable("abs_lw", "f8", ("lw_band", "rh_idx"))
            abs_lw.units = "meter^2 kilogram^-1"
            ext = nc.createVariable("ext_sw", "f8", ("sw_band", "rh_idx"))
            ssa = nc.createVariable("ssa_sw", "f8", ("sw_band", "rh_idx"))
            asm = nc.createVariable("asm_sw", "f8", ("sw_band", "rh_idx"))
            ssa.units, ext.units, asm.units = "fraction", "meter^2 kilogram^-1", "-"
            ri = {}
            for nm in ("refindex_real_aer_sw", "refindex_im_aer_sw"):
                ri[nm] = nc.createVariable(nm, "f8", ("sw_band",))
            for nm in ("refindex_real_aer_lw", "refindex_im_aer_lw"):
                ri[nm] = nc.createVariable(nm, "f8", ("lw_band",))
            for nm in ri:
                ri[nm].units = "-"
            ri["refindex_real_aer_sw"].long_name = "real refractive index of aerosol - shortwave"
            ri["refindex_im_aer_sw"].long_name = "imaginary refractive index of aerosol - shortwave"
            ri["refindex_real_aer_lw"].long_name = "real refractive index of aerosol - longwave"
            ri["refindex_im_aer_lw"].long_name = "imaginary refractive index of aerosol - longwave"
            nc.createDimension("opticsmethod_len", 32)
            om = nc.createVariable("opticsmethod", "S1", ("opticsmethod_len",))
            nc.createDimension("namelength", 20)
            an = nc.createVariable("aername", "S1", ("namelength",))
            nc.createDimension("name_len", 32)
            nmv = nc.createVariable("name", "S1", ("name_len",))
            sc = {}
            for nm in ("density", "sigma_logr", "dryrad", "radmin_aer", "radmax_aer",
                       "hygroscopicity", "num_to_mass_ratio"):
                sc[nm] = nc.createVariable(nm, "f8", ())
            units = {"density": "kg m^-3", "sigma_logr": "-", "dryrad": "m", "radmin_aer": "m",
                     "radmax_aer": "m", "hygroscopicity": "-", "num_to_mass_ratio": "kg^-1"}
            lnames = {"density": "aerosol material density",
                      "sigma_logr": "geometric standard deviation of aerosol",
                      "dryrad": "dry number mode radius of aerosol",
                      "radmin_aer": "minimum dry radius of aerosol for bin",
                      "radmax_aer": "maximum dry radius of aerosol for bin",
                      "hygroscopicity": "hygroscopicity of aerosol",
                      "num_to_mass_ratio": "ratio of number to mass of aerosol"}
            for nm in sc:
                sc[nm].units = units[nm]
            for nm in sc:
                sc[nm].long_name = lnames[nm]

            nc.title = f"CARMA RRTMG optics for {model} bin {ib + 1} (pre-built offline)"
            nc.source = ("tools/carma_optics_offline: CARMA base/ Mie code via carma_optics_driver, "
                         "same calculation as CARMA_CreateOpticsFile")
            nc.carma_model = model
            nc.mie_routine = {"bohren": "I_MIERTN_BOHREN1983 (bhmie, patched: d(n) recurrence starts at "
                                        "2*max(xstop,|m x|)+15; see tools/carma_optics_offline/bhmie_patched.F90)",
                              "toon": "I_MIERTN_TOON1981 (miess)"}[mie]

            rhv[:] = [p["rh"]]
            lwv[:] = wave_m[:NLW]
            swv[:] = wave_m[NLW:]
            ri["refindex_real_aer_sw"][:] = p["n"][NLW:]
            ri["refindex_im_aer_sw"][:] = p["k"][NLW:]
            ri["refindex_real_aer_lw"][:] = p["n"][:NLW]
            ri["refindex_im_aer_lw"][:] = p["k"][:NLW]
            nmv[:] = np.array(list(c_name.ljust(32)), "S1")
            an[:] = np.array(list(c_name.ljust(20)[:20]), "S1")
            om[:] = np.array(list("insoluble".ljust(32)), "S1")
            sc["density"][...] = rho * 1e-3 / 1e-6
            sc["sigma_logr"][...] = 0.0
            sc["dryrad"][...] = r * 1e-2
            sc["radmin_aer"][...] = rlow * 1e-2
            sc["radmax_aer"][...] = rup * 1e-2
            sc["hygroscopicity"][...] = 0.0
            sc["num_to_mass_ratio"][...] = 1.0 / rmass / 1e-3
            abs_lw[:, 0] = o["abs_lw"][ib]
            ext[:, 0] = o["ext_sw"][ib]
            ssa[:, 0] = o["ssa_sw"][ib]
            asm[:, 0] = o["asm_sw"][ib]
        paths.append(path)
    return paths


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("models", nargs="*")
    ap.add_argument("--all-sizes", action="store_true", help="the 10 alumina_d<D>um models")
    ap.add_argument("--mie", choices=["bohren", "toon"], default="bohren")
    ap.add_argument("--outdir", type=Path, default=DIN / "atm/cam/physprops/alumina")
    args = ap.parse_args()
    models = list(args.models) + ([f"alumina_d{d}um" for d in SIZES] if args.all_sizes else [])
    if not models:
        ap.error("give model names or --all-sizes")

    for model in models:
        p = parse_model(model)
        o = run_driver(p, args.mie)
        paths = write_files(model, p, o, args.mie, args.outdir)
        print(f"{model}: nbin={p['nbin']} rmin={p['rmin']:.4e} cm vmrat={p['vmrat']} rho={p['rho']} "
              f"mie={args.mie} -> {len(paths)} files in {args.outdir}")


if __name__ == "__main__":
    main()
