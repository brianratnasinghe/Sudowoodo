#!/usr/bin/env python3
# sweep_eps.py
import argparse, shutil, re, subprocess, textwrap
from pathlib import Path

def parse_args():
    p = argparse.ArgumentParser(description="Sweep epsilon scaling; make run dirs and submit.")
    p.add_argument("--scales", type=float, nargs="+", required=True,
                   help="Epsilon scale factors, e.g. 0.5 1.0 1.5")
    p.add_argument("--out", type=Path, default=Path("sweep_eps"),
                   help="Output parent folder (default: sweep_eps)")
    p.add_argument("--template-top", type=Path, default=Path("afm_system.top"),
                   help="Reference topology to copy (default: afm_system.top)")
    p.add_argument("--template-gro", type=Path, default=Path("afm_system.gro"),
                   help="Reference structure to copy (default: afm_system.gro)")
    p.add_argument("--toppar-dir", type=Path, default=Path("toppar_custom"),
                   help="Directory with ITPs (default: toppar_custom)")
    p.add_argument("--base-itp", type=str, default="sudowoodo_base.itp",
                   help="The base ITP file (inside toppar_custom) where epsilons live")
    p.add_argument("--gmx", type=str, default="gmx",
                   help="GROMACS launcher (default: gmx)")
    p.add_argument("--ntomp", type=int, default=24,
                   help="Number of OpenMP threads (default: 24)")
    p.add_argument("--ntmpi", type=int, default=1,
                   help="Number of MPI ranks (default: 1)")
    p.add_argument("--cluster", choices=["none","slurm"], default="none",
                   help="Run mode: local shell or SLURM (default: none)")
    p.add_argument("--submit", action="store_true",
                   help="Submit/launch runs immediately")
    return p.parse_args()

def ensure_inputs(args):
    missing = []
    for p in [args.template_top, args.template_gro, args.toppar_dir]:
        if not Path(p).exists():
            missing.append(str(p))
    if missing:
        raise FileNotFoundError("Missing required input(s): " + ", ".join(missing))

def copy_tree(src: Path, dst: Path):
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)

def scale_eps_in_itp(itp_path: Path, scale: float):
    txt = itp_path.read_text()
    out_lines = []
    in_nb = False
    line_re = re.compile(r"""
        ^(\s*
        [^;\s]+
        \s+[^;\s]+
        \s+\d+
        \s+([0-9eE\.\+\-]+)
        \s+([0-9eE\.\+\-]+)
        (.*))$
    """, re.VERBOSE)

    for raw in txt.splitlines():
        line = raw
        if re.match(r"^\s*\[\s*nonbond_params\s*\]", line):
            in_nb = True
            out_lines.append(line)
            continue
        if in_nb and re.match(r"^\s*\[", line):
            in_nb = False
            out_lines.append(line)
            continue
        if in_nb and line.strip().startswith(";"):
            out_lines.append(line); continue
        if in_nb and line.strip():
            m = line_re.match(line)
            if m:
                sigma = float(m.group(2))
                epsilon = float(m.group(3))
                rest = m.group(4)
                eps_scaled = epsilon * scale
                new_line = re.sub(r"[ \t]+", " ", f"{m.group(0)}")
                fields = new_line.split()
                fields[4] = f"{sigma:.6f}"
                fields[5] = f"{eps_scaled:.6f}"
                comment = ""
                if ";" in line:
                    comment = " " + line[line.index(";"):]
                out_lines.append(f"{fields[0]:>6} {fields[1]:>6} {fields[2]:>3} {fields[3]:>10} {fields[4]:>12}{comment}")
            else:
                out_lines.append(line)
        else:
            out_lines.append(line)
    itp_path.write_text("\n".join(out_lines) + "\n")

def write_em_mdp(path: Path):
    em_mdp = textwrap.dedent("""\
        integrator  = steep
        emtol       = 100.0
        emstep      = 0.01
        nsteps      = 50000

        nstlist         = 1
        cutoff-scheme   = Verlet
        ns_type         = grid
        coulombtype     = reaction-field
        vdw-type        = cutoff
        vdw-modifier    = Potential-shift-verlet
        rcoulomb        = 2.5
        rvdw            = 2.5
        pbc             = xyz
        verlet-buffer-tolerance  = 0.005
    """)
    path.write_text(em_mdp)

def write_eq_mdp(path: Path):
    eq_mdp = textwrap.dedent("""\
        integrator               = sd
        dt                       = 0.1

        nsteps                   = 400000
        nstcomm                  = 100
        comm-grps                = system

        nstxout                  = 0
        nstvout                  = 0
        nstfout                  = 0
        nstlog                   = 50000
        nstenergy                = 50000
        nstxout-compressed       = 5000
        compressed-x-precision   = 1000
        compressed-x-grps        = system
        energygrps               = system

        cutoff-scheme            = Verlet
        nstlist                  = 20
        ns_type                  = grid
        pbc                      = xyz
        verlet-buffer-tolerance  = 0.005

        coulombtype              = reaction-field
        rcoulomb                 = 6.0
        epsilon_r                = 15
        epsilon_rf               = 0
        vdw_type                 = cutoff
        vdw-modifier             = Potential-shift-verlet
        rvdw                     = 6.0
        rlist                    = 6.0

        tcoupl                   = v-rescale
        tc-grps                  = system
        tau_t                    = 0.05
        ref_t                    = 300
        Pcoupl                   = no; parrinello-rahman
        Pcoupltype               = semiisotropic
        tau_p                    = 12.0         ;parrinello-rahman is more stable with larger tau-p, DdJ, 20130422
        compressibility          = 3e-4  0
        ref_p                    = 1.0   1.0

        gen_vel                  = no
    """)
    path.write_text(eq_mdp)

def write_prod_mdp(path: Path):
    prod_mdp = textwrap.dedent("""\
        integrator               = sd
        dt                       = 0.5
        nsteps                   = 200000000
        nstcomm                  = 100

        nstxout                  = 0
        nstvout                  = 0
        nstfout                  = 0
        nstlog                   = 100000
        nstenergy                = 100000
        nstxout-compressed       = 5000
        compressed-x-precision   = 100

        cutoff-scheme            = Verlet
        nstlist                  = 20
        rlist                    = 20
        rvdw                     = 7
        pbc                      = xyz
        verlet-buffer-tolerance  = 0.005

        coulombtype              = Reaction-Field
        rcoulomb                 = 7
        epsilon_r                = 15
        epsilon_rf               = 0
        vdw_type                 = cutoff
        vdw-modifier             = Potential-shift-verlet

        tcoupl                   = v-rescale
        tc-grps                  = system
        tau_t                    = 0.05
        ref_t                    = 300
        Pcoupl                   = no; parrinello-rahman
        Pcoupltype               = anisotropic
        tau_p                    = 12.0         ;parrinello-rahman is more stable with larger tau-p, DdJ, 20130422
        compressibility          = 3e-4 3e-4 3e-4 0 0 0
        ref_p                    = 1 1 1 0 0 0

        gen_vel                  = no
        gen_temp                 = 300
    """)
    path.write_text(prod_mdp)

def write_run_sh(case_dir: Path, args, case_name: str):
    gmx = args.gmx
    ntmpi = args.ntmpi
    ntomp = args.ntomp
    if args.cluster == "none":
        script = textwrap.dedent(f"""\
            #!/bin/bash
            set -euo pipefail

            {gmx} grompp -f EM.mdp -c afm_system.gro -p afm_system.top -o EM.tpr
            {gmx} mdrun -deffnm EM -ddcheck -ntmpi {ntmpi} -ntomp {ntomp} -dlb no

            {gmx} grompp -f EQ.mdp -c EM.gro -p afm_system.top -o EQ.tpr -maxwarn 2
            {gmx} mdrun -deffnm EQ -ddcheck -ntmpi {ntmpi} -ntomp {ntomp} -dlb no -v

            {gmx} grompp -f production.mdp -c EQ.gro -p afm_system.top -o production.tpr
            {gmx} mdrun -deffnm production -ddcheck -ntmpi {ntmpi} -ntomp {ntomp} -dlb no -v
        """)
    else:
        script = textwrap.dedent(f"""\
            #!/bin/bash
            #SBATCH -J {case_name}
            #SBATCH -o run.out
            #SBATCH -e run.err
            #SBATCH -N 1
            #SBATCH -c {ntomp}
            #SBATCH -t 24:00:00

            set -euo pipefail
            # module load gromacs   # <-- uncomment/adjust for your site

            {gmx} grompp -f EM.mdp -c afm_system.gro -p afm_system.top -o EM.tpr
            srun {gmx} mdrun -deffnm EM -ddcheck -ntmpi {ntmpi} -ntomp {ntomp} -dlb no

            {gmx} grompp -f EQ.mdp -c EM.gro -p afm_system.top -o EQ.tpr -maxwarn 2
            srun {gmx} mdrun -deffnm EQ -ddcheck -ntmpi {ntmpi} -ntomp {ntomp} -dlb no -v

            {gmx} grompp -f production.mdp -c EQ.gro -p afm_system.top -o production.tpr
            srun {gmx} mdrun -deffnm production -ddcheck -ntmpi {ntmpi} -ntomp {ntomp} -dlb no -v
        """)
    sh = case_dir / "run.sh"
    sh.write_text(script)
    sh.chmod(0o755)
    return sh

def main():
    args = parse_args()
    ensure_inputs(args)
    args.out.mkdir(parents=True, exist_ok=True)

    top_src = args.template_top.resolve()
    gro_src = args.template_gro.resolve()
    toppar_src = args.toppar_dir.resolve()

    for scale in args.scales:
        tag = f"{scale:.2f}".rstrip("0").rstrip(".")
        case_name = f"eps_{tag}"
        case_dir = (args.out / case_name).resolve()
        if case_dir.exists():
            shutil.rmtree(case_dir)
        case_dir.mkdir(parents=True)

        shutil.copy2(top_src, case_dir / "afm_system.top")
        shutil.copy2(gro_src, case_dir / "afm_system.gro")
        toppar_dst = case_dir / "toppar_custom"
        copy_tree(toppar_src, toppar_dst)

        base_itp_path = toppar_dst / args.base_itp
        if not base_itp_path.exists():
            raise FileNotFoundError(f"Could not find {args.base_itp} inside {toppar_dst}")
        scale_eps_in_itp(base_itp_path, scale)

        # Write all three MDP files
        write_em_mdp(case_dir / "EM.mdp")
        write_eq_mdp(case_dir / "EQ.mdp")
        write_prod_mdp(case_dir / "production.mdp")

        # Write run script for all three steps
        run_sh = write_run_sh(case_dir, args, case_name)

        print(f"[ok] prepared {case_dir}")

        if args.submit:
            if args.cluster == "none":
                subprocess.Popen(["bash", str(run_sh)], cwd=case_dir)
                print(f"[launched] {case_name} locally")
            else:
                subprocess.run(["sbatch", str(run_sh)], cwd=case_dir, check=True)
                print(f"[submitted] {case_name} via SLURM")

if __name__ == "__main__":
    main()