#!/usr/bin/env python3
"""
AFM cell wall builder tool with custom epsilon mapping.
- Copies template .gro and .itp files.
- Modifies sudowoodo_base.itp with user-specified epsilon for bead pairs.
- Generates topology, .mdp files, run.sh, and afm_build.log.
- Calls build_afm_system.py to create afm_system.gro in the output folder.

Usage:
  python afm_build_sweep.py --out run_$(date +%s) --epsilon CC=1.0,CX=0.8,CP=0.7,XX=0.6,XP=0.5,PP=0.4
Optional:
  --seed 123456
"""

import argparse, shutil, os, re, random, textwrap, subprocess
from pathlib import Path

def get_args():
    p = argparse.ArgumentParser(description="AFM cell wall builder and sweep tool (custom epsilon mapping)")
    p.add_argument('--out', type=Path, required=True, help="Output folder")
    p.add_argument('--epsilon', type=str, required=True,
                   help="Comma-separated epsilon mapping, e.g. CC=1.0,CX=0.8,CP=0.7,XX=0.6,XP=0.5,PP=0.4")
    p.add_argument('--seed', type=int, help="Random seed (int). If not set, random seed is chosen and logged.")
    p.add_argument('--nxylo', type=int, default=458)
    p.add_argument('--npctn', type=int, default=5501)
    p.add_argument('--ncell', type=int, default=146)
    p.add_argument('--gmx', type=str, default="gmx")
    p.add_argument('--ntomp', type=int, default=24)
    p.add_argument('--ntmpi', type=int, default=1)
    return p.parse_args() 

def ensure_dir(path):
    path.mkdir(parents=True, exist_ok=True)

def copy_file(src, dst):
    shutil.copy2(src, dst)

def write_text(path, txt):
    path.write_text(txt)

def parse_epsilon_map(epsilon_str):
    mapping = {}
    for item in epsilon_str.split(','):
        key, val = item.split('=')
        key = key.strip().upper()
        val = float(val)
        if key in ['CC','CX','CP','XX','XP','PP']:
            mapping[(key[0], key[1])] = val
            mapping[(key[1], key[0])] = val
    return mapping

def scale_epsilon_in_itp(itp_path, new_path, epsilon_map):
    re_lj = re.compile(r'^(\s*\w+\s+\w+\s+\d+\s+([0-9eE\.\+\-]+)\s+([0-9eE\.\+\-]+))')
    lines = itp_path.read_text().splitlines()
    out = []
    in_nb = False
    for line in lines:
        if '[ nonbond_params' in line:
            in_nb = True
            out.append(line)
            continue
        if in_nb and line.strip().startswith('['):
            in_nb = False
        if in_nb and re_lj.match(line):
            parts = line.split()
            i, j = parts[0], parts[1]
            sigma = float(parts[3])
            epsilon = float(parts[4])
            new_epsilon = epsilon_map.get((i, j), epsilon)
            parts[4] = f"{new_epsilon:.6f}"
            out.append(' '.join(parts))
        else:
            out.append(line)
    new_path.write_text('\n'.join(out) + '\n')

def randomize_structures(seed, out_dir):
    for src, dst in [
        ('X.gro', out_dir / 'X.gro'),
        ('P.gro', out_dir / 'P.gro'),
        ('C.gro', out_dir / 'C.gro')
    ]:
        copy_file(src, dst)
    # Optionally: call your actual randomization logic here, passing seed as needed
    # Example: subprocess.run(["python", "build_afm_system.py", "--seed", str(seed), ...])
    return

def generate_topology(args, out_dir):
    top_txt = textwrap.dedent(f"""\
        ;;;;;; AFM-Based Combined Topology
        #include "toppar_custom/sudowoodo_base.itp"
        #include "toppar_custom/sudowoodo_xyloglucan.itp"
        #include "toppar_custom/sudowoodo_pectin.itp"
        #include "toppar_custom/sudowoodo_cellulose.itp"

        [ system ]
        AFM-Based Cell Wall System

        [ molecules ]
        Cell      {args.ncell}
        Xylo      {args.nxylo}
        Pctn      {args.npctn}
    """)
    write_text(out_dir / "afm_system.top", top_txt)

def generate_itps(args, out_dir, epsilon_map):
    toppar_dir = out_dir / "toppar_custom"
    ensure_dir(toppar_dir)
    scale_epsilon_in_itp(Path('toppar_custom/sudowoodo_base.itp'), toppar_dir / "sudowoodo_base.itp", epsilon_map)
    for src, dst in [
        ('toppar_custom/sudowoodo_xyloglucan.itp', toppar_dir / "sudowoodo_xyloglucan.itp"),
        ('toppar_custom/sudowoodo_pectin.itp', toppar_dir / "sudowoodo_pectin.itp"),
        ('toppar_custom/sudowoodo_cellulose.itp', toppar_dir / "sudowoodo_cellulose.itp")
    ]:
        copy_file(src, dst)

def write_mdp_files(args, out_dir):
    def mdp_default_em():
        return textwrap.dedent("""\
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
    def mdp_default_eq():
        return textwrap.dedent("""\
            integrator               = sd
            dt                       = 0.05
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
            tau_p                    = 12.0
            compressibility          = 3e-4  0
            ref_p                    = 1.0   1.0
            gen_vel                  = no
        """)
    def mdp_default_prod():
        return textwrap.dedent("""\
            integrator               = sd
            dt                       = 0.1
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
            tau_p                    = 12.0
            compressibility          = 3e-4 3e-4 3e-4 0 0 0
            ref_p                    = 1 1 1 0 0 0
            gen_vel                  = no
            gen_temp                 = 300
        """)
    write_text(out_dir / "EM.mdp", mdp_default_em())
    write_text(out_dir / "EQ.mdp", mdp_default_eq())
    write_text(out_dir / "production.mdp", mdp_default_prod())

def write_run_sh(args, out_dir):
    sh_txt = textwrap.dedent(f"""\
        #!/bin/bash
        set -euo pipefail

        {args.gmx} grompp -f EM.mdp -c afm_system.gro -p afm_system.top -o EM.tpr
        {args.gmx} mdrun -deffnm EM -ddcheck -ntmpi {args.ntmpi} -ntomp {args.ntomp} -dlb no

        {args.gmx} grompp -f EQ.mdp -c EM.gro -p afm_system.top -o EQ.tpr -maxwarn 2
        {args.gmx} mdrun -deffnm EQ -ddcheck -ntmpi {args.ntmpi} -ntomp {args.ntomp} -dlb no -v

        {args.gmx} grompp -f production.mdp -c EQ.gro -p afm_system.top -o production.tpr
        {args.gmx} mdrun -deffnm production -ddcheck -ntmpi {args.ntmpi} -ntomp {args.ntomp} -dlb no -v
    """) 
    sh_path = out_dir / "run.sh"
    write_text(sh_path, sh_txt)
    os.chmod(sh_path, 0o755)

def write_log(out_dir, seed, args, epsilon_map):
    eps_map_str = ', '.join([f"{k[0]}{k[1]}={v}" for k,v in epsilon_map.items() if k[0]<=k[1]])
    log_txt = textwrap.dedent(f"""\
        AFM-Build Sweep Run Log
        ======================
        Output directory: {out_dir}
        Epsilon mapping: {eps_map_str}
        Polymer counts: Xylo={args.nxylo}  Pctn={args.npctn}  Cell={args.ncell}
        Seed used: {seed}
    """)
    write_text(out_dir / "afm_build.log", log_txt)

def build_afm_system(seed, out_dir=None):
    """
    Call build_afm_system.py with the given seed inside the output folder.
    """
    print(f"[info] Building afm_system.gro using build_afm_system.py ...")
    builder = Path(__file__).parent / "build_afm_system.py"
    if not builder.exists():
        raise FileNotFoundError(f"Could not find build_afm_system.py in {builder.parent}")

    cmd = ["python", str(builder), "--seed", str(seed)]
    subprocess.run(cmd, cwd=out_dir, check=True)

def main():
    args = get_args()
    ensure_dir(args.out)
    seed = args.seed if args.seed is not None else random.randint(1,99999999)
    random.seed(seed)
    epsilon_map = parse_epsilon_map(args.epsilon)
    write_log(args.out, seed, args, epsilon_map)
    randomize_structures(seed, args.out)
    generate_topology(args, args.out)
    generate_itps(args, args.out, epsilon_map)
    write_mdp_files(args, args.out)
    write_run_sh(args, args.out)
    build_afm_system(seed, args.out)
    print(f"[ok] Setup complete in {args.out} (seed={seed})")

if __name__ == "__main__":
    main()