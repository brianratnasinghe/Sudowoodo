#!/usr/bin/env python3
"""
AFM cell wall builder tool with custom epsilon mapping and ktheta modifications.
- Copies template .gro and .itp files.
- Modifies sudowoodo_base.itp with user-specified epsilon for bead pairs.
- Modifies polymer .itp files with user-specified ktheta values.
- Generates topology, .mdp files, run.sh, and afm_build.log.
- Calls build_afm_system.py to create afm_system.gro in the output folder.

Usage:
  python afm_build_sweep.py --out run_$(date +%s) --epsilon CC=1.0,CX=0.8,CP=0.7,XX=0.6,XP=0.5,PP=0.4
Optional:
  --seed 123456
  --ktheta "120,150,180"    # pectin,cellulose,xyloglucan
  --ktheta ",150,180"       # use default pectin, custom cellulose and xyloglucan
  --ktheta "120,,"          # custom pectin, default cellulose and xyloglucan
"""

import argparse, shutil, os, re, random, textwrap, subprocess
from pathlib import Path

def get_args():
    p = argparse.ArgumentParser(description="AFM cell wall builder and sweep tool (custom epsilon mapping)")
    p.add_argument('--out', type=Path, required=True, help="Output folder")
    p.add_argument('--epsilon', type=str, required=True,
                   help="Comma-separated epsilon mapping, e.g. CC=1.0,CX=0.8,CP=0.7,XX=0.6,XP=0.5,PP=0.4")
    p.add_argument('--ktheta', type=str, 
                   help="Comma-separated ktheta values for pectin,cellulose,xyloglucan. "
                        "Use empty values to keep defaults, e.g. '120,150,180' or ',150,180' or '120,,'")
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

def parse_ktheta_values(ktheta_str):
    """
    Parse ktheta string into values for pectin, cellulose, xyloglucan.
    Returns dict with keys 'pectin', 'cellulose', 'xyloglucan' and float values or None for defaults.
    """
    if not ktheta_str:
        return {}
    
    # Split and clean values, handling up to 3 values
    parts = ktheta_str.split(',')
    if len(parts) > 3:
        raise ValueError("ktheta accepts at most 3 values: pectin,cellulose,xyloglucan")
    
    # Pad with empty strings if fewer than 3 values provided
    while len(parts) < 3:
        parts.append('')
    
    result = {}
    polymer_names = ['pectin', 'cellulose', 'xyloglucan']
    
    for i, (name, value_str) in enumerate(zip(polymer_names, parts)):
        value_str = value_str.strip()
        if value_str:  # Not empty
            try:
                result[name] = float(value_str)
            except ValueError:
                raise ValueError(f"Invalid ktheta value for {name}: '{value_str}'. Must be a number.")
    
    return result

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

def modify_ktheta_in_itp(itp_path, new_path, ktheta_value=None):
    """
    Modify the k_theta value in an .itp file.
    If ktheta_value is None, copy the file unchanged.
    """
    lines = itp_path.read_text().splitlines()
    out = []
    
    for line in lines:
        if line.startswith('#define k_theta') and ktheta_value is not None:
            # Replace the k_theta value
            out.append(f'#define k_theta {ktheta_value}')
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

def generate_itps(args, out_dir, epsilon_map, ktheta_values=None):
    toppar_dir = out_dir / "toppar_custom"
    ensure_dir(toppar_dir)
    
    # Base file (LJ parameters only)
    scale_epsilon_in_itp(Path('toppar_custom/sudowoodo_base.itp'), toppar_dir / "sudowoodo_base.itp", epsilon_map)
    
    # Polymer-specific files with potential ktheta modification
    itp_files = [
        ('toppar_custom/sudowoodo_xyloglucan.itp', toppar_dir / "sudowoodo_xyloglucan.itp", 'xyloglucan'),
        ('toppar_custom/sudowoodo_pectin.itp', toppar_dir / "sudowoodo_pectin.itp", 'pectin'),
        ('toppar_custom/sudowoodo_cellulose.itp', toppar_dir / "sudowoodo_cellulose.itp", 'cellulose')
    ]
    
    ktheta_values = ktheta_values or {}
    
    for src, dst, polymer_name in itp_files:
        ktheta_value = ktheta_values.get(polymer_name)
        modify_ktheta_in_itp(Path(src), dst, ktheta_value)

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

def write_log(out_dir, seed, args, epsilon_map, ktheta_values=None):
    eps_map_str = ', '.join([f"{k[0]}{k[1]}={v}" for k,v in epsilon_map.items() if k[0]<=k[1]])
    
    log_txt = textwrap.dedent(f"""\
        AFM-Build Sweep Run Log
        ======================
        Output directory: {out_dir}
        Epsilon mapping: {eps_map_str}
        Polymer counts: Xylo={args.nxylo}  Pctn={args.npctn}  Cell={args.ncell}
        Seed used: {seed}
    """)
    
    if ktheta_values:
        ktheta_str = ', '.join([f"{k}={v}" for k, v in ktheta_values.items()])
        log_txt += f"        Custom ktheta values: {ktheta_str}\n"
    
    write_text(out_dir / "afm_build.log", log_txt)

def build_afm_system(seed, out_dir=None, ktheta_str=None):
    """
    Call build_afm_system.py with the given seed inside the output folder.
    """
    print(f"[info] Building afm_system.gro using build_afm_system.py ...")
    builder = Path(__file__).parent / "build_afm_system.py"
    if not builder.exists():
        raise FileNotFoundError(f"Could not find build_afm_system.py in {builder.parent}")

    cmd = ["python", str(builder), "--seed", str(seed)]
    if ktheta_str:
        cmd.extend(["--ktheta", ktheta_str])
    
    subprocess.run(cmd, cwd=out_dir, check=True)

def main():
    args = get_args()
    ensure_dir(args.out)
    seed = args.seed if args.seed is not None else random.randint(1,99999999)
    random.seed(seed)
    
    epsilon_map = parse_epsilon_map(args.epsilon)
    ktheta_values = parse_ktheta_values(args.ktheta) if args.ktheta else {}
    
    write_log(args.out, seed, args, epsilon_map, ktheta_values)
    randomize_structures(seed, args.out)
    generate_topology(args, args.out)
    generate_itps(args, args.out, epsilon_map, ktheta_values)
    write_mdp_files(args, args.out)
    write_run_sh(args, args.out)
    build_afm_system(seed, args.out, args.ktheta)
    print(f"[ok] Setup complete in {args.out} (seed={seed})")

    if ktheta_values:
        print(f"[info] Custom ktheta values used: {ktheta_values}")

if __name__ == "__main__":
    main()