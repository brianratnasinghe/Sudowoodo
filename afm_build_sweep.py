#!/usr/bin/env python3
"""
AFM cell wall builder tool with custom epsilon mapping and ktheta modifications.
- Copies template .gro and .itp files.
- Modifies sudowoodo_base.itp with user-specified epsilon for bead pairs.
- Modifies polymer .itp files with user-specified ktheta values.
- Generates topology, .mdp files, run.sh, and afm_build.log.
- Calls build_afm_system.py to create afm_system.gro in the output folder.

Usage:
  python afm_build_sweep.py --out run_$(date +%s) --epsilon CC=1.0,CX=0.8,CP=0.7,XX=0.6,XP=0.5,PP=0.4 --epsilon-pr -0.5 --epsilon-pn 2.0 --epsilon-pc 5.0
Optional:
  --seed 123456
  --ktheta "120,150,180"    # pectin,cellulose,xyloglucan
  --ktheta ",150,180"       # use default pectin, custom cellulose and xyloglucan
  --ktheta "120,,"          # custom pectin, default cellulose and xyloglucan
  --multilayer              # Generate 4-layer fiber system with 180° rotation between layers
"""

import argparse, shutil, os, re, random, textwrap, subprocess
from pathlib import Path

PECTIN_NEUTRAL_TYPE = "PN"
PECTIN_REPULSIVE_TYPE = "PR"
PECTIN_CROSSLINK_TYPE = "PC"
DEFAULT_PECTIN_NEUTRAL_EPSILON = 2.0
DEFAULT_PECTIN_REPULSIVE_EPSILON = -0.5
DEFAULT_PECTIN_CROSSLINK_EPSILON = 5.0
DEFAULT_PECTIN_OTHER_PAIR_EPSILON = 2.0
# GRO files use fixed-width fields; these 0-based slice bounds parse the residue number at 0:5 and residue name at 5:10.
GRO_RESIDUE_NUMBER_START = 0
GRO_RESIDUE_NUMBER_END = 5
GRO_RESIDUE_NAME_START = 5
GRO_RESIDUE_NAME_END = 10
MIN_GRO_ATOM_LINE_LENGTH = GRO_RESIDUE_NAME_END
# Pectin variant atomtypes reuse the base pectin atomtype metadata: atomic number, mass (amu), charge (e), particle type, sigma, and epsilon.
PECTIN_ATOMTYPE_ATOMIC_NUMBER = 1
PECTIN_ATOMTYPE_MASS = 26.6
PECTIN_ATOMTYPE_CHARGE = 0.000
PECTIN_ATOMTYPE_PARTICLE_TYPE = "A"
PECTIN_ATOMTYPE_SIGMA = 0.0
PECTIN_ATOMTYPE_EPSILON = 0.0

def get_args():
    p = argparse.ArgumentParser(description="AFM cell wall builder and sweep tool (custom epsilon mapping)")
    p.add_argument('--out', type=Path, required=True, help="Output folder")
    p.add_argument('--epsilon', type=str, required=True,
                   help="Comma-separated epsilon mapping, e.g. CC=1.0,CX=0.8,CP=0.7,XX=0.6,XP=0.5,PP=0.4")
    p.add_argument('--epsilon-pr', type=float, default=DEFAULT_PECTIN_REPULSIVE_EPSILON,
                   help="Exact epsilon for the PN/PR pectin pair.")
    p.add_argument('--epsilon-pn', type=float, default=DEFAULT_PECTIN_NEUTRAL_EPSILON,
                   help="Exact epsilon for the PN/PN pectin pair.")
    p.add_argument('--epsilon-pc', type=float, default=DEFAULT_PECTIN_CROSSLINK_EPSILON,
                   help="Exact epsilon for the PN/PC pectin pair.")
    p.add_argument('--ktheta', type=str, 
                   help="Comma-separated ktheta values for pectin,cellulose,xyloglucan. "
                        "Use empty values to keep defaults, e.g. '120,150,180' or ',150,180' or '120,,'")
    p.add_argument('--seed', type=int, help="Random seed (int). If not set, random seed is chosen and logged.")
    p.add_argument('--multilayer', action='store_true',
                   help="Generate 4-layer fiber system. Each layer is rotated 180° relative to the previous layer.")
    p.add_argument('--nxylo', type=int, default=458)
    p.add_argument('--npctn', type=int, default=5501)
    p.add_argument('--ncell', type=int, default=146)
    p.add_argument('--gmx', type=str, default="gmx")
    p.add_argument('--ntomp', type=int, default=24)
    p.add_argument('--ntmpi', type=int, default=1)
    p.add_argument('--deform', type=str, default=None,
                   help="Deform tensor for production.mdp (e.g. '0 0 0.0001 0 0 0'). "
                        "Written as 'deform = <value>' in production.mdp when provided.")
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
    Gracefully handles extra values by ignoring them.
    """
    if not ktheta_str or not ktheta_str.strip():
        return {}
    
    # Split and clean values, handling up to 3 values
    parts = ktheta_str.split(',')
    
    # Warn about extra values but don't error
    if len(parts) > 3:
        print(f"[warning] ktheta has {len(parts)} values but only first 3 (pectin,cellulose,xyloglucan) will be used")
    
    # Limit to first 3 parts, ignore any excess
    parts = parts[:3]
    
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

def build_pectin_variant_epsilon_map(epsilon_pn, epsilon_pr, epsilon_pc):
    mapping = {
        (PECTIN_NEUTRAL_TYPE, PECTIN_NEUTRAL_TYPE): epsilon_pn,
        (PECTIN_REPULSIVE_TYPE, PECTIN_REPULSIVE_TYPE): DEFAULT_PECTIN_OTHER_PAIR_EPSILON,
        (PECTIN_CROSSLINK_TYPE, PECTIN_CROSSLINK_TYPE): DEFAULT_PECTIN_OTHER_PAIR_EPSILON,
        (PECTIN_NEUTRAL_TYPE, PECTIN_REPULSIVE_TYPE): epsilon_pr,
        (PECTIN_NEUTRAL_TYPE, PECTIN_CROSSLINK_TYPE): epsilon_pc,
        (PECTIN_REPULSIVE_TYPE, PECTIN_CROSSLINK_TYPE): DEFAULT_PECTIN_OTHER_PAIR_EPSILON,
    }
    for left, right in list(mapping):
        mapping[(right, left)] = mapping[(left, right)]
    return mapping

def build_pectin_nonbond_lines(cp_sigma, cp_epsilon, xp_sigma, xp_epsilon, pp_sigma, pectin_variant_epsilons):
    if None in (cp_sigma, cp_epsilon, xp_sigma, xp_epsilon, pp_sigma):
        raise ValueError("Could not derive all required C/X/P nonbond parameters from the base ITP file")
    out = []
    for base_type, sigma, epsilon in [("C", cp_sigma, cp_epsilon), ("X", xp_sigma, xp_epsilon)]:
        for pectin_type in [PECTIN_NEUTRAL_TYPE, PECTIN_REPULSIVE_TYPE, PECTIN_CROSSLINK_TYPE]:
            out.append(f"{base_type} {pectin_type} 1 {sigma:.6f} {epsilon:.6f}")
    for left, right in [
        (PECTIN_NEUTRAL_TYPE, PECTIN_NEUTRAL_TYPE),
        (PECTIN_REPULSIVE_TYPE, PECTIN_REPULSIVE_TYPE),
        (PECTIN_CROSSLINK_TYPE, PECTIN_CROSSLINK_TYPE),
        (PECTIN_NEUTRAL_TYPE, PECTIN_REPULSIVE_TYPE),
        (PECTIN_NEUTRAL_TYPE, PECTIN_CROSSLINK_TYPE),
        (PECTIN_REPULSIVE_TYPE, PECTIN_CROSSLINK_TYPE),
    ]:
        out.append(f"{left} {right} 1 {pp_sigma:.6f} {pectin_variant_epsilons[(left, right)]:.6f}")
    return out

def scale_epsilon_in_itp(itp_path, new_path, epsilon_map, pectin_variant_epsilons):
    re_lj = re.compile(r'^(\s*\w+\s+\w+\s+\d+\s+([0-9eE\.\+\-]+)\s+([0-9eE\.\+\-]+))')
    lines = itp_path.read_text().splitlines()
    out = []
    in_nb = False
    in_atomtypes = False
    cp_sigma = cp_epsilon = None
    xp_sigma = xp_epsilon = None
    pp_sigma = None
    added_pectin_atomtypes = False
    added_pectin_nonbond = False

    def _add_pectin_nonbond_params_once():
        nonlocal added_pectin_nonbond
        if added_pectin_nonbond:
            return
        out.extend(build_pectin_nonbond_lines(
            cp_sigma,
            cp_epsilon,
            xp_sigma,
            xp_epsilon,
            pp_sigma,
            pectin_variant_epsilons,
        ))
        added_pectin_nonbond = True
    for line in lines:
        if line.strip().startswith("[ atomtypes"):
            in_atomtypes = True
            out.append(line)
            continue
        if in_atomtypes and line.strip().startswith("[") and not line.strip().startswith("[ atomtypes"):
            in_atomtypes = False
        parts = line.split()
        if in_atomtypes and parts and parts[0] == "P" and not added_pectin_atomtypes:
            out.append(line)
            for pectin_type in [PECTIN_NEUTRAL_TYPE, PECTIN_REPULSIVE_TYPE, PECTIN_CROSSLINK_TYPE]:
                out.append(f"{pectin_type} {PECTIN_ATOMTYPE_ATOMIC_NUMBER} {PECTIN_ATOMTYPE_MASS} {PECTIN_ATOMTYPE_CHARGE:.3f} {PECTIN_ATOMTYPE_PARTICLE_TYPE} {PECTIN_ATOMTYPE_SIGMA:.1f} {PECTIN_ATOMTYPE_EPSILON:.1f}")
            added_pectin_atomtypes = True
            continue
        if '[ nonbond_params' in line:
            in_nb = True
            out.append(line)
            continue
        if in_nb and line.strip().startswith('['):
            in_nb = False
            _add_pectin_nonbond_params_once()
        if in_nb and re_lj.match(line):
            parts = line.split()
            i, j = parts[0], parts[1]
            sigma = float(parts[3])
            epsilon = float(parts[4])
            new_epsilon = epsilon_map.get((i, j), epsilon)
            parts[4] = f"{new_epsilon:.6f}"
            if (i, j) in [("C", "P"), ("P", "C")]:
                cp_sigma, cp_epsilon = sigma, new_epsilon
            elif (i, j) in [("X", "P"), ("P", "X")]:
                xp_sigma, xp_epsilon = sigma, new_epsilon
            elif i == "P" and j == "P":
                pp_sigma = sigma
            out.append(' '.join(parts))
        else:
            out.append(line)
    if in_nb:
        _add_pectin_nonbond_params_once()
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

def _choose_distributed_positions(total_beads):
    if total_beads < 4:
        raise ValueError(f"Pectin template must have at least 4 beads; found {total_beads}")
    positions = []
    for quarter in range(4):
        start = (quarter * total_beads) // 4
        end = ((quarter + 1) * total_beads) // 4
        positions.append(random.randrange(start, end))
    neg_positions = set(random.sample(positions, 2))
    pos_positions = set(positions) - neg_positions
    return neg_positions, pos_positions

def _write_randomized_pectin_itp(src_path, dst_path, molecule_name, ktheta_value=None):
    lines = src_path.read_text().splitlines()
    out = []
    in_moleculetype = False
    in_atoms = False
    atom_line_indices = []
    atom_ids = []

    for line in lines:
        if line.startswith('#define k_theta') and ktheta_value is not None:
            out.append(f'#define k_theta {ktheta_value}')
            continue
        if line.strip().startswith('[moleculetype]'):
            in_moleculetype = True
            out.append(line)
            continue
        if in_moleculetype and line.strip() and not line.strip().startswith(';'):
            parts = line.split()
            parts[0] = molecule_name
            out.append("  " + " ".join(parts))
            in_moleculetype = False
            continue
        if line.strip().startswith('[atoms]'):
            in_atoms = True
            out.append(line)
            continue
        if in_atoms and line.strip().startswith('['):
            in_atoms = False
        if in_atoms and line.strip() and not line.strip().startswith(';'):
            atom_line_indices.append(len(out))
            atom_ids.append(line.split()[0])
        out.append(line)

    neg_positions, pos_positions = _choose_distributed_positions(len(atom_ids))
    for idx, out_idx in enumerate(atom_line_indices):
        parts = out[out_idx].split()
        if idx in neg_positions:
            bead_type = PECTIN_REPULSIVE_TYPE
        elif idx in pos_positions:
            bead_type = PECTIN_CROSSLINK_TYPE
        else:
            bead_type = PECTIN_NEUTRAL_TYPE
        parts[1] = bead_type
        parts[4] = f"{bead_type}{atom_ids[idx]}"
        out[out_idx] = "  " + " ".join(parts)

    dst_path.write_text('\n'.join(out) + '\n')

def count_pectin_fibers_from_gro(gro_path):
    pectin_residues = set()
    lines = gro_path.read_text().splitlines()
    if len(lines) < 3:
        return 0
    for line in lines[2:-1]:
        if len(line) < MIN_GRO_ATOM_LINE_LENGTH:
            continue
        try:
            residue_number = int(line[GRO_RESIDUE_NUMBER_START:GRO_RESIDUE_NUMBER_END])
        except ValueError:
            continue
        residue_name = line[GRO_RESIDUE_NAME_START:GRO_RESIDUE_NAME_END].strip()
        if residue_name == "Pctn":
            pectin_residues.add(residue_number)
    return len(pectin_residues)

def _get_atom_types_from_itp(itp_path):
    """
    Parse the [atoms] section of a per-fiber pectin ITP and return a dict
    mapping 1-based atom id (int) to atom type string (e.g., 'PN', 'PR', 'PC').
    """
    atom_types = {}
    in_atoms = False
    for line in Path(itp_path).read_text().splitlines():
        if line.strip().startswith('[atoms]'):
            in_atoms = True
            continue
        if in_atoms and line.strip().startswith('['):
            in_atoms = False
        if in_atoms and line.strip() and not line.strip().startswith(';'):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    atom_types[int(parts[0])] = parts[1]
                except ValueError:
                    pass
    return atom_types


def update_gro_pectin_atomnames(gro_path, toppar_dir):
    """
    Rewrite pectin atom names in the GRO file to reflect bead types from the
    per-fiber ITP files generated by generate_itps.  For each Pctn residue N
    the atom name is changed from 'Pi' to '{type}i', e.g. 'P1' -> 'PN1',
    'P7' -> 'PR7', so that visualization and analysis tools can distinguish
    between neutral (PN), repulsive (PR), and crosslink (PC) beads.
    """
    gro_path = Path(gro_path)
    toppar_dir = Path(toppar_dir)
    lines = gro_path.read_text().splitlines()
    if len(lines) < 3:
        return

    itp_cache = {}

    def _types_for_residue(res_num):
        if res_num not in itp_cache:
            itp_path = toppar_dir / f"sudowoodo_pectin_{res_num}.itp"
            itp_cache[res_num] = _get_atom_types_from_itp(itp_path) if itp_path.exists() else {}
        return itp_cache[res_num]

    out = [lines[0], lines[1]]
    for line in lines[2:-1]:
        if len(line) < MIN_GRO_ATOM_LINE_LENGTH:
            out.append(line)
            continue
        residue_name = line[GRO_RESIDUE_NAME_START:GRO_RESIDUE_NAME_END].strip()
        if residue_name != "Pctn":
            out.append(line)
            continue
        try:
            residue_number = int(line[GRO_RESIDUE_NUMBER_START:GRO_RESIDUE_NUMBER_END])
        except ValueError:
            out.append(line)
            continue
        # GRO atom name occupies positions 10:15 (5 chars, right-justified)
        atom_name = line[10:15].strip()
        if not atom_name.startswith('P'):
            out.append(line)
            continue
        try:
            atom_idx = int(atom_name[1:])  # strip leading 'P', parse index
        except (ValueError, IndexError):
            out.append(line)
            continue
        bead_type = _types_for_residue(residue_number).get(atom_idx)
        if bead_type:
            new_name = f"{bead_type}{atom_idx}"
            line = line[:10] + f"{new_name:>5}" + line[15:]
        out.append(line)
    out.append(lines[-1])
    gro_path.write_text('\n'.join(out) + '\n')


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

def generate_itps(args, out_dir, epsilon_map, pectin_count, ktheta_values=None):
    if not isinstance(pectin_count, int) or pectin_count <= 0:
        raise ValueError(f"pectin_count must be a positive integer representing the number of pectin fibers, got {pectin_count!r}")

    toppar_dir = out_dir / "toppar_custom"
    ensure_dir(toppar_dir)
    
    # Base file (LJ parameters only)
    pectin_variant_epsilons = build_pectin_variant_epsilon_map(
        epsilon_pn=args.epsilon_pn,
        epsilon_pr=args.epsilon_pr,
        epsilon_pc=args.epsilon_pc,
    )
    scale_epsilon_in_itp(
        Path('toppar_custom/sudowoodo_base.itp'),
        toppar_dir / "sudowoodo_base.itp",
        epsilon_map,
        pectin_variant_epsilons,
    )
    
    # Polymer-specific files with potential ktheta modification
    itp_files = [
        ('toppar_custom/sudowoodo_xyloglucan.itp', toppar_dir / "sudowoodo_xyloglucan.itp", 'xyloglucan'),
        ('toppar_custom/sudowoodo_cellulose.itp', toppar_dir / "sudowoodo_cellulose.itp", 'cellulose')
    ]
    
    ktheta_values = ktheta_values or {}
    
    for src, dst, polymer_name in itp_files:
        ktheta_value = ktheta_values.get(polymer_name)
        modify_ktheta_in_itp(Path(src), dst, ktheta_value)

    pectin_template = Path('toppar_custom/sudowoodo_pectin.itp')
    pectin_ktheta = ktheta_values.get('pectin')
    for pectin_idx in range(1, pectin_count + 1):
        _write_randomized_pectin_itp(
            pectin_template,
            toppar_dir / f"sudowoodo_pectin_{pectin_idx}.itp",
            f"Pctn_{pectin_idx}",
            pectin_ktheta
        )

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
    def mdp_default_prod(deform=None):
        txt = textwrap.dedent("""\
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
        if deform is not None:
            txt += f"deform                   = {deform}\n"
        return txt
    write_text(out_dir / "EM.mdp", mdp_default_em())
    write_text(out_dir / "EQ.mdp", mdp_default_eq())
    deform = getattr(args, 'deform', None)
    write_text(out_dir / "production.mdp", mdp_default_prod(deform))

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
        Pectin variant epsilons: PN/PR={args.epsilon_pr} PN/PN={args.epsilon_pn} PN/PC={args.epsilon_pc} other pairs={DEFAULT_PECTIN_OTHER_PAIR_EPSILON}
        Polymer counts: Xylo={args.nxylo}  Pctn={args.npctn}  Cell={args.ncell}
        Seed used: {seed}
    """)
    
    if ktheta_values:
        ktheta_str = ', '.join([f"{k}={v}" for k, v in ktheta_values.items()])
        log_txt += f"        Custom ktheta values: {ktheta_str}\n"
    
    if hasattr(args, 'multilayer') and args.multilayer:
        log_txt += "        Multi-layer mode: enabled (4 layers)\n"
    
    if hasattr(args, 'deform') and args.deform is not None:
        log_txt += f"        Deform settings: {args.deform}\n"
    
    write_text(out_dir / "afm_build.log", log_txt)

def build_afm_system(seed, out_dir=None, ktheta_str=None, multilayer=False):
    """
    Call build_afm_system.py with the given seed inside the output folder.
    The builder writes both afm_system.gro and afm_system.top.
    """
    print(f"[info] Building afm_system.gro using build_afm_system.py ...")
    builder = Path(__file__).parent / "build_afm_system.py"
    if not builder.exists():
        raise FileNotFoundError(f"Could not find build_afm_system.py in {builder.parent}")

    cmd = ["python", str(builder), "--seed", str(seed)]
    if ktheta_str:
        cmd.extend(["--ktheta", ktheta_str])
    if multilayer:
        cmd.append("--multilayer")
    
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
    write_mdp_files(args, args.out)
    write_run_sh(args, args.out)
    build_afm_system(seed, args.out, args.ktheta, args.multilayer)
    # build_afm_system.py writes afm_system.gro, which determines how many per-fiber pectin ITPs are needed.
    pectin_count = count_pectin_fibers_from_gro(args.out / "afm_system.gro")
    generate_itps(args, args.out, epsilon_map, pectin_count, ktheta_values)
    update_gro_pectin_atomnames(args.out / "afm_system.gro", args.out / "toppar_custom")
    print(f"[ok] Setup complete in {args.out} (seed={seed})")

    if ktheta_values:
        print(f"[info] Custom ktheta values used: {ktheta_values}")
    
    if args.multilayer:
        print("[info] Multi-layer mode enabled (4 layers)")

if __name__ == "__main__":
    main()