from __future__ import print_function
import math
import time
import argparse
import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm

# ----------------------------
# Config (halved box & counts)
# ----------------------------
BOX = np.array([900.0, 450.0, 40.0])  # nm
NUM_CELL = 73
NUM_XYLO = 229
NUM_PCTN = 2750

# LJ-bead "sizes" (used only for overlap checks)
sigma_lookup = {
    "Cell": 3.0,
    "Xylo": 1.5,
    "Pctn": 0.8
}

SAFETY = 1.2            # min-distance factor vs. sigma
TRY_MAIN_CELL = 600     # placement attempts for cellulose (main/off-plane)
TRY_OTHER = 500         # placement attempts for Xylo/Pctn

# ----------------------------
# Utilities
# ----------------------------
def rotate_chain(coords, angle_deg):
    a = math.radians(angle_deg)
    R = np.array([[math.cos(a), -math.sin(a), 0.0],
                  [math.sin(a),  math.cos(a), 0.0],
                  [0.0,          0.0,         1.0]])
    return np.dot(coords, R.T)

def translate_chain(coords, shift):
    return coords + shift

def load_chain_gro(gro_path):
    atoms = []
    with open(gro_path) as f:
        lines = f.readlines()
    for line in lines[2:-1]:
        x = float(line[20:28]); y = float(line[28:36]); z = float(line[36:44])
        atoms.append([x, y, z])
    # last line is box (we ignore template box, we write our own BOX later)
    return np.array(atoms, dtype=float)

def write_combined_gro(output_path, systems, box):
    n_atoms = sum(len(coords) for (_t, coords) in systems)
    with open(output_path, 'w') as f:
        f.write("AFM-Based Multi-Polymer System\n")
        f.write("%5d\n" % n_atoms)
        res_counters = {"Cell": 1, "Xylo": 1, "Pctn": 1}
        prefix = {"Cell": "C", "Xylo": "X", "Pctn": "P"}
        for (chain_type, coords) in systems:
            resname = chain_type[:4].capitalize()
            res_index = res_counters[chain_type]
            for i, (x, y, z) in enumerate(coords):
                atomname = "%s%d" % (prefix[chain_type], i + 1)
                f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %
                        (res_index, resname, atomname, i + 1, x, y, z))
            res_counters[chain_type] += 1
        f.write("%10.5f%10.5f%10.5f\n" % (box[0], box[1], box[2]))

def write_combined_top(output_path, chain_specs):
    with open(output_path, "w") as top:
        top.write(";;;;;; AFM-Based Combined Topology\n;\n")
        top.write("#include \"toppar_custom/sudowoodo_base.itp\"\n")
        included = set()
        for entry in chain_specs:
            chain_type = entry[0]
            if chain_type not in included:
                top.write("#include \"toppar_custom/sudowoodo_%s.itp\"\n" % chain_type.lower())
                included.add(chain_type)
        top.write("\n[ system ]\nAFM-Based Cell Wall System\n\n[molecules]\n")
        counts = {}
        for entry in chain_specs:
            t = entry[0]
            counts[t] = counts.get(t, 0) + 1
        # Keep a stable order
        for t in ["Xylo", "Pctn", "Cell"]:
            if t in counts:
                top.write("%-10s %d\n" % (t, counts[t]))

class SpatialIndex(object):
    def __init__(self):
        self.coords = []
        self.types = []
        self.tree = None
    def _rebuild(self):
        if self.coords:
            self.tree = cKDTree(np.array(self.coords, dtype=float))
        else:
            self.tree = None
    def add_chain(self, coords, ctype):
        self.coords.extend(coords.tolist())
        self.types.extend([ctype] * len(coords))
        self._rebuild()
    def has_overlap(self, coords, ctype):
        if self.tree is None:
            return False
        new_sigma = sigma_lookup.get(ctype, 3.0)
        search_r = max(5.0, SAFETY * new_sigma)
        for atom in coords:
            idxs = self.tree.query_ball_point(atom, r=search_r)
            for idx in idxs:
                other_sigma = sigma_lookup.get(self.types[idx], 3.0)
                min_dist = SAFETY * max(new_sigma, other_sigma)
                if np.linalg.norm(atom - self.tree.data[idx]) < min_dist:
                    return True
        return False

# ----------------------------
# Build
# ----------------------------
def build(seed):
    np.random.seed(seed)
    print("[INFO] Using random seed: %d" % seed)
    print("[INFO] Target box: %.1f x %.1f x %.1f nm" % (BOX[0], BOX[1], BOX[2]))
    print("[INFO] Target counts: Cell=%d, Xylo=%d, Pctn=%d" % (NUM_CELL, NUM_XYLO, NUM_PCTN))

    t0 = time.time()
    print("[STEP] Loading templates: C.gro, X.gro, P.gro ...")
    C_coords = load_chain_gro("C.gro")
    X_coords = load_chain_gro("X.gro")
    P_coords = load_chain_gro("P.gro")

    systems = []
    chain_specs = []
    indexer = SpatialIndex()

    # --- Place cellulose ---
    mid_z = BOX[2] * 0.5
    num_main = int(NUM_CELL * 0.8)
    num_off = NUM_CELL - num_main

    print("[STEP] Placing cellulose (main lamella): %d" % num_main)
    placed_ok = 0
    for _ in tqdm(range(num_main), desc="Cell(main)", ncols=80):
        base = rotate_chain(C_coords.copy(), 0.0)
        success = False
        for _try in range(TRY_MAIN_CELL):
            x = np.random.uniform(0, BOX[0])
            y = np.random.uniform(0, BOX[1])
            z = mid_z
            cand = translate_chain(base, [x, y, z] - base.mean(axis=0))
            if np.any(cand.min(0) < 0) or np.any(cand.max(0) > BOX):
                continue
            if indexer.has_overlap(cand, "Cell"):
                continue
            systems.append(("Cell", cand))
            chain_specs.append(("Cell", "C.gro", x, y, 0.0))
            indexer.add_chain(cand, "Cell")
            success = True
            placed_ok += 1
            break
        if not success:
            print("[WARN] Failed to place a main cellulose fibril after %d tries" % TRY_MAIN_CELL)
    print("[INFO] Placed cellulose (main): %d/%d" % (placed_ok, num_main))

    print("[STEP] Placing cellulose (off-plane): %d" % num_off)
    placed_off = 0
    for _ in tqdm(range(num_off), desc="Cell(off)", ncols=80):
        angle = np.random.choice([-30.0, 30.0])
        zoff = 4.5 if angle > 0 else -4.5
        base = rotate_chain(C_coords.copy(), angle)
        success = False
        for _try in range(TRY_MAIN_CELL):
            x = np.random.uniform(0, BOX[0])
            y = np.random.uniform(0, BOX[1])
            z = mid_z + zoff
            cand = translate_chain(base, [x, y, z] - base.mean(axis=0))
            if np.any(cand.min(0) < 0) or np.any(cand.max(0) > BOX):
                continue
            if indexer.has_overlap(cand, "Cell"):
                continue
            systems.append(("Cell", cand))
            chain_specs.append(("Cell", "C.gro", x, y, angle))
            indexer.add_chain(cand, "Cell")
            success = True
            placed_off += 1
            break
        if not success:
            print("[WARN] Failed to place an off-plane cellulose fibril after %d tries" % TRY_MAIN_CELL)
    print("[INFO] Placed cellulose (off-plane): %d/%d" % (placed_off, num_off))

    # --- Place xylo + pectin ---
    print("[STEP] Placing Xyloglucan (%d) & Pectin (%d) ..." % (NUM_XYLO, NUM_PCTN))
    tail = [("Xylo", "X.gro")] * NUM_XYLO + [("Pctn", "P.gro")] * NUM_PCTN
    # Place longer chains first: assume Xylo generally larger than Pctn
    def vol_like(ctype):
        return sigma_lookup.get(ctype, 3.0)  # simple proxy; templates already sized
    tail.sort(key=lambda t: vol_like(t[0]), reverse=True)

    placed_tail = 0
    for ctype, tpl in tqdm(tail, desc="Xylo+Pctn", ncols=80):
        base = X_coords if ctype == "Xylo" else P_coords
        success = False
        for _try in range(TRY_OTHER):
            ang = np.random.uniform(-45.0, 45.0)
            rot = rotate_chain(base.copy(), ang)
            x = np.random.uniform(0, BOX[0])
            y = np.random.uniform(0, BOX[1])
            z = np.random.uniform(0, BOX[2])
            cand = translate_chain(rot, [x, y, z] - rot.mean(axis=0))
            if np.any(cand.min(0) < 0) or np.any(cand.max(0) > BOX):
                continue
            if indexer.has_overlap(cand, ctype):
                continue
            systems.append((ctype, cand))
            chain_specs.append((ctype, tpl, x, y, ang))
            indexer.add_chain(cand, ctype)
            success = True
            placed_tail += 1
            break
        if not success:
            print("[WARN] Failed to place %s after %d tries" % (ctype, TRY_OTHER))

    # --- Write outputs ---
    print("[STEP] Writing outputs: afm_system.gro / afm_system.top")
    write_combined_gro("afm_system.gro", systems, BOX)
    write_combined_top("afm_system.top", chain_specs)

    # --- Summary ---
    dt = time.time() - t0
    counts = {"Cell":0, "Xylo":0, "Pctn":0}
    for t, _ in systems: counts[t] += 1
    print("[DONE] Build completed in %.1f s" % dt)
    print("[INFO] Final counts: Cell=%d, Xylo=%d, Pctn=%d, TOTAL chains=%d" %
          (counts["Cell"], counts["Xylo"], counts["Pctn"], len(systems)))
    total_atoms = sum(len(c) for _t, c in systems)
    print("[INFO] Total atoms (beads): %d" % total_atoms)

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build AFM-based multi-polymer system")
    parser.add_argument("--seed", type=int, default=None, help="Random seed; default: current time")
    parser.add_argument("--ktheta", type=float, default=None, help="Override ktheta value in polymer .itp files. If provided, updates k_theta values in pectin, cellulose, and xyloglucan .itp files in current directory.")
    args = parser.parse_args()
    seed = args.seed if args.seed is not None else int(time.time())
    
    # Update ktheta in .itp files if requested
    if args.ktheta is not None:
        import re
        itp_files = [
            "toppar_custom/sudowoodo_pectin.itp",
            "toppar_custom/sudowoodo_cellulose.itp", 
            "toppar_custom/sudowoodo_xyloglucan.itp"
        ]
        
        for itp_file in itp_files:
            try:
                # Read the file
                with open(itp_file, 'r') as f:
                    lines = f.readlines()
                
                # Update the ktheta line
                for i, line in enumerate(lines):
                    if line.strip().startswith('#define k_theta'):
                        lines[i] = f"#define k_theta {args.ktheta}\n"
                        break
                
                # Write back to file
                with open(itp_file, 'w') as f:
                    f.writelines(lines)
                    
                print(f"[INFO] Updated ktheta to {args.ktheta} in {itp_file}")
            except FileNotFoundError:
                print(f"[WARN] Could not find {itp_file}, skipping ktheta update")
            except Exception as e:
                print(f"[WARN] Error updating ktheta in {itp_file}: {e}")
    
    build(seed)