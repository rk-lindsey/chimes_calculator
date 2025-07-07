import re
from lammps_logfile import File
import sys

def parse_log_file(log_path):
    stress = {
        'sxx': None,
        'syy': None,
        'szz': None,
        'sxy': None,
        'sxz': None,
        'syz': None
    }

    log = File(log_path)
    energy = float(log.get("PotEng")[0])
    volume = float(log.get("Volume")[0])
    temp = float(log.get("Temp")[0])
    n_atoms = float(log.get("Atoms")[0])
    kB = (1.3806488 * 10 **(-23)) * 0.000239006 # kcal/K
    sub_constant = (kB*temp*n_atoms/volume) / 0.14393 * (6.022 * 10 **(23))
    print(sub_constant)

    try:
        stress['sxx'] = float(log.get("Pxx")[0]) * 0.000101325 - sub_constant
        stress['syy'] = float(log.get("Pyy")[0]) * 0.000101325 - sub_constant
        stress['szz'] = float(log.get("Pzz")[0]) * 0.000101325 - sub_constant
        stress['sxy'] = float(log.get("Pxy")[0]) * 0.000101325 - sub_constant
        stress['sxz'] = float(log.get("Pxz")[0]) * 0.000101325 - sub_constant
        stress['syz'] = float(log.get("Pyz")[0]) * 0.000101325 - sub_constant
    except (IndexError, ValueError):
        pass

    return energy, stress

def parse_traj_file(traj_path):
    forces = []
    with open(traj_path, 'r') as f:
        lines = f.readlines()
        reading_atoms = False
        for i, line in enumerate(lines):
            if "ITEM: ATOMS" in line:
                headers = line.strip().split()[2:]
                fx_idx = headers.index('fx') if 'fx' in headers else None
                fy_idx = headers.index('fy') if 'fy' in headers else None
                fz_idx = headers.index('fz') if 'fz' in headers else None
                reading_atoms = True
                continue

            if reading_atoms and fx_idx is not None:
                parts = line.strip().split()
                try:
                    fx = float(parts[fx_idx]) #/ 627.509 * 0.529177249
                    fy = float(parts[fy_idx]) #/ 627.509 * 0.529177249
                    fz = float(parts[fz_idx]) #/ 627.509 * 0.529177249
                    forces.append((fx, fy, fz))
                except (ValueError, IndexError):
                    pass
    return forces

def write_output(filename, energy, stress, forces):
    with open(filename, 'w') as f:
        # f.write("Energy (kcal/mol)\n")
        f.write(f"{energy:.6f}\n" if energy is not None else "N/A\n")

        # f.write("sxx (kcal/mol/A^3)\n")
        for key in ['sxx', 'syy', 'szz', 'sxy', 'sxz', 'syz']:
            f.write(f"{stress[key]:.6f}\n" if stress[key] is not None else "N/A\n")

        for i, (fx, fy, fz) in enumerate(forces):
            f.write(f"{fx:.6f}\n")
            f.write(f"{fy:.6f}\n")
            f.write(f"{fz:.6f}\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 gen_compare.py <lammps_output.log> <tray.lammpstrj> <output.txt>")
        sys.exit(1)

    log_file = sys.argv[1]
    traj_file = sys.argv[2]
    output_file = sys.argv[3]

    energy, stress = parse_log_file(log_file)
    forces = parse_traj_file(traj_file)
    write_output(output_file, energy, stress, forces)

    print(f"Output written to {output_file}")

if __name__ == "__main__":
    main()