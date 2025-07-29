import re
from subprocess import check_output 
from subprocess import CalledProcessError
import sys

# Unit converters 
Kb      = 0.0019872067          # Boltzmann constant in kcal/mol-K. (Copied from LAMMPS).
GPa     = 6.9476955    # Unit conversion factor... kcal/mol/A^3 * (this constant) ==> GPa
atm     = GPa*9869.2326671      # stress conversion to atm (for LAMMPS).
GPa2atm = 9869.23266716       # x_GPa * GPa2atm = x_atm
NA = 6.02213676e23              # Avogadro's number (1/mol)


def run_bash_cmnd(cmnd_str):

    """ 

    Runs a (bash) shell command - captures and returns any resulting output. 

    Usage: run_bash_cmnd("my command string")

    Notes: Linux wildcards will not work as expected. Use the glob if needed.

    """

    msg = ""

    try:
        msg = check_output(cmnd_str.split())
    except CalledProcessError as err_msg:
        msg = err_msg.output

    return msg.decode('utf-8')

def parse_log_file(log_path):
    stress = {
        'sxx': None,
        'syy': None,
        'szz': None,
        'sxy': None,
        'sxz': None,
        'syz': None
    }
    # LAMMPS tests have units = real

    dat = run_bash_cmnd("awk /Step/{getline;print}" + " " + log_path).split()

    energy = float(dat[4]) # kcal/mol
    volume = float(dat[14]) # A^3
    n_atoms = float(run_bash_cmnd("awk /Loop/{print($(NF-1))}" + " " + log_path)) 


    try:
        # LAMMPS output pressure in units of atm
        stress['sxx'] = float(dat[ 8]) / GPa2atm  # GPa
        stress['syy'] = float(dat[ 9]) / GPa2atm  # GPa
        stress['szz'] = float(dat[10]) / GPa2atm  # GPa
        stress['sxy'] = float(dat[11]) / GPa2atm  # GPa
        stress['sxz'] = float(dat[12]) / GPa2atm  # GPa
        stress['syz'] = float(dat[13]) / GPa2atm  # GPa

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
        f.write(f"{energy:15.6f}\n" if energy is not None else "N/A\n")

        # f.write("sxx (kcal/mol/A^3)\n")
        for key in ['sxx', 'syy', 'szz', 'sxy', 'sxz', 'syz']:
            f.write(f"{stress[key]:15.6f}\n" if stress[key] is not None else "N/A\n")

        # Write forces in scientific notation, full precision
        for fx, fy, fz in forces:
            f.write(f"{fx:15.6f}\n")
            f.write(f"{fy:15.6f}\n")
            f.write(f"{fz:15.6f}\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 generate_compare.py <lammps_output.log> <tray.lammpstrj> <output.txt>")
        sys.exit(1)

    log_file = sys.argv[1]
    traj_file = sys.argv[2]
    output_file = sys.argv[3]

    energy, stress = parse_log_file(log_file)
    forces = parse_traj_file(traj_file)
    write_output(output_file, energy, stress, forces)

    print(f"Output written to {output_file}")
    print(f"Moved {output_file} to the generated_output folder")

if __name__ == "__main__":
    main()
 