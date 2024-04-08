import argparse

parser = argparse.ArgumentParser(description='Copy best pdb from one cycle to the next one in PMPNN_AF2_Loop')
parser.add_argument('old_rank_output_csv_file', type=str)
parser.add_argument('previous_pdb_file', type=str, help="path to previous best pdb")
parser.add_argument('pdb_file', type=str,help = "path to where the best of the previous cycle will be")

# Parse the arguments
args = parser.parse_args()

import pymol
from pymol import cmd

try:
    # Initialize PyMOL in headless mode (no GUI)
    pymol.pymol_argv = ['pymol', '-c']  # -q for quiet, -c for no GUI
    pymol.finish_launching()
    cmd.load(args.previous_pdb_file,"prot")
    atom_iterator = cmd.get_model("prot and name CA")
    residues_inspected = []
    pLDDT_sum = 0
    for atom in atom_iterator.atom:
        resi = int(atom.resi)
        if resi in residues_inspected: continue
        pLDDT_sum += atom.b
        residues_inspected.append(resi)
    cmd.delete("prot")
    previous_pLDDT = pLDDT_sum * 1.0 / len(residues_inspected)
except Exception as e:
    print("Error while determining pLDDT of last input")
    print(str(e))
    previous_pLDDT = 0
print(f"Previous pLDDT: {previous_pLDDT:.1f}")

import shutil

#retrieve path of best one
with open(args.old_rank_output_csv_file,"r") as ranked_csv:
    data_keys = ranked_csv.readline().strip().split(',')
    pLDDT_index = data_keys.index("pLDDT")
    seq_index = data_keys.index("sequence")
    path_index = data_keys.index("path")
    best_data = ranked_csv.readline().strip().split(',')
    best_pLDDT = float(best_data[pLDDT_index])
    if best_pLDDT > previous_pLDDT:
        best_path = best_data[path_index]
        best_seq = best_data[seq_index]
        print(f"Best PDB: {best_path}")
        print(f"Sequence: {best_seq}")
        print(f"pLDDT: {best_pLDDT:.1f}")
        shutil.copy(best_path,args.pdb_file)  
    else:
        print(f"No improvement in pLDDT. Best is {best_pLDDT:.1f}")
        shutil.copy(args.previous_pdb_file,args.pdb_file)  

cmd.quit() #MUST ALWAYS BE AT THE END