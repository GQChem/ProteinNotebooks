#Copyright © 2024 LOCBP @ University of Zürich
#Distributed under MIT license
"""
This script makes a PSE with plddt and a fasta file containing the output of the series of cycles
"""

import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('original_pdb_file', type=str, help = "")
parser.add_argument('folder', type=str, help = "OmegaFold or AlphaFold. The indeces of residues are different between the two")
parser.add_argument('best_pdb_file', type=str, help = "")
parser.add_argument('best_pse_file', type=str, help = "Create a pymol session contaning the N best models aligned with the original protein and colored by pLDDT")
parser.add_argument('best_fasta_file', type=str, help = "Create a fasta file with the sequences of the best N best models")
# Parse the arguments
args = parser.parse_args()

original_pdb_file = args.original_pdb_file
folder = args.folder
best_pdb_file = args.best_pdb_file
best_pse_file = args.best_pse_file
best_fasta_file = args.best_fasta_file

"""Initialize Pymol"""
import pymol
from pymol import cmd
try:
    # Initialize PyMOL in headless mode (no GUI)
    pymol.pymol_argv = ['pymol', '-c']  # -q for quiet, -c for no GUI
    pymol.finish_launching()
    #Some settings for the session to have good pictures
    cmd.do("show cartoon")
    cmd.set("seq_view", 1)
    cmd.set("cartoon_gap_cutoff", 0)
    cmd.set("sphere_scale", 0.2)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("spec_reflect", 0)
    cmd.set("ray_trace_frames", 1)
    cmd.set("ray_trace_color", "gray20")
except Exception as e:
    print("Error while initializing pymol")
    print(str(e))

"""
PYMOL PSE CREATION
"""
def plddt(selection="all"):     
    blue_rgb = [0,76,202]
    blue = []
    for c in blue_rgb:
        blue.append(c/255.0)
    lightblue_rgb = [73, 196, 238]
    lightblue = []
    for c in lightblue_rgb:
        lightblue.append(c/255.0)
    yellow_rgb = [255, 213, 57]
    yellow = []
    for c in yellow_rgb:
        yellow.append(c/255.0)
    orange_rgb = [255, 113, 67]
    orange = []
    for c in orange_rgb:
        orange.append(c/255.0)     
    cmd.set_color('blue_plddt', blue)
    cmd.set_color('lightblue_plddt', lightblue)
    cmd.set_color('yellow_plddt', yellow)
    cmd.set_color('orange_plddt', orange)
    #select and color blue
    blue_upper = 100.0
    blue_lower = 90.0
    blue_sel_str = selection + " & ! b < " + str(blue_lower) + " & ! b > " + str(blue_upper)
    cmd.color('blue_plddt', blue_sel_str)
    #select and color lightblue
    lightblue_upper = 90.0
    lightblue_lower = 70.0
    lightblue_sel_str = selection + " & ! b < " + str(lightblue_lower) + " & ! b > " + str(lightblue_upper)
    cmd.color('lightblue_plddt', lightblue_sel_str)
    #select and color yellow
    yellow_upper = 70.0
    yellow_lower = 50.0
    yellow_sel_str = selection + " & ! b < " + str(yellow_lower) + " & ! b > " + str(yellow_upper)
    cmd.color('yellow_plddt', yellow_sel_str)
    #select and color orange
    orange_upper = 50.0
    orange_lower = 0.0
    orange_sel_str = selection + " & ! b < " + str(orange_lower) + " & ! b > " + str(orange_upper)
    cmd.color('orange_plddt', orange_sel_str)

    prots = []
    if selection != "all": prots = [selection]
    else:
        for name in cmd.get_names():
            prots.append(name)
    for prot in prots:
        atom_iterator = cmd.get_model(f"{prot} and name CA")
        residues_inspected = []
        pLDDT_sum = 0
        for atom in atom_iterator.atom:
            resi = int(atom.resi)
            if resi in residues_inspected: continue
            pLDDT_sum += atom.b
            residues_inspected.append(resi)
        pLDDT_avg = pLDDT_sum * 1.0 / len(residues_inspected)
        print(f"Object: {prot}\tpLDDT:{pLDDT_avg}")
cmd.extend('rank_plddt', plddt)

#Converts a pymol selection into an array
def sele_to_list(s):
    a = []
    if s == "": return a
    elif '+' in s:
        plus_parts = s.split('+')
        for pp in plus_parts:
            if '-' in pp:
                min,max = pp.split('-')
                for ri in range(int(min),int(max)+1):
                    a.append(ri)
            else:
                a.append(int(pp))
    else:
        if '-' in s:
            min,max = s.split('-')
            for ri in range(int(min),int(max)+1):
                a.append(ri)
        else:
            a.append(int(s))     
    return a
def list_to_sele(a):
    s = ""
    i = 0
    while i < len(a):   
        if i > 0: s += "+"
        s += f"{a[i]}"   
        #represent consecutive indeces with a dash
        if i < len(a) - 1:
            if int(a[i])+1 == int(a[i+1]):
                s += "-"
                j = i + 2
                while j < len(a):
                    if int(a[j]) != int(a[j-1])+1: break
                    j += 1
                i = j - 1
                s += f"{a[i]}" 
        i += 1        
    return s

import os
from pymol import stored

original_name = os.path.basename(original_pdb_file)[:-4]
best_name = os.path.basename(best_pdb_file)[:-4]
if original_name == best_name: original_name += "_original"
best_name_mutations = f"{best_name}_mutations"
best_name_plddt = f"{best_name}_plddt"

cmd.load(best_pdb_file,best_name)
cmd.save(best_fasta_file)
cmd.delete(best_name)

cmd.load(original_pdb_file,original_name)
cmd.load(best_pdb_file,best_name_mutations)
cmd.load(best_pdb_file,best_name_plddt)

original_sequence = []
best_sequence = []

atom_iterator = cmd.get_model(f"{original_name} and name CA")
residues_inspected = []
for atom in atom_iterator.atom:
    if atom.resi in residues_inspected: continue
    residues_inspected.append(atom.resi)
    original_sequence.append(atom.resn)
atom_iterator = cmd.get_model(f"{best_name_mutations} and name CA")
residues_inspected = []
for atom in atom_iterator.atom:
    if atom.resi in residues_inspected: continue
    residues_inspected.append(atom.resi)
    best_sequence.append(atom.resn)

unmutated, mutated = [],[]
shift = 1 if folder == "AlphaFold" else 0
for resi in range(shift,shift+len(best_sequence)):
    best_resn = best_sequence[resi-shift]
    original_resn = original_sequence[resi-shift]
    if best_resn == original_resn: unmutated.append(resi)
    else: mutated.append(resi)

unmutated_sele = list_to_sele(unmutated)
mutated_sele = list_to_sele(mutated)

cmd.do(f"rank_plddt {original_name}")
cmd.do(f"rank_plddt {best_name_plddt}")
cmd.select("unmutated",f"{best_name_mutations} and resi {unmutated_sele}")
cmd.select("mutated",f"{best_name_mutations} and resi {mutated_sele}")
cmd.color("gray70",best_name_mutations)
cmd.color("hotpink","mutated")
cmd.alignto(original_name)

cmd.save(best_pse_file)

cmd.quit()