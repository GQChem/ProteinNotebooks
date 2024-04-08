import argparse

parser = argparse.ArgumentParser(description='Makes pse file from af2 output')
parser.add_argument('af2_out_folder', type=str, help = "path to af2 generated pdbs")
parser.add_argument('pymol_pse_file', type=str, help = "Path to pymol session to be created")
parser.add_argument('only_first', type=bool, help = "Only compare the best folding of each sequence generated")

# Parse the arguments
args = parser.parse_args()

#Check prody is part of your environment, otherwise install it beforehand
import os

import pymol
from pymol import cmd


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
        print(f"Object: {prot}\tpLDDT:{pLDDT_avg:.2f}")
cmd.extend('rank_plddt', plddt)


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

pdb_files = [f for f in os.listdir(args.af2_out_folder) if f.endswith(".pdb")]
allow_shortname = True
if len(pdb_files) > 9: #two give 10
    allow_shortname = False

prot_names = []
for pdb in pdb_files:
    longname = pdb.split("_unrelaxed_")[0] if "_unrelaxed_" in pdb else pdb.split("_relaxed_")[0]
    rank = pdb.split("_rank_")[1][2]
    if args.only_first and rank != "1": 
        continue
    model = pdb.split("_model_")[1][0]
    short_name = f"r{rank}m{model}" if allow_shortname else f"{longname}_r{rank}m{model}"
    pdb_file = os.path.join(args.af2_out_folder,pdb)
    cmd.load(pdb_file, short_name)
    cmd.do(f"rank_plddt {short_name}")
    prot_names.append(short_name)

if len(prot_names) > 1:
    for sn in prot_names[1:]:
        cmd.align(sn, prot_names[0])

cmd.save(args.pymol_pse_file)

cmd.quit()