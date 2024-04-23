#Copyright © 2024 LOCBP @ University of Zürich
#Distributed under MIT license
import argparse

parser = argparse.ArgumentParser(description='Makes csv files from ProteinMPNN.fa output ')
parser.add_argument('of_log_file', type=str, help = "Path to of.log file")
parser.add_argument('queries_csv_file', type=str, help = "Path to queries file containing id, sequence, ...")
parser.add_argument('num_pmpnn_seq', type=int, help = "Number of sequences per design")
parser.add_argument('sele_csv_file', type=str, help = "Path to of.log file")
parser.add_argument('of_out_folder', type=str, help = "path to of generated pdbs")
parser.add_argument('pdb_file', type=str, help = "RMSD will be included as a metrics. Write - otherwise, don't leave empty!")
parser.add_argument('rank_output_csv_file', type=str, help = "where to save")
parser.add_argument('metric', type=str, help = "pLDDT or RMSD")
parser.add_argument('alignment', type=str, help = "align or cealign or super")
parser.add_argument('pymol_pse_file', type=str, help = "Path to pymol session to be created")
parser.add_argument('pymol_best_pse', type=int, help = "Create a pymol session contaning the N best models aligned with the original protein and colored by pLDDT")
parser.add_argument('rank_best_fasta_file', type=str, help = "Create a fasta file with the sequences of the best N best models")
parser.add_argument('only_first', type=bool, help = "Only compare the best folding of each sequence generated")

"""
metric cannot be pTM with omegafold because it is not given in output (calculation could be implemented later)
Only_first is not used here because OmegaFold only gives one pdb in output

General warning: In OmegaFold output, first residue is indexed as 0, whereas in AlphaFold it is indexed as 1
"""

# Parse the arguments
args = parser.parse_args()

"""Initialize Pymol"""
#Check prody is part of your environment, otherwise install it beforehand
import os

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

if args.pdb_file.endswith(".pdb"):
    try:    
        cmd.load(args.pdb_file, "original")
    except Exception as e:
        print("Error while loading original")
        print(str(e))

def calculate_plddt(folded_path):
    try:    
        # Load the two protein structures
        cmd.load(folded_path, "folded")
        # Align the proteins and calculate RMSD from backbone
        atom_iterator = cmd.get_model(f"folded and name CA")
        residues_inspected = []
        pLDDT_sum = 0
        for atom in atom_iterator.atom:
            resi = int(atom.resi)
            if resi in residues_inspected: continue
            pLDDT_sum += atom.b
            residues_inspected.append(resi)
        cmd.delete("folded")  
        pLDDT_avg = pLDDT_sum * 1.0 / len(residues_inspected)      
        return pLDDT_avg
    except Exception as e:
        print("Error while calculating plddt")
        print(str(e))
    return -1

def calculate_rmsd(folded_path):
    global args
    try:    
        # Load the two protein structures
        cmd.load(folded_path, "folded")
        # Align the proteins and calculate RMSD
        if args.alignment == 'align':
            rmsd = cmd.align("folded and name CA+C+N+O","original and name CA+C+N+O")[0]  
        elif args.alignment == 'cealign':
            cealign = cmd.cealign("original and name CA+C+N+O", "folded and name CA+C+N+O")
            rmsd = cealign["RMSD"]
        elif args.alignment == 'super':            
            rmsd = cmd.super("folded and name CA+C+N+O","original and name CA+C+N+O")[0]

        cmd.delete("folded")
        return rmsd
    except Exception as e:
        print("Error while calculating rmsd")
        print(str(e))
    return -1

"""
Queries CSV contains the following data: (From Protein MPNN)
id,sequence,T,sample,score,global_score,seq_recovery

PMPNN_sele.csv contains this:
id,fixed,mobile

here we merge the two data tables based on id
"""

#queries_csv_file
#sele_csv_file
pmpnn_data = {}
pmpnn_keys = []
with open(args.sele_csv_file,"r") as sele_csv:
    keys = sele_csv.readline().strip().split(",") #id,fixed,mobile
    non_id_keys = keys[1:]
    for k in non_id_keys:
        pmpnn_keys.append(k)
    for line in sele_csv:
        values = line.strip().split(',')
        id = values[0]
        non_id_values = values[1:]
        for n in range(args.num_pmpnn_seq):
            name = f"{id}_{n+1}"
            pmpnn_data[name] = dict()
            for i in range(len(non_id_keys)):
                pmpnn_data[name][non_id_keys[i]] = non_id_values[i]
with open(args.queries_csv_file,"r") as queries_csv:
    keys = queries_csv.readline().strip().split(",") #id,sequence,T,sample,score,global_score,seq_recovery
    non_id_keys = keys[1:]
    for k in non_id_keys:
        pmpnn_keys.append(k)
    for line in queries_csv:
        values = line.strip().split(',')
        name = values[0]
        non_id_values = values[1:]
        for i in range(len(non_id_keys)):
            pmpnn_data[name][non_id_keys[i]] = non_id_values[i]

"""
Here we read the log file of OmegaFold

This is a sample
INFO:root:Loading weights from /data/gquarg/OmegaFold/release1.pt
INFO:root:Constructing OmegaFold
INFO:root:Reading /scratch/gquarg/RFD_output/test_004/of_queries.fasta
...
INFO:root:Predicting 2th chain in /scratch/gquarg/RFD_output/test_004/of_queries.fasta
INFO:root:190 residues in this chain.
INFO:root:Finished prediction in 19.09 seconds.
INFO:root:Saving prediction to /scratch/gquarg/RFD_output/test_004/test_004_design_1_2.pdb
INFO:root:Saved
...
INFO:root:Done!

Each string starts with INFO:root: so we can discard the first 10 characters

Each prediction gives 5 lines from which we can derive:
Line 2: Length of the protein
Line 3: Time for prediction
Line 4: PDB path and name of folded protein

"""

data = []
ranked_data = []
with open(args.of_log_file,"r") as oflog:
    scores = dict()
    for logline in oflog:
        line = logline.strip()
        line = line[10:] #remove INFO:root:
        if line.startswith("Predicting"):
            scores = dict()
        elif line.endswith("residues in this chain."):
            scores["length"]=line.split(' ')[0]
        elif line.startswith("Finished prediction"):
            scores["duration [s]"] = line.split(' ')[-2]
        elif line.startswith("Saving prediction"):
            path = line.split(' ')[-1]
            name=os.path.basename(path)[:-4] #remove .pdb extension
            scores["name"]=name
            scores["pLDDT"] = "{:.2f}".format(calculate_plddt(path))
            if args.pdb_file.endswith(".pdb"): #a pdb is given as input
                rmsd = calculate_rmsd(path)
                scores["RMSD"] = "{:.4f}".format(rmsd)
            #Add data from MPNN
            for k in pmpnn_keys:
                scores[k] = pmpnn_data[name][k]
            #path at last
            if os.path.exists(path):
                scores["path"] = path
            else:
                scores["path"] = "-"
            data.append(scores)

if len(data) == 0: 
    print("Ranking script couldn't parse output correctly")
else:
    if args.metric in list(data[0].keys()):
        small_to_big = args.metric == "RMSD"
        ranked_data = sorted(data,key=lambda x: float(x[args.metric]),reverse=not small_to_big)
    else:
        print(f"Error while ranking, key {args.metric} not found")
        ranked_data = data
    with open(args.rank_output_csv_file,"w") as csv_file:
        keys = list(ranked_data[0].keys())
        csv_file.write("name")
        for key in keys:
            if key != "name":
                csv_file.write(",")
                csv_file.write(key)
        for ranked_scores in ranked_data:
            csv_file.write("\n")
            csv_file.write(ranked_scores["name"])
            for key in keys:
                if key != "name": 
                    csv_file.write(",")
                    csv_file.write(str(ranked_scores[key]))

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
        print(f"Object: {prot}\tpLDDT: {pLDDT_avg:.2f}")
cmd.extend('rank_plddt', plddt)

if len(ranked_data) > 0 and args.pymol_best_pse > 0:
    N_best = args.pymol_best_pse if args.pymol_best_pse < len(ranked_data) else len(ranked_data)
    num_d = 0 #Number of best designs considered
    best_designs = []
    for i in range(N_best):
        scores = ranked_data[i]
        name = scores["name"]
        undscore_split = name.split("_")
        design = undscore_split[-2]
        sequence = undscore_split[-1]
        if design.isnumeric():
            short_name = f"d{design}s{sequence}"
            if design in best_designs:
                continue
            else:
                best_designs.append(design)
                num_d += 1
        else:
            short_name = f"s{sequence}"
            num_d += 1
        pLDDT_obj = f"{short_name}_pLDDT"
        cmd.load(scores["path"], pLDDT_obj)
        cmd.do(f"rank_plddt {pLDDT_obj}")
        mobile_obj = f"{short_name}_mobile"
        cmd.load(scores["path"], mobile_obj)
        cmd.color("hotpink",mobile_obj)
        fixed_sele = scores["fixed"]
        """
        OmegaFold indeces start from 0, not 1, so here we are shifting the selection accordingly
        """
        plus_parts = fixed_sele.split('+')
        new_pluses = []
        for plus in plus_parts:
            if '-' in plus:
                start, end = plus.split('-')
                adjusted = f"{int(start)-1}-{int(end)-1}"
            else:
                adjusted = str(int(plus) - 1)
            new_pluses.append(adjusted)
        fixed_sele_shifted = '+'.join(new_pluses)
        cmd.color("gray80",f"{mobile_obj} and resi {fixed_sele_shifted}")
        segments_obj = f"{short_name}_segments"
        cmd.load(scores["path"], segments_obj)
        cmd.color("gray80",segments_obj)
        cmd.spectrum(selection=f"{segments_obj} and resi {fixed_sele_shifted}")
    cmd.save(args.rank_best_fasta_file)
    cmd.load(args.pdb_file, "original")
    cmd.color("bluewhite","original")
    cmd.alignto("original",args.alignment)
    cmd.save(args.pymol_pse_file)

cmd.quit()