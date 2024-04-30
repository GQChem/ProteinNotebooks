#Copyright © 2024 LOCBP @ University of Zürich
#Distributed under MIT license
"""
This is used to rank output from the log file of alphafold2
The reason we use it is that in the log file you can also see pTM, not stored in the pdb file
"""

import argparse

parser = argparse.ArgumentParser(description='Makes csv files from ProteinMPNN.fa output ')
parser.add_argument('af2_log_file', type=str, help = "Path to af2.log file")
parser.add_argument('queries_csv_file', type=str, help = "Path to queries file containing id, sequence, ...")
parser.add_argument('num_pmpnn_seq', type=int, help = "Number of sequences per design")
parser.add_argument('sele_csv_file', type=str, help = "Path to af2.log file")
parser.add_argument('af2_out_folder', type=str, help = "path to af2 generated pdbs")
parser.add_argument('pdb_file', type=str, help = "")
parser.add_argument('rank_output_csv_file', type=str, help = "where to save")
parser.add_argument('metric', type=str, help = "pLDDT or pTM or RMSD or pLDDT/RMSD")
parser.add_argument('alignment', type=str, help = "align or cealign or super")
parser.add_argument('pymol_pse_file', type=str, help = "Path to pymol session to be created")
parser.add_argument('pymol_best_pse', type=int, help = "Create a pymol session contaning the N best models aligned with the original protein and colored by pLDDT")
parser.add_argument('rank_best_fasta_file', type=str, help = "Create a fasta file with the sequences of the best N best models")
parser.add_argument('only_first', type=bool, help = "Only compare the best folding of each sequence generated")

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

here we merged the two data tables based on id
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
Here we read the log file of AlphaFold

This is a sample relative to one query

...
2024-02-26 21:19:08,458 Query 1/200: cp_2_1C_2N_001_design_29_1 (length 190)
2024-02-26 21:19:09,403 Sleeping for 9s. Reason: PENDING
...
2024-02-26 21:19:53,885 Sleeping for 10s. Reason: RUNNING
2024-02-26 21:20:24,803 Padding length to 200
2024-02-26 21:21:00,991 alphafold2_ptm_model_1_seed_000 recycle=0 pLDDT=64.9 pTM=0.598
2024-02-26 21:21:17,170 alphafold2_ptm_model_1_seed_000 recycle=1 pLDDT=63.3 pTM=0.586 tol=2.17
2024-02-26 21:21:33,528 alphafold2_ptm_model_1_seed_000 recycle=2 pLDDT=70.8 pTM=0.701 tol=2.24
2024-02-26 21:21:50,062 alphafold2_ptm_model_1_seed_000 recycle=3 pLDDT=69.7 pTM=0.694 tol=0.807
2024-02-26 21:21:50,063 alphafold2_ptm_model_1_seed_000 took 85.3s (3 recycles)
...
2024-02-26 21:25:34,990 alphafold2_ptm_model_5_seed_000 recycle=0 pLDDT=63.2 pTM=0.607
2024-02-26 21:25:52,490 alphafold2_ptm_model_5_seed_000 recycle=1 pLDDT=63.2 pTM=0.591 tol=2.44
2024-02-26 21:26:09,967 alphafold2_ptm_model_5_seed_000 recycle=2 pLDDT=66.9 pTM=0.667 tol=1.55
2024-02-26 21:26:27,473 alphafold2_ptm_model_5_seed_000 recycle=3 pLDDT=66.9 pTM=0.678 tol=1.01
2024-02-26 21:26:27,474 alphafold2_ptm_model_5_seed_000 took 70.0s (3 recycles)
2024-02-26 21:26:27,524 reranking models by 'plddt' metric
2024-02-26 21:26:27,525 rank_001_alphafold2_ptm_model_4_seed_000 pLDDT=69.8 pTM=0.686
2024-02-26 21:26:27,526 rank_002_alphafold2_ptm_model_1_seed_000 pLDDT=69.7 pTM=0.694
2024-02-26 21:26:27,528 rank_003_alphafold2_ptm_model_2_seed_000 pLDDT=67.7 pTM=0.678
2024-02-26 21:26:27,529 rank_004_alphafold2_ptm_model_5_seed_000 pLDDT=66.9 pTM=0.678
2024-02-26 21:26:27,530 rank_005_alphafold2_ptm_model_3_seed_000 pLDDT=61.9 pTM=0.513
...

This is the format of each line:
#date time ...
The first two can be used to calculate the time requiered for the folding

Each new protein start with the word "Query" and with the format:
#date time Query n/N NAME (length LENGTH)
from which we can derive the NAME and length of the folded protein

Then after all the calculations it gives a reranking line, and the next line contains data of the best model
#date time rank_00N_..._model_M_seed_x pLDDT=... pTM=..
"""

data = []
ranked_data = []
with open(args.af2_log_file,"r") as af2log:
    scores = dict()
    query_found = False
    reranking_found = False
    time_zero = ""
    for logline in af2log:
        line = logline.strip()
        if line.startswith("Exception"):
            query_found = False
            reranking_found = False
            continue
        if not query_found:
            if "Query" in line:
                comps = line.split(' ') #date time Query n/N NAME (length LENGTH)
                scores["name"] = comps[4]
                scores["length"] = comps[-1][:-1]
                time_zero = comps[1]
                query_found = True
        else:
            if not reranking_found:
                if "reranking" in line:
                    reranking_found = True
            else:
                comps = line.split(' ') #date time rank_00N_..._model_M_seed_x pLDDT=... pTM=..
                folded_pdb = comps[2]
                rank = folded_pdb.split("rank_")[1][2]
                scores["model"] = folded_pdb.split("model_")[1][0]
                scores["pLDDT"] = comps[3].split("=")[1]
                scores["pTM"] = comps[4].split("=")[1]    
                name = scores["name"]
                folded_pdb_file = os.path.join(args.af2_out_folder,f"{name}_unrelaxed_{folded_pdb}"+".pdb")
                if args.pdb_file.endswith(".pdb"):
                    if os.path.exists(folded_pdb_file):
                        rmsd = calculate_rmsd(folded_pdb_file)
                        scores["RMSD"] = "{:.4f}".format(rmsd)
                    else:
                        scores["RMSD"] = "-"                
                #Add data from MPNN
                for k in pmpnn_keys:
                    scores[k] = pmpnn_data[name][k]
                #path at last
                if os.path.exists(folded_pdb_file):
                    scores["path"] = folded_pdb_file
                else:
                    scores["path"] = "-"
                time_end = comps[1]
                zero = sum([int(x)*(60**(2-i)) for i,x in enumerate(time_zero.split(',')[0].split(':'))])
                end = sum([int(x)*(60**(2-i)) for i,x in enumerate(time_end.split(',')[0].split(':'))])
                scores["duration [s]"] = str(end-zero)
                data.append(scores)
                scores = dict()
                if args.only_first or rank == "5":
                    query_found = False
                    reranking_found = False

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
        for i,key in enumerate(keys):
            if i > 0: csv_file.write(",")
            csv_file.write(key)
        for ranked_scores in ranked_data:
            csv_file.write("\n")
            for i,key in enumerate(keys):
                if i > 0: csv_file.write(",")
                csv_file.write(ranked_scores[key])

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
        print(f"Object: {prot}\tpLDDT:{pLDDT_avg:.1f}")
cmd.extend('rank_plddt', plddt)

if len(ranked_data) > 0 and args.pymol_best_pse > 0:
    N_best = args.pymol_best_pse if args.pymol_best_pse < len(ranked_data) else len(ranked_data)
    #Find bests
    best_designs = []
    num_d = 0 #Number of best designs considered
    designs_considered = []
    for scores in ranked_data:
        name = scores["name"]
        undscore_split = name.split("_")
        design = undscore_split[-2]
        sequence = undscore_split[-1]
        if design.isnumeric():
            short_name = f"d{design}s{sequence}"
            if design in designs_considered:
                continue
            else:
                designs_considered.append(design)
                num_d += 1
        else:
            short_name = f"s{sequence}"
            num_d += 1
        best_designs.append((short_name,scores))
    """Saving .fasta file"""
    if args.rank_best_fasta_file != '-':
        for short_name,scores in best_designs:
            cmd.load(scores["path"], short_name)
        cmd.save(args.rank_best_fasta_file)
        cmd.delete("all")
    """Saving PSE file"""
    for short_name,scores in best_designs:
        pLDDT_obj = f"{short_name}_pLDDT"
        cmd.load(scores["path"], pLDDT_obj)
        cmd.do(f"rank_plddt {pLDDT_obj}")
        mobile_obj = f"{short_name}_mobile"
        cmd.load(scores["path"], mobile_obj)
        cmd.color("hotpink",mobile_obj)
        fixed_sele = scores["fixed"]
        if fixed_sele != "":
            cmd.color("gray80",f"{mobile_obj} and resi {fixed_sele}")
            segments_obj = f"{short_name}_segments"
            cmd.load(scores["path"], segments_obj)
            cmd.color("gray80",segments_obj)
            cmd.spectrum(selection=f"{segments_obj} and resi {fixed_sele}")
    cmd.load(args.pdb_file, "original")
    cmd.color("bluewhite","original")
    cmd.alignto("original")
    cmd.save(args.pymol_pse_file,args.alignment)

cmd.quit()