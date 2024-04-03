import argparse

parser = argparse.ArgumentParser(description='Makes csv files from ProteinMPNN.fa output ')
parser.add_argument('af2_log_file', type=str, help = "Path to af2.log file")
parser.add_argument('queries_csv_file', type=str, help = "Path to af2.log file")
parser.add_argument('num_pmpnn_seq', type=int, help = "Number of sequences per design")
parser.add_argument('sele_csv_file', type=str, help = "Path to af2.log file")
parser.add_argument('af2_out_folder', type=str, help = "path to af2 generated pdbs")
parser.add_argument('pdb_file', type=str, help = "RMSD will be included as a metrics. Write - otherwise, don't leave empty!")
parser.add_argument('rank_output_csv_file', type=str, help = "where to save")
parser.add_argument('metric', type=str, help = "pLDDT or pTM or RMSD")
parser.add_argument('pymol_pse_file', type=str, help = "Path to pymol session to be created")
parser.add_argument('pymol_best_pse', type=int, help = "Create a pymol session contaning the N best models aligned with the original protein and colored by pLDDT")
parser.add_argument('only_first', type=bool, help = "Only compare the best folding of each sequence generated")

# Parse the arguments
args = parser.parse_args()

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
    try:    
        # Load the two protein structures
        cmd.load(folded_path, "folded")
        # Align the proteins and calculate RMSD from backbone
        rmsd = cmd.align("original and name CA+C+N+O", "folded and name CA+C+N+O")[0]  # cmd.align returns a tuple, RMSD is the first element
        cmd.delete("folded")
        return rmsd
    except Exception as e:
        print("Error while calculating rmsd")
        print(str(e))
    return -1

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

data = []
ranked_data = []
with open(args.af2_log_file,"r") as af2log:
    scores = dict()
    query_found = False
    reranking_found = False
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
                        scores["RMSD"] = "{:.4f}".format(calculate_rmsd(folded_pdb_file))
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
        print(f"Object: {prot}\tpLDDT:{pLDDT_avg}")
cmd.extend('rank_plddt', plddt)

if len(ranked_data) > 0 and args.pymol_best_pse > 0:
    N_best = args.pymol_best_pse if args.pymol_best_pse < len(ranked_data) else len(ranked_data)
    cmd.load(args.pdb_file, "original")
    cmd.color("bluewhite","original")
    for i in range(N_best):
        scores = ranked_data[i]
        name = scores["name"]
        model = scores["model"]
        if "_design_" in name:
            design,sequence = name.split("_design_")[1].split('_')
            short_name = f"d{design}s{sequence}m{model}"
        else:
            sequence = name.split("_unrelaxed_")[0].split("_")[-1]
            short_name = f"s{sequence}m{model}"
        pLDDT_obj = f"{short_name}_pLDDT"
        cmd.load(scores["path"], pLDDT_obj)
        cmd.do(f"rank_plddt {pLDDT_obj}")
        cmd.align(pLDDT_obj, "original")
        mobile_obj = f"{short_name}_mobile"
        cmd.load(scores["path"], mobile_obj)
        cmd.color("hotpink",mobile_obj)
        fixed_sele = scores["fixed"]
        cmd.color("gray80",f"{mobile_obj} and resi {fixed_sele}")
        cmd.align(mobile_obj, "original")
    cmd.save(args.pymol_pse_file)

cmd.quit()