#Copyright © 2024 LOCBP @ University of Zürich
#Distributed under MIT license
import argparse

parser = argparse.ArgumentParser(description='Reads RFD log to create fixed positions')
parser.add_argument('INPAINT_FOLDER', type=str, help="")
parser.add_argument('rfd_log_file', type=str, help = "Path to rfd.log file")
parser.add_argument('rfd_sh_file', type=str, help = "Path to rfd.sh file")
parser.add_argument('CONTIGS', type=str, help = "")
parser.add_argument('INPAINT', type=str, help = "")
parser.add_argument('INPAINT_AUTO_NUM_DESIGNS', type=int, help = "")
parser.add_argument('INPAINT_AUTO_DISTANCE', type=float, help = "")
parser.add_argument('INPAINT_AUTO_MIN_OCCURENCY', type=int, help = "")
parser.add_argument('INPAINT_AUTO_EXCLUDE', type=str, help = "")
parser.add_argument('pdb_file', type=str, help = "Path to input file")
parser.add_argument('inpaint_csv_file', type=str, help = "Path to inpaint.csv file")
parser.add_argument('inpaint_pse_file', type=str, help = "Path to inpaint.pse file")

# Parse the arguments
args = parser.parse_args()

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
import pymol
from pymol import cmd
import numpy as np

design_names = [d[:-4] for d in os.listdir(args.INPAINT_FOLDER) if d.endswith(".pdb")] #exclude pdb extension

fixed = dict()
mobile = dict()
#derive fixed from logfile PKILFEDP-----PLSEDWQ----...
with open(args.rfd_log_file,"r") as rfdlog:
    log_found = False
    seq = ""
    found_name = ""
    for line in rfdlog:
        line = line.strip()
        if not log_found:
            if "Making design" in line:
                found_name = os.path.basename(line.split(" ")[-1])
                log_found = True
                #print("Design: "+found_name)
        if log_found:
            if "Sequence init" in line:
                seq = line.split(" ")[-1]
                fixed[found_name] = [i+1 for i, x in enumerate(seq) if x != "-"]
                mobile[found_name] = [i+1 for i, x in enumerate(seq) if x == "-"]
                log_found = False
                #print("Sequence: "+seq)

#derives all fixed aminoacids in the original structure form contigs
#example 5-10/A22-25/1/A11-18/9 will result in 22,23,24,25,11,12,13,14,15,16,17,18
fixed_original = []
contig_chains = args.CONTIGS.split('/')
for chain in contig_chains:
    if chain[0].isnumeric(): continue
    min,max=chain[1:].split('-')
    for resi in range(int(min),int(max)+1):
        fixed_original.append(resi)
#now we can map fixed in a given design and in the original

#this variable will store in how many designs each residue in the fixed chains is proximal to a designed residue
fixed_occurancy = [0 for r in fixed_original]
try:
    # Initialize PyMOL in headless mode (no GUI)
    pymol.pymol_argv = ['pymol', '-c']  # -q for quiet, -c for no GUI
    pymol.finish_launching()
except Exception as e:
    print("Error while initializing pymol")
    print(str(e))

#Find proximity
print("Processing data..")
for name in design_names:
    print(name)
    design_pdb = os.path.join(args.INPAINT_FOLDER,name+".pdb")
    cmd.load(design_pdb,"design")
    for f_r in fixed[name]:
        fA = f"design and resi {f_r} and name CA"
        fC = f"design and resi {f_r} and name C"
        fN = f"design and resi {f_r} and name N"
        occ_index = fixed[name].index(f_r)
        for m_r in mobile[name]:
            mA = f"design and resi {m_r} and name CA"
            if cmd.get_distance(fA,mA) <= args.INPAINT_AUTO_DISTANCE:
                #check that the side chain of the fixed residue is pointing toward the designed ones
                #for this, considers the angle between CAf-CAm and CAf-CBf
                fA_c = cmd.get_coords(fA)[0]
                mA_c = cmd.get_coords(mA)[0]

                fC_c = cmd.get_coords(fC)[0]
                fN_c = cmd.get_coords(fN)[0]
                fCN_c = 0.5 * (fC_c + fN_c)

                fNC_fA = fA_c - fCN_c
                mA_fA = mA_c - fCN_c

                angle = 180.0 / np.pi * np.arccos(np.dot(fNC_fA,mA_fA)/(np.linalg.norm(fNC_fA)*np.linalg.norm(mA_fA)))

                if angle < 90: ##COULD BE A PARAMETER
                    fixed_occurancy[occ_index] += 1
                    break
    cmd.delete("design")

#Write rfd inpaint flag
inpaint_flag = set()
for i in range(len(fixed_original)):
    if fixed_occurancy[i] >= args.INPAINT_AUTO_MIN_OCCURENCY:
        inpaint_flag.add(fixed_original[i])

if args.INPAINT != "-":
    inpaint_input_chains = args.INPAINT.split('/')
    for chain in inpaint_input_chains:
        if chain[0].isnumeric(): continue
        min,max=chain[1:].split('-')
        for resi in range(int(min),int(max)+1):
            inpaint_flag.add(resi)

if args.INPAINT_AUTO_EXCLUDE != "-":
    to_exclude = sele_to_list(args.INPAINT_AUTO_EXCLUDE)
    for resi in to_exclude: 
        try:
            inpaint_flag.remove()
        except Exception:
            pass

sorted_inpaint_sele = ""
if len(inpaint_flag) > 0: 
    sorted_inpaint_flag = sorted(list(inpaint_flag))
    sorted_inpaint_sele = list_to_sele(sorted_inpaint_flag)
    inpaint_seq = "A"+sorted_inpaint_sele.replace('+','/A')
    print(f"Inpainting: {inpaint_seq}")
    rfd_cmd = ""
    with open(args.rfd_sh_file,'r') as rfd_sh:
        for line in rfd_sh:
            line=line.strip()
            if line != "":
                rfd_cmd = line
    rfd_cmd += f" 'contigmap.inpaint_seq=[{inpaint_seq}]'"
    with open(args.rfd_sh_file,'w') as rfd_sh:
        rfd_sh.write(rfd_cmd)

#Save data
#CSV file
with open(args.inpaint_csv_file,'w') as inpaint_csv:
    inpaint_csv.write("resi,occurancy,inpainted")
    #first, order base on resi
    fixed_touples = [[fixed_original[i],fixed_occurancy[i]] for i in range(len(fixed_original))]
    fixed_ordered = sorted(fixed_touples,key=lambda x: x[0])
    for resi,occurancy in fixed_ordered:
        inpaint_csv.write('\n')
        inpaint_csv.write(str(resi))
        inpaint_csv.write(',')
        inpaint_csv.write(str(occurancy))
        inpaint_csv.write(',')
        if occurancy >= args.INPAINT_AUTO_MIN_OCCURENCY:
            inpaint_csv.write('1')
        else:
            inpaint_csv.write('0')

#PyMol Session


pdb_name = os.path.basename(args.pdb_file)
cmd.load(args.pdb_file, pdb_name)
cmd.load(args.pdb_file, "Occurancy")

fixed_original_selection = list_to_sele(fixed_original)
cmd.select("fixed_residues",f"Occurancy and resi {fixed_original_selection}")
cmd.select("non_fixed_residues",f"Occurancy and not resi {fixed_original_selection}")
cmd.remove("non_fixed_residues")
cmd.delete("non_fixed_residues")
cmd.delete("fixed_residues")

q1 = int(args.INPAINT_AUTO_NUM_DESIGNS*1.0/3)
q2 = int(args.INPAINT_AUTO_NUM_DESIGNS*2.0/3)
quartile1 = [fixed_original[i] for i in range(len(fixed_original)) if fixed_occurancy[i] > 0 and fixed_occurancy[i] <= q1]
quartile2 = [fixed_original[i] for i in range(len(fixed_original)) if fixed_occurancy[i] > q1 and fixed_occurancy[i] <= q2]
quartile3 = [fixed_original[i] for i in range(len(fixed_original)) if fixed_occurancy[i] > q2]
quartile1_sele = list_to_sele(quartile1)
quartile2_sele = list_to_sele(quartile2)
quartile3_sele = list_to_sele(quartile3)
name1 = f"rare_1-{q1}"
name2 = f"moderate_{q1+1}-{q2}"
name3 = f"common_{q2+1}-{args.INPAINT_AUTO_NUM_DESIGNS}"
if len(quartile1) != 0: cmd.select(name1,f"Occurancy and resi {quartile1_sele}")
if len(quartile2) != 0: cmd.select(name2,f"Occurancy and resi {quartile2_sele}")
if len(quartile3) != 0: cmd.select(name3,f"Occurancy and resi {quartile3_sele}")

cmd.color("gray80")
if len(quartile1) != 0: cmd.color("yelloworange",name1)
if len(quartile2) != 0: cmd.color("orange",name2)
if len(quartile3) != 0: cmd.color("red",name3)
if sorted_inpaint_sele != "": 
    cmd.select("inpainted",f"{pdb_name} and resi {sorted_inpaint_sele}")
    cmd.color("gray20","inpainted")


print(f"Saving inpaint to {args.inpaint_pse_file}")
cmd.save(args.inpaint_pse_file)

cmd.quit() #MUST ALWAYS BE AT THE END OF THE SCRIPT