import argparse

parser = argparse.ArgumentParser(description='Reads RFD log to create fixed positions')
parser.add_argument('JOB_FOLDER', type=str, help="")
parser.add_argument('INPAINT_FOLDER', type=str, help="")
parser.add_argument('rfd_log_file', type=str, help = "Path to rfd.log file")
parser.add_argument('rfd_sh_file', type=str, help = "Path to rfd.sh file")
parser.add_argument('CONTIGS', type=str, help = "")
parser.add_argument('INPAINT', type=str, help = "")
parser.add_argument('INPAINT_AUTO_NUM_DESIGNS', type=str, help = "")
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
                print("Design: "+found_name)
        if log_found:
            if "Sequence init" in line:
                seq = line.split(" ")[-1]
                fixed[found_name] = [i+1 for i, x in enumerate(seq) if x != "-"]
                mobile[found_name] = [i+1 for i, x in enumerate(seq) if x == "-"]
                log_found = False
                print("Sequence: "+seq)

#derives all fixed aminoacids in the original structure form contigs
#example 5-10/A11-18/1/A22-25/9 will result in 11,12,13,14,15,16,17,18,22,23,24,25
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
for name in design_names:
    design_pdb = os.path.join(args.INPAINT_FOLDER,name+".pdb")
    cmd.load(design_pdb,"design")
    
    fixed_selection = list_to_sele(fixed[name])
    mobile_selection = list_to_sele(mobile[name])

    cmd.select('fixed_residues', f'design and resi {fixed_selection}')
    cmd.select('mobile_residues', f'design and resi {mobile_selection}')

    # Calculate distances and create a named selection 'close_pairs'
    cmd.distance('close_pairs', 'fixed_residues', 'mobile_residues', cutoff=args.INPAINT_AUTO_DISTANCE)

    # Retrieve list of distances and associated residues
    pairs = cmd.get_session()['names']['close_pairs']['measure'][0]['measurements']
    
    # Identify unique fixed residues that are within the distance threshold
    close_residues = set()
    for pair in pairs:
        index1, index2 = pair[0][0]-1, pair[1][0]-1  # PyMOL indices are 1-based, convert to 0-based
        atom1, atom2 = cmd.index('fixed_residues')[index1], cmd.index('mobile_residues')[index2]
        resi1, resi2 = atom1[0][1], atom2[0][1]  # Extract residue indices
        if resi1 in fixed:
            close_residues.add(resi1)
        if resi2 in fixed:
            close_residues.add(resi2)

    #Update occurancy
    for close_res in close_residues:
        fixed_occurancy[fixed[name].index(close_res)] += 1

    cmd.delete("design")
    cmd.delete("fixed_residues")
    cmd.delete("mobile_residues")

#Write rfd inpaint flag
inpaint_flag = set()
for i in range(len(fixed_original)):
    if fixed_occurancy[i] >= args.INPAINT_AUTO_MIN_OCCURENCY:
        inpaint_flag.append(fixed_original[i])

if args.INPAINT != "-":
    inpaint_input_chains = args.INPAINT.split('/')
    for chain in inpaint_input_chains:
        if chain[0].isnumeric(): continue
        min,max=chain[1:].split('-')
        for resi in range(int(min),int(max)+1):
            inpaint_flag.append(resi)

if args.INPAINT_AUTO_EXCLUDE != "-":
    to_exclude = sele_to_list(args.INPAINT_AUTO_EXCLUDE)
    for resi in to_exclude: 
        try:
            inpaint_flag.remove()
        except Exception:
            pass

sorted_inpaint_flag = list(sorted(inpaint_flag))
inpaint_seq = "A"+list_to_sele(sorted_inpaint_flag).replace('+','/A')
        
rfd_cmd = ""
with open(args.rfd_sh_file,'r') as rfd_sh:
    rfd_cmd = rfd_sh.readline.strip()
rfd_cmd += f" 'contigmap.inpaint_seq=[{inpaint_seq}]'"
with open(args.rfd_sh_file,'w') as rfd_sh:
    rfd_sh.write(rfd_cmd)

#Save data
#CSV file
with open(args.inpaint_csv_file,'w') as inpaint_csv:
    inpaint_csv.write("resi,occurancy")
    for i in len(fixed_original):
        inpaint_csv.write('\n')
        inpaint_csv.write(str(fixed_original[i]))
        inpaint_csv.write(',')
        inpaint_csv.write(str(fixed_occurancy[i]))

#PyMol Session
q1 = int(args.INPAINT_AUTO_NUM_DESIGNS/3.0)
q2 = int(args.INPAINT_AUTO_NUM_DESIGNS*2.0/3)
quartile1 = [fixed_original[i] for i in range(len(fixed_original)) if fixed_occurancy[i] <= q1]
quartile2 = [fixed_original[i] for i in range(len(fixed_original)) if fixed_occurancy[i] > q1 and fixed_occurancy[i] <= q2]
quartile3 = [fixed_original[i] for i in range(len(fixed_original)) if fixed_occurancy[i] > q2]
quartile1_sele = list_to_sele(quartile1)
quartile2_sele = list_to_sele(quartile2)
quartile3_sele = list_to_sele(quartile3)
name1 = f"rare (1-{q1})"
name2 = f"moderate ({q1+1}-{q2})"
name3 = f"common ({q2+1}-{args.INPAINT_AUTO_NUM_DESIGNS})"
cmd.load(args.pdb_file)

fixed_original_selection = list_to_sele(fixed_original)
cmd.select("fixed",f"resi {fixed_original_selection}")
cmd.select("non_fixed",f"not resi {fixed_original_selection}")
cmd.remove("non_fixed")
cmd.delete("non_fixed")
cmd.delete("fixed")

cmd.select(name1,f"resi {quartile1_sele}")
cmd.select(name2,f"resi {quartile2_sele}")
cmd.select(name3,f"resi {quartile3_sele}")

cmd.color("gray80")
cmd.color("yelloworange",name1)
cmd.color("orange",name2)
cmd.color("red",name3)

cmd.save(args.inpaint_pse_file)

cmd.quit() #MUST ALWAYS BE AT THE END OF THE SCRIPT