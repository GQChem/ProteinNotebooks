import argparse

parser = argparse.ArgumentParser(description='Makes csv files from ProteinMPNN.fa output ')
parser.add_argument('PMPNN_FA_FOLDER', type=str)
parser.add_argument('queries_csv_file', type=str)
parser.add_argument('all_queries_csv_file', type=str)
parser.add_argument('cycle', type=str)

# Parse the arguments
args = parser.parse_args()

import os

all_keys_str = ""
all_keys = []
all_data = []
all_sequences = []
if os.path.exists(args.all_queries_csv_file):
    with open(args.all_queries_csv_file,"r") as all_queries:
        all_lines = all_queries.readlines()
        all_keys_str = all_lines[0].strip()
        all_keys = list(all_keys_str.split(','))
        all_data = all_lines[1:]
        seq_index = all_keys.index("sequence")
        for data in all_data:
            all_sequences.append(data[seq_index])

fa_files = os.listdir(args.PMPNN_FA_FOLDER)
for fa in fa_files:
    if fa.endswith(".fa"):
        data = []
        with open(os.path.join(args.PMPNN_FA_FOLDER,fa),"r") as file:
            lines = [line.strip() for line in file.readlines()]
            #Starts from 2 to skip the original protein
            for i in range(2,len(lines),2):
                seq_data = dict()
                seq_data["id"] = "" #put here just to preserve the order
                seq_data["sequence"] = lines[i+1]
                params = lines[i][1:].split(", ")
                for p in params:
                    key,value = p.split("=")
                    seq_data[key] = value
                seq_data["id"] = fa[:-3] + "_" + seq_data["sample"]
                if seq_data["sequence"] in all_sequences:
                    id = seq_data["id"]
                    print(f"Skipped duplicate: {id}")
                else:
                    data.append(seq_data)
                    all_sequences.append(seq_data["sequence"])
        if len(data) == 0: continue
        old_data = []
        if os.path.exists(args.queries_csv_file):
            with open(args.queries_csv_file,"r") as old_queries:
                old_data = old_queries.readlines()
        with open(args.queries_csv_file,"w") as queries:
            keys = list(data[0].keys())
            if len(old_data) != 0:
                queries.writelines(old_data)
            #Header is already part of old_data
            else:
                for i,k in enumerate(keys):
                    if i>0: queries.write(',')
                    queries.write(k)                    
            for seq_data in data:
                queries.write("\n")
                for i,k in enumerate(keys):
                    if i>0: queries.write(',')
                    queries.write(seq_data[k])

if os.path.exists(args.queries_csv_file): #there were some queries with different sequence
    with open(args.queries_csv_file,"r") as new_queries:
        new_data = [l.strip() for l in new_queries.readlines()] #contains header
        if os.path.exists(args.all_queries_csv_file):
            with open(args.all_queries_csv_file,"w") as all_queries: #there were some data already before
                all_queries.write(all_keys_str)
                all_queries.write("\n")
                all_queries.writelines(all_data)
                for nd in new_data[1:]:         
                    all_queries.write("\n")           
                    all_queries.write(args.cycle)  
                    all_queries.write(',')
                    all_queries.write(nd)
        else:
            with open(args.all_queries_csv_file,"w") as all_queries: #there were some data already before
                all_queries.write("cycle,")
                all_queries.write(new_data[0])
                for nd in new_data[1:]:         
                    all_queries.write("\n")           
                    all_queries.write(args.cycle)  
                    all_queries.write(',')
                    all_queries.write(nd)
else:
    if os.path.exists(args.all_queries_csv_file):
        with open(args.all_queries_csv_file,"w") as all_queries: #there were some data already before
            all_queries.writelines(all_keys_str)
            all_queries.write("\n")
            all_queries.writelines(all_data)