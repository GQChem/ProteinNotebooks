import argparse

parser = argparse.ArgumentParser(description='Slices the queries for alphafold refolding of the best ones')
parser.add_argument('queries_csv_file', type=str)
parser.add_argument('rank_output_csv_file', type=str)
parser.add_argument('FOLD_BEST_WITH_ALPHAFOLD', type=int)
parser.add_argument('best_queries_csv_file', type=str)

# Parse the arguments
args = parser.parse_args()

import os

ranked_keys = []
best_ids = []
with open(args.rank_output_csv_file,"r") as ranked_csv:
    ranked_keys = ranked_csv.readline().strip().split(',')
    name_i = ranked_keys.index('name')
    best_i = 0
    for line in ranked_csv:
        if best_i >= args.FOLD_BEST_WITH_ALPHAFOLD:
            break
        best_ids.append(line.strip().split(',')[name_i])
        best_i += 1

queries_keys = []
best_queries_rows = []
with open(args.queries_csv_file,"r") as queries_csv:
    queries_keys = queries_csv.readline().strip().split(',')
    id_i = queries_keys.index("id")
    for line in queries_csv:
        line = line.strip()
        queries_rows = line.split(',')
        if queries_rows[id_i] in best_ids:
            best_queries_rows.append(line)

with open(args.best_queries_csv_file,"w") as best_queries_csv:
    for i,k in enumerate(queries_keys):
        if i > 0: best_queries_csv.write(',')
        best_queries_csv.write(k)
    for best_query in best_queries_rows:
        best_queries_csv.write('\n')
        best_queries_csv.write(best_query)

