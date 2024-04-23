#Copyright © 2024 LOCBP @ University of Zürich
#Distributed under MIT license
"""
Script used by all models to reverse translate to DNA and optimize the sequence
Usage: python dna_encoder.py path/to/input.fasta path/to/output/prefix
Requires: biopython
"""
print("DNA Encoder is based on CoCoPUTs (HIVE, tables updated April 2024). Please cite accordingly.")

import argparse

parser = argparse.ArgumentParser(description='Slices the queries for alphafold refolding of the best ones')
parser.add_argument('proteins_fasta_file', type=str)
parser.add_argument('output_prefix', type=str)
args = parser.parse_args()
proteins_fasta_file=args.proteins_fasta_file
output_prefix=args.output_prefix

from Bio import SeqIO

protein_sequences = dict()
for seq_record in SeqIO.parse(proteins_fasta_file, "fasta"):
    protein_sequences[seq_record.id] = str(seq_record.seq)

#Organism-specific codon usage tables (Source: https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=cocoputs)
CoCoPUTs_tables = {
"Homo sapiens":"""TTT	17.14	(1385301)	TCT	16.93	(1368632)	TAT	12.11	( 978774)	TGT	10.40	( 841042)
TTC	17.48	(1413268)	TCC	17.32	(1399962)	TAC	13.49	(1090514)	TGC	10.81	( 873765)
TTA	 8.71	( 703680)	TCA	14.14	(1142684)	TAA	 0.44	(  35218)	TGA	 0.79	(  63801)
TTG	13.44	(1086777)	TCG	 4.03	( 325925)	TAG	 0.35	(  28499)	TGG	11.60	( 937286)
 	 	 	 	 	 	 	 	 	 	 	 
CTT	14.08	(1138433)	CCT	19.31	(1560898)	CAT	11.83	( 956479)	CGT	 4.55	( 367659)
CTC	17.81	(1439345)	CCC	19.11	(1544626)	CAC	14.65	(1184041)	CGC	 8.71	( 704401)
CTA	 7.44	( 601662)	CCA	18.92	(1529004)	CAA	14.06	(1136523)	CGA	 6.42	( 518818)
CTG	36.10	(2918400)	CCG	 6.22	( 503096)	CAG	35.53	(2872161)	CGG	10.79	( 871786)
 	 	 	 	 	 	 	 	 	 	 	 
ATT	16.48	(1331901)	ACT	14.26	(1152700)	AAT	18.43	(1489775)	AGT	14.05	(1135376)
ATC	18.67	(1508988)	ACC	17.85	(1442511)	AAC	18.30	(1478832)	AGC	19.69	(1591829)
ATA	 8.08	( 652939)	ACA	16.52	(1335468)	AAA	27.48	(2221062)	AGA	13.28	(1073213)
ATG	21.53	(1739992)	ACG	 5.59	( 452037)	AAG	31.77	(2567940)	AGG	12.13	( 980476)
 	 	 	 	 	 	 	 	 	 	 	 
GTT	11.74	( 949137)	GCT	18.99	(1534685)	GAT	24.03	(1942185)	GGT	10.83	( 875715)
GTC	13.44	(1086717)	GCC	25.84	(2088762)	GAC	24.27	(1961667)	GGC	19.79	(1599325)
GTA	 7.66	( 618960)	GCA	17.04	(1377145)	GAA	33.65	(2719693)	GGA	17.12	(1384137)
GTG	25.87	(2090923)	GCG	 5.91	( 477758)	GAG	39.67	(3206546)	GGG	15.35	(1240793)""",
"Escherichia coli":"""TTT	22.28	( 884859304)	TCT	 8.59	( 341230176)	TAT	16.28	( 646775431)	TGT	 5.15	( 204369628)
TTC	16.26	( 645822743)	TCC	 8.82	( 350462264)	TAC	12.18	( 483602182)	TGC	 6.31	( 250671692)
TTA	13.72	( 544875130)	TCA	 7.50	( 297755617)	TAA	 2.01	(  79698775)	TGA	 1.00	(  39709081)
TTG	13.33	( 529575852)	TCG	 8.79	( 349261974)	TAG	 0.23	(   9277733)	TGG	15.22	( 604449078)
 	 	 	 	 	 	 	 	 	 	 	 
CTT	11.39	( 452555895)	CCT	 7.19	( 285670480)	CAT	12.77	( 507363784)	CGT	20.72	( 822947687)
CTC	10.95	( 434982409)	CCC	 5.54	( 219881409)	CAC	 9.44	( 374998113)	CGC	21.55	( 855993095)
CTA	 3.88	( 154295599)	CCA	 8.42	( 334296050)	CAA	15.03	( 596935496)	CGA	 3.64	( 144705009)
CTG	52.32	(2078211248)	CCG	22.80	( 905697146)	CAG	29.26	(1162220193)	CGG	 5.68	( 225733345)
 	 	 	 	 	 	 	 	 	 	 	 
ATT	30.19	(1199278030)	ACT	 8.98	( 356533130)	AAT	18.08	( 717967207)	AGT	 9.02	( 358422205)
ATC	24.68	( 980290510)	ACC	23.00	( 913523340)	AAC	21.50	( 853806657)	AGC	15.92	( 632395429)
ATA	 4.72	( 187294284)	ACA	 7.53	( 299007362)	AAA	33.81	(1342811432)	AGA	 2.35	(  93340146)
ATG	27.65	(1098018887)	ACG	14.52	( 576519827)	AAG	10.67	( 423892201)	AGG	 1.40	(  55660373)
 	 	 	 	 	 	 	 	 	 	 	 
GTT	18.38	( 729993910)	GCT	15.54	( 617228810)	GAT	32.42	(1287642772)	GGT	24.54	( 974821805)
GTC	15.12	( 600554770)	GCC	25.62	(1017548717)	GAC	19.23	( 763756668)	GGC	28.83	(1145183663)
GTA	10.96	( 435276257)	GCA	20.59	( 817719553)	GAA	39.62	(1573751166)	GGA	 8.37	( 332414382)
GTG	25.97	(1031522759)	GCG	32.96	(1309274393)	GAG	18.29	( 726478043)	GGG	11.27	( 447491989)"""
}

frequency_tables = dict() #Contains CoCoPUTs tables in format "Codon": Frequency
ATGME_tables = dict() #Contains CoCoPUTs tables converted in a format compatible with atgme.org (also T->U)
for organism in CoCoPUTs_tables.keys():
    table = CoCoPUTs_tables[organism]
    frequency_table = dict()
    atgme_table = ""
    lines = table.split('\n')
    for line in lines:
        line = line.strip()
        if line == "": continue
        tab_split = line.split('\t')
        codon_frequency = [(tab_split[i],tab_split[i+1]) for i in [0,3,6,9]]
        for i, x in enumerate(codon_frequency):
            codon,frequency_s=x
            #frequency table
            frequency=float(frequency_s.replace('\xa0',' '))
            frequency_table[codon]=frequency
            #atmge table
            if i > 0: atgme_table += "  "
            codon_u = codon.replace("T","U") if "T" in codon else codon
            atgme_table += codon_u
            frequency_s_f="{:.1f}".format(frequency) #Formatted with one decimal for atmge.org
            for i in range(5-len(frequency_s_f)):
                atgme_table += " "
            atgme_table += frequency_s_f
            atgme_table += "(      )" #the website doesn't care about the number within the brackets, onnly the number of spaces
        atgme_table += "\n"
    frequency_tables[organism]=frequency_table
    ATGME_tables[organism]=atgme_table

#Now define a dictionary containing, for each amino acid, a vector called "codons" and another called "frequencies"
from Bio.Data import CodonTable
std_table=CodonTable.standard_dna_table.forward_table #dictionary of the form "codon":"1-letter AA"

aa_codon_frequency_tables = dict() #Contains for each amino acid, a vector with codons and another with frequencies
for organism in CoCoPUTs_tables.keys():
    aa_codon_frequency_table=dict()
    for codon in std_table.keys():
        aa = std_table[codon]
        if not aa in aa_codon_frequency_table.keys():
            aa_codon_frequency_table[aa]=[]
            aa_codon_frequency_table[aa].append([])
            aa_codon_frequency_table[aa].append([])
        aa_codon_frequency_table[aa][0].append(codon)
        aa_codon_frequency_table[aa][1].append(frequency_tables[organism][codon])
    aa_codon_frequency_tables[organism]=aa_codon_frequency_table
#print(aa_codon_frequency_tables)
#Make DNA sequences based on relative probabilities of amino acids
import random
weighted_DNAs = dict()
most_frequent_codons_DNAs = dict()
half_max_DNAs = dict()
threshold_weighted_DNAs = dict()
for protein in protein_sequences.keys():
    weighted_DNAs[protein]=dict()
    most_frequent_codons_DNAs[protein]=dict()
    half_max_DNAs[protein]=dict()
    threshold_weighted_DNAs[protein]=dict()

    sequence = protein_sequences[protein]
    for organism in CoCoPUTs_tables.keys():
        weighted_DNAseq=""
        most_frequent_codons_DNAseq=""
        half_max_DNAseq=""
        threshold_weighted_DNAseq=""
        aa_count=0
        for aa in sequence:
            if aa_count>0: 
                weighted_DNAseq+='\t'
                most_frequent_codons_DNAseq+='\t'
                half_max_DNAseq+='\t'
                threshold_weighted_DNAseq+='\t'
            codons = aa_codon_frequency_tables[organism][aa][0]
            frequencies = aa_codon_frequency_tables[organism][aa][1]

            weighted_DNAseq += random.choices(codons, weights=frequencies, k=1) [0]

            max_frequency = max(frequencies)
            most_frequent_codons_DNAseq += codons[frequencies.index(max_frequency)]

            half_max_indeces = [hmi for hmi in range(len(frequencies)) if frequencies[hmi] >= max_frequency/2.0]
            codons_half_max = [codons[i] for i in half_max_indeces]
            frequencies_half_max = [frequencies[i] for i in half_max_indeces]
            half_max_DNAseq += random.choices(codons_half_max,weights=frequencies_half_max,k=1)[0]

            threshold_indeces = [hmi for hmi in range(len(frequencies)) if frequencies[hmi] >= 10.0]
            if len(threshold_indeces)>0:
                codons_threshold = [codons[i] for i in threshold_indeces]
                frequencies_threshold = [frequencies[i] for i in threshold_indeces]
                threshold_weighted_DNAseq += random.choices(codons_threshold,weights=frequencies_threshold,k=1)[0]
            else:
                threshold_weighted_DNAseq += codons[frequencies.index(max_frequency)]

            aa_count+=1

        weighted_DNAs[protein][organism]=weighted_DNAseq
        most_frequent_codons_DNAs[protein][organism]=most_frequent_codons_DNAseq
        half_max_DNAs[protein][organism]=half_max_DNAseq
        threshold_weighted_DNAs[protein][organism]=threshold_weighted_DNAseq
#Save
for organism in CoCoPUTs_tables.keys():
    organism_u = organism.strip().replace(' ','_')
    with open(output_prefix+"_DNA_"+organism_u+"_weighted.fasta",'w') as file_fasta:
        for i, protein in enumerate(list(protein_sequences.keys())):
            if i > 0: file_fasta.write('\n')
            file_fasta.write(f">{protein}_{organism_u}_w")
            file_fasta.write('\n')
            file_fasta.write(weighted_DNAs[protein][organism])
    with open(output_prefix+"_DNA_"+organism_u+"_most_frequent_codons.fasta",'w') as file_fasta:
        for i, protein in enumerate(list(protein_sequences.keys())):
            if i > 0: file_fasta.write('\n')
            file_fasta.write(f">{protein}_{organism_u}_mfc")
            file_fasta.write('\n')
            file_fasta.write(most_frequent_codons_DNAs[protein][organism])
    with open(output_prefix+"_DNA_"+organism_u+"_more_than_half_maximum_frequency.fasta",'w') as file_fasta:
        for i, protein in enumerate(list(protein_sequences.keys())):
            if i > 0: file_fasta.write('\n')
            file_fasta.write(f">{protein}_{organism_u}_mthmf")
            file_fasta.write('\n')
            file_fasta.write(half_max_DNAs[protein][organism])
    with open(output_prefix+"_DNA_"+organism_u+"_threshold_10.fasta",'w') as file_fasta:
        for i, protein in enumerate(list(protein_sequences.keys())):
            if i > 0: file_fasta.write('\n')
            file_fasta.write(f">{protein}_{organism_u}_t10")
            file_fasta.write('\n')
            file_fasta.write(threshold_weighted_DNAs[protein][organism])