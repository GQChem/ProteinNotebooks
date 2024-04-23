#Copyright © 2024 LOCBP @ University of Zürich
#Distributed under MIT license

#It is assumed that the most important files don't start with an underscore and are not in subfolders
import argparse
import os
import zipfile

# Set up argument parser
parser = argparse.ArgumentParser(description='Creates a zip containing the most important files')
parser.add_argument('JOB_FOLDER', type=str)

# Parse the arguments
args = parser.parse_args()

# Define the zip file name based on the job folder
zip_name = os.path.basename(args.JOB_FOLDER)
zip_file = os.path.join(args.JOB_FOLDER, f"{zip_name}.zip")

# Create a list to hold all files and directories to be zipped
result_files = []

# Walk through the directory, adding files and the DNA folder if present
for root, dirs, files in os.walk(args.JOB_FOLDER):
    # Filter out hidden files and directories
    files = [f for f in files if not f.startswith('_') and '.' in f]
    dirs[:] = [d for d in dirs if d != 'DNA']  # Exclude all but DNA folder from further walking

    for file in files:
        result_files.append(os.path.join(root, file))

    # If the 'DNA' folder is found, include it and its contents
    if 'DNA' in dirs:
        dna_path = os.path.join(root, 'DNA')
        for dna_root, dna_dirs, dna_files in os.walk(dna_path):
            for file in dna_files:
                result_files.append(os.path.join(dna_root, file))

# Write collected files and folders to the zip file
with zipfile.ZipFile(zip_file, 'w') as zip:
    for file in result_files:
        arcname = os.path.relpath(file, args.JOB_FOLDER)  # Set the archive name to relative path for organization
        zip.write(file, arcname=arcname)



