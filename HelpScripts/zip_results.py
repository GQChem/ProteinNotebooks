#It is assumed that the most important files don't start with an underscore and are not in subfolders

import argparse

parser = argparse.ArgumentParser(description='Creates a zip containing the most important files')
parser.add_argument('JOB_FOLDER', type=str)

# Parse the arguments
args = parser.parse_args()

import os
import zipfile

zip_name = os.path.basename(args.JOB_FOLDER)
result_files = [file for file in os.listdir(args.JOB_FOLDER) if not file.startswith('_')]
zip_file = os.path.join(args.JOB_FOLDER, f"{zip_name}.zip")
with zipfile.ZipFile(zip_file, 'w') as zip:
    for file in result_files:
        zip.write(os.path.join(args.JOB_FOLDER,file),arcname=file)


