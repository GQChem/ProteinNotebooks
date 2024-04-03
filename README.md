# ProteinNotebooks
INSTALLATION OF RFDIFFUSION, PROTEINMPNN AND ALPHAFOLD2 ON SCIENCECLUSTER

1. Get access to science cluster project.
S3IT-responsible (Gianluca) or PI should write an email to help@s3it.uzh.ch ccing and providing the contact details of the person to be added, and explicitly asking to add that person to the laboratory project on Science Cluster named locbp.chem.uzh.
Open a terminal from your computer and connect to the cluster with your UZH short name. Connect even if there are fingerprint issues. You will be prompted to insert your UZH password.  
ssh -l shortname cluster.s3it.uzh.ch
This is required to create a home folder for your account. More info here: 
	https://docs.s3it.uzh.ch/cluster/overview/
2. Connect to ScienceApps:
https://apps.s3it.uzh.ch/ 

3. INSTALLATION
Notes:
	All the models will be installed in your personal data folder, which allows permanent (though not backed-up) storage of up to 200 GB. The output as well will be located there by default. Make sure you download the output of your programs to be safe. 
	I like to have my output in the scratch folder, especially for RFdiffusion (∼10 GB every 1000 seq), which allows storage of up to 20 TB and then download the most relevant files (PyMol session and ranked datasheet). However, files older than one month in scratch are automatically deleted. You can choose where to save from the notebooks.
	Lines in green are prompts to be copy-pasted and run in the shell as they are. In orange should be modified before running. In red are only for comment purposes or additional features.
	You can copy paste blocks of lines instead of running them one by one. However, if things go wrong it is better to check the output of each line.
	To open a shell in a given folder from ScienceApps, you can either use Clusters > Shell Access and then use the common shell commands (cd foldername changes directory, cd .. goes to enclosing folder, ls lists folders and files in the current directory) or Files > Home, navigate to the desired folder, then Open in terminal.
🌀 Combined installer
You can install all the models at once. Open a cluster shell and run:
module load mamba
mamba init
Close the shell and open a new one in the installation folder (data) and run:
wget https://raw.githubusercontent.com/GQChem/ProteinNotebooks/main/InstallationFiles/protein_models.sh
chmod +x protein_models.sh
bash protein_models.sh
Make sure you are in the ProteinEnv environment: when the code is terminated: you should see (ProteinEnv) in the shell input. If it is not the case, installation failed and you should follow the manual installation (Separate file). 
For later update of the notebooks and scripts only, run this from the installation folder.
bash protein_models.sh


4. USAGE
Run Jupyter (ScienceApps > Interactive Apps > Jupyter server) and set your parameters. 
A. Single job or Small batch calculations (e.g. 1-4 h): choose appropriate time and GPU. Performance is T4 < V100 < A100. Do NOT run models on CPU - they will crash!
B. Batch calculations (<24 h): No need for GPU nor extensive time. We only use the Jupyter interface to generate jobs to be submitted via Slurm.
Launch > Wait for the process to start (usually 1 min) > Connect to Jupyter.

Open the notebook you need. The usage is the same for all. There are three cells:
Cell 1. The first lines define the model input and should be modified according to your needs. Running the cell will create the appropriate bash files for execution. 
Cell 2. This cell runs inside Jupyter notebook (A) whatever was configured in before. By default, it runs MODEL_batch.sh, which contains the whole queue. Don't remove the %%bash command.
Cell 3. Generates required input for Slurm (B). Run it and follow the instructions in the output. You'll find the console output in the job.out file (see JobsComposer interface), and the script output in the usual data directory. You can examine the output using ScienceApps > Files to avoid running a new Jupyter server session.
Notes: 
	PDB files should be uploaded in the home folder when needed as input. The scripts will make a copy of them in the output folder.
	Consider building a slurm file with a few small pipelines instead of a big one. It helps if anything goes wrong on the cluster side (e.g. your resources have been reallocated). E.g. instead of running 50 designs x 3 sequences each, run 5 pipelines with 10 designs x 3 sequences each. Log files also become easier to analyze.
	There are four kind of queues on the cluster: < 24 h, < 48 h, < 7 days, < 28 days. 
For more information check:
	https://docs.s3it.uzh.ch/cluster/job_submission/
	https://slurm.schedmd.com/quickstart.html

5. Troubleshooting
If you encounter any issue while installing or using the model from this procedure (or with S3IT in general) even if you solved it, get in touch with the S3IT-responsible labmate or simply write the solution in the following section of the shared protocol.
a. USAGE
i. It happened to not being able to connect to ScienceApps after reloading the page due to "Bad request: Header size issue". Solution: Delete the cache and cookies from the browser (relative to app use).
ii. You've run a 24 h calculation with 1000 proteins in output but guess what? The job stopped when AF2 was folding protein 996 due to time limit. You can still rank these proteins! Open the last pipeline_NNN.sh file in data/bash_files, then simply copy-paste and run the rank instruction in the ProteinEnv conda environment (either ScienceClusters > Shell access or in a Jupyter server terminal). It will look like this:
	python /home/USER-NAME/data/RFDutils/rank.py ... True
Max 50-150 sequences (num designs X seq per target) are recommended for RFdiffusion in 24 h, depending on the protein length.
a. INSTALLATION
i. Got this message when activating NVIDIA SE3 Transformer: Solving environment: \ Killed. Yo can avoid it by having more RAM, but to solve it, do NOT install RFdiffusion from a Jupyter Server session terminal. Rather, use the shell as outlined in this protocol.
ii. When running RFdiffusion: AttributeError: module 'jax' has no attribute 'linear_util'. or RFdiffusion: ModuleNotFoundError: no module named 'rfdiffusion'. This is a due to package versioning of cuda (NVIDIA GPU) and pytorch:
	https://github.com/YoshitakaMo/localcolabfold/issues/212
This was encountered when attempting to adapt the colab notebook to the server. Solution: install the adapted environment specs (yml file) and use the provided notebooks.
6. Various
i. You can untar the scaffolded protein binder design (PPI) examples of rfdiffusion if you need them, open a shell in RFdiffusion folder then run:
tar -xvf examples/ppi_scaffolds_subset.tar.gz -C examples/
ii. Jupyter suggests creating a kernel for installation of packages. However, it was found that with this method the notebooks don't work either and the only way to have RFdiffusion to work is to use bash commands in the SE3nv-cuda117 environment, as implemented already. In case you're still interested, these are the instructions to add any environment to Jupyter Notebook:
source activate ENVIRONMENT_NAME 
ipython kernel install --user --name ENVIRONMENT_NAME
When starting a notebook, you can now choose ENVIRONMENT_NAME.
iii. This whole protocol was written with the idea of using the user-friendly ScienceApps website only. If you're a shell person, you can connect with these instructions instead: 
	https://docs.s3it.uzh.ch/cluster/connecting/
iv. When installing localcolabfold, you should run the following if you plan to use alphafold from console (not needed with provided notebooks). When running the installer, the output will automatically suggest what command to run for this purpose. It will look like this:
export PATH="/path/to/your/localcolabfold/colabfold-conda/bin:$PATH"
v. If you want to implement MMSEQ locally: Create a new Job with ScienceApps > Jobs > JobComposer > From DefaultTemplate. Edit the job name (Job Options) with something like MMSeqServer. Open the job.sh in editor, then copy-paste the following and substitute your shortname:
#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --time=23:59:00
cd /scratch/SHORTNAME
wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
tar xvfz mmseqs-linux-avx2.tar.gz
export PATH=$(pwd)/mmseqs/bin/:$PATH
wget https://raw.githubusercontent.com/sokrypton/ColabFold/main/setup_databases.sh
chmod +x setup_databases.sh
mkdir MMSeqsDB
MMSEQS_NO_INDEX=1 ./setup_databases.sh MMSeqsDB
rm setup_databases.sh
You will need to resubmit the same script once every month. It takes approx. 1 day.
7. Additional resources
It is recommended to familiarize with the models before starting:
AlphaFold:
	https://www.nature.com/articles/s41586-021-03819-2
	https://github.com/google-deepmind/alphafold
	https://github.com/sokrypton/ColabFold
	https://github.com/YoshitakaMo/localcolabfold
RFdiffusion:
	https://www.science.org/doi/10.1126/science.abj8754
	https://www.nature.com/articles/s41586-023-06415-8
	https://github.com/RosettaCommons/RFdiffusion
ProteinMPNN:
	https://www.science.org/doi/10.1126/science.add2187
	https://github.com/dauparas/ProteinMPNN	https://huggingface.co/spaces/simonduerr/ProteinMPNN/blob/main/ProteinMPNN/README.md
	https://meilerlab.org/wp-content/uploads/2022/12/protein_mpnn_tutorial_Nov2022.pdf
![image](https://github.com/GQChem/ProteinNotebooks/assets/161027511/9cdcbd12-1656-4d32-a99f-dff441e28319)
