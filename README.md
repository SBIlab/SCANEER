# SCANEER
This package provides an implementation of a pipeline of SCANEER. SCANEER scans evolution of protein sequences and direct mutation strategy to improve enzyme activity.

## Requirements
+ Python (v 2.7.13)
+ Biopython (v 1.72)
+ NumPy (v 1.15.4)
+ Java (v 1.7.0)

## Installation instruction
+ All python packages can be installed via pip (https://pypi.org/project/pip/).
+ Installations would take few minutes for each python package.

## Running SCANEER
1. Clone this repository and ```cd``` into it.
2. Check input files in a directory ```input/```.
    + The input file is a multiple sequence alignment of desired enzyme to be engineered.
    + Input files should be CLUSTAL formatted (*.aln).
    + You can run SCANEER using the provided multiple sequence alignments in ```input/```.
    + You can change the path containing input files by modifying ```input_path``` in a script ```run_SCANEER.py```.
3. Run ```run_SCANEER.py```.
    ```bash
    python run_SCANEER.py
    ```
4. The outputs will be in a directory ```output/```.
    + The directory containing output files will be created automatically.
    + You can change the path containing output files by modifying ```output_path``` in a script ```run_SCANEER.py```.
    + The contents of each output file are as follows:
        + ```*.aln_cn``` - A CLUSTAL format file containing a multiple sequence alignment to calculate covarying strength.
        + ```*.coe_out_mcbasc``` - A text file containing the calculated covarying strengths of all combination of residue pairs.
        + ```*.cn``` - A text file containing the number of co-evolutionary relationships of residues.
        + ```*.coenet``` - A text file containing a residue-residue co-evolutionary network of the enzyme
        + ```*.txt``` - A text file containing final SCI scores of mutations.
