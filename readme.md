[PROJECT TITLE]
Author: Dionysus Morris-Evans


This script is designed to ...


USAGE: bash example_script.sh [ARG1] [ARG2]


######################


INPUTS:

REQUIRED FLAG INPUTS:
File1 - PRJNA Accession ID number 
File2 - ...


OPTIONAL FLAG INPUTS:
Arg1 - Input File ...
Arg2 - Output File ... 


######################


OUTPUTS:

Directories {6}:

Directory 1: raw_data - contains raw data
Directory 2: meta - contains...
Directory 3: logs - contains...
Directory 4: results - contains...
Directory 5: scripts - contains...
Directory 6: reference_data - contains...

Files {xxx00}:
File1: Returns a text file output containing listed SRA IDs and. Outputs to SRA_Run
File2:
File3:
File(s)4: aligned_bam files - Returns n number of files relating to the index of provided motifs in FASTA format within Directory1. Outputs to Directory1 of relative path (./motifs/)
File5: Output Summary File - Returns Summary text file report.

######################

FILE NAMING CONVENTION:

[chronology#]_[YYYY-MM-DD]_[taxonomy]_[prjctID]_[expname-smplname-repl#]_[###]_[researcherinitials]_[v##].xxx

example:
01_2025-03-24_spombe_prjctA_expA-sample1-replicate1_DME_v3.xlsx


BUG REPORTS:
- File(n) outputs occassionally contain -- output in FASTA files
- Needs rechecked to ensure relative paths function to always execute smoothly
- User inputs should provide relative paths as arguments.
- Inputs must exist otherwise the script will not execute successfully

V0.0.0-beta
YYYY-MM-DD
