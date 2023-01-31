# LeaderTrailerAntiterminationComplexIn30SBiogenesis
This directory contains the custom scripts used for the data analysis in the manuscript "Roles of the leader-trailer helix and antitermination complex in biogenesis of the 30S ribosomal subunit" by Benjamin R. Warner, Ralf Bundschuh, and Kurt Fredrick.

# Requirements
The following packages are required to run this pipeline:
  - Python3
  - Biopython
  - Pandas
  - locaRNA (including RNAclust)
  - Newick utils

# Quick user guide
The first script to run is get_16S.py. It collects all the 16S leader-trailer sequences for the organisms defined in the file AllGTDBEnterobacteriaceaeWithNCBIType.csv.
Then, main.sh clusters the RNAs by structure and predicts consensus structures for each of the identified subgroups. Finally, GetSupplementalTableInfo.py generates
Supplemental Table S1 of the manuscript.

