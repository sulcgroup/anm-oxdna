# Oxdna with DNA/RNA - Protein Hybrid Models
An extension of the oxDNA/RNA models that also includes Anisotropic Network Model for representation of proteins. 
Included in the repository is the simulation code in the oxDNA directory as well as the scripts for generating peptides/proteins.

# Models 
Below are the Interaction Types for protein and DNA/RNA - protein hybrid simulations along with a brief description of their purpose:

  #### AC -> For Simulation of Proteins exclusively, Treats the system as the classic ANM
  #### ACT -> For Simulation of Proteins exclusively, Uses classic ANM and bending/torsional modulation
  #### DNANM -> For Simulation of Proteins and DNA, Uses classic ANM for protein and oxDNA2 for DNA
  #### DNACT -> For Simulation of Proteins and DNA, Uses ANMT model for protein and oxDNA2 for DNA
  #### RNANM -> For Simulation of Proteins and RNA, Uses classic ANM for protein and oxRNA2 for RNA
  #### RNACT -> For Simulation of Proteins and RNA, Uses ANMT model for protein and oxRNA2 for RNA
  #### CUDADNANM -> CUDA version of DNANM
  #### CUDADNACT -> CUDA version of DNACT
  #### CUDARNANM -> CUDA version of RNANM
  #### CUDARNACT -> CUDA version of RNACT

# Documentation
Basic Documentation relating to file formats (.top, .par, .dat, and input) can be found in BLANK directory

Example Usage of the model is shown in a series of Jupyter Notebooks availabe in the ANM_EXAMPLES directory
