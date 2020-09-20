# Oxdna with DNA/RNA - Protein Hybrid Models
An extension of the oxDNA/RNA models that also includes Anisotropic Network Model for representation of proteins. 
Included in the repository is the simulation code in the oxDNA directory as well as the scripts for generating peptides/proteins with examples.
Examples of setting up a simulation of protein (represented using ANM model, or ANMT model) with oxDNA or oxRNA, along with the necessary conversion scripts, are provided in [ANMUtils directory](https://github.com/sulcgroup/anm-oxdna/tree/master/oxDNA/ANMUtils).

An overview of the models is provided in the manuscript as well as in [documentation](https://github.com/sulcgroup/anm-oxdna/blob/master/oxDNA/ANMUtils/anm-oxdna-overview.pdf) online.

# Models 
Below are the Interaction Types for protein and DNA/RNA - protein hybrid simulations along with a brief description of their purpose. The models
are specified as interaction_type in the input file. The respective interaction types are discussed in the [documentation](https://github.com/sulcgroup/anm-oxdna/blob/master/oxDNA/ANMUtils/anm-oxdna-overview.pdf) and include:

  * AC -> For Simulation of Proteins exclusively, Treats the system as the classic ANM
  * ACT -> For Simulation of Proteins exclusively, Uses classic ANM and bending/torsional modulation
  * DNANM -> For Simulation of Proteins and DNA, Uses classic ANM for protein and oxDNA2 for DNA
  * DNACT -> For Simulation of Proteins and DNA, Uses ANMT model for protein and oxDNA2 for DNA
  * RNANM -> For Simulation of Proteins and RNA, Uses classic ANM for protein and oxRNA2 for RNA
  * RNACT -> For Simulation of Proteins and RNA, Uses ANMT model for protein and oxRNA2 for RNA

# Documentation
Basic Documentation relating to file formats (.top, .par, .dat, and input) can be found in [documentation](https://github.com/sulcgroup/anm-oxdna/blob/master/oxDNA/ANMUtils/anm-oxdna-overview.pdf) and their genration from PDB file are provided in [ANMUtils directory](https://github.com/sulcgroup/anm-oxdna/tree/master/oxDNA/ANMUtils). 

Example Usage of the model with the provided scripts is shown in a series of example [Jupyter Notebooks](https://github.com/sulcgroup/anm-oxdna/tree/master/oxDNA/ANMUtils) that are designed as walk-through of settign up simulation files and model parameters, as well as running simulations of the respective systems.
