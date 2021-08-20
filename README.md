This repository contains the data and scripts for the manuscript: "Nanoconfinement of Microvilli Alters Gene Expression and Boosts T cell Activation" by Aramesh et al. 
It contains :

# 1_ Gene Transcriptomics 
a) Raw data for gene transcriptomics (porous/non-porous +/- samples with 3 replicates)
b) R code used for the analysis of the diffrentially expressed genes

# 2_ Image analysis

Detecting and tracing actin protrusions in T-cells. 

The following Python package versions are used:<br>
<pre>
python:       3.7.5

numpy:        1.16.4 
h5py:         2.10.0 
skimage:      0.16.2 
scipy:        1.4.1  
bokeh:        2.0.0  
ipywidgets:   7.5.1  
matplotlib:   3.1.1  
tifffile:     2020.2.16  
</pre>

### Tracking protrusions in a Z-stack

<code>depth_tracking_of_protrusions.ipynb</code>

### Tracking protrusions over time

<code>temporal_tracking_of_protrusions.ipynb</code>

