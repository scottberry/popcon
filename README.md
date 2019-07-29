# Installation
Works on the cluster, can have issues on Macs depending on the installed C++ libraries. Deleting and installing a fresh version of Xcode did the trick for me (Joel).
Create & activate your virtualenv, then:

```
pip install pybind11
cd popcon
pip install .
```

# Quick start
The popcon module runs based on feature-data & metadata downloaded from TissueMaps. Ensure that you run a current version of measure morphology on these objects, such that they have Morphology_Local_Centroid_y & Morphology_Local_Centroid_x measurements. Depending on the size of your wells (number of sites per well), the density measurements can get quite memory heavy. Also, when running in parallel mode, ensure enough memory is given for the edge detection (scale memory with # of cores).
Then call it the following way:
```
popcon.py PATH/TO/TISSUEMAPS/CSVs /OUTPUT/PATH -m OBJECT_TO_BE_MEASURED -c OBJECT_THAT_FEATURES_WERE_DOWNLOADED_FOR --calculate_density -r RADI_FOR_DENSITY MEASUREMENT --parallel -x X_DIM_OF_IMAGE -y Y_DIM_OF_IMAGE -d DOWNSAMPLING_FOR_DENSITY --find_edge -exp AMOUNT_OF_EXPANSION_FOR_EDGE_DETECTION -er RADIUS_OF_SITES_FOR_EDGE_DETECTION
```
Leave away the `--calculate_density` or `--find_edge` if you don't need the density measurement or the edge detection.


# Examples for running on cluster

Submit on cluster using the following examples. The run_python_default.sh script activates the virtual environment and calls the popcon python script. Go to the popcon folder to call this:

```
$ sbatch -c 16 ~/batch/run_python_default.sh popcon.py ~/NascentRNA/20190114-NascentRNA-Multiplex/FEATURES/Cells/ ~/NascentRNA/20190114-NascentRNA-Multiplex/FEATURES/POPCON -r 100 200 300 400 500 --parallel -x 2048 -y 2048 -d 4

$ sbatch -c 8 --mem 29000 /data/homes/jluethi/popcon/run_python_default.sh popcon.py /data/homes/jluethi/20190702-StainingTest1/feature-values /data/homes/jluethi/20190702-StainingTest1/popcon_features -m Nuclei -c Nuclei --calculate_density -r 300 400 500 --parallel -d 2 --find_edge -exp 250 -er 1
```
