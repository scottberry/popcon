# Running on cluster

Submit on cluster using (e.g.):
```
$ sbatch -c 16 ~/batch/run_python_default.sh popcon.py ~/NascentRNA/20190114-NascentRNA-Multiplex/FEATURES/Cells/ ~/NascentRNA/20190114-NascentRNA-Multiplex/FEATURES/POPCON -r 100 200 300 400 500 --parallel -x 2048 -y 2048 -d 4

$ sbatch -c 8 --mem 29000 /data/homes/jluethi/popcon/run_python_default.sh popcon.py /data/homes/jluethi/20190702-StainingTest1/feature-values /data/homes/jluethi/20190702-StainingTest1/popcon_features -m Nuclei -c Nuclei --calculate_density -r 300 400 500 --parallel -d 2 --find_edge -exp 250 -er 1
```
