# Running on cluster

Submit on cluster using (e.g.):
```
$ sbatch -c 16 ~/batch/run_python_default.sh popcon.py ~/NascentRNA/20190114-NascentRNA-Multiplex/FEATURES/Cells/ ~/NascentRNA/20190114-NascentRNA-Multiplex/FEATURES/POPCON -r 100 200 300 400 500 --parallel -x 2048 -y 2048 -d 4
```


Joel's example:
sbatch -c 4 --mem 15000 /data/homes/jluethi/popcon/run_python_default.sh popcon.py ~/20190505-WTC-PermeabilizationTest5 /data/homes/jluethi/20190505-WTC-PermeabilizationTest5/popcon_features -m Nuclei -c Nuclei -r 100 200 300 400 500 --parallel -d 2
