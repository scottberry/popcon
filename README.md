# Running on cluster

Submit on cluster using (e.g.):
```
$ sbatch -c 16 ~/batch/run_python_default.sh popcon.py ~/NascentRNA/20190114-NascentRNA-Multiplex/FEATURES/Cells/ ~/NascentRNA/20190114-NascentRNA-Multiplex/FEATURES/POPCON -r 100 200 300 400 500 --parallel -x 2048 -y 2048 -d 4
```
