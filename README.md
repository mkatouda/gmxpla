# gmxpla

Gromax protein-ligand MD trajectory analysis tools to obtain csv and png file summary of 
protein-ligand interaction energy, ligand RMSD, and ligand RMSF.

## Licence

This package is distributed under the MIT License.

## Tutorial

You can try tutorial of meekovina to know how to install and run:  
https://colab.research.google.com/github/mkatouda/gmxpla/blob/main/gmxpla_tutorial_jp.ipynb

## Required softwares

1. python: 3.7 or later
2. pyyaml (https://pyyaml.org/)
3. matplotlib (https://matplotlib.org/)
4. gromacs (https://www.gromacs.org/)

## Installation

### Install gromacs

```
See install guide: https://manual.gromacs.org/current/install-guide/index.html
```

### Install from github

```
pip install git+https://github.com/mkatouda/gmxpla.git
```

### Local install

```
git clone https://github.com/mkatouda/gmxpla.git
cd gmxpla
pip install .
```

## Command usage

```
usage: gmxpla [-h] [-i INP] [-e EDR] [-t TPR] [-x XTC] [-n NDX] [-oc OUTCSV] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i INP, --inp INP     yaml style input file, overwriting argument values (default: None)
  -e EDR, --edr EDR     Gromacs energy file (edr file) (default: None)
  -t TPR, --tpr TPR     Gromacs topology file (tpr or gro file) (default: None)
  -x XTC, --xtc XTC     Gromacs trajectory file (xtc file) (default: None)
  -n NDX, --ndx NDX     Gromacs index file (ndx file) (default: None)
  -oc OUTCSV, --outcsv OUTCSV
                        docking score output (csv file) (default: docking_score.csv)
  -v, --verbose         Verbose output. (default: False)
```

## Exmaples of command line usage

```
gmxpla -e prod.edr -t prod.tpr -x prod.xtc -n index.ndx -oc prod_docking_score.csv
```

## Exmaples of yaml input usage

Prepare input yaml file input.yml:

```
edr: './prod.edr'
tpr: './prod.tpr'
xtc: './prod.xtc'
ndx: './index.ndx'
outcsv: './prod_docking_score.csv'
verbose: False
```

Then, run gmxpla in command line:

```
gmxpla -i input.yml
```

Keywards of yaml file are the same in the name of command line options.  
See above explation of command line options.  

## Author

Michio Katouda (katouda@rist.or.jp)  
