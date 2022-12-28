#!/usr/bin/env python

"""
# gmxpla

Gromax protein-ligand MD trajectory analysis tools

## Licence

This package is distributed under the MIT License.

## Required softwares

1. python: 3.7 or later
2. pyyaml (https://pyyaml.org/)
3. matplotlib (https://matplotlib.org/)
4. gromacs (https://www.gromacs.org/)

## Optional softwares
5. oddt (https://github.com/oddt/oddt)

## Installation

- Install gromacs
See install guide: https://manual.gromacs.org/current/install-guide/index.html

- Install from github
pip install git+https://github.com/mkatouda/gmxpla.git

- Local install
git clone https://github.com/mkatouda/gmxpla.git
cd gmxpla
pip install .

"""

import os
import re
import argparse
import subprocess

import yaml
import matplotlib.pyplot as plt

def is_exe(fpath):
    if not fpath: return False
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

def isInPath(program):
    fpath = os.path.split(program)[0]
    if fpath:
        if is_exe(program): return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file): return True

    return False

if not isInPath('gmx'):
    if isInPath('gmx_seq'):
        gmxSuffix = '_seq'
    elif isInPath('gmx_mpi'):
        gmxSuffix = '_mpi'
else:
    gmxSuffix = ''
gmxbin = 'gmx' + gmxSuffix

def get_parser():
    class customHelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                              argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        formatter_class=customHelpFormatter,
        description="gromax protein-ligand MD trajectory analysis tools"
    )
    parser.add_argument(
        '-i', '--inp', type=str,
        help = "yaml style input file, overwriting argument values"
    )
    parser.add_argument(
        '-e', '--edr', type=str,
        help = "Gromacs energy file (edr file)"
    )
    parser.add_argument(
        '-t', '--tpr', type=str,
        help = "Gromacs topology file (tpr or gro file)"
    )
    parser.add_argument(
        '-x', '--xtc', type=str,
        help = "Gromacs trajectory file (xtc file)"
    )
    parser.add_argument(
        '-n', '--ndx', type=str,
        help = "Gromacs index file (ndx file)"
    )
    parser.add_argument(
        '-oc', '--outcsv', type=str, default='docking_score.csv',
        help = "docking score output (csv file)"
    )
    parser.add_argument(
        '--ifpcos', action='store_true',
        help = "Whether to calc IFPcos score"
    )
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help = 'Verbose output.'
    )
    args = parser.parse_args()

    print(args)

    return args 

def set_config(args):
    # Read config yaml file
    if args.inp is not None and os.path.isfile(args.inp):
        with open(args.inp, 'r') as f:
            conf = yaml.safe_load(f)
    else:
        conf = {}

    # Set up default config values from program arguments
    conf_def = vars(args).copy()
    [conf.setdefault(k, v) for k, v in conf_def.items()]

    return conf

def run_cmd(cmd, debug=False):
    print(' '.join(cmd).replace('\n'))
    print(cmd)
    results = subprocess.run(cmd, capture_output=True, check=True, text=True)
    print(results.stdout, results.stderr)

def run_twocmd_pype(cmd1, cmd2, debug=False):
    print(' '.join(cmd1) + ' | ' + ' '.join(cmd2))
    p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    output = p2.communicate()[0] 
    p2.returncode

def plot2d(x, y, xlabel, ylabel, label, pngfile):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111)
    ax.plot(x, y, label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig(pngfile)
    plt.clf()
    plt.close()

def xvg2csvpng(xvgfile):
    with open(xvgfile) as fin:
        lines = fin.readlines()

        x = []; y = []
        p = '"([^"]*)"'
        for line in lines:
            if line[0] == '@':
                if 'xaxis' in line:
                    xlabel = re.search(p, line).group().strip('"')
                elif 'yaxis' in line:
                    ylabel = re.search(p, line).group().strip('"')
                elif 'title' in line:
                    title = re.search(p, line).group().strip('"')
            
            elif line[0] != '#' and line[0] != '@':
                list = line.split()
                x.append(float(list[0]))
                y.append(float(list[1]))

    #print(xlabel, ylabel)
    #for i in range(len(x)):
    #    print(x[i], y[i])

    s = '{},{}\n'.format(xlabel, ylabel)
    for i in range(len(x)):
        s += '{},{}\n'.format(x[i], y[i])

    basename = os.path.splitext(os.path.basename(xvgfile))[0]
    csvfile = basename + '.csv'
    with open(csvfile, 'w') as f:
        f.write(s)

    pngfile = basename + '.png'
    plot2d(x, y, xlabel, ylabel, title, pngfile)

    yavg = sum(y) / float(len(y))
    print('yavg: ', yavg)

    return yavg

def gmx_rms(tpr_path, xtc_path, ndx_path, xvg_path, selection, tu='ps', debug=False):
    cmd1 = ['echo', selection]
    cmd2 = [gmxbin, 'rms', '-s', tpr_path, '-f', xtc_path, '-n', ndx_path, '-o', xvg_path, '-tu', tu]

    run_twocmd_pype(cmd1, cmd2, debug=debug)
    return xvg2csvpng(xvg_path)

def gmx_rmsf(tpr_path, xtc_path, ndx_path, xvg_path, selection, debug=False):
    cmd1 = ['echo', selection]
    cmd2 = [gmxbin, 'rmsf', '-s', tpr_path, '-f', xtc_path, '-n', ndx_path, '-o', xvg_path]

    run_twocmd_pype(cmd1, cmd2, debug=debug)
    return xvg2csvpng(xvg_path)

def gmx_trjconv_prolig(tpr_path, xtc_path, ndx_path, debug=False):

    xtc_in_path = xtc_path
    xtc_out_path = os.path.splitext(os.path.basename(xtc_path))[0] + '_nowat.xtc'
    cmd1 = ['echo', 'TempCouplingA\n']
    cmd2 = [gmxbin, 'trjconv', '-s', tpr_path, '-f', xtc_in_path, '-n', ndx_path, '-o', xtc_out_path, '-pbc', 'nojump']
    run_twocmd_pype(cmd1, cmd2, debug=debug)

    xtc_in_path = xtc_out_path 
    xtc_out_path = os.path.splitext(os.path.basename(xtc_path))[0] + '_nowat_center.xtc'
    cmd1 = ['echo', 'LIG\n TempCouplingA\n']
    cmd2 = [gmxbin, 'trjconv', '-s', tpr_path, '-f', xtc_in_path, '-n', ndx_path, '-o', xtc_out_path, '-center', '-pbc', 'mol', '-ur', 'compact']
    run_twocmd_pype(cmd1, cmd2, debug=debug)

    xtc_in_path = xtc_out_path 
    xtc_out_path = os.path.splitext(os.path.basename(xtc_path))[0] + '_nowat_fit.xtc'
    cmd1 = ['echo', 'Backbone\n TempCouplingA\n']
    cmd2 = [gmxbin, 'trjconv', '-s', tpr_path, '-f', xtc_in_path, '-n', ndx_path, '-o', xtc_out_path, '-fit', 'rot+trans']
    run_twocmd_pype(cmd1, cmd2, debug=debug)

    xtc_in_path = xtc_out_path 
    gro_out_path = os.path.splitext(os.path.basename(xtc_path))[0] + '_nowat_fit_start.gro'
    cmd1 = ['echo', 'TempCouplingA\n']
    cmd2 = [gmxbin, 'trjconv', '-s', tpr_path, '-f', xtc_in_path, '-n', ndx_path, '-o', gro_out_path, '-dump', '0']
    run_twocmd_pype(cmd1, cmd2, debug=debug)

    pdb_out_path = os.path.splitext(os.path.basename(xtc_path))[0] + '_nowat_fit_protein.pdb'
    cmd1 = ['echo', 'Protein\n']
    cmd2 = [gmxbin, 'trjconv', '-s', tpr_path, '-f', xtc_in_path, '-n', ndx_path, '-o', pdb_out_path]
    run_twocmd_pype(cmd1, cmd2, debug=debug)

    pdb_out_path = os.path.splitext(os.path.basename(xtc_path))[0] + '_nowat_fit_ligand.pdb'
    cmd1 = ['echo', 'LIG\n']
    cmd2 = [gmxbin, 'trjconv', '-s', tpr_path, '-f', xtc_in_path, '-n', ndx_path, '-o', pdb_out_path]
    run_twocmd_pype(cmd1, cmd2, debug=debug)

    return xtc_out_path

def gmx_energy_intr(edr_path, debug=False):

    basename = os.path.splitext(os.path.basename(edr_path))[0]
    xvg_comp_path = basename + '_ie_component.xvg'
    cmd1 = ['echo', 'Coul-SR:Protein-LIG\n LJ-SR:Protein-LIG\n 0\n']
    cmd2 = [gmxbin, 'energy', '-f', edr_path, '-o', xvg_comp_path]
    run_twocmd_pype(cmd1, cmd2, debug=debug)

    xvg_sum_path = basename + '_ie_sum.xvg'
    cmd1 = ['echo', 'Coul-SR:Protein-LIG\n LJ-SR:Protein-LIG\n 0\n']
    cmd2 = [gmxbin, 'energy', '-f', edr_path, '-o', xvg_sum_path, '-sum']
    run_twocmd_pype(cmd1, cmd2, debug=debug)

    return xvg2csvpng(xvg_sum_path)

def gmx_ifpcos(protein_pdb_path, ligand_pdb_path, outbasename):
    from oddt import toolkit
    from oddt import fingerprints
    import numpy as np

    protein = []
    for pro in toolkit.readfile('pdb', protein_pdb_path):
        pro.protein = True
        protein.append(pro)

    ligand = list(toolkit.readfile('pdb', ligand_pdb_path))

    IFP0 = fingerprints.SimpleInteractionFingerprint(ligand[0], protein[0])
    IFPCos = []
    for i in range(0, len(protein)):
        IFP1 = fingerprints.SimpleInteractionFingerprint(ligand[i], protein[i])
        cos_sim = - (np.dot(IFP0, IFP1) / (np.linalg.norm(IFP0) * np.linalg.norm(IFP1)))
        IFPCos.append(cos_sim)

    IFPCos = np.array(IFPCos)
    score_IFPCos = np.mean(IFPCos)

    csvfile = outbasename + '.csv'
    np.savetxt(csvfile, IFPCos)

    time = np.array([float(i) for i in range(len(IFPCos))])
    xlabel = 'Time (ps)'
    ylabel = 'IFP_cos'
    title = 'IFP_cos'
    pngfile = outbasename + '.png'
    plot2d(time, IFPCos, xlabel, ylabel, title, pngfile)

    return score_IFPCos

def gmxpla_prolig_run(edr_path, tpr_path, xtc_path, ndx_path, score_csv_path, ifpcos=False, debug=False):

    score_ie = gmx_energy_intr(edr_path, debug=debug)

    xtc_nowat_fit_path = gmx_trjconv_prolig(tpr_path, xtc_path, ndx_path, debug=debug)

    xvg_nowat_fit_rms_path = os.path.splitext(os.path.basename(xtc_nowat_fit_path))[0] + '_rms.xvg'
    selection = 'Backbone\n LIG_Heavy\n'
    score_rms = gmx_rms(tpr_path, xtc_nowat_fit_path, ndx_path, xvg_nowat_fit_rms_path, selection, debug=debug)

    xvg_nowat_fit_rmsf_path = os.path.splitext(os.path.basename(xtc_nowat_fit_path))[0] + '_rmsf.xvg'
    selection = 'LIG_Heavy\n'
    score_rmsf = gmx_rmsf(tpr_path, xtc_nowat_fit_path, ndx_path, xvg_nowat_fit_rmsf_path, selection, debug=debug)

    if ifpcos:
        outbasename_ifpcos = os.path.splitext(os.path.basename(xtc_nowat_fit_path))[0] + '_ifpcos'
        protein_pdb_path = os.path.splitext(os.path.basename(xtc_path))[0] + '_nowat_fit_protein.pdb'
        ligand_pdb_path = os.path.splitext(os.path.basename(xtc_path))[0] + '_nowat_fit_ligand.pdb'
        score_ifpcos = gmx_ifpcos(protein_pdb_path, ligand_pdb_path, outbasename_ifpcos)
        with open(score_csv_path, mode='w') as fout:
            fout.write('IE_score,RMSD_score,RMSF_score,IFPcos_score\n')
            fout.write('{},{},{},{}\n'.format(score_ie, score_rms, score_rmsf, score_ifpcos))
    else:
        with open(score_csv_path, mode='w') as fout:
            fout.write('IE_score,RMSD_score,RMSF_score\n')
            fout.write('{},{},{}\n'.format(score_ie, score_rms, score_rmsf))

def gmxpla_main(conf):
    edr_path = conf['edr']
    tpr_path = conf['tpr']
    xtc_path = conf['xtc']
    ndx_path = conf['ndx']
    score_csv_path = conf['outcsv']
    ifpcos = conf['ifpcos']
    debug = conf['verbose']
    gmxpla_prolig_run(edr_path, tpr_path, xtc_path, ndx_path, score_csv_path, ifpcos=ifpcos, debug=debug)

def main():
    args = get_parser()
    if args.verbose: print(args)

    conf = set_config(args)

    print('======= Input configulations =======')
    for k, v in conf.items():
        print('{}: {}'.format(k, v))
    print('====================================')

    gmxpla_main(conf)

if __name__ == '__main__':
    main()
