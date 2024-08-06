import argparse
import dataclasses
import logging
import math
import os

R = 1.9872036e-3
T = 298.15


@dataclasses.dataclass
class Conformer:
    fname: str
    energy_Eh: float
    energy_relative_kcalpermol: float = 0
    population_relative: float = 0


logging.basicConfig(level=logging.DEBUG)

argparser = argparse.ArgumentParser(prog='crest2chimerax',
                                    description='A script that takes CREST conformer ensemble and produces ChimeraX '
                                                'script to visualize', )
argparser.add_argument('infile', help='CREST file')
argparser.add_argument('outfile')
argparser.add_argument('-r', '--refatoms', default='*', help='refatoms string as used in ChimeraX align command')
args = argparser.parse_args()

infile = os.path.abspath(args.infile)
outfile = os.path.abspath(args.outfile)

workdir = os.path.dirname(infile)
os.chdir(workdir)
logging.info(f'Working in {workdir}')

try:
    os.mkdir('./conformers')
except FileExistsError:
    logging.error('Directory "conformers" already exists and will be overwritten!')

conformers = []

with open(infile, 'r') as infile:
    lines = infile.readlines()
    n_atoms = int(lines[0])
    n_confromers = int(len(lines) / (n_atoms + 2))
    logging.debug(f'{n_atoms} atoms in molecule')
    logging.debug(f'{len(lines)} lines in file => {n_confromers} conformers')

    for i in range(n_confromers):
        conformers.append(Conformer(energy_Eh=float(lines[i * (n_atoms + 2) + 1]), fname=f'./conformers/c{i}.xyz'))
        logging.debug(f'Conformer {i} with energy {conformers[i].energy_Eh}')
        with open(conformers[i].fname, 'w') as c:
            logging.debug(f'Writing conformer {i} in file {os.path.abspath(f'./conformers/c{i}.xyz')}')
            c.writelines(lines[i * (n_atoms + 2):(i + 1) * (n_atoms + 2)])

e_min = min([c.energy_Eh for c in conformers])
for c in conformers:
    c.energy_relative_kcalpermol = (c.energy_Eh - e_min) * 625.5
    c.population_relative = math.exp(-c.energy_relative_kcalpermol / R / T)

with open(outfile, 'w') as outfile:
    outfile.writelines(f'cd {workdir}\n')
    for i, c in enumerate(conformers):
        outfile.writelines(f'open {c.fname}\n')
        outfile.writelines(f'transparency #{i + 1} {100 - c.population_relative * 100} target ab\n')

        if i == 0:
            continue
        outfile.writelines(f'align #{i+1}@{args.refatoms} toAtoms #1@{args.refatoms}\n')
