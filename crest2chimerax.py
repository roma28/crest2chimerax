import argparse
import dataclasses
import logging
import math
import os

R = 1.9872036e-3  # kcal/mol
T = 298.15  # K
Eh_to_kcalpermol = 625.5


@dataclasses.dataclass
class Conformer:
    fname: str
    energy_Eh: float
    energy_relative_kcalpermol: float = 0
    population_relative: float = 0


logging.basicConfig(level=logging.DEBUG)

argparser = argparse.ArgumentParser(prog='crest2chimerax',
                                    description='A script that takes CREST .xyz conformer ensemble and produces '
                                                'ChimeraX script to visualize', )
argparser.add_argument('infile', help='CREST file')
argparser.add_argument('outfile', help='ChimeraX .cxc file')
argparser.add_argument('-r', '--refatoms', default='*',
                       help='Refatoms string as used in ChimeraX align command (* by default)')
argparser.add_argument('-e', '--energy', type=float, default=6.0,
                       help='Cutoff energy in kcal/mol (6 kcla/mol by default)')
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
    answer = input("Continue? ")
    if not answer.upper() in ["Y", "YES"]:
        exit(-1)

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

del lines

e_min = min(conformers, key=lambda x: x.energy_Eh).energy_Eh
for c in conformers:
    c.energy_relative_kcalpermol = (c.energy_Eh - e_min) * Eh_to_kcalpermol
    c.population_relative = math.exp(-c.energy_relative_kcalpermol / R / T)

conformers = filter(lambda x: x.energy_Eh < e_min, conformers)

with open(outfile, 'w') as outfile:
    outfile.writelines(f'cd {workdir}\n')
    for i, c in enumerate(conformers):
        outfile.writelines(f'open {c.fname}\n')
        outfile.writelines(f'transparency #{i + 1} {100 - c.population_relative * 100} target ab\n')
        if i == 0:
            continue
        outfile.writelines(f'align #{i + 1}@{args.refatoms} toAtoms #1@{args.refatoms}\n')
