# crest2chimerax

This a script to create .cxc command file for ChimeraX
to align the conformers from CREST calculation at selected
atoms and sets transparency proportionally to the conformer
population.

# Usage
The basic usage is 
```
python crest2chimerax.py crest_conformers.xyz out.cxc
```
it will import all the conformers with relative
energy < 6 kcal/mol, align them using `*` as `refatoms` and set
their transparency equal to their population.

You can modify energy cutoff by using `-e` or `--energy` option
and reference atoms for alignment using `-r` or `--reference`
option.

You can use `-p` or `--population` option to limit displayed
conformers to a fraction of total ensemble (useful for
an ensemble with a couple of major conformers and lots of 
barely populated ones)