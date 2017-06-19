#!/usr/bin/env python
# reads output from gen_t_skel and creates pngs of compounds with different names
# but the same t_skel.

import sys
from openeye.oechem import *
from openeye.oedepict import *

t_skels = {}
in_smis = {}

with open(sys.argv[1], "r") as fp:
    all_lines = fp.readlines()
    print("Number of lines {}".format(len(all_lines)))
    for nextline in all_lines:
        t_skel, in_smi, name = nextline.split()
        if t_skel not in t_skels:
            t_skels[t_skel] = []
        t_skels[t_skel].append(name)
        if t_skel not in in_smis:
            in_smis[t_skel] = []
        in_smis[t_skel].append(in_smi)

dups = []
for t_skel, names in t_skels.items():
    if len(names) > 1:
        # print("{} : {} : {}".format(t_skel, names, in_smis[t_skel]))
        dups.append((t_skel,(names[0], in_smis[t_skel][0]),(names[1], in_smis[t_skel][1])))


image = OEImage(1000, 2000)
grid = OEImageGrid(image, 5, 3)
opts = OE2DMolDisplayOptions(grid.GetCellWidth(), grid.GetCellHeight(), OEScale_AutoScale)

file_num = 1
file_row = 1

for i in range(len(dups)):
    print(i, len(dups))
    dup = dups[i]
    
    print(dup)
    ts_mol = OEGraphMol()
    OESmilesToMol(ts_mol, dup[0])
    OEPrepareDepiction(ts_mol)
    disp = OE2DMolDisplay(ts_mol, opts)
    cell = grid.GetCell(file_row, 1)
    OERenderMolecule(cell, disp)

    mol1 = OEGraphMol()
    OESmilesToMol(mol1, dup[1][1])
    mol1.SetTitle(dup[1][0])
    OEPrepareDepiction(mol1)
    disp = OE2DMolDisplay(mol1, opts)
    cell = grid.GetCell(file_row, 2)
    OERenderMolecule(cell, disp)

    mol2 = OEGraphMol()
    OESmilesToMol(mol2, dup[2][1])
    mol2.SetTitle(dup[2][0])
    OEPrepareDepiction(mol2)
    disp = OE2DMolDisplay(mol2, opts)
    cell = grid.GetCell(file_row, 3)
    OERenderMolecule(cell, disp)

    file_row += 1
    if file_row == 6:
        fn = "t_skel_dups_{}.png".format(file_num)
        OEWriteImage(fn, image)
        file_num += 1
        file_row = 1
        image = OEImage(1000, 2000)
        grid = OEImageGrid(image, 5, 3)
