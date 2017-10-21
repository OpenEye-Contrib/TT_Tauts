#!/usr/bin/env python
# reads output from gen_t_skel and creates pngs of compounds with different names
# but the same t_skel.

import sys
from openeye.oechem import *
from openeye.oedepict import *

################################################################################
def write_same_t_skels_to_pdf(t_skels, pdf_name):
    # take a set of molecules with the same t_skel and write them to a pdf
    rows = 6
    cols = 3
    reportopts = OEReportOptions(rows, cols)
    reportopts.SetPageOrientation(OEPageOrientation_Portrait)
    reportopts.SetCellGap(20)
    reportopts.SetPageMargins(20)
    report = OEReport(reportopts)

    opts = OE2DMolDisplayOptions(report.GetCellWidth(), report.GetCellHeight(), OEScale_AutoScale)
    bond_pen = opts.GetDefaultBondPen()
    # print('Current bond pen width : {}'.format(bond_pen.GetLineWidth()))
    # the default width in PDF is too fat for me.
    bond_pen.SetLineWidth(1)

    for t_skel in t_skels:
        for mol in t_skel:
            cell = report.NewCell()
            disp = OE2DMolDisplay(mol, opts)
            OERenderMolecule(cell, disp)
        
    OEWriteReport(pdf_name, report)
    

################################################################################
def total_h_count(mol):
    tot_h = 0
    for atom in mol.GetAtoms():
        tot_h += atom.GetTotalHCount()

    return tot_h

################################################################################
t_skels = {}
in_smis = {}

with open(sys.argv[1], "r") as fp:
    all_lines = fp.readlines()
    print("Number of lines {}".format(len(all_lines)))
    for nextline in all_lines:
        if nextline.startswith('T_Skel') :
            bits = nextline.split()
            name = bits[2]
            in_smi = bits[4]
            t_skel = bits[6]
            # print('{} {} {}'.format(name, in_smi, t_skel))

            if t_skel not in t_skels:
                t_skels[t_skel] = []
            t_skels[t_skel].append(name)
            if t_skel not in in_smis:
                in_smis[t_skel] = []
            in_smis[t_skel].append(in_smi)

print('Number of t_skels : {}'.format(len(t_skels)))

dups = []
for t_skel, names in t_skels.items():
    if len(names) > 1:
        # print("{} : {} : {}".format(t_skel, names, in_smis[t_skel]))
        dups.append((t_skel,(names[0], in_smis[t_skel][0]),(names[1], in_smis[t_skel][1])))

print('Number of duplicates : {}'.format(len(dups)))

same_smis = []
diff_smis = []

for i in range(len(dups)):
    # print(i, len(dups))
    dup = dups[i]
    # print(dup)
    ts_mol = OEGraphMol()
    OESmilesToMol(ts_mol, dup[0])
    OEPrepareDepiction(ts_mol)

    mol1 = OEGraphMol()
    OESmilesToMol(mol1, dup[1][1])
    mol1.SetTitle(dup[1][0])
    OEPrepareDepiction(mol1)

    mol2 = OEGraphMol()
    OESmilesToMol(mol2, dup[2][1])
    mol2.SetTitle(dup[2][0])
    OEPrepareDepiction(mol2)

    smi1 = OECreateSmiString(mol1)
    smi2 = OECreateSmiString(mol2)
    if smi1 != smi2:
        tot_h1 = total_h_count(mol1)
        tot_h2 = total_h_count(mol2)
        if tot_h1 != tot_h2:
            mol1.SetTitle(mol1.GetTitle() + ' : {}H'.format(tot_h1))
            mol2.SetTitle(mol2.GetTitle() + ' : {}H'.format(tot_h2))
        diff_smis.append((ts_mol, mol1, mol2))
    else:
        same_smis.append((ts_mol, mol1, mol2))
    
print('Number of duplicates with same non-isomeric SMILES : {}'.format(len(diff_smis)))
fn = sys.argv[1] + "_t_skel_dups_same_smi.pdf"
write_same_t_skels_to_pdf(same_smis, fn)

print('Number of duplicates with different non-isomeric SMILES : {}'.format(len(same_smis)))
fn = sys.argv[1] + "_t_skel_dups_diff_smi.pdf"
write_same_t_skels_to_pdf(diff_smis, fn)

