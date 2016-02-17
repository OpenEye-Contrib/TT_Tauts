#!/bin/bash

for smi_file in `ls [0-9]*.smi agl[12].smi chebi2295_small.smi` ; do
	echo "XXXXXXXXXXXXX" , $smi_file
	/home/cosgrove/SVN_WORK/TT_Tauts/exe_DEBUG/tt_tauts_batch $smi_file | grep -e "AWOOGA" -e "dodgy SMILES"
	echo "YYYYYYYYYYYYY"
done
