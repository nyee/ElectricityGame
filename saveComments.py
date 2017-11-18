"""
This script saves all comments from a python style file and saves it in output.
When we save database files, we often lose all the comments that were originally there.

Technically its much better if all comments are written out in an entries longDesc, but
you never know what somebody will do
"""
import os.path
from os import walk
import re

databaseDirectory="/Users/Nate/code/RMG-database/input"

families=['Intra_R_Add_Exocyclic',
          'Intra_RH_Add_Endocyclic',
          'Intra_R_Add_Endocyclic',
          'Intra_R_Add_ExoTetCyclic',
          'Intra_Disproportionation',
          'intra_substitutionS_cyclization',
          'intra_H_migration',
          'Intra_RH_Add_Exocyclic',
          'Cyclic_Ether_Formation',
          'Birad_recombination',
          'intra_substitutionCS_cyclization',
          '1,2_shiftS',
          'intra_OH_migration',
          'intra_substitutionCS_isomerization',
          'intra_substitutionS_isomerization',
          'Diels_alder_addition'
          ]

output = "/Users/Nate/code/RMG-database/comments.txt"

comments = []
#comments to ignore
ignore = ["#!/usr/bin/env python",
          "# encoding: utf-8"]

for family in families:
    familyDirectory = os.path.join(databaseDirectory, 'kinetics/families/{0}'.format(family))
    files = []
    for (dirpath, dirnames, filenames) in walk(familyDirectory):
        files.extend([os.path.join(dirpath, filename) for filename in filenames])
        # print (dirpath, dirnames, filenames)
    # print files
    for filename in files:
        if os.path.splitext(filename)[1] == '.py':
            with open(filename, 'rb') as inputfile:
                for linenumber, line in enumerate(inputfile):
                    if '#' in line:
                        for ignoredLine in ignore:
                            if re.match(re.escape(ignoredLine), line):
                                break
                        else:
                            comments.append([filename, linenumber+1, line])

with open(output, 'wb') as outputFile:
    for index, comment in enumerate(comments):
        if index > 0 and not comment[0] == comments[index-1][0]:
            pass
        else:
            outputFile.write("In {0}:\n".format(comment[0]))
        outputFile.write("line {0}: {1}".format(comment[1],comment[2]))