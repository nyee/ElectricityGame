from rmgpy.data.rmg import RMGDatabase
from rmgpy import settings
from rmgpy.molecule.group import Group
from rmgpy.data.thermo import ThermoDatabase, ThermoGroups
from collections import OrderedDict
import os
import re
import copy
from rmgpy.molecule.atomtype import atomTypes

"""
Working script to make database wide changes to groups
"""


database = RMGDatabase()
database.load(settings['database.directory'], thermoLibraries = [], kineticsFamilies='all', kineticsDepositories=[], reactionLibraries=[])
path=settings['database.directory']
groupName="group"

# targetLabel=['S', 'Sd']
# replace1=atomTypes['CS']
# targetLabel=['O', 'Od']
# replace1=atomTypes['CO']
targetLabel=['N1d']
targetAtomTypes=[atomTypes[x] for x in targetLabel]
requiredLonePairs = [1]

ignore=[]
#fix thermo
# for databaseName, thermoGroup in database.kinetics.families.iteritems():
for databaseName, thermoGroup in database.thermo.groups.iteritems():
    databaseModify=0
    print databaseName
    if hasattr(thermoGroup, 'ownReverse'):
        if not thermoGroup.ownReverse:
            for product in thermoGroup.forwardTemplate.products:
                ignore.append(product)
                ignore.extend(product.children)
    for entryName, entry in thermoGroup.entries.iteritems():
    # for entryName, entry in thermoGroup.groups.entries.iteritems():
        if entry in ignore: continue
        matched=False
        if isinstance(entry.item, Group):
            targetList=[]
            for groupAtom in entry.item.atoms:
                for atomtype1 in groupAtom.atomType:
                    if atomtype1 in targetAtomTypes:
                        matched=True
                        targetList.append(groupAtom)
        if matched:
            databaseModify+=1
            # print entry.label
            for targetAtom in targetList:
                targetAtom.lonePairs = requiredLonePairs
                # for atom2, bonds12 in targetAtom.bonds.iteritems():
                #     if atomTypes['Cb'] in atom2.atomType or atomTypes['Cb'] in atom2.atomType:
                #         if not 'B' in bonds12.order:
                #             print entryName
                # modify=''
                # removeLigands=[]
                # num_of_Dbonds=sum([1 if x.order[0] is 'D' and len(x.order)==1 else 0 for x in cdAtom.bonds.values()])
                # if num_of_Dbonds == 2:
                #     print entry.label
                #     # print entry.item.toAdjacencyList()
                #     modify='Cdd'
                #     databaseModify+=1
                #     continue
                # elif num_of_Dbonds == 1:
                #     for ligand, bond in cdAtom.bonds.iteritems():
                #         #Ignore ligands that are not double bonded
                #         if 'D' in bond.order and len(bond.order) == 1:
                #             #Ignore if there are several atomTypes
                #             if len(ligand.atomType) == 1:
                #                 for atomtype1 in ligand.atomType:
                #                     if atomtype1 is atomTypes['O'] or atomtype1 in atomTypes['O'].specific:
                #                         print entry.label
                #                         # print entry.item.toAdjacencyList()
                #                         modify='CO'
                #                         removeLigands.append(ligand)
                #                         databaseModify+=1
                #
                #                     elif atomtype1 is atomTypes['S'] or atomtype1 in atomTypes['S'].specific:
                #                         print entry.label
                #                         # print entry.item.toAdjacencyList()
                #                         modify='CS'
                #                         cdAtom.atomType=[atomTypes['CS']]
                #                         removeLigands.append(ligand)
                #                         databaseModify+=1
                #
                #     if modify:
                #         for ligand in removeLigands:
                #             #Don't remove the ligand if it is a starred atom
                #             if ligand.label == '':
                #                 del cdAtom.bonds[ligand]
                #                 entry.item.atoms.remove(ligand)
                #         #if other atomTypes currently on the ligand, just remove
                #         if len(cdAtom.atomType)>1:
                #             cdAtom.atomType.remove(atomTypes['Cd'])
                #         else:
                #             cdAtom.atomType=[atomTypes[modify]]

    if databaseModify>0:
        print databaseName, databaseModify
        thermoGroup.save(os.path.join(path, "thermo/groups/"+databaseName+'.py'))
        # kineticsPath="kinetics/families/"+databaseName
        # thermoGroup.save(os.path.join(path, kineticsPath))
