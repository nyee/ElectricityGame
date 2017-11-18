"""
This script is written to make forbidden groups for intra-H-migration of substituted benzenes
"""
import re
from rmgpy.molecule.group import Group, GroupAtom, GroupBond
from rmgpy.molecule.atomtype import atomTypes
from rmgpy.molecule.pathfinder import find_shortest_path

def addAtomsToBases(bases, aromaticPosition, site1Extension, site2Extension):
    """
    site1 is related to the radical
    site2 is related to the H
    """

    names = []
    baseGroups = [Group().fromAdjacencyList(adjlist) for adjlist in bases]
    for baseIndex, group in enumerate(baseGroups):
        site1 = 0 #atomIndex
        if aromaticPosition == 'ortho':
            site2 = 1
            distance = 0
            name = "benzene_ortho"
        elif aromaticPosition == 'meta':
            site2 = 2
            distance = 1
            name = "benzene_meta"
        elif aromaticPosition == 'para':
            site2 = 3
            distance = 2
            name = "bezene_para"
        else:
            raise Exception("Invalid aromatic Position")

        distance = distance + site1Extension + site2Extension
        name+=str(site1Extension)+"_"
        name+=str(site2Extension)
        name+="_r"+str(baseIndex)
        names.append(name)

        #extend site1
        currentAtomIndex = site1
        for x in range(site1Extension):
            newAtom =GroupAtom(atomType=[atomTypes['R!H']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
            newBond = GroupBond(newAtom, group.atoms[currentAtomIndex], ['S'])
            group.addAtom(newAtom)
            group.addBond(newBond)
            currentAtomIndex = group.atoms.index(newAtom)
        #at end change the last added atom to, and label
        else:
            site1Atom = group.atoms[currentAtomIndex]
            site1Atom.radicalElectrons = [1]
            site1Atom.label = "*1"

        currentAtomIndex = site2
        #extend site2
        for x in range(site2Extension):
            newAtom =GroupAtom(atomType=[atomTypes['R!H']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
            newBond = GroupBond(newAtom, group.atoms[currentAtomIndex], ['S'])
            group.addAtom(newAtom)
            group.addBond(newBond)
            currentAtomIndex = group.atoms.index(newAtom)
        #at end, add one more H atom
        else:
            site2Atom = group.atoms[currentAtomIndex]
            site2Atom.label = "*2"
            newAtom =GroupAtom(atomType=[atomTypes['H']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
            newBond = GroupBond(newAtom, site2Atom, ['S'])
            group.addAtom(newAtom)
            group.addBond(newBond)
            currentAtomIndex = group.atoms.index(newAtom)
            group.atoms[currentAtomIndex].label = "*3"

        #Number the in between atoms based on distance
        path = find_shortest_path(site1Atom, site2Atom)
        assert len(path) -2 == distance
        if distance == 1:
            distanceLabel = ["*4"]
        elif distance == 2:
            distanceLabel = ["*4", "*5"]
        elif distance == 3:
            distanceLabel = ["*4", "*6", "*5"]
        elif distance == 4:
            distanceLabel = ["*4", "*6", "*7", "*5"]
        elif distance == 5:
            distanceLabel = ["*4", "*6", "*7", "*8", "*5"]
        elif distance == 0:
            distanceLabel=[]
        else:
            raise Exception("Distance too large for {0} with site extensions {1} and {2}".format(aromaticPosition, site1Extension, site2Extension))
        for atomIndex, atom in enumerate(path[1:-1]):
            atom.label = distanceLabel[atomIndex]


    return names, baseGroups
####################################################################################
if __name__ == "__main__":

    adjlist0 = """
    1   C u0 {2,B} {6,B}
    2   C u0 {1,B} {3,B}
    3   C u0 {2,B} {4,B}
    4   C u0 {3,B} {5,B}
    5   C u0 {4,B} {6,B}
    6   C u0 {1,B} {5,B}
"""

    adjlist1 = """
    1   C u0 {2,D} {6,S}
    2   C u0 {1,D} {3,S}
    3   C u0 {2,S} {4,D}
    4   C u0 {3,D} {5,S}
    5   C u0 {4,S} {6,D}
    6   C u0 {1,S} {5,D}
"""

    adjlist2 = """
    1   C u0 {2,S} {6,D}
    2   C u0 {1,S} {3,D}
    3   C u0 {2,D} {4,S}
    4   C u0 {3,S} {5,D}
    5   C u0 {4,D} {6,S}
    6   C u0 {1,D} {5,S}
"""

    bases = [adjlist0, adjlist1, adjlist2]
    names = []
    comment = "Forbidden because TS would be too strained. May be resonance structure of other groups."

    sitesToForbid={'ortho': [(0,0), (0,1), (1,1)],
                   'meta': [(0,0), (0,1), (1,1), (2,0), (2,1), (2,2), (3,0)],
                   'para':[(0,0), (0,1), (1,1), (2,0), (2,1), (3,0)]}

    # test = addAtomsToBases(bases, 'ortho', 0, 1)
    forbiddenGroups = []
    for aromaticPosition, extensions in sitesToForbid.iteritems():
        for x in extensions:
            newNames, newGroups = addAtomsToBases(bases, aromaticPosition, x[0], x[1])
            forbiddenGroups.extend(newGroups)
            names.extend(newNames)
            #then do for reverse
            if not x[0] == x[1]:
                newNames, newGroups = addAtomsToBases(bases, aromaticPosition, x[1], x[0])
                forbiddenGroups.extend(newGroups)
                names.extend(newNames)

    with open("/Users/Nate/code/newTest/newForbiddenGroups.py", 'wb') as outFile:
        for index, group in enumerate(forbiddenGroups):
            outFile.write("forbidden(\n")
            outFile.write("    label = \"{0}\",\n".format(names[index]))
            outFile.write("    group = \n")
            outFile.write("\"\"\"\n")
            outFile.write(group.toAdjacencyList())
            outFile.write("\"\"\",\n")
            outFile.write("    shortDesc = u\"\"\"\"\"\",\n")
            outFile.write("    longDesc = \n")
            outFile.write("\"\"\"\n")
            outFile.write(comment)
            outFile.write("\"\"\",\n")
            outFile.write(")\n")
            outFile.write("\n")

