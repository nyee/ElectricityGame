from rmgpy.data.rmg import RMGDatabase
from rmgpy import settings
from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.atomtype import atomTypes
from rmgpy.molecule.group import Group, ImplicitBenzeneError, GroupBond, GroupAtom, UnexpectedChargeError
from rmgpy.species import Species
import logging
from rmgpy.data.base import DatabaseError



if __name__ == "__main__":
#     family_name = "Disproportionation"
#     # carbonMonosulfide = Molecule(SMILES="[C-]#[S+]")
#     test =Group().fromAdjacencyList("""
# 1 * Ct u0 p1 c-1 {2,T}
# 2   Ot u0 p1 c+1 {1,T}""")
#
#     test =Molecule().fromAdjacencyList("""
# 1  C u0 p0 c0 {2,B} {6,B} {7,S}
# 2  C u0 p0 c0 {1,B} {3,B} {8,S}
# 3  C u0 p0 c0 {2,B} {4,B} {9,S}
# 4  C u0 p0 c0 {3,B} {5,B} {10,S}
# 5  C u0 p0 c0 {4,B} {6,B} {11,S}
# 6  C u0 p0 c0 {1,B} {5,B} {12,S}
# 7  H u0 p0 c0 {1,S}
# 8  H u0 p0 c0 {2,S}
# 9  H u0 p0 c0 {3,S}
# 10 H u0 p0 c0 {4,S}
# 11 H u0 p0 c0 {5,S}
# 12 H u0 p0 c0 {6,S}""")

    # test = Species().fromAdjacencyList(test.toAdjacencyList())

    # print test.toAdjacencyList()
    # print test.toSMILES()

    # try:
    #     sample = test.makeSampleMolecule()
    #     print sample.toAdjacencyList()
    # except UnexpectedChargeError, e:
    #     print "UnExpected Charge"
    #     print e.graph.toAdjacencyList()
#
#     # print matchNodeToStructure(test, sample, sample.getLabeledAtoms())
#
#
#     resonance = sample.generateResonanceIsomers()
#     print resonance
#     # test4=test3.addImplicitBenzene()
#     # print test4.toAdjacencyList()


    family_name = "R_Addition_MultipleBond"
    databaseDirectory = settings['database.directory']
    database = RMGDatabase()
    database.load(databaseDirectory, kineticsFamilies=[family_name], thermoLibraries=['primaryThermoLibrary'])
    print "database loaded"

    # print database.kinetics.families[family_name].getReactionTemplate(database.kinetics.families[family_name].depositories[0].entries[73].item)
    family = database.kinetics.families[family_name]


