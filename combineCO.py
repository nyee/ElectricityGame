"""
This script takes a chemkin file and calculates the reverse reaction coefficient at inputted temperatures.
Then, it outputs reverse rates which seem abnormally high, usually indicative that thermo of the some
species is incorrect.

Last edited by nyee on 10/30/2016
"""

from rmgpy.chemkin import loadChemkinFile
from rmgpy.kinetics.arrhenius import Arrhenius
from rmgpy.species import Species

def printReactionLine(reaction):
    newLine = [x.label for x in reaction.reactants + reaction.products]
    return newLine

############################################
if __name__ == "__main__":
    chemkinPath = "/Users/Nate/Dropbox (MIT)/Research/butanolSurfaces/RMG models/fullSeed_noPdep3/chem_annotated.inp"
    dictionaryPath = "/Users/Nate/Dropbox (MIT)/Research/butanolSurfaces/RMG models/fullSeed_noPdep3/species_dictionary.txt"
    # chemkinPath = "/Users/Nate/code/RMG-Py/examples/rmg/minimal/chemkin/chem_annotated.inp"
    # dictionaryPath = "/Users/Nate/code/RMG-Py/examples/rmg/minimal/chemkin/species_dictionary.txt"

    adjlist1 = """
    multiplicity 3
    1 C u2 p0 c0 {2,D}
    2 O u0 p2 c0 {1,D}
    """
    adjlist2 = """
    1 C u0 p1 c-1 {2,T}
    2 O u0 p1 c+1 {1,T}
    """

    CO1 = Species().fromAdjacencyList(adjlist1)
    CO2 = Species().fromAdjacencyList(adjlist2)


    speciesList, reactionList = loadChemkinFile(chemkinPath, dictionaryPath, readComments = True)
    print "loaded chemkin file"

    newReactionList1 = []
    newReactionList2 = []
    hasboth = []
    for reaction in reactionList:
        for reactionSpecies in reaction.products+reaction.reactants:
            if reactionSpecies.isIsomorphic(CO2):
                newReactionList2.append(reaction)
            if reactionSpecies.isIsomorphic(CO1):
                newReactionList1.append(reaction)

    #find both
    for reaction in newReactionList1:
        if reaction in newReactionList2:
            hasboth.append(reaction)

    print "The has boths"
    print len(hasboth)
    for reaction in hasboth:
        print printReactionLine(reaction)
    print "length r1:", len(newReactionList1)
    print "length r2:", len(newReactionList2)

    duplicates =[]
    for reaction1 in newReactionList1:
        copy1 = reaction1.reactants+reaction1.products
        speciesToRemove = None
        for species1 in copy1:
            if species1.isIsomorphic(CO1):
                speciesToRemove = species1
                break
        if speciesToRemove: copy1.remove(species1)
        else: raise Exception("some reaction doesn't have CO1")
        print [species1.label for species1 in copy1]
        for reaction2 in newReactionList2:
            #check forward
            copy2 = reaction2.reactants+reaction2.products
            matches = []
            for species1 in copy1:
                if not species1 in copy2:
                    break
            #we hit a match
            else: duplicates.append((reaction1, reaction2))

    print "\n"
    print
    #
    print "duplicates"
    for duplicate in duplicates:
        printReactionLine(duplicate[0])
        printReactionLine(duplicate[1])
