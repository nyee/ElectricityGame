"""
purpose: 
helps with inputing large data sets into training reactions.  the reactions must be in the form A*T^n*exp(-Ea/RT).

details: 
You give it a reaction family, path to the directory holding RMG-database, reactant SMILES (as an array), 
and rate information (as an array).  Optional arguments include 
references, short description, long description, and rank.  The script will append 
all the reactions to dictionary.txt, and reactions.py in the corresponding family's training reaction folder 

to use: 
import the method addTraining in python, and call the method using the arguments in the function.  
An example code exists with the name addTraining-run.py (watch out: the example will change your dictionary)

limitations:
1. current units supported are : A is in cm3/mols (or s-1 for unimolecular rxns) and E is in kcal/mol.
2. written to work with bimolecular reactions but never tested. let me know if it doesn't work

last edited by Mark on July 29th
last editted by Nate on Aug 30th
"""

from rmgpy.data.kinetics.database import *
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule.molecule import *
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.base import ForbiddenStructures
from rmgpy.data.reference import *
from rmgpy.data.kinetics.family import KineticsFamily
import os.path


def addTraining(reactionFamily, A, n, Ea, smiles1, degeneracy,
                smiles2=None,
                productSmiles1=None,
                productSmiles2 = None,
                rank=3,
                shortDesc='',
                longDesc='',
                reference=None,
                directory='',
                startingIndex=1):
    """
this method inputs a reaction into training reactions using just the reactants and reaction family.
A is in cm3/mols and E is in kcal/mol.  A, n, Ea, smiles1, and smiles2 are lists containting multiple reactions
from the family.

productSmiles are

forward is whether the training reaction is writen in the direction of the reaction family or not

This currently does not if multiple reaction paths found, or with bimolecular reactions.

Fails if resonance structures
    """

    path = os.path.join(directory, 'RMG-database/input/kinetics/families')
    # throw errors
    if len(A) != len(n) or len(n) != len(Ea) or len(n) != len(smiles1):
        # 		error('lengths of matrixes not equal')
        raise Exception('lengths of matrixes not equal')
    # load databases
    # 	kdatabase=KineticsDatabase()
    rdatabase = RMGDatabase()
    # 	rdatabase.load(os.path.join(directory,'RMG-database/input'),kineticsFamilies=[reactionFamily])
    rdatabase.load(os.path.join(directory, 'RMG-database/input'),
                   thermoLibraries=None,
                   transportLibraries=None,
                   reactionLibraries=[],
                   seedMechanisms=None,
                   kineticsFamilies=[reactionFamily],
                   kineticsDepositories=[],
                   statmechLibraries=None,
                   depository=True,
                   solvation=False,
                   )
    # 	kdatabase.loadFamilies(path,families=[reactionFamily])
    # convert the smiles into molecule objects
    productsSpecified = 0
    reactant2 = []
    products1 = []
    products2 = []
    reactant1 = [Molecule(SMILES=smiles) for smiles in smiles1]
    if smiles2: reactant2 = [Molecule(SMILES=smiles) for smiles in smiles2]
    if productSmiles1:
        products1 = [Molecule(SMILES=smiles) for smiles in productSmiles1]
        productsSpecified +=1
    if productSmiles2:
        products2 = [Molecule(SMILES=smiles) for smiles in productSmiles2]
        productsSpecified +=1

    # find possible reaction paths for reactants
    reactionList = []
    for index, smile in enumerate(smiles1):
        if smiles2:
            reactantList = [reactant1[index], reactant2[index]]
        else:
            reactantList = [reactant1[index]]
        reactions = rdatabase.kinetics.families[reactionFamily].generateReactions(reactantList)

        # throw error if multiple paths exist
        reaction = None
        if productsSpecified ==2:
            #match products
            for reaction1 in reactions:
                if products1[index].isIsomorphic(reaction1.products[0]) and products2[index].isIsomorphic(reaction1.products[1]):
                    reaction = reaction1
                    break
                elif products1[index].isIsomorphic(reaction1.products[1]) and products2[index].isIsomorphic(reaction1.products[0]):
                    reaction = reaction1
                    break
        elif productsSpecified == 1:
            for reaction1 in reactions:
                if products1[index].isIsomorphic(reaction1.products[0]):
                    reaction = reaction1
                    break
        elif len(reactions) > 1:
            raise Exception("multiple reactions found for index" + index)
        #If only one reaction path choose this one
        else: reaction = reactions[0]

        if not reaction: raise Exception("No reaction matching for index {0} reactants and products found!".format(index))

        #####now we have to do some workaround to obtain the graph values for the reactants
        ####note for BIMOLECULAR reactions: this part may mix up molecule and species lists. fix
        # convert to list of molecules from list of species to use 'KineticsFamily.applyRecipe'
        products_with_molecules = []
        for product in reaction.products:
            products_with_molecules.append(product.molecule[0])
        # find reactants with the atom labels using 'applyRecipe'
        newreactants = rdatabase.kinetics.families[reactionFamily].applyRecipe(products_with_molecules, forward = False)
        # 		print newreactants
        # convert back to species list from molecule list
        react_with_species = []
        for reactant in newreactants:
            react_with_species.append(Species(molecule=[reactant]))
        reaction.reactants = react_with_species

        # make labels for species...the = in species labels leads to erros parsing the reaction, so it is replaced with !
        reaction.reactants[0].label = reaction.reactants[0].molecule[0].toSMILES().replace('=', '!') + '' + str(
            index) + 'r1'
        reaction.products[0].label = reaction.products[0].molecule[0].toSMILES().replace('=', '!') + '' + str(
            index) + 'p1'
        if len(reaction.reactants) > 1:
            reaction.reactants[1].label = reaction.reactants[1].molecule[0].toSMILES().replace('=', '!') + '' + str(
                index) + 'r2'
        if len(reactions[0].products) > 1:
            reaction.products[1].label = reaction.products[1].molecule[0].toSMILES().replace('=', '!') + '' + str(
                index) + 'p2'

        print reaction.__str__()
        reactionList.append(reaction)
    # print reactionList

    # write to dictionary file
    with open(os.path.join(path, reactionFamily, 'training', 'dictionary.txt'), 'ab') as f:
        f.write('\n')
        for index, reaction in enumerate(reactionList):
            for species in reaction.reactants:
                f.write(species.toAdjacencyList())
                f.write('\n')
            for species in reaction.products:
                f.write(species.toAdjacencyList())
                f.write('\n')
    print 'finished writing dictionary to file'
    # create strings for putting together a reactions.py entry
    startstr = """
entry(
    index = """
    afterindex = """,
    label = \""""
    afterlabel = """\",
    degeneracy = """
    afterdegneracy = """
    kinetics = Arrhenius(
        A = ("""
    if smiles2:  # if using bimolecular reactions, need to change the units for prefactor
        afterprefactor = """, 'cm^3/(mol*s), '*|/', 2.51189),
        n = """
    else:
        afterprefactor = """, 's^-1', '*|/', 2.51189),
        n = """
    afterexponential = """,
        Ea = ("""
    afterEa = """, 'kcal/mol', '+|-', 1.5),
        T0 = (1, 'K'),
    ),
    reference = """
    afterreference = """,
    rank = """
    afterrank = """,
    referenceType = \"theory\",
    shortDesc =u\' """
    aftershortdesc = """\',
    longDesc = u\"\"\""""
    endstr = """\"\"\",
)"""

    # write to reactions.py file - 'ab' appends to bottom of the file
    with open(os.path.join(path, reactionFamily, 'training', 'reactions.py'), 'ab') as f:
        for index, reaction in enumerate(reactionList):
            f.write('\n\n')
            f.write(startstr)
            f.write(str(startingIndex + index))
            f.write(afterindex)
            f.write(reaction.__str__())
            f.write(afterlabel)
            f.write(str(degeneracy[index]))
            f.write(afterdegneracy)
            f.write(str(A[index]))
            f.write(afterprefactor)
            f.write(str(n[index]))
            f.write(afterexponential)
            f.write(str(Ea[index]))
            f.write(afterEa)
            if reference:
                f.write(reference.toPrettyRepr())
            else:
                f.write('None')
            f.write(afterreference)
            f.write(str(rank))
            f.write(afterrank)
            f.write(shortDesc)
            f.write(aftershortdesc)
            f.write(longDesc)
            f.write(endstr)

    print 'finished writing reactions.py to file'


####
if __name__ == "__main__":
    # input optional parameters to put in RMG
    shortDesc0 = u"""CBS-QB3 calculation with 1-d rotor treatment at B3LYP/631G(d)"""
    longDesc0 = u"""
    Quantum chemistry calculations CBS-QB3 calculation with 1-d rotor treatment at
    B3LYP/631G(d)" using Gaussian 03 and Gaussian 09. High-pressure-limit rate
    coefficient computed TST with Eckart Tunnelling"
    """
    reference0 = Article(
        authors=["K. Wang", "S. Villano", "A. Dean"],
        title=u'Reactions of allylic radicals that impact molecular weight growth kinetics',
        journal="Phys. Chem. Chem. Phys.",
        volume="17",
        pages="""6255-6273""",
        year="2015",
    )
    # required conditions from paper
    family = 'H_Abstraction'
    A = [18.06, 26.3, 18.4, 64.3, 29.7, 19.17, 16.9, 23.0, 120.0, 80.50, 84.7, 275.00, 333.0, 1110.0, 429.5, 427.5,
         945.0, 890.0, 21.0, 299.7, 10.5, 8.32, 667.5, 570.0, 1500, 590.0, 1.52e3, 6.07]
    n = [3.27, 3.24, 3.27, 3.13, 3.22, 3.28, 3.30, 3.22, 3.09, 3.12, 3.12, 3.08, 2.92, 2.73, 2.95, 2.93, 2.84, 2.86,
         3.21, 2.84, 3.31, 3.29, 2.90, 2.91, 2.80, 2.90, 2.93, 3.42]
    Ea = [6.85, 7.03, 7.15, 7.32, 7.04, 6.73, 6.48, 5.48, 5.73, 5.82, 5.78, 4.37, 4.91, 5.32, 5.24, 4.01, 3.72, 3.77,
          5.77, 5.63, 5.46, 5.72, 3.45, 4.67, 3.38, 5.06, 2.93, 8.66]
    reactant1smiles = ['CC=CC', 'CC=CCC', 'C=CC', 'C=C(C)CC', 'C=C(C)C', 'CC=C(C)C', 'CC=C(C)C',
                       'C=CCCC=C', 'CC=CCC', 'C=CCC', 'C=CCCC', 'C=C(C)CC',
                       'C=CC(C)C', 'C=CC(C)CC',
                       'C1C=CCC1', 'C1C=CCCC1',
                       'CC1=CCCC1', 'CC1=CCCCC1',
                       'C=CCC=C', 'C=CC=CC', 'C=CC(C)=CC', 'C=CC=C(C)C',
                       'C1=CCC=CC1', 'C1=CC=CCC1', 'CC1=CC=CC1', 'C1=CC=C1', 'c1ccccc1', 'c1cccc1']
    reactant2smiles = ['[CH3]'] * 28
    degeneracy = [6, 3, 3, 3, 6, 3, 3, 4, 2, 2, 2, 2, 1, 1, 4, 4, 1, 1, 2, 3, 3, 6, 4, 4, 1, 2, 6, 3]
    forward = True

    # print len(A)
    # print len(n)
    # print len(Ea)
    # print len(reactant1smiles)
    # print len(reactant2smiles)
    # print len(degeneracy)

    addTraining(family, A, n, Ea, reactant1smiles,
                degeneracy,
                smiles2=reactant2smiles,  # for bimolecular reactions
                rank=3,
                reference=reference0,
                shortDesc=shortDesc0,
                longDesc=longDesc0,
                directory=r'/Users/Nate/code')
