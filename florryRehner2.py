import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize as opt

#global variables
colors = ['-k', '-r', '-g', '-b', '--r', '--g', '--b', '-k', '-k']

savePath = "/Users/Nate/Dropbox (MIT)/Research/Peng/biooils/Figures/FlorryRehner.png"


def getPureChi(species1, species2, VDict, solDict, R, T):
    diffSquared = (solDict[species1] - solDict[species2])**2
    chi = diffSquared * VDict[species1]/ (R*T)
    return chi

def getElasticChemPotI(phiE, Vi, pE,):
    """
    Returns the elastic chemical potential for species i, normalized by RT

    :param phiE: guess at elatomer volume fraction
    :param Vi: molar volume of solvent species in cm^3/mol
    :param pE: crosslink density of elatomer
    :return:
    """
    return pE*Vi*(phiE**(1.0/3.0)-phiE/2.0)

def getMixingChemPotI(x, species_i, allSpecies, VDict, chiDict):
    """
    Returns the mixing chemical potential for species i, normalized by RT

    :param x: list of floats, guesses for volume fractions of species in the phase
    :param species_i: str defining the species which we are calculating the chemical potential for
    :param allSpecies: list of all species, defines the indexing of the species, and must match x's indexing
    :param VDict: Dictionary of molar volumes
    :param chiDict: Dictionary of chi values. keys: tuple in indexed order of two species, value: chi
    """
    print "start iteration"
    #determine species i index
    i = allSpecies.index(species_i)

    #species_i non-cross terms
    term1 = math.log(x[i]) + (1- x[i])
    print "term1 {0}".format(term1)

    #cross terms of i not equal j
    term2 = 0
    for j, species_j in enumerate(allSpecies[0:-1]): #Don't need to do it for elatomer as Vi/Velatomer <<<<1
        if i==j: continue
        term2 += -VDict[species_i]/VDict[species_j] * x[j] #Make sure each term is negative!
    print "term2 {0}".format(term2)

    #cross terms i<j
    term3 = 0
    for j, species_j in enumerate(allSpecies):
        if i < j:
            term3 += chiDict[(species_i, species_j)] * x[j] * (1 - x[i])
    print "term3 {0}".format(term3)

    #cross terms h not equal i and j not equal i
    term4 = 0
    for h, species_h in enumerate(allSpecies):
        if h == i: continue
        for j, species_j in enumerate(allSpecies[h+1:]):
            if j == i: continue
            term4 += -chiDict[(species_h, species_j)] * VDict[species_i]/VDict[species_h] * x[h] * x[j] #Make sure each term is negative!
    print "term4 {0}".format(term4)

    #cross terms h < i
    term5 = 0
    for h, species_h in enumerate(allSpecies[0:i]):
        term5 += chiDict[species_h, species_i] * VDict[species_i]/VDict[species_h] * x[h] * (1 - x[i])
    print "term5 {0}".format(term5)

    return term1 + term2 + term3 + term4 + term5

def getEquilibriumRoots(phiE, newSpecies, VDict, chiDict, phiS, pE):
    """

    ai3 = chi_i3
    bi3 = chi_12 - chi_i3
    chi_Ei3 = ai3 + bi3*phi_jm
    """
    #unpack args
    species = args[0]
    VDict = args[1]
    chiDict = args[2]
    phiS = args[3]
    pE = args[4]

    #calculate liquid phase chemical potential for each species i
    liquidChemPotList = []
    for species_i in species[:-1]:
        #don't input elastomer in speciesList because it does not exist in liquid phase
        # liquidChemPotList.append(getMixingChemPotI(phiS, species_i, species[0:-1], VDict, chiDict))
        liquidChemPotList.append(getLiquidChemPotI(phiS, species_i, species[0:-1], VDict, chiDict))


    #calculate elatomer phase chemical potential for each species i
    elastomerChemPotList = []
    for i, species_i in enumerate(species[:-1]):
        chemPot = getMixingChemPotI(phiE, species_i, species, VDict, chiDict)
        chemPot += getElasticChemPotI(phiE[-1], VDict[species_i], pE,)
        elastomerChemPotList.append(chemPot)

    roots = []
    for i, species_i in enumerate(species[:-1]):
        roots.append(liquidChemPotList[i] - elastomerChemPotList[i])

    #One more equation for total volume fraction = 1
    ResidualFraction = 1
    for i in range(len(phiE)):
        ResidualFraction += -phiE[i]
    roots.append(ResidualFraction)

    return sum([root**2 for root in roots]) #try as an optimization of least squares
    # return roots

def constructArgs(species, VDict, solDict, phiSDict, pE, R, T):
    """
    creates the bundle of args to send to getEquilibriumRoots

    species: list of species, the last one must be the elastomer
    phiS: list of volume fractions in solvent phase, should have one less entry than species

    """

    #Sort solvent species by molar volume:
    newSpecies = species[0:-1]
    newSpecies = sorted(newSpecies, key = lambda spc: VDict[spc] )
    print newSpecies
    newSpecies.append(species[-1])

    phiS = [phiSDict[species_i] for species_i in newSpecies[:-1]]

    #construct chiDict
    chiDict={}
    for i, species_i in enumerate(newSpecies[0:-1]):
        for species_j in newSpecies[i+1:]:
            chiDict[(species_i, species_j)]=getPureChi(species_i, species_j, VDict, solDict, R, T)
    print chiDict
    args = [newSpecies, VDict, chiDict, phiS, pE]
    return args

def getLiquidChemPotI(x, species_i, allSpecies, VDict, chiDict):
    """
    Returns the mixing chemical potential for species i, normalized by RT

    :param x: list of floats, guesses for volume fractions of species in the phase
    :param species_i: str defining the species which we are calculating the chemical potential for
    :param allSpecies: list of all species, defines the indexing of the species, and must match x's indexing
    :param VDict: Dictionary of molar volumes
    :param chiDict: Dictionary of chi values. keys: tuple in indexed order of two species, value: chi
    """
    print "start iteration"
    #determine species i index
    i = allSpecies.index(species_i)

    #species_i non-cross terms
    term1 = math.log(x[i]) + (1- x[i])
    print "term1 {0}".format(term1)

    #cross terms of i not equal j
    term2 = 0
    for j, species_j in enumerate(allSpecies): #Don't need to do it for elatomer as Vi/Velatomer <<<<1
        if i==j: continue
        term2 += -VDict[species_i]/VDict[species_j] * x[j] #Make sure each term is negative!
    print "term2 {0}".format(term2)

    #cross terms i<j
    term3 = 0
    for j, species_j in enumerate(allSpecies):
        if i < j:
            term3 += chiDict[(species_i, species_j)] * x[j] * (1 - x[i])
    print "term3 {0}".format(term3)

    #cross terms h not equal i and j not equal i
    term4 = 0
    for h, species_h in enumerate(allSpecies):
        if h == i: continue
        for j, species_j in enumerate(allSpecies[h+1:]):
            if j == i: continue
            term4 += -chiDict[(species_h, species_j)] * VDict[species_i]/VDict[species_h] * x[h] * x[j] #Make sure each term is negative!
    print "term4 {0}".format(term4)

    #cross terms h < i
    term5 = 0
    for h, species_h in enumerate(allSpecies[0:i]):
        term5 += chiDict[species_h, species_i] * VDict[species_i]/VDict[species_h] * x[h] * (1 - x[i])
    print "term5 {0}".format(term5)

    return term1 + term2 + term3 + term4 + term5
if __name__ == '__main__':


    #dictionary of physical parameters
    #cm^3/mol
    VDict = {'isoOctane': 162.48,
             'toluene': 106.27,
             'anisole': 108.62,
             'phenol': 90.49,
             'furan': 72.64,
             'furfural': 82.83}
    #MPa^0.5
    solDict = {'isoOctane': 14.169,
             'toluene': 18.028,
             'anisole': 19.552,
             'phenol': 24.948,
             'furan': 18.674,
             'furfural': 23.970,
               'bunaN': 21.24}

    # phiSDict = {'toluene': 0.45,
    #             'isoOctane': 0.52,
    #             'phenol': 0.03}
    pE = 1.17E-3 #mol/cm^3
    R = 8.314 #J/mok-K
    T = 298.15 #K
    #guess

    #test fuelC
    volFractionC = 0.465 #volume fraction toluene in C
    additiveFraction = 0.35

    phiSDict = {'toluene': (1-additiveFraction)*volFractionC,
                'isoOctane': (1-additiveFraction)*(1.0 - volFractionC),
                'phenol': additiveFraction}
    species = ['isoOctane', 'toluene', 'phenol', 'bunaN']

    # phiSDict = {'toluene': volFractionC,
    #             'isoOctane':(1.0 - volFractionC),}
    # species = ['isoOctane', 'toluene', 'bunaN']

    # phiSDict = {'anisole': 1}
    # species = ['anisole', 'bunaN']

    # phiSDict ={'': 1}
    # species = ['', 'bunaN']
    args = tuple(constructArgs(species, VDict, solDict, phiSDict, pE, R, T))
    phiS = args[3]
    # print args[0]
    # print phiS
    x0 = [0.5 * x for x in phiS]
    x0.append(0.5)

    bounds = [(1E-6, 1.0-1E-6) for spc in species]
    # bounds = [(1E-6, 1.0-1E-6) for spc in species[0:-1]]
    # bounds.append((0.3, 1-1E-6))
    # print bounds
    # print x0
    # print args,

    res = opt.minimize(getEquilibriumRoots, x0, args=args, method = "SLSQP", bounds = bounds, options={'disp': False})
    print args[0]
    print res['x']

    # species1 = 'toluene'
    # species2 = 'isoOctane'
    # print getPureChi(species1, species2, VDict, solDict, R, T)
    #
    # species1 = 'isoOctane'
    # species2 = 'bunaN'
    # print getPureChi(species1, species2, VDict, solDict, R, T)