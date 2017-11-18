import math
import scipy
import numpy
from scipy import optimize as opt

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

def getsolubility(x, solvent, V1Dict, v2Dict, R, T):

    n = x[0]
    d2 = x[1]

    v2 = v2Dict[solvent]
    V1 = V1Dict[solvent]

    # print "V1*n", V1*n
    # print "term1", (v2**(1.0/3.0)-v2/2)
    # print "term2", math.log(1.0-v2)
    # print "term3", v2
    # print "sum", (V1*n*(v2**(1.0/3.0)-v2/2.0) + math.log(1.0-v2) + v2)

    chi = -(V1*n*(v2**(1.0/3.0)-v2/2.0) + math.log(1.0-v2) + v2)/v2**2
    # print "chi", chi
    diffSquared = (R*T/V1 * chi)
    # print "diff squared", diffSquared

    if solvent == 'furfural':
        d1 = diffSquared**(0.5) +d2
    else:
        d1 = d2 - (diffSquared**(0.5))

    return d1


def getEquilibriumRoots(x, solvents, V1Dict, v2Dict, R, T):

    resids = []
    for solvent in solvents:
        resids.append(solDict[solvent] - getsolubility(x, solvent, V1Dict, v2Dict, R, T))

    return sum([resid**2 for resid in resids])



if __name__ == '__main__':

    v2Dict = {'toluene': 0.426,
              'anisole': 0.341,
              'furan': 0.448,
              'furfural': 0.421}

    V1Dict = {'isoOctane': 162.48,
             'toluene': 106.27,
             'anisole': 108.62,
             'phenol': 90.49,
             'furan': 72.64,
             'furfural': 82.83}

    solDict = {'isoOctane': 14.169,
             'toluene': 18.028,
             'anisole': 19.552,
             'phenol': 24.948,
             'furan': 18.674,
             'furfural': 23.970,
               'bunaN': 21.2}

    R = 8.314
    T = 298.15

    solvents = ["toluene", 'anisole', 'furan', 'furfural']
    args = tuple([solvents, V1Dict, v2Dict, R, T])
    # x = [  1.17079564e-03,   2.12445711e+01]
    # for solvent in solvents:
    #     print solvent, getsolubility(x, solvent, V1Dict, v2Dict, R, T), solDict[solvent] - getsolubility(x, solvent, V1Dict, v2Dict, R, T)


    x0 = [1E-4, 21.2]
    bounds = [(1E-6, 1), (18,28) ]

    res = opt.minimize(getEquilibriumRoots, x0, args=args, method = "SLSQP", bounds = bounds, options={'disp': False})
    print res