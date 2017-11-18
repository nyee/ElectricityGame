import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize as opt

#global variables
colors = ['-k', '-r', '-g', '-b', '--r', '--g', '--b', '-k', '-k']

savePath = "/Users/Nate/Dropbox (MIT)/Research/Peng/biooils/Figures/FlorryRehner.png"
# savePath1 = "/Users/Nate/Dropbox (MIT)/Research/Peng/biooils/Figures/Swelling.png"

# def getSolutionsForPhi2(n, V1):
#     chiList = np.linspace(0.1, 5, (5-0.1)/0.1)
#     phi2List = []
#     for chi in chiList:
#         func = lambda phi2: math.log(1-phi2, math.e) + phi2 + chi*phi2**2 +V1*n*(phi2**(1/3) - phi2/2)
#         phi2_solution = brenth(func, 1E-9, 1.0-1E-9)
#         phi2List.append(phi2_solution)
#     return phi2List
#
# def getChiFromPhi2(n,V1, phi2):
#     func = lambda chi: math.log(1-phi2, math.e) + phi2 + chi*phi2**2 +V1*n*(phi2**(1/3) - phi2/2)
#     guess = 1
#     return fsolve(func, guess)
#
# def getPhi2FromChi(n, V1, chi):
#     func = lambda phi2: math.log(1-phi2, math.e) + phi2 + chi*phi2**2 +V1*n*(phi2**(1/3) - phi2/2)
#     return brenth(func, 1E-9, 1.0-1E-9)


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
    return pE*Vi*(phiE**(1/3)-phiE/2)

def getMixingChemPotI(x, species_i, allSpecies, VDict, chiDict):
    """
    Returns the mixing chemical potential for species i, normalized by RT

    :param x: list of floats, guesses for volume fractions of species in the phase
    :param species_i: str defining the species which we are calculating the chemical potential for
    :param allSpecies: list of all species, defines the indexing of the species, and must match x's indexing
    :param VDict: Dictionary of molar volumes
    :param chiDict: Dictionary of chi values. keys: tuple in indexed order of two species, value: chi
    """

    #determine species i index
    i = allSpecies.index(species_i)

    #species_i non-cross terms
    term1 = math.log(x[i]) + (1- x[i])

    #cross terms of i not equal j
    term2 = 0
    for j, species_j in enumerate(allSpecies[0:-1]): #Don't need to do it for elatomer as Vi/Velatomer <<<<1
        if i==j: continue
        term2 += -VDict[species_i]/VDict[species_j] * x[j] #Make sure each term is negative!

    #cross terms i<j
    term3 = 0
    for j, species_j in enumerate(allSpecies):
        if i < j:
            term3 += chiDict[(species_i, species_j)] * x[j] * (1 - x[i])

    #cross terms h < i < j
    term4 = 0
    for h, species_h in enumerate(allSpecies[0:i]):
        for j, species_j in enumerate(allSpecies):
            if i < j:
                term4 += -chiDict[(species_h, species_j)] * VDict[species_i]/VDict[species_h] * x[h] * x[j] #Make sure each term is negative!

    #cross terms h < i
    term5 = 0
    for h, species_h in enumerate(allSpecies[0:i]):
        term5 += chiDict[species_h, species_i] * VDict[species_i]/VDict[species_h] * x[h] * (1 - x[i])

    return term1 + term2 + term3 + term4 + term5

def constrainedFunction(x, f, lower, upper, args, minIncr=0.001):
    x = np.asarray(x)
    lower = np.asarray(lower)
    upper = np.asarray(upper)
    xBorder = np.where(x<lower, lower, x)
    xBorder = np.where(x>upper, upper, xBorder)
    fBorder = np.asarray(f(xBorder, args))

    distFromBorder =[]

    for i, guess in enumerate(fBorder):
        if guess <= lower[i]:
            distFromBorder.append(guess - lower[i] - minIncr)
        elif guess >= upper[i]:
            distFromBorder.append(guess - upper[i] + minIncr)
        else:
            distFromBorder.append(0)

     # distFromBorder = [(np.sum(np.where(fBorder<lower, lower-fBorder, 0.))+np.sum(np.where(f>upper, x-upper, 0.)))]
     # distFromBorder = [0,1,1]
    show = fBorder - distFromBorder
    print fBorder, distFromBorder, show
     # show = fBorder + (fBorder+np.where(fBorder>0, minIncr, -minIncr))*distFromBorder
    return show

def getEquilibriumRoots(phiE, args):
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
        liquidChemPotList.append(getMixingChemPotI(phiS, species_i, species[0:-1], VDict, chiDict))

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

    return roots

    #varaibles to solve for
    # phi1m = phi[0]
    # phi2m = phi[1]
    # phi3m = phi[2]
    #
    # mu1s = math.log(phi1s) + 1 - phi1s - V1/V2*phi2s + chi12*(1-phi1s)**2
    # mu2s = math.log(1-phi1s) + phi2s - V2/V1*phi1s + chi12*V2/V1*phi1s**2
    #
    # mu1m_1 = math.log(phi1m) + 1 - phi1m - V1/V2*phi2m
    # mu1m_2 = (V1/V2*chi12*phi2m + (a13 + b13*phi2m)*phi3m)*(1-phi1m) - phi1m*phi2m*phi3m*b13
    # mu1m_3 = -V1/V2*(a23 + b23*phi1m)*phi2m*phi3m + V1/V2*phi2m*phi3m*(phi2m + phi3m)*b23
    # mu1m_4 = n*V1*(phi3m**(1/3) - 0.5*phi3m)
    # mu1m = mu1m_1 + mu1m_2 + mu1m_3 + mu1m_4
    #
    # mu2m_1 = math.log(phi2m) + 1 - phi2m + V2/V1*phi1m
    # mu2m_2 = (V2/V1*chi12*phi1m + (a23 + b23*phi1m)*phi3m)*(1-phi2m) - phi1m*phi2m*phi3m*b23
    # mu2m_3 = -V2/V1*(a13 + b13*phi2m)*phi1m*phi3m + V2/V1*phi1m*phi3m*(phi1m + phi3m)*b13
    # mu2m_4 = n*V2*(phi3m**(1/3)- 0.5*phi3m)
    # mu2m = mu2m_1 + mu2m_2 + mu2m_3 +mu2m_4

def constructArgs(species, VDict, solDict, phiSDict, pE, R, T):
    """
    creates the bundle of args to send to getEquilibriumRoots

    species: list of species, the last one must be the elastomer
    phiS: list of volume fractions in solvent phase, should have one less entry than species

    """

    #Sort solvent species by molar volume:
    newSpecies = species[0:-1]
    newSpecies = sorted(newSpecies, key = lambda spc: VDict[spc] )
    newSpecies.append(species[-1])

    phiS = [phiSDict[species_i] for species_i in newSpecies[:-1]]

    #construct chiDict
    chiDict={}
    for i, species_i in enumerate(newSpecies[0:-1]):
        for species_j in newSpecies[i+1:]:
            chiDict[(species_i, species_j)]=getPureChi(species_i, species_j, VDict, solDict, R, T)

    args = [newSpecies, VDict, chiDict, phiS, pE]
    return args

if __name__ == '__main__':


    #dictionary of physical parameters
    #cm^3/mol
    VDict = {'isoOctane': 162.48,
             'toluene': 106.27,
             'anisole': 108.62,
             'phenol': 90.49,
             'furan': 72.639,
             'furfural': 82.83}
    #MPa^0.5
    solDict = {'isoOctane': 14.169,
             'toluene': 18.028,
             'anisole': 19.552,
             'phenol': 24.948,
             'furan': 18.674,
             'furfural': 23.970,
               'bunaN': 21.466}

    #test fuelC
    volFractionC = 0.465 #volume fraction toluene in C
    phiSDict = {'toluene': volFractionC,
                'isoOctane': 1.0 - volFractionC}
    pE = 6.3E-4 #mol/cm^3
    R = 8.314 #J/mok-K
    T = 298.15 #K
    #guess

    species = ['isoOctane', 'toluene', 'bunaN']

    args = constructArgs(species, VDict, solDict, phiSDict, pE, R, T)
    phiS = args[3]
    chiDict = args[2]
    species = args[0]
    x0 = [0.5 * x for x in phiS]
    x0.append(0.5)

    phim1 = variable(0,1)
    phim2 = variable(0,1)
    phim3 = variable(0,1)

    phiE = [phim1, phim2, phim3]

    constraint(math.log(phim1) ==0.5)

    # constraint(phim1 + phim2 + phim3 == 1)
    # for species_i in species[0:-1]:
    #     constraint(getMixingChemPotI(phiS, species_i, species[0:-1], VDict, chiDict) == getMixingChemPotI(phiE, species_i, species, VDict, chiDict)
    #                + getElasticChemPotI(phiE[-1], VDict[species_i], pE,))
    #
    # if solve([phim1, phim2, phim3]):
    #     print ("Solution found a=%d, b=%d and c=%d" %
    #            (phim1.value(), phim2.value(), phim3.value()))

    # solution = opt.root(constrainedFunction, x0=x0,
    #                     args=(getEquilibriumRoots, [0.0 for species_i in species], [1.0 for species_i in species], args))
    #
    # print solution

#constrainedFunction(x, f, lower, upper, args, minIncr=0.001)
#getEquilibriumRoots(phiE, args)



    # chiList = np.linspace(0.1, 5, (5-0.1)/0.1)
    # V1list = np.linspace(50, 200, 4)
    #
    # print V1list
    #check different V1:
    # for index, V1 in enumerate(V1list):
    #     phi2List = getSolutionsForPhi2(n, V1)
    #     swellList = [phi2**(-1) for phi2 in phi2List]
    #     plt.plot(chiList, swellList, colors[index], label = "V1 = {0} mL/mol".format(V1))
    # plt.show()

    # nList = [1E-5, 5E-5, 1E-4, 5E-4, 1E-3]
    # V1 = 120
    # for index, n in enumerate(nList):
    #     phi2List = getSolutionsForPhi2(n, V1)
    #     swellList = [phi2**(-1) for phi2 in phi2List]
    #     plt.plot(chiList, swellList, colors[index], label = "n1 = {0} mL/mol".format(n))
    # plt.show()

    # chiList=[]
    # swellList=[]
    # for run1 in range(300):
    #     chi=0.5*run1+0.1
    #     bnds=[(1,None)]
    #     start_pos=np.ones(1)*1e2
    #     x =  minimize(equations, start_pos, jac=jacobian, bounds=bnds, method="L-BFGS-B")
    #     if x.success:
    #         n1=x.x[0]
    #         phi1=n1/(n1+r*n2)
    #         phi2=r*n2/(n1+r*n2)
    #         print chi, n1, phi1, phi2, phi1/phi2
    #     chiList.append(chi)
    #     swellList.append(phi1/phi2)
    #
    # fig = plt.figure()
    # ax = plt.axes()
    #
    # print "start here"
    # for chi in chiList:
    #     print chi
    # for swell in swellList:
    #     print swell

    # plt.semilogy(chiList, swellList)
    # plt.xlabel("$\chi$", fontsize = 16)
    # plt.ylabel("Percent Swelling")
    #
    # fig.set_size_inches(6,5)
    # # plt.show()
    # fig.savefig(savePath1, dpi = 140)

    # ##plot volume of phenol vs chi
    # R=8.314
    # T=298
    # V=112.63
    #
    # dgas=16
    # dphenol=25.3
    #
    # dElastomerList=[15, 16.25, 17, 23.0]
    #
    # volumePercentList=[]
    # dtotalList=[]
    # for run1 in range(100):
    #     volumePercent=run1*0.001
    #     dsolvent=dgas*(1-volumePercent)+dphenol*volumePercent
    #     volumePercentList.append(volumePercent)
    #     dtotalList.append(dsolvent)
    #
    # results=[]
    # print dtotalList
    # for dElast in dElastomerList:
    #         chiList=[]
    #         for dSolvent in dtotalList:
    #             chiList.append(V*(dElast-dSolvent)**2/(R*T))
    #         results.append(chiList)
    #
    # fig = plt.figure()
    # ax = plt.axes()
    #
    # for index, dElast in enumerate(dElastomerList):
    #     plt.semilogy(volumePercentList, results[index], label="$\delta_{elast}=$"+str(dElast))
    #
    # plt.legend(loc="bottom left")
    # plt.ylabel("$\chi$")
    # plt.xlabel("Volume Percent Phenol",)
    # fig.set_size_inches(6,5)
    # # plt.show()
    #
    # fig.savefig(savePath, dpi = 140)

    # for n1 in range(100):
    #
    # plt.semilogx()

    #
    # # for index, data in enumerate(chiResults):
    # #     print data, swellingResults[index]
    # for chi, data in zip(chiList, dataList):
    #     plt.semilogx(data[0], data[1], label="chi="+str(chi))
    # plt.legend().draggable()
    # plt.ylabel("Percent swelling")
    # plt.xlabel("n1_tot")
    # plt.show()