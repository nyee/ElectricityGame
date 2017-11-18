# from rmgpy.data.rmg import RMGDatabase
# from rmgpy.chemkin import loadSpeciesDictionary, loadChemkinFile, writeKineticsEntry
# import numpy as np
# from rmgpy.kinetics.arrhenius import Arrhenius
# from rmgpy.kinetics import ArrheniusEP, MultiArrhenius
# from rmgpy import settings
import re



if __name__ == "__main__":

    # supportedTypes = [Arrhenius, ArrheniusEP]
    chemkinPath = "/Users/Nate/Dropbox (MIT)/Research/Peng/cresols/merge_v2/merge_v201.inp"
    # dictPath = "/Users/Nate/Dropbox (MIT)/Research/Peng/cresols/RMG models/v2_gasoline/v2_gasoline_sd.txt"
    outputPath = "/Users/Nate/Dropbox (MIT)/Research/Peng/cresols/merge_v2/merge_v2011.inp"
    #
    # loadChemkinFile(chemkinPath, dictPath)
    # speciesList, reactionList = loadChemkinFile(chemkinPath, dictionaryPath=dictPath, transportPath=None, readComments = True)


    linesList = []
    with open(chemkinPath, 'rb') as inputFile:
        for line in inputFile:
            linesList.append(line)
    #
    modLines = []
    for line in linesList:
        if re.search(re.escape('2.510e-11 6.770     -8.600'), line ):
            newLine = re.sub(re.escape('2.510e-11 6.770     -8.600'), '1.700e+13 0.000     20.460', line)
            modLines.append(newLine)
        else:
            modLines.append(line)


    # for line in modLines:
    #     newLine = re.sub(re.escape('2.510e-11 6.770     -8.600'), '1.700e+13 0.000     20.460', line)
    #     newModLines.append(newLine)
    #     print newLine


    #
    # Tlist = np.array([298, 600,700,800,900,1000,1100,1200,1300,1400,1500,2000]) #K
    #
    # irreversible = []
    # for reaction in reactionList:
    #     Kc = reaction.getEquilibriumConstants(Tlist, type='Kc')
    #     Kclarge = [x > 1e8 for x in Kc]
    #     if all(Kclarge) and type(reaction.kinetics) in supportedTypes:
    #         irreversible.append(reaction)
    #
    # print len(irreversible)
    # irChemkin = []
    # irChemkin = [writeKineticsEntry(reaction, speciesList, verbose = False) for reaction in irreversible]
    # # for reaction in irreversible:
    # #     irChemkin.append(writeKineticsEntry(reaction, speciesList, verbose = False))
    #
    # print len(irChemkin)
    #
    # originalLines = []
    # with open(chemkinPath, 'rb') as originalChemkin:
    #     for line in originalChemkin:
    #         originalLines.append(line)
    #
    # newLines = []
    # found = 0
    # for line in originalLines:
    #     if line in irChemkin:
    #         try:
    #             newChunk = []
    #             found+=1
    #             newLine = re.sub("\=", "<=>", line)
    #             newChunk.append(newLine)
    #             reactionIndex = irChemkin.index(line)
    #             #This chunk checks for multiArrhenius (as well as try block), which may not be formatted the same
    #             ######################################################################
    #             reverseKinetics = reactionList[reactionIndex].generateReverseRateCoefficient()
    #
    #             # print reverseKinetics.A
    #             # print reverseKinetics.Ea
    #
    #             string = '{0:<9.3e} {1:<9.3f} {2:<9.3f}'.format(
    #                 reverseKinetics.A.value_si/ (reverseKinetics.T0.value_si ** reverseKinetics.n.value_si) * 1.0e6 ** (len(reactionList[reactionIndex].products) - 1),
    #                 reverseKinetics.n.value_si,
    #                 reverseKinetics.Ea.value_si / 4184.
    #             )
    #             string = "REV/ {0} /".format(string)
    #             ######################################################################
    #             string = "REV/ 0 0 0 /"
    #
    #             newChunk.append(string)
    #             newChunk.append("\n")
    #             newLines.extend(newChunk)
    #         except AttributeError:
    #             found+=-1
    #             newLines.append(line)
    #     else: newLines.append(line)
    #
    # print found

    with open(outputPath, 'wb') as outFile:
        outFile.writelines(modLines)
