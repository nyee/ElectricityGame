from rmgpy.chemkin import loadChemkinFile
from numpy import std
import os.path
import csv

model1 = "/Users/Nate/Dropbox (MIT)/Research/Peng/chem_annotated_1.inp"
dict1 = "/Users/Nate/Dropbox (MIT)/Research/Peng/species_dictionary_1.txt"
model2 = "/Users/Nate/Dropbox (MIT)/Research/Peng/2_chem.inp"
dict2 = "/Users/Nate/Dropbox (MIT)/Research/Peng/2_species_dictionary.txt"
model3 = "/Users/Nate/Dropbox (MIT)/Research/Peng/3_chem.inp"
dict3 = "/Users/Nate/Dropbox (MIT)/Research/Peng/3_species_dictionary.txt"

temperature = 700
outPath = os.path.join("/Users/Nate/Dropbox (MIT)/Research/Peng/", "G_"+str(temperature)+".csv")

speciesList1, reactionList1 = loadChemkinFile(model1, dict1, readComments = False)
speciesList2, reactionList2 = loadChemkinFile(model2, dict2, readComments = False)
speciesList3, reactionList3 = loadChemkinFile(model3, dict3, readComments = False)

allSpeciesList = [speciesList1, speciesList2, speciesList3]
speciesDict1 = {}
for species in speciesList1:
    speciesDict1[species.label] = species


resultsDict = {}
for index, speciesList in enumerate(allSpeciesList):
    for species in speciesList:
        if species.label in resultsDict:
            resultsDict[specieslabel2][0].append(species.label)
            resultsDict[specieslabel2][1].append(species.thermo.getFreeEnergy(temperature)*0.239/1000.0)
        else:
            if index > 0:
                for specieslabel2, species2 in speciesDict1.iteritems():
                    if species2.isIsomorphic(species):
                        resultsDict[specieslabel2][0].append(species.label)
                        resultsDict[specieslabel2][1].append(species.thermo.getFreeEnergy(temperature)*0.239/1000.0)
                        break
                else:
                    resultsDict[species.label] = [[species.label], [species.molecule[0].toSMILES(), species.thermo.getFreeEnergy(temperature)*0.239/1000.0]]

            else:
                resultsDict[species.label] = [species.label], [species.molecule[0].toSMILES(), species.thermo.getFreeEnergy(temperature)*0.239/1000.0]

print len(resultsDict)

sortedList=[]
for label in resultsDict:
    stddev = std(resultsDict[label][1][1:])
    newLine = resultsDict[label][0]+resultsDict[label][1]+[stddev]
    if len(sortedList) == 0:
        sortedList.append(newLine)
    else:
        for index2, savedLine in enumerate(sortedList):
            if stddev == savedLine[-1] or stddev > savedLine[-1]:
                break
        sortedList.insert(index2, newLine)


heading = ["label1", "label2", "label3", 'SMILES', 'model1 G (kcal/mol)', 'model2 G (kcal/mol)', "model3 G (kcal/mol)", "std dev"]
with open(outPath, 'wb') as csvfile:
    spamwriter = csv.writer(csvfile)
    spamwriter.writerow(heading)
    for line in sortedList:
        spamwriter.writerow(line)