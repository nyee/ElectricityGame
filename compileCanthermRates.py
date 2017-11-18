import os
import re
import csv
from rmgpy.species import Species

d1="/Users/Nate/Dropbox (MIT)/Research/butanolSurfaces/alpha/canthermResults/f12-tz"
d2="/Users/Nate/Dropbox (MIT)/Research/butanolSurfaces/iC4H8O_beta/canthermResults/ibutanol_beta_f12tz"
d3="/Users/Nate/Dropbox (MIT)/Research/butanolSurfaces/gamma_iBuOH/canthermresults/f12-vtz"

dic1='/Users/Nate/Dropbox (MIT)/Research/butanolSurfaces/alpha/alpha_Jobs.csv'
dic2='/Users/Nate/Dropbox (MIT)/Research/butanolSurfaces/iC4H8O_beta/beta_Jobs.csv'
dic3='/Users/Nate/Dropbox (MIT)/Research/butanolSurfaces/gamma_iBuOH/gamma_Job.csv'

folderList=[d1,d2,d3]
dicList=[dic1, dic2, dic3]
initialList=['A','B','G']

#walk through folder looking for chem.inp

#start with small molecules
smilesDict={'HO2': 'O[O]',
            'OH': '[OH]',
            'H_bmk': '[H]',
            'CH3_bmk': '[CH3]',
            'H2O': 'O',
            'CH2O': 'CO'}

reactionList=[]
reactionPieces=[]
speciesList=[]
speciesDict={}
thermoDict={}
isomorphicDict={}
smilesIndex=0
outputDirectory="/Users/Nate/Dropbox (MIT)/Research/butanolSurfaces/f12_library2"

#Read in dictionary csv
for dic, initial in zip(dicList, initialList):
    fileList=[]
    with open(dic, 'rb') as inputFile:
        spamreader=csv.reader(inputFile)
        for line in spamreader:
            fileList.append(line)
    #get index of SMILES

    for columnIndex, column in enumerate(fileList[0]):
        if re.search("SMILES", column):
            smilesIndex=columnIndex
            break
    for line in fileList:
        if re.search('p[0-9]*', line[0]) and line[smilesIndex]:
            #fix smiles:
            fixedSmiles=re.sub("_bmk", '', line[smilesIndex])
            newName=re.sub('p', initial+'p', line[0])
            # print newName, fixedSmiles
            smilesDict[newName]=fixedSmiles

#copy reactions from chem.inp
for topdir, initial in zip(folderList, initialList):
    for dirName, subdirList, fileList in os.walk(os.path.join(topdir, "reactions"), topdown=False):
        if "chem.inp" in fileList:
            with open(os.path.join(dirName, "chem.inp"), 'rb') as inputFile:
                for line in inputFile:
                    if "=" in line:
                        newLine=re.sub('p', initial+'p', line.strip())
                        reactionList.append(newLine)

#break down reaction into pieces
for reaction in reactionList:
    piece=re.split("\s\s+", reaction)
    reactionPieces.append(piece)

#parse out thermo from output.py
# for topdir, initial in zip(folderList, initialList):
#     for dirName, subdirList, fileList in os.walk(os.path.join(topdir, "species"), topdown=False):
#         if "output.py" in fileList:
#             with open(os.path.join(dirName, "output.py"), 'rb') as inputFile:
#                 for line in inputFile:
#                     if re.search("\# Thermodynamics for p", line):
#                         newName=re.sub("\# Thermodynamics for p", initial+"p", line)
#                         newName=re.sub('\:', '', newName)
#                         newName=newName.strip()
#                         thermoDict[newName]=""
#
# for thermoName in thermoDict:
#     if not thermoName in smilesDict:
#         print "cant find thermoName {0}".format(thermoName)

#parse out species from reactions
for reaction in reactionList:
    split1=re.split('\s', reaction)
    for piece in split1:
        if re.match('[A-Za-z]', piece):
            if not piece in speciesList:
                speciesList.append(piece)

# print speciesList
# print len(reactionList)

for speciesName1 in speciesList:
    if speciesName1 in smilesDict:
        #Check if duplicate
        speciesObject1=Species().fromSMILES(smilesDict[speciesName1])
        for speciesName2, speciesObject2 in speciesDict.iteritems():
            if speciesObject1.isIsomorphic(speciesObject2):
                print "found duplicate {0} and {1}".format(speciesName1, speciesName2)
                isomorphicDict[speciesName1]=speciesName2
                break
        else:
            speciesDict[speciesName1]=speciesObject1
    else:
        print "name not found yet {0}".format(speciesName1)

#Rewrite reaction pieces with duplicates
correctedReactionPieces=[]

for reaction in reactionPieces:
    correctedReaction=[]
    for piece in reaction:
        correctedPiece=piece
        for duplicate, correctName in isomorphicDict.iteritems():
            if duplicate in piece:
                correctedPiece=re.sub(r'\b{0}\b'.format(duplicate), correctName, piece)
        correctedReaction.append(correctedPiece)
    correctedReactionPieces.append(correctedReaction)

reactionPieces=correctedReactionPieces


#write out dictionary:
with open(os.path.join(outputDirectory, 'dictionary.txt'), 'wb') as outFile:
    for name, species in speciesDict.iteritems():
        outFile.write(name+'\n')
        outFile.write(species.molecule[0].toAdjacencyList())
        outFile.write("\n")

#write out thermoLibrary:


#write out reactions:
with open(os.path.join(outputDirectory, 'reactions.py'), 'wb') as outFile:
    outFile.write("#!/usr/bin/env python\n")
    outFile.write("# encoding: utf-8\n")
    outFile.write("\n")
    outFile.write("name = 'iButanol_peroxy'\n")
    outFile.write("shortDesc = u'Peroxy chemistry for isobutanol'\n")
    outFile.write("longDesc = ''\n")
    outFile.write("\n")

    for index, reaction in enumerate(reactionPieces):
        outFile.write("entry(\n")
        outFile.write("    index = "+str(index+1)+",\n")
        outFile.write("    label = "+"'"+reaction[0]+"',\n")
        outFile.write("    degeneracy = 1,\n")
        outFile.write("    kinetics = Arrhenius(\n")
        #test if unimolecular or bimolecular:
        uniMolecularTest=re.split('\<\=\>', reaction[0])
        if re.search('\+', uniMolecularTest[0]):
            outFile.write("        A = ("+str(reaction[1]+", 'cm^3/(mol*s)'),\n"))
        else:
            outFile.write("        A = "+"("+str(reaction[1]+", 's^-1'),\n"))
        outFile.write("        n = "+str(reaction[2]+",\n"))
        outFile.write("        Ea = ("+str(reaction[3]+", 'kcal/mol'),\n"))
        outFile.write("        T0 = (1, 'K'),\n")
        outFile.write("    ),\n")
        outFile.write(")\n")
        outFile.write("\n")
