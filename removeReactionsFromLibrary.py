"""
purpose:
remove all reactions from a library that have species in species dictionary
"""

from rmgpy.data.rmg import RMGDatabase
from rmgpy import settings
from os import mkdir
from shutil import copyfile
import os.path
from rmgpy.chemkin import loadSpeciesDictionary

dictionaryPath = "/Users/Nate/code/RMG-database/input/kinetics/libraries/isobutanol/dictionary.txt"
libraryName = "publishedIBuOH"
# libraryName = "test"
newName = "publishedIBuOH_peroxyRemoved"
# libraryNameToRemove = "isobutanol"

speciesDictionary = loadSpeciesDictionary(dictionaryPath)
FullDatabase=RMGDatabase()
FullDatabase.load(settings['database.directory'],
                  kineticsFamilies='all',
                  kineticsDepositories='all',
                  thermoLibraries=['primaryThermoLibrary'],   # Use just the primary thermo library, which contains necessary small molecular thermo
                  # reactionLibraries=[libraryName, libraryNameToRemove],
                  reactionLibraries=[libraryName],
                  )

library = FullDatabase.kinetics.libraries[libraryName]
# libraryToRemove = FullDatabase.kinetics.libraries[libraryNameToRemove]

reactionIndexToRemove = []

for entryIndex, entry in library.entries.iteritems():
    reaction = entry.item
    for species1 in reaction.reactants + reaction.products:
        for label, species2, in speciesDictionary.iteritems():
            if species1.isIsomorphic(species2): break
        else: break
    else: reactionIndexToRemove.append(entryIndex)

# for entryIndex1, entry1 in libraryToRemove.entries.iteritems():
#     reaction1 = entry1.item
#     for entryIndex2, entry2 in library.entries.iteritems():
#         if reaction1.isIsomorphic(entry2.item):
#             reactionIndexToRemove.append(entryIndex2)
#             break

originalNumber = len(library.entries)

for index in reactionIndexToRemove:
    print "Removing:", library.entries[index]
    print "    kinetics:", library.entries[index].data
    del library.entries[index]

# print "Potential reactions to remove:", len(libraryToRemove.entries)
print "Original number of entries:", originalNumber
print "New number of entries:", len(library.entries)
oldDirectory = os.path.join(settings['database.directory'], '/Users/Nate/code/RMG-database/input/kinetics/libraries/'+ libraryName)
newDirectory = os.path.join(settings['database.directory'], '/Users/Nate/code/RMG-database/input/kinetics/libraries/'+ newName)
# if not os.path.exists(newDirectory):
#     mkdir(newDirectory)
# copyfile(os.path.join(oldDirectory, "dictionary.txt"), os.path.join(newDirectory, "dictionary.txt"))
# library.save(os.path.join(newDirectory, "reactions.py"))
