from rmgpy.data.thermo import *
from rmgpy.data.rmg import RMGDatabase
from rmgpy import settings
import re
import copy
import os.path
from rmgpy.data.kinetics.family import KineticsFamily
from rmgpy.thermo import *

def checkOverlappingChildren(database, ignore = []):
    #checks that overlapping are in order from most specific to least specific in terms of parent relationships:

    for name, entry in database.entries.iteritems():
        #ignore this check for to product
        if entry in ignore: continue
        numberOfChildren=len(entry.children)
        #scan from top to bottom
        for index, child1 in enumerate(entry.children):
            if index==numberOfChildren-1:
                break
            for child2 in entry.children[index+1:]:
                if database.matchNodeToChild(child1, child2):
                    return (child1, child2)
                elif database.matchNodeToChild(child2, child1):
                    return (child2, child1)
    return None

if __name__ == "__main__":

    database = RMGDatabase()
    database.load(settings['database.directory'], thermoLibraries = [], kineticsFamilies='all', kineticsDepositories=['training'], reactionLibraries=[])

    #solvation
    #statmech
    #transport

    # for databaseName, specificThermoDatabase in database.thermo.groups.iteritems():
    for databaseName, specificThermoDatabase in database.transport.groups.iteritems():
        numberChanged=0
    # for databaseName, specificThermoDatabase in test.iteritems():
        saveDatabase=False
        # saveDatabase=True
        error=checkOverlappingChildren(specificThermoDatabase)
        if not error is None:
            print databaseName
            saveDatabase=True
        while not error is None:
            numberChanged+=1
            print error
            (upperChild, lowerChild)= error
            lowerChild.parent=upperChild
            upperChild.parent.children.remove(lowerChild)
            upperChild.children.append(lowerChild)
            error=checkOverlappingChildren(specificThermoDatabase)
        if saveDatabase:
            specificThermoDatabase.save(os.path.join(settings['database.directory'], "/groups/"+ databaseName+".py"))
            print databaseName, "had {0} nodes moved".format(numberChanged)

    # database.thermo.save(os.path.join(settings['database.directory'], "thermo"))

    # for familyName, family in database.kinetics.families.iteritems():
    #     numberChanged=0
    # # for databaseName, specificThermoDatabase in test.iteritems():
    #     saveDatabase=False
    #     # saveDatabase=True
    #     if not family.ownReverse:
    #             ignore = family.forwardTemplate.products
    #     else: ignore =[]
    #     error=checkOverlappingChildren(family.groups, ignore)
    #     if not error is None:
    #         print familyName
    #         saveDatabase=True
    #     while not error is None:
    #         numberChanged+=1
    #         print error[0].label, error[1].label
    #         (upperChild, lowerChild)= error
    #         lowerChild.parent=upperChild
    #         upperChild.parent.children.remove(lowerChild)
    #         upperChild.children.append(lowerChild)
    #         error=checkOverlappingChildren(family.groups)
    #     if saveDatabase:
    #         family.save(os.path.join(settings['database.directory'], "kinetics/families/" + familyName))
    #         print familyName, "had {0} nodes moved".format(numberChanged)