from rmgpy.data.rmg import RMGDatabase
from rmgpy import settings
from rmgpy.molecule.group import Group
from rmgpy.data.thermo import ThermoDatabase, ThermoGroups
from collections import OrderedDict
import os

database = RMGDatabase()
database.load(settings['database.directory'], thermoLibraries = [], kineticsFamilies='none', kineticsDepositories='none', reactionLibraries=[])
path=settings['database.directory']
CSpath=os.path.join(path, "CSThermoGroups")

#dictionary of databasename to dictionary of CsEntries
newDatabase=ThermoDatabase()
for databaseName, thermoGroup in database.thermo.groups.iteritems():
    CSEntries={}
    for entryName, entry in thermoGroup.entries.iteritems():
        matched=False
        if isinstance(entry.item, Group):
            for groupAtom in entry.item.atoms:
                """
                We only want to remove groups with explicitly defined CS, not and "Or" atom Type with CS an option
                """
                if len(groupAtom.atomType)>1: continue
                for atom in groupAtom.atomType:
                    if atom.label=="CS": matched=True
        if matched:
            CSEntries[entryName]=entry
    if len(CSEntries)>0:
        newGroups=ThermoGroups(label=thermoGroup.label, name=thermoGroup.name, shortDesc=thermoGroup.shortDesc,
                               longDesc=thermoGroup.longDesc)
        newGroups.entries=OrderedDict(CSEntries)
        newDatabase.groups[databaseName]=newGroups
#write out CS groups to save
newDatabase.save(os.path.join(CSpath))
print "Originally there are", str(len(database.thermo.groups["group"].entries)), "groups in groups.py"
print "There are", str(len(newDatabase.groups["group"].entries)), "CS groups in groups.py"

#remove entries
for databaseName, thermoGroup in database.thermo.groups.iteritems():
    if databaseName in newDatabase.groups:
        for entryName, entry in newDatabase.groups[databaseName].entries.iteritems():
            thermoGroup.removeGroup(entry)
        databasePath=os.path.join(path,"thermo/groups/"+databaseName+".py")
        thermoGroup.save(databasePath)

print "After removing CS groups, there are", str(len(database.thermo.groups["group"].entries)), "groups in groups.py"
