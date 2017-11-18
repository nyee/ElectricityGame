from rmgpy.data.rmg import RMGDatabase
from rmgpy import settings
from rmgpy.molecule.group import Group
import os
from rmgpy.molecule.atomtype import atomTypes

database = RMGDatabase()
database.load(settings['database.directory'], thermoLibraries = [], kineticsFamilies='none', kineticsDepositories='none', reactionLibraries=[])
path=settings['database.directory']
groupName="radical"
# groupName = 'group'

thermoGroup = database.thermo.groups[groupName]

#fix thermo
for entryName, entry in thermoGroup.entries.iteritems():
    if isinstance(entry.item, Group):
        print entryName
        if entry.item.standardizeGroup(): print entryName, "was modified"
thermoGroup.save(os.path.join(path, "thermo/groups/"+groupName+'.py'))