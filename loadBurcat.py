from rmgpy.chemkin import readThermoEntry
from rmgpy.species import Species
from rmgpy.molecule import Molecule
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.data.base import Entry
import re

path="/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/input_burcat_Nick.txt"
fileList=[]
outpath = "/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/burcat/library.py"
with open(path, 'rb') as inputfile:
    for line in inputfile:
        fileList.append(line)

entryList=[]
newEntry=[]
#separate by entry
for line in fileList:
    if line.strip() =='':
        entryList.append(newEntry)
        newEntry=[]
    else:   newEntry.append(line)

letterList=['A', 'B', 'C', 'D']

#parse into species
speciesThermo={}
for entry in entryList:
    for index, line in enumerate(entry):
        if re.search("InChI", line):
            newSpecies=Species()
            newSpecies.label = re.sub('\!', '', entry[index-1].strip())
            inchiStr=line.strip()
            if re.search('/b.*', inchiStr):
                inchiStr = re.sub('/b.*', '', inchiStr)
            # elif re.search(r'/b.*?/', inchiStr):
            #     re.sub(r'/b.*?/', '/', inchiStr)
            # print inchiStr
            newSpecies.molecule= [Molecule().fromInChI(inchiStr)]
            if newSpecies.label in speciesThermo:
                if not newSpecies.label[-1] in letterList: newSpecies.label += 'A'
                elif newSpecies.label[-1] == "A": newSpecies.label +='B'
                elif newSpecies.label[-1] == 'B': newSpecies.label += 'C'
                else: raise Exception("Need more if statements")
            speciesThermo[newSpecies.label] = newSpecies
            thermoString=''
            for thermoLine in entry[index+2:]:
                thermoString += thermoLine
            (name1, thermo, formula) = readThermoEntry(thermoString)
            newSpecies.thermo = thermo


index=0
library=ThermoLibrary(label='Burcat_hydrocarbons')
for label, species in speciesThermo.iteritems():
    index+=1
    print label, species.molecule[0].InChI
    library.loadEntry(index = index,
                      label = label,
                      molecule = species.molecule[0].toAdjacencyList(),
                      thermo = species.thermo.toThermoData())

library.save(outpath)