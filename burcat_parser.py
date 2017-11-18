import re
import csv

path="/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/burcat_therm.dat"
fileList=[]
speciesPath="/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/burcat_species.csv"
with open(path, 'rb') as inputfile:
    for line in inputfile:
        fileList.append(line)

#list out all the species

canStartWith=["Q",  #code filler
           "t-",
           "s",
           "p",
           "o",
           "m",
           "n",
           "Tetryl",
           "T-",
           "P-",
           "O",
           "N-",
           "H",  #requires further sorting
            "C",  #requires further sorting
           "BENZO"
              ]
# blackList=[]
# blackList=["beta",
#            "cy",
#            "Zr", #zirconium
#            "Zn", #zinc
#            "ZN", #not sure
#            "Xe", #xenon
#            "W", #tungston
#            "T", #Technetium, telenium, tritium
#            "Sn", #tin
#            "Si", #silicon
#            "SI", #disilane
#            "Sb", #antimony
#            "S", #sulfur
#            "R", #radon, radium
#            "Pt", #platinium
#            "Po", #polonium
#            "Pd", #paladium
#            "Pb", #lead
#            "P", #phosphorous and other stuff manually curated
#            "Ni", #nickel
#            "Na", #sodium
#            "Ne", #neon
#            "NT", #nitrogen with tritrium
#            "N", #bad stuffs
#            "M", #molybdenum, magnesium, manganese
#            "K", #kryton, potassium
#            "J", #jet fuels
#            "I", #iodine, iridium
#            "Hg", #mercury
#            "Ge", #germanium
#            "F", #iron, fluorine
#            "E", #electron gas
#            "D", #deuterium
#            "Cu", #copper
#            "Cr", #chromium
#            "Cl", #chlorine
#            "Ca", #calcium
#            "B", #boron, berylium
#            "A", #argon, silver, aliuminum, air
# ]

filterList=["Br",
            "Cl",
            "Cu",
            "Cr",
            "Hg",
            "\+",
            "\(L",
            "F",
            "I",
            "CHLORO",
            "N[0-9]",
            "Ca",
            "\(l",
            "N",
            "CL",
            "BR",
            "D",
            "P",
            "AL",
            'B',
            "S",
            "\cr",
            "T",
        ]

superWhiteList=["BENZOTRIFUROXAN",
                "H2CO-Formald",
                "N-C10H22",
                "N-UNDECANE",
                "N-DODECANE",
                "T-C3H5",
                "T-C4H10O",
                "Tetryl",
                ]

examineList=[]
examined=[]

speciesList=[]

#print lines for burcat
for line in fileList[20:]:
    if re.match("[A-Za-z]", line):
        for starter in canStartWith:
            if re.match(starter, line):
                species=re.split('\s', line)[0]
                for filter in filterList+examineList:
                    if re.search(filter, species):
                        if filter in examineList: examined.append(species)
                        break
                else: speciesList.append(line)



speciesList.sort()

# for species in examined:
with open(speciesPath, 'wb') as outputFile:
    spamwriter=csv.writer(outputFile)

    for line in speciesList:
        species=re.split('\s', line)[0:4]
        spamwriter.writerow(species)

print "burcat has this many species left", len(speciesList)