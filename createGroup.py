import csv
import re

def createGroup(line, index, rad = False):
    #get adjList
    #Radical HBI corrections to complicated to write script for, just use normal name and write manually
    if rad: adjGroup = ''
    else: adjGroup=createSaturatedAdjGroup(line[0])

    #getLabel
    if rad: label = line[0]
    else: label=createSatLabel(line[0])
    label="'"+label+"'"

    #parse numbers
    if re.search('\+\/\-', line[1]):
        [H298, Herror] = re.split('\+\/\-', line[1], maxsplit=1)
        H298 = str(round(float(H298),1))
        Herror=",'+|-'," +Herror
    else:
        H298=line[1]
        Herror=''
    H298 = str(round(float(H298), 1))
    if re.search('\+\/\-', line[2]):
        [S298, Serror] = re.split('\+\/\-', line[2], maxsplit=1)
        Serror=",'+|-'," +Serror
    else:
        S298=line[2]
        Serror=''


    cpData=""
    cpError=""
    for cp in line[3:]:
        if re.search('\+\/\-', cp):
            [newCp, newCpError] = re.split('\+\/\-', cp, maxsplit=1)
            cpData=cpData+newCp+","
            cpError=cpError+newCpError+","
        else:
            cpData=cpData+cp+","

    #remove last comma
    cpData=cpData[0:-1]

    #more manipulation of cpError based on fact that confidence interval decreases at increasing temperatures
    if not cpError=="":
        split=re.split(',', cpError[0:-1])
        splitLength=len(split)
        additions=7-len(split)

        for x in range(additions):
            cpError=cpError+split[-1]+","

        #remove last comma
        cpError=cpError[0:-1]
        cpError=",'+|-',[" +cpError+']'


    entryLines=["entry(\n"]
    entryLines.append("    index = "+str(index)+",\n")
    entryLines.append("    label = "+label+",\n")
    entryLines.append("    group = \n")
    entryLines.append("\"\"\"\n")
    entryLines.extend(adjGroup)
    entryLines.append("\"\"\",\n")

    entryLines.append("    thermo = ThermoData(\n")
    entryLines.append("        Tdata = ([300,400,500,600,800,1000,1500],'K'),\n")
    entryLines.append("        Cpdata = (["+cpData+"],'J/(mol*K)'"+cpError+"),\n")
    entryLines.append("        H298 = ("+H298+",'kJ/mol'"+Herror+"),\n")
    entryLines.append("        S298 = ("+S298+",'J/(mol*K)'"+Serror+"),\n")
    entryLines.append("    ),\n")
    entryLines.append("    shortDesc = u\"\Derived from CBS-QB3 calculation with 1DHR treatment\",\n")
    entryLines.append("    longDesc = \n")
    entryLines.append("u\"\"\"\n")
    entryLines.append("Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors\n")
    entryLines.append("optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,\n")
    entryLines.append("DOI: 10.1002/chem.201301381 \n")
    entryLines.append("\"\"\",\n")
    entryLines.append(")\n")
    entryLines.append("\n")
    return entryLines

def createSaturatedAdjGroup(label):
    #assumes that X is a C=O group
    group=[]
    filler=" * "
    bonds={}
    bonds[2]="{2,S}\n"
    bonds[3]="{2,S} {3,S}\n"
    bonds[4]="{2,S} {3,S} {4,S}\n"
    bonds[5]="{2,S} {3,S} {4,S} {5,S}\n"

    #stores index of ketene carbons
    ketene=[]
    double=[]
    for index, letter in enumerate(label):
        labelLength=len(label)
        bondStr=bonds[labelLength]
        if index>0:
            filler="   "
            bondStr= "{1,S}\n"
        if letter=="X":
            group.append(str(index+1)+filler+"CO"+" u0 "+bondStr)
        elif letter=="K":
            ketene.append(index+1)
            group.append(str(index+1)+filler+"Cd"+" u0 "+bondStr)
        elif letter=="D":
            group.append(str(index+1)+filler+"Cd"+" u0 "+bondStr)
            double.append(index+1)
        #This needs to be turned into Os or Cs
        elif letter=="C" or letter=="O":
            group.append(str(index+1)+filler+letter+"s"+" u0 "+bondStr)
        elif letter=="H":
            group.append(str(index+1)+filler+letter+"  u0 "+bondStr)

        #radicals
        elif letter=="x":
            group.append(str(index+1)+filler+"CO"+" u1 "+bondStr)
        elif letter=="k":
            ketene.append(index+1)
            group.append(str(index+1)+filler+"Cd"+" u1 "+bondStr)
        elif letter=="d":
            group.append(str(index+1)+filler+"Cd"+" u1 "+bondStr)
        elif letter=="c" or letter=="o":
            group.append(str(index+1)+filler+letter.upper()+"s"+" u1 "+bondStr)
    #The following adds the ketene group
    #catches up index to correct number
    index=index+2
    for keteneIndex in ketene:
        group[keteneIndex-1]=group[keteneIndex-1].strip()+" {"+str(index)+",D}\n"
        group.append(str(index)+filler+"Cdd"+" u0 "+"{"+str(keteneIndex)+",D} {"+str(index+1)+",D}\n")
        index+=1
        group.append(str(index)+filler+"Od"+" u0 "+"{"+str(index-1)+",D}\n")
        index+=1
    for doubleIndex in double:
        group[doubleIndex-1]=group[doubleIndex-1].strip()+" {"+str(index)+",D}\n"
        group.append(str(index)+filler+"C"+" u0 "+"{"+str(doubleIndex)+",D}\n")
        index+=1

    return group

def createSatLabel(label):
    finalLabel=""
    newLetter=""
    firstCd=False
    for index, letter in enumerate(label):
        if index==0:
            if letter=="X":
                newLetter="CO"
            elif letter=="K":
                newLetter="CCO"
            elif letter=="D":
                newLetter="Cd"
                firstCd=True
            else: newLetter=letter.upper()+"s"

            if letter=="x":
                newLetter="CJO"
            elif letter=="k":
                newLetter="CJCO"
            elif letter=="d":
                newLetter="CdJ"
            elif letter=='o' or letter=='c':
                newLetter=newLetter+"J"
            newLetter=newLetter+"-"

        else:
            if letter=="X":
                newLetter="(CO)"
            elif letter=="K":
                newLetter="(CCO)"
            elif letter=="D":
                newLetter="Cd"
            elif letter=="C" or letter=="O":
                newLetter=letter+"s"
            elif letter=="H":
                newLetter=letter

            if letter=="x":
                newLetter="(CJO)"
            elif letter=="k":
                newLetter="(CJCO)"
            elif letter=="d":
                newLetter="CdJ"
            elif letter=='o' or letter=='c':
                newLetter=letter.upper()+"sJ"

        finalLabel+=newLetter
    if firstCd:finalLabel=finalLabel[0:3]+"Cd"+finalLabel[3:]
    return finalLabel



if __name__ == "__main__":
    fileList=[]
    outputName="/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/oxygenates/NNI.py"
    # outputName="/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/oxygenates/oxy_species2.py"
    nextGroupIndex=14
    # nextGroupIndex = 3009

    with open("/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/oxygenates/NNI.csv", 'rb') as inputFile:
    # with open("/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/oxygenates/oxy_species2.csv", 'rb') as inputFile:
        spamreader=csv.reader(inputFile)
        for row in inputFile:
            newRow=row.strip()
            rowList=re.split(',', newRow)
            # print rowList
            fileList.append(rowList)

    allEntries=[]

    for index, line in enumerate(fileList):
        print line
        newEntry=createGroup(line, index+nextGroupIndex, rad = True)
        # newEntry=createGroup(line, index+nextGroupIndex, rad = False)
        allEntries.extend(newEntry)
        # print newEntry

    with open(outputName, 'wb') as outFile:
        for line in allEntries:
            print line
            outFile.write(line)

    print createSatLabel("cCCK")
    # print createSaturatedAdjGroup("CCCKH")
