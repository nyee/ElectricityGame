"""
This script goes through a library and rewrites the Delocalized carbons to include CS
"""

import re
path='/Users/Nate/code/RMG-database/input/kinetics/families/R_Addition_MultipleBond/groups.py'
outPath=path

fileList=[]
with open(path, 'rb') as inputFile:
    for line in inputFile:
        fileList.append(line)

newFileList=[]
for line in fileList:
    matchObject=re.search('\[.*CO', line)
    if matchObject:
        if re.search('CS', line):
            newFileList.append(line)
        else:
            newLine=re.sub('CO', 'CO,CS', line)
            newFileList.append(newLine)
    else: newFileList.append(line)

with open(outPath, 'wb') as outFile:
    for line in newFileList:
        outFile.write(line)


