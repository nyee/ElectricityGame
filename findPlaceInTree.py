from rmgpy.data.thermo import *
from rmgpy.data.rmg import RMGDatabase
from rmgpy import settings
import re
import copy
import os.path
from rmgpy.thermo import *

def getEntryDepth(entry, topGroupName):
    if entry.label==topGroupName:
        return 0
    else:
        depth=1
        parent=entry.parent
        while not parent.label==topGroupName:
            parent=parent.parent
            depth+=1
        return depth


if __name__ == "__main__":
    path="/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/oxy_species2.py"
    newGroups = ThermoGroups()
    newGroups.local_context['ThermoData']=ThermoData
    newGroups.load(path)
    outputPath="/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/parents.txt"
    groupName='group'
    topGroupName="R"

    database = RMGDatabase()
    database.load(settings['database.directory'], thermoLibraries = [], kineticsFamilies='none', kineticsDepositories='none', reactionLibraries=[])

    thermoDatabase = database.thermo

    #list of all parents
    results={}
    for name1, entry1 in newGroups.entries.iteritems():
        results[name1]=[]
        for name2, entry2 in thermoDatabase.groups[groupName].entries.iteritems():
            if newGroups.matchNodeToChild(entry2, entry1):
                results[name1].append(name2)

    #keeps just the deepest result, I think
    newResults={}

    for name1, parents in results.iteritems():
        newResults[name1]=[]
        depths=[]
        for ancestor in parents:
            if ancestor==name1: continue
            depths.append(getEntryDepth(thermoDatabase.groups[groupName].entries[ancestor], topGroupName))
        print name1, depths
        maximum=max(depths)
        for index, depth in enumerate(depths):
            if depth==maximum: newResults[name1].append(parents[index])

    #This code is to help me decide if tree is structured correctly after adding things in
    duplicateCheck={}
    for name1, parents in results.iteritems():
        depths=[]
        for ancestor in parents:
            depths.append(getEntryDepth(thermoDatabase.groups[groupName].entries[ancestor], topGroupName))
        if len(depths)!=len(set(depths)):
            duplicateCheck[name1]=[parents,depths]
    #Print out results to check for proper tree structuring:
    with open(outputPath, 'wb') as outputFile:
        for name1, parents in duplicateCheck.iteritems():
            print name1, parents
            outputFile.write(name1+": "+ str(parents)+"\n")
            outputFile.write("\n")

    #modify groups
    groupsFile=[]
    with open(os.path.join("/Users/Nate/code/RMG-database/input/thermo/groups", groupName+".py"), 'rb') as inputFile:
        for line in inputFile:
            groupsFile.append(line)

    #get treeIndex:
    treeIndex=None
    for index, line in enumerate(groupsFile):
        if re.match("tree", line):
            treeIndex=index
            break
    #Add new group names into the tree

    firstTime=True
    for name1, parents in newResults.iteritems():
        groupFound=False
        #initialize file
        if firstTime:
            firstTime=False
        else:
            groupsFile=copy.copy(newGroupsFile)
        newGroupsFile=[]
        treePast=False
        print name1, parents, len(groupsFile)
        parentDepth=getEntryDepth(thermoDatabase.groups['group'].entries[parents[0]], topGroupName)
        #search through file for place in tree
        for index, line in enumerate(groupsFile):
            newGroupsFile.append(line)
            if treePast and not groupFound:
                if re.sub("L[0-9]\:", "", line).strip()==parents[0]:
                    print "found!", index
                    newGroupsFile.append(" "*4*(parentDepth+1)+"L" +str(parentDepth+2)+": "+name1 +"\n")
                    groupFound=True
            elif index>treeIndex:
                treePast=True

    #new new entries imported from a file
    newGroupText=[]
    with open(path, 'rb') as inputFile:
        for line in inputFile:
            newGroupText.append(line)
    #add in the new groups
    newGroupsFile=newGroupsFile[0:treeIndex]+newGroupText+newGroupsFile[treeIndex:]

    with open(os.path.join("/Users/Nate/code/RMG-database/input/thermo/groups", groupName+".py"), 'wb') as outFile:
        for line in newGroupsFile:
            outFile.write(line)

