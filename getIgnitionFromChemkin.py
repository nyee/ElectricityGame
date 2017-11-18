import csv
import re
import os.path
from os import getcwd

temperatures = range(650, 1225, 25)
print len(temperatures)
# temperatures = [650,675]
ignitionDelays = []

name = os.path.basename(os.path.normpath(getcwd()))
print name
currentDirectory = getcwd()

print currentDirectory
for y in range(len(temperatures)):
    solnNumber = y + 1
    fileName = os.path.join(currentDirectory, "CKSoln_solution_no_{0}.csv".format(solnNumber))
    csvList = []
    with open(fileName, 'rb') as csvfile:
        spamreader = csv.reader(csvfile)
        for row in spamreader:
            csvList.append(row)

    #get indexes
    for index, heading in enumerate(csvList[0]):
        if re.search('Time', heading): timeColumn = index
        elif re.search('Temperature', heading): tempColumn = index

    for index, line in enumerate(csvList[1:]):
        if float(line[tempColumn]) > 1500:
            ignitionDelays.append(line[timeColumn])
            break
    else:
        ignitionDelays.append(0)

    for index, time in enumerate(ignitionDelays):
        if time == 0:
            print "No ignition delay for solution ".format(index)

with open("{0}_ignitionDelays.csv".format(name), 'wb') as csvfile:
    spamwriter = csv.writer(csvfile)
    spamwriter.writerow(['T', 'IDT_MAX_DP/DT'])
    for x in range(len(temperatures)):
        row = [temperatures[x], ignitionDelays[x]]
        print row
        spamwriter.writerow(row)

