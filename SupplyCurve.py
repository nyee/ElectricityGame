import csv
import matplotlib.pyplot as plt
import copy
import os.path

filename="/Users/Nate/Dropbox (MIT)/MIT_homework/EnergyMarkets/electricity_game/totalMarketMC.csv"
basepath="/Users/Nate/Dropbox (MIT)/MIT_homework/EnergyMarkets/electricity_game/supply"

fileList=[]
# portfolioIndex=1
# indCapacity=2
# aggrCapacity=3
# totalMC=6

with open(filename, 'rb') as csvFile:
    spamreader=csv.reader(csvFile)
    for line in spamreader:
        fileList.append(line)

xList=[[],[],[],[],[],[],[]]
yList=[[],[],[],[],[],[],[]]
legend=["Ostrom", "Modigliani", "Vickrey", "Phelps", "Tobin", "Nash", "Samuelson"]


prevAgg=0
for index, line in enumerate(fileList[3:]):
    if line [1]=='': continue
    portfolio=int(line[1])-1
    agg=int(line[3])
    totalMC=float(line[6])


    newX=range(prevAgg,agg)
    newY=[totalMC]*len(newX)

    #To draw line add more points
    nextMC=fileList[3+index][6]
    if not nextMC=='':
        nextMC=float(nextMC)
        extraY=range(totalMC, nextMC)
        extraX=[agg]*len(extraY)
        newX.extend(extraX)
        newY.extend(extraY)


    xList[portfolio].extend(newX)
    yList[portfolio].extend(newY)

    prevAgg=copy.copy(agg)

for xData, yData in zip(xList, yList):
    plt.plot(xData,yData, '.')

plt.xlabel("Capacity (MW)")
plt.ylabel("Marginal Cost (dollars)")
plt.legend(legend, loc=2)
# plt.axes([0,23000, 0, 65])

totalX=xList[0]+xList[1]+xList[2]+xList[3]+xList[4]+xList[5]+xList[6]
totalX.sort()
totalY=yList[0]+yList[1]+yList[2]+yList[3]+yList[4]+yList[5]+yList[6]
totalY.sort()

plt.savefig(os.path.join(basepath, "total.png"))
with open(os.path.join(basepath, "total.csv"), 'wb') as output:
    spamwriter=csv.writer(output)
    spamwriter.writerow(["Capacity (MW)", "Margical Cost ($)"])
    for index in range(len(totalX)):
        spamwriter.writerow([totalX[index],totalY[index]])


indXList=[[],[],[],[],[],[],[]]
indYList=[[],[],[],[],[],[],[]]
for z in range(6):
    prevAgg=0
    for line in fileList[3:]:
        if line [1]=='': continue
        portfolio=int(line[1])-1
        if not z==portfolio: continue
        agg=int(line[3])
        totalMC=float(line[6])

        newX=range(prevAgg,agg)
        newY=[totalMC]*len(newX)
        indXList[portfolio].extend(newX)
        indYList[portfolio].extend(newY)

        prevAgg=copy.copy(agg)

for z in range(6):
    plt.clf()
    plt.plot(indXList[z], indYList[z])
    plt.title=(legend[z] + " Individual Supply Curve")
    plt.xlabel("Capacity (MW)")
    plt.ylabel("Margical Cost (dollars)")
    plt.savefig(os.path.join(basepath, legend[z] + ".png"))


    with open(os.path.join(basepath, legend[z]+".csv"), 'wb') as output:
        spamwriter=csv.writer(output)
        spamwriter.writerow(["Capacity (MW)", "Margical Cost ($)"])
        for index in range(len(indXList[z])):
            spamwriter.writerow([indXList[z][index],indYList[z][index]])