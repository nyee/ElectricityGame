import csv
import copy

#returns supply data with (portfolio, aggregate supply, mc), eliminating duplicates
def getSupplyCurve(filename):
    fileList=[]
    with open(filename, 'rb') as csvFile:
        spamreader=csv.reader(csvFile)
        for line in spamreader:
            fileList.append(line)

    mcList=[]

    legend=["Ostrom", "Modigliani", "Vickrey", "Phelps", "Tobin", "Nash", "Samuelson"]

    prevAgg=0
    for index, line in enumerate(fileList[3:]):
        if line [1]=='': continue
        portfolio=legend[int(line[1])-1]
        agg=int(line[3])
        totalMC=float(line[6])

        #combine if same portolio:
        if index>0:
            if totalMC==mcList[-1][2] and portfolio==mcList[-1][0]:
                # print "found similar", line
                mcList[-1]=(portfolio, agg, totalMC)
            else: mcList.append((portfolio, agg, totalMC))
        else: mcList.append((portfolio, agg, totalMC))

        prevAgg=copy.copy(agg)
    return mcList

def getSupplyNoMerge(filename):

    fileList=[]
    with open(filename, 'rb') as csvFile:
        spamreader=csv.reader(csvFile)
        for line in spamreader:
            fileList.append(line)

    mcList=[]

    legend=["Ostrom", "Modigliani", "Vickrey", "Phelps", "Tobin", "Nash", "Samuelson"]
    prevAgg=0
    for index, line in enumerate(fileList[3:]):
        if line [1]=='': continue
        portfolio=legend[int(line[1])-1]
        agg=int(line[3])
        totalMC=float(line[6])

        mcList.append((portfolio, agg, totalMC))
        prevAgg=copy.copy(agg)
    return mcList

#gets mc of the very last quantity produced
def getMc(quant, mcN):
    for index, data in enumerate(mcN):
        #less or equal okay
        if quant < data[1] or quant==data[1]:
            return mcN[index][2]
    else:
        return data[2]

def getMcFromPrice(price, mcN):
    for index, data in enumerate(mcN):
        if price < data[2]:
            return mcN[index-1][2]
        #Is different from above because of boundaries! if equal, upgrade mc
        elif price==data[2]:
            return mcN[index][2]
    else:
        return data[2]

def getMargCombCap(price, mcN):
    mc=getMcFromPrice(price, mcN)

    #calculate if another dataPoint is sharing in at that price
    margQuant=0.0
    for index,data in enumerate(mcN):
        if data[2]==mc:
            margQuant+=data[1]-mcN[index-1][1]
    return margQuant

def getMargIndvCap(price, mcN, identity=0):
    mc=getMcFromPrice(price, mcN)

    matchList=[]
    for index,data in enumerate(mcN):
        if mc==data[2]:
            matchList.append((index, data))

    #highest mc
    if identity==0:
        port=matchList[0][1][0]
    elif identity==1:
        port=matchList[-1][1][0]

    total=0.0
    for match in matchList:
        if port==match[1][0]:
            cap=match[1][1]-mcN[match[0]-1][1]
            total+=cap

    return total

def getNextHighestMC(mc, mcN):
    for index, data in enumerate(mcN):
        if mc<data[2]:
            return data[2]

#gets the amount to be produced at the marginal price
def getAvailable(quant, mcN):
    mc=getMc(quant, mcN)
    mcLower=getMcFromPrice(mc-0.01, mcN)

    for index, data in enumerate(mcN):
        if mcLower==data[2]:
            result=quant-data[1]
    return result


#identity is whether you are the first firm or the second firm in an IS1, IS2
def getPrice(quant, mcN, increment, identity=0):
    mc=getMc(quant, mcN)

    #two ways to share have same MC or increment up to somebody's MC
    #see how much we share if we have same MC
    nextHighestMC=getNextHighestMC(mc, mcN)
    price1=mc+increment
    price2=nextHighestMC-increment

    if price1>nextHighestMC:
        return price2
    elif price2<mc:
        return price1
    else:
        margIndvCap=getMargIndvCap(mc, mcN, identity)
        margCombCap=getMargCombCap(mc, mcN)

        #we have to share...
        if margCombCap>margIndvCap:
            # agg=getAggregate(mc, mcN)

            available=getAvailable(quant, mcN)

            myCap=margIndvCap
            yourCap=margCombCap-margIndvCap

            leftover=0
            #Calculate splitMine
            splitYours=min(yourCap, available/2)
            splitMine=min(myCap, available/2)
            if not splitMine+splitYours==available:
                leftover=available-splitMine-splitYours
                if splitYours<splitMine: splitMine+=leftover
                else: splitYours+=leftover

            #I go price2, you go price1:
            smallMine=max(0, available-yourCap)

            #I go price1 you go price2:
            bigMine=min(myCap, available)

            #Check options, case where split is dominant strategy:
            if price1*splitMine>price2*smallMine: return price1
            #Case with mixed equliibrium, I guess average...
            else:
                if smallMine*price2>splitMine*price1:
                    return price2
                elif yourCap==splitMine:
                    return price2
                elif bigMine*price2+price1*splitMine>smallMine*price2+splitMine*price2:
                    return price1
                else: return price2
        #we don't have to share and price

        else: return max(price1, price2)


def getIdentity(quant, mcN, identity=0):
    mc=getMc(quant, mcN)

    matchList=[]
    for index,data in enumerate(mcN):
        if mc==data[2]:
            matchList.append((index, data))

    #highest mc
    if identity==0:
        return matchList[0][1][0]
    elif identity==1:
        return matchList[-1][1][0]

if __name__ == '__main__':
    filename="/Users/Nate/Dropbox (MIT)/MIT_homework/EnergyMarkets/electricity_game/totalMarketMC.csv"
    # mcList=getSupplyCurve(filename)
    mcN=getSupplyNoMerge(filename)
    # print getAvailable(16250, mcN)
    # print getMargIndvCap(27.6, mcList)
    # print getNextHighestMC(27.5, mcList)
    # print getPrice(11051, mcList, increment=0.01, identity=0)
    # print getPrice(11349, mcList, increment=0.01, identity=1)
    # print getPrice(11750, mcList, increment=0.01, identity=1)
    #
    # firstPrice=getPrice(11050, mcList, increment=0.01, identity=1)
    # for quant in range(11050, 11750):
    #     newPrice=getPrice(quant, mcList, increment=0.01, identity=1)
    #     if newPrice>firstPrice:
    #         print quant,newPrice
    #         break


    firstResultList=[[],[]]

    for identity in [0,1]:
        currentPrice=0.0
        # print identity 10400
        for quant in range(10400, 20100):
            newPrice=getPrice(quant, mcN, 0.01, identity, )
            portfolio=getIdentity(quant, mcN, identity)

            if newPrice>currentPrice:
                print newPrice
                currentPrice=copy.copy(newPrice)
                firstResultList[identity].append((quant, portfolio, newPrice))

    # print firstResultList
    resultList=firstResultList[0]+firstResultList[1]




    with open("/Users/Nate/Dropbox (MIT)/MIT_homework/EnergyMarkets/electricity_game/price.csv", 'wb') as outFile:
        spamwriter=csv.writer(outFile)
        for line in resultList:
            spamwriter.writerow(line)
