from scipy.optimize import minimize
import math
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl

#global variables
r=float(100)
chi=1e-5
n2=float(1e5)

savePath = "/Users/Nate/Dropbox (MIT)/Research/Peng/biooils/Figures/FlorryHuggins.png"
savePath1 = "/Users/Nate/Dropbox (MIT)/Research/Peng/biooils/Figures/Swelling.png"

def equations(n1):
    phi1=n1/(n1+r*n2)
    phi2=r*n2/(n1+r*n2)

    obj=n1*math.log(phi1)+chi*n1*phi2+n2*math.log(phi2)
    return obj

def jacobian(n1):
    phi1=n1/(n1+r*n2)
    phi2=r*n2/(n1+r*n2)

    dGdn1=math.log(phi1)+1+chi-phi1-2*chi*phi1+chi*phi1**2-phi2**2/r

    return np.array([dGdn1])

results={}

chiList=[]
swellList=[]
for run1 in range(300):
    chi=0.02*run1+0.1
    bnds=[(1,None)]
    start_pos=np.ones(1)*1e2
    x =  minimize(equations, start_pos, jac=jacobian, bounds=bnds, method="L-BFGS-B")
    if x.success:
        n1=x.x[0]
        phi1=n1/(n1+r*n2)
        phi2=r*n2/(n1+r*n2)
        print chi, n1, phi1, phi2, phi1/phi2
    chiList.append(chi)
    swellList.append(phi1/phi2)

fig = plt.figure()
ax = plt.axes()

print "start here"
for chi in chiList:
    print chi
for swell in swellList:
    print swell

# plt.semilogy(chiList, swellList)
# plt.xlabel("$\chi$", fontsize = 16)
# plt.ylabel("Percent Swelling")
#
# fig.set_size_inches(6,5)
# # plt.show()
# fig.savefig(savePath1, dpi = 140)

# ##plot volume of phenol vs chi
# R=8.314
# T=298
# V=112.63
#
# dgas=16
# dphenol=25.3
#
# dElastomerList=[15, 16.25, 17, 23.0]
#
# volumePercentList=[]
# dtotalList=[]
# for run1 in range(100):
#     volumePercent=run1*0.001
#     dsolvent=dgas*(1-volumePercent)+dphenol*volumePercent
#     volumePercentList.append(volumePercent)
#     dtotalList.append(dsolvent)
#
# results=[]
# print dtotalList
# for dElast in dElastomerList:
#         chiList=[]
#         for dSolvent in dtotalList:
#             chiList.append(V*(dElast-dSolvent)**2/(R*T))
#         results.append(chiList)
#
# fig = plt.figure()
# ax = plt.axes()
#
# for index, dElast in enumerate(dElastomerList):
#     plt.semilogy(volumePercentList, results[index], label="$\delta_{elast}=$"+str(dElast))
#
# plt.legend(loc="bottom left")
# plt.ylabel("$\chi$")
# plt.xlabel("Volume Percent Phenol",)
# fig.set_size_inches(6,5)
# # plt.show()
#
# fig.savefig(savePath, dpi = 140)

# for n1 in range(100):
#
# plt.semilogx()

# chi=1e-5
# chiList=[]
# dataList=[]
# for run2 in range(10)[1:]:
#     chi=1e-5*10**(run2-4)
#     chiResults=[]
#     swellingResults=[]
#     chiList.append(chi)
#     for run in range(100)[1:]:
#         print n1tot
#         n1tot=math.exp(float(run)/5.0)
#         # chi=1e-10*math.exp(float(run)/5)
#         x0=np.array([5e4, 5e2])
#         bnds=((1e-12,n1tot), (1e-12, n2tot))
#         start_pos = np.ones(2)*1e2
#         try:
#             #for scipy.optimize.minimize, equation should be a function that returns value of objective (delta G)
#             #start_pos is initial guess
#             #jacobian should return the dG/dn1
#             #bounds are limits for each x paramter
#             x =  minimize(equations, start_pos, jac=jacobian, bounds=bnds, method="L-BFGS-B")
#             if x.success:
#                 n1= x.x[0]
#                 n2= x.x[1]
#
#                 phi2B= r*(n2tot - n2) / ((n1tot - n1) + r*(n2tot - n2))
#                 phi1B= (n1tot-n1)/ ((n1tot - n1) + r*(n2tot - n2))
#
#                 # chiResults.append(chi)
#                 chiResults.append(n1tot)
#                 swellingResults.append((r*(n2tot - n2)+(n1tot-n1))/(r*n2tot))
#         except ValueError: pass
#     dataList.append([chiResults, swellingResults])
#
# # for index, data in enumerate(chiResults):
# #     print data, swellingResults[index]
# for chi, data in zip(chiList, dataList):
#     plt.semilogx(data[0], data[1], label="chi="+str(chi))
# plt.legend().draggable()
# plt.ylabel("Percent swelling")
# plt.xlabel("n1_tot")
# plt.show()