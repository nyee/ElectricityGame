from scipy.optimize import minimize
import math
import numpy as np
import matplotlib.pyplot as plt

#global variables
r=float(100)
n1tot=float(1e6)
n2tot=float(1e5)

def equations(p):
    n1, n2 = p

    phi2A= r*n2 / (n1 + r*n2)
    phi2B= r*(n2tot - n2) / ((n1tot - n1) + r*(n2tot - n2))
    obj1=n1*math.log(1-phi2A)+n2*math.log(phi2A)+chi*n1*phi2A
    obj2=(n1tot-n1)*math.log(1-phi2B)+(n2tot-n2)*math.log(phi2B)+chi*(n1tot-n1)*phi2B

    obj=obj1+obj2

    return obj

def jacobian(p):
    n1, n2 = p

    phi2A= r*n2 / (n1 + r*n2)
    phi2B= r*(n2tot - n2) / ((n1tot - n1) + r*(n2tot - n2))

    dp2A1=-r*n2 / (n1 + r*n2)**2
    dp2B1= r*(n2tot - n2) / ((n1tot - n1) + r * (n2tot - n2))**2
    dp2A2= n1 / ((n1 + r*n2)**2)
    dp2B2=-(n1tot - n1)/ ((n1tot - n1) + r * (n2tot - n2))**2

    dGdn1I= dp2A1*(n1/(1 - phi2A) + n2/phi2A + chi*n1)
    dGdn1II= dp2B1*((n1tot - n1)/(1 - phi2B)+ (n2tot - n2)/phi2B + chi*(n1tot - n1))
    dGdn1III= chi*((n1tot - n1) + (phi2A - phi2B)) + math.log(1 - phi2A) - math.log (1-phi2B)
    dGdn1= dGdn1I+dGdn1II+dGdn1III


    dGdn2I= dp2A2*(n1/(1 - phi2A) + n2/phi2A + chi*n1)
    dGdn2II= dp2B2*((n1tot - n1)/(1 - phi2B)+ (n2tot - n2)/phi2B + chi*(n1tot - n1))
    dGdn2III= math.log(phi2A) - math.log(phi2B)
    dGdn2= dGdn2I+dGdn2II+dGdn2III

    return np.array([dGdn1, dGdn2])

results={}

chi=1e-5
chiList=[]
dataList=[]
for run2 in range(10)[1:]:
    chi=1e-5*10**(run2-4)
    chiResults=[]
    swellingResults=[]
    chiList.append(chi)
    for run in range(100)[1:]:
        print n1tot
        n1tot=math.exp(float(run)/5.0)
        # chi=1e-10*math.exp(float(run)/5)
        x0=np.array([5e4, 5e2])
        bnds=((1e-12,n1tot), (1e-12, n2tot))
        start_pos = np.ones(2)*1e2
        try:
            #for scipy.optimize.minimize, equation should be a function that returns value of objective (delta G)
            #start_pos is initial guess
            #jacobian should return the dG/dn1
            #bounds are limits for each x paramter
            x =  minimize(equations, start_pos, jac=jacobian, bounds=bnds, method="L-BFGS-B")
            if x.success:
                n1= x.x[0]
                n2= x.x[1]

                phi2B= r*(n2tot - n2) / ((n1tot - n1) + r*(n2tot - n2))
                phi1B= (n1tot-n1)/ ((n1tot - n1) + r*(n2tot - n2))

                # chiResults.append(chi)
                chiResults.append(n1tot)
                swellingResults.append((r*(n2tot - n2)+(n1tot-n1))/(r*n2tot))
        except ValueError: pass
    dataList.append([chiResults, swellingResults])

# for index, data in enumerate(chiResults):
#     print data, swellingResults[index]
for chi, data in zip(chiList, dataList):
    plt.semilogx(data[0], data[1], label="chi="+str(chi))
plt.legend().draggable()
plt.ylabel("Percent swelling")
plt.xlabel("n1_tot")
plt.show()