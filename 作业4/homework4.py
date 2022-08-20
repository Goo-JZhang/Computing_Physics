import numpy as np
import random
import math
import pandas as pd
import time
import numba
from numba import jit
import vegas
import matplotlib.pyplot as plt

#@jit(nopython=True)
def Pc(x,y):
    if abs(x-y)<=12:
        return (x-y)**2
    else:
        return (24-abs(x-y))**2

def Pt(x,y):
    if abs(x-y)<=32:
        return (x-y)**2
    else:
        return (64-abs(x-y))**2

def distance(r1,r2):
    s=Pc(r1[0],r2[0])+Pc(r1[1],r2[1])+Pc(r1[2],r2[2])+Pt(r1[3],r2[3])
    return math.sqrt(s)

def getx(r1,r2,r3):
    dis=sorted([distance(r1,r2)/dl,distance(r1,r3)/dl,distance(r2,r3)/dl])
    return math.ceil(dis[2]),math.ceil(dis[1]),math.ceil(dis[0])

def Rpoint(R):
    Result=[]
    for i in range(-R,R+1):
        s1=R**2-i**2
        for j in range(-math.floor(pow(s1,0.5)),math.floor(pow(s1,0.5))+1):
            s2=s1-j**2
            for k in range(-math.floor(pow(s2,0.5)),math.floor(pow(s2,0.5))+1):
                s3=s2-k**2
                for t in range(-math.floor(pow(s3,0.5)),math.floor(pow(s3,0.5))+1):
                    Result.append((i,j,k,t))
    return Result

lmax=pow(3*12**2+32**2,0.5)

SN=5000
pN=1024
N=20
tN=10000
dl=lmax/N
rN=4
R=6
Rp=Rpoint(R)
#print(Rp)
Lp=len(Rp)
#print(Lp)

def cnorm(a):
    if a<0:
        return a+24
    elif a>23:
        return a-24
    else:
        return a

def tnorm(t):
    if t<0:
        return t+64
    elif t>63:
        return t-64
    else:
        return t

def randpoint(c):
    result=[]
    for i in range(rN):
        d=Rp[random.randint(0,Lp-1)]
        result.append((cnorm(c[0]+d[0]),cnorm(c[1]+d[1]),cnorm(c[2]+d[2]),tnorm(c[3]+d[3])))
    return result


def homework4_1_a():
    #创建字典存储三角形出现个数
    #start=time.time()
    print("1A start")
    D={}
    for x1 in range(0,N+1):
        for x2 in range(0,x1+1):
            for x3 in range(max(0,x1-x2-1),x2+1):
                D[str((x1,x2,x3))]=np.zeros(SN,dtype=int)
    #print("D complete",time.time()-start)
    start=time.time()
    for i in range(SN):
        #print(i,'complete')
        #dstart=time.time()
        S=list(zip(np.random.randint(low=0,high=24,size=pN),\
        np.random.randint(low=0,high=24,size=pN),np.random.randint(low=0,high=24,size=pN),\
        np.random.randint(low=0,high=64,size=pN)))
        #样本中取三角形
        t1=np.random.randint(low=0,high=pN,size=tN)
        t2=np.random.randint(low=0,high=pN,size=tN)
        t3=np.random.randint(low=0,high=pN,size=tN)
        #print("randtime:",time.time()-dstart)
        #dstart=time.time()
        for j in range(tN):
            #t1=random.randint(0,pN-1)
            #t2=random.randint(0,pN-1)
            #t3=random.randint(0,pN-1)
            #D[str(getx(S[t1],S[t2],S[t3]))][i]+=1
            D[str(getx(S[t1[j]],S[t2[j]],S[t3[j]]))][i]+=1
        #print("caltime:",time.time()-dstart)
        if not i%100:
            print('process:',100*i/SN,' time used:',time.time()-start)
    #print(time.time()-start)
    D=pd.DataFrame(D)
    D.to_excel("H4a.xlsx",index=False)
    avr=D.sum()
    var=(D**2).sum()-(avr**2)/SN
    avr=list(avr/(SN*tN))
    var=list(np.sqrt(var/(SN*(SN-1)*(tN**2))))
    P=pd.DataFrame(dict(zip(D.keys(),list(zip(avr,var)))))
    P.to_excel('Possibility-a.xlsx',index=False)
    #print(time.time()-start)
    #print(P.apply(lambda x:x.sum(),axis=1))

def homework4_1_b():
    print("1B start")
    D={}
    for x1 in range(0,N+1):
        for x2 in range(0,x1+1):
            for x3 in range(max(0,x1-x2-1),x2+1):
                D[str((x1,x2,x3))]=np.zeros(SN,dtype=int)
    start=time.time()
    for i in range(SN):
        #print(i,'complete')
        Sc=list(zip(np.random.randint(low=0,high=24,size=pN//rN),\
        np.random.randint(low=0,high=24,size=pN),np.random.randint(low=0,high=24,size=pN//rN),\
        np.random.randint(low=0,high=64,size=pN//rN)))
        S=[]
        for c in Sc:
            S+=randpoint(c)
        #print(S)
        #样本中取三角形
        t1=np.random.randint(low=0,high=pN,size=tN)
        t2=np.random.randint(low=0,high=pN,size=tN)
        t3=np.random.randint(low=0,high=pN,size=tN)
        for j in range(tN):
            #t1=random.randint(0,pN-1)
            #t2=random.randint(0,pN-1)
            #t3=random.randint(0,pN-1)
            #D[str(getx(S[t1],S[t2],S[t3]))][i]+=1
            D[str(getx(S[t1[j]],S[t2[j]],S[t3[j]]))][i]+=1
        if not i%100:
            print('process:',100*i/SN,' time used:',time.time()-start)
    #print(time.time()-start)
    D=pd.DataFrame(D)
    D.to_excel("H4b.xlsx",index=False)
    avr=D.sum()
    var=(D**2).sum()-(avr**2)/SN
    avr=list(avr/(SN*tN))
    var=list(np.sqrt(var/(SN*(SN-1)*(tN**2))))
    P=pd.DataFrame(dict(zip(D.keys(),list(zip(avr,var)))))
    P.to_excel('Possibility-b.xlsx',index=False)

def homework4_1_c():
    Pa=pd.read_excel('Possibility-a.xlsx')
    Pb=pd.read_excel('Possibility-b.xlsx')
    Dr={}
    for cl in Pa.columns.values:
        if Pb[cl][0]:
            r=Pa[cl][0]/Pb[cl][0]
            if r<0.95 or r>1.05:
                Dr[cl]=[r]
            else:
                pass
        elif Pa[cl][0]:
            Dr[cl]='inf'
        else:
            pass
    Dr=pd.DataFrame(Dr)
    Dr.T.to_excel('Ratio.xlsx')

##第二题

def expSlat(x,xlist):
    S=(x-xlist[0])**2+(x-xlist[-1])**2+0.25*(x**2)+0.25*(xlist[-1])**2
    for i in range(len(xlist)-1):
        S+=(xlist[i+1]-xlist[i])**2+0.25*(xlist[i]**2)
    return np.exp(-S)

No=8

Ao=pow(1/math.pi,No/2)

def homework4_2_a():
    Vol=[[-10,10] for i in range(No-1)]
    integ=vegas.Integrator(Vol)
    nstrat=(No-1)*[4]
    def PathInt(x):
        integ(lambda y:expSlat(x,y),nitn=5,nstrat=nstrat)
        result=integ(lambda y:expSlat(x,y),nitn=5,nstrat=nstrat)
        #print(result.summary())
        return (Ao*result.mean,Ao*result.sdev)
    xcal=[]
    xstan=[]
    xerr=[]
    for x in np.linspace(0,2,9):
        y,ye=PathInt(x)
        xcal.append(y)
        xerr.append(ye)
        ys=np.exp(-x**2-2)/pow(math.pi,0.5)
        xstan.append(ys)
        print('x={}, PathInt={},error={}, Standard={}'.format(x,y,ye,ys))
    #print(xcal)
    #print(xstan)
    x=np.linspace(0,2,9)
    plt.errorbar(x,xcal,yerr=xerr,color='red',fmt='o',label='PathInt')
    plt.scatter(x,xstan,label='Standard')
    plt.grid()
    plt.legend()
    plt.show()
    
ep=1.4
Ncor=20
Npar=20

def S(j,x):
    Nx=len(x)
    return (x[(j+1)%Nx]-x[j])**2+(x[j]-x[j-1])**2+0.25*(x[j]**2)

def Metropolis(xarray):
    for i in range(len(xarray)):
        oldS=S(i,xarray)
        oldx=xarray[i]
        xarray[i]=xarray[i]+random.uniform(-ep,ep)
        newS=S(i,xarray)-oldS
        if newS<0:
            pass
        else:
            eta=random.uniform(0,1)
            if np.exp(-newS)>eta:
                pass
            else:
                xarray[i]=oldx
    return xarray

def Jacknife(ConG,N):
    JK=pd.DataFrame()
    G=np.zeros(Npar)
    for t in range(Npar):
        for i in range(N):
            G[t]+=np.sum(ConG[i]*np.roll(ConG[i],t))
    for key in range(N):
        tempG=np.zeros(Npar)
        for t in range(Npar):
            tempG[t]+=np.sum(ConG[key]*np.roll(ConG[key],t))
        JK[key]=(G-tempG)/((N-1)*Npar)
    bar=JK.sum(axis=1)/N
    Err=np.zeros(JK.shape[0])
    for key in JK.columns:
        Err+=(JK[key]-bar)**2
    Err=np.sqrt(Err/N)
    return JK,bar,Err

def JackE(G,N):
    E=pd.DataFrame()
    for i in range(N):
        tempE=np.zeros(Npar-1)
        for j in range(Npar-1):
            tempE[j]=G.iloc[j,i]/G.iloc[j+1,i]
        E[i]=tempE
    bar=E.sum(axis=1)/N
    Err=np.zeros(E.shape[0])
    for key in E.columns:
        Err+=(E[key]-bar)**2
    Err=np.sqrt(Err/N)
    return Err



def homework4_2_bc():
    ##b
    xarray=np.zeros(Npar)
    #预热
    for i in range(10*Ncor):
        xarray=Metropolis(xarray)
    ConG=[]
    for i in range(1000):
        for j in range(Ncor):
            xarray=Metropolis(xarray)
        ConG.append(xarray.copy())
    D=pd.DataFrame(dict(zip(range(1000),ConG)))
    N1,N2,N3,N4=25,100,400,1000
    ##c
    #G=np.zeros(Npar)
    JK1,bar1,err1=Jacknife(ConG,N1)#Nconf=25
    JK2,bar2,err2=Jacknife(ConG,N2)#Nconf=100
    JK3,bar3,err3=Jacknife(ConG,N3)#Nconf=400
    JK4,bar4,err4=Jacknife(ConG,N4)#Nconf=1000
    
    #Nconf=25
    '''
    for t in range(Npar):
        for i in range(N1):
            G[t]+=np.sum(ConG[i]*np.roll(ConG[i],t))
        G1[t]=G[t]/(N1*Npar)
        for i in range(N1,N2):
            G[t]+=np.sum(ConG[i]*np.roll(ConG[i],t))
        G2[t]=G[t]/(N2*Npar)
        for i in range(N2,N3):
            G[t]+=np.sum(ConG[i]*np.roll(ConG[i],t))
        G3[t]=G[t]/(N3*Npar)
        for i in range(N3,N4):
            G[t]+=np.sum(ConG[i]*np.roll(ConG[i],t))
        G4[t]=G[t]/(N4*Npar)
    '''
    f=open("texG.txt",'w+')
    for i in range(Npar):
        f.write("{}&{:.6}&{:.6}&{:.6}&{:.6}\\\\\n".format(i,bar1[i],bar2[i],bar3[i],bar4[i]))
    f.close()
    E1=np.zeros(Npar-1)
    for i in range(Npar-1):
        E1[i]=np.log(bar1[i]/bar1[i+1])
    E2=np.zeros(Npar-1)
    for i in range(Npar-1):
        E2[i]=np.log(bar2[i]/bar2[i+1])
    E3=np.zeros(Npar-1)
    for i in range(Npar-1):
        E3[i]=np.log(bar3[i]/bar3[i+1])
    E4=np.zeros(Npar-1)
    for i in range(Npar-1):
        E4[i]=np.log(bar4[i]/bar4[i+1])
    '''
    Err1=np.zeros(Npar-1)
    for i in range(Npar-1):
        Err1[i]=np.sqrt((err1[i]/bar1[i])**2+(err1[i+1]/bar1[i+1])**2)
    Err2=np.zeros(Npar-1)
    for i in range(Npar-1):
        Err2[i]=np.sqrt((err2[i]/bar2[i])**2+(err2[i+1]/bar2[i+1])**2)
    Err3=np.zeros(Npar-1)
    for i in range(Npar-1):
        Err3[i]=np.sqrt((err3[i]/bar3[i])**2+(err3[i+1]/bar3[i+1])**2)
    Err4=np.zeros(Npar-1)
    for i in range(Npar-1):
        Err4[i]=np.sqrt((err4[i]/bar4[i])**2+(err4[i+1]/bar4[i+1])**2)
    '''
    Err1=JackE(JK4,N1)
    Err2=JackE(JK4,N2)
    Err3=JackE(JK4,N3)
    Err4=JackE(JK4,N4)
    #print(err1)
    #print(err2)
    #print(err3)
    #print(err4)
    xn=list(range(20))
    plt.figure('G1')
    plt.grid()
    plt.xlabel("n")
    plt.ylabel("Gn")
    plt.errorbar(xn,bar1,yerr=err1,color='lightcoral',ecolor='lightcoral',label='Nconf=25')
    plt.legend()
    plt.figure('G2')
    plt.grid()
    plt.xlabel("n")
    plt.ylabel("Gn")
    plt.errorbar(xn,bar2,yerr=err2,color='sienna',ecolor='sienna',label='Nconf=100')
    plt.legend()
    plt.figure('G3')
    plt.grid()
    plt.xlabel("n")
    plt.ylabel("Gn")
    plt.errorbar(xn,bar3,yerr=err3,color='orange',ecolor='orange',label='Nconf=400')
    plt.legend()
    plt.figure('G4')
    plt.grid()
    plt.xlabel("n")
    plt.ylabel("Gn")
    plt.errorbar(xn,bar4,yerr=err4,color='indigo',ecolor='indigo',label='Nconf=1000')
    plt.legend()
    plt.figure('Npar=25')
    plt.grid()
    plt.xlabel("n")
    plt.ylabel("En")
    plt.errorbar(list(range(19)),E1,yerr=Err1,color='lightcoral',ecolor='lightcoral',label='Nconf=25')
    plt.legend()
    plt.figure('Npar=100')
    plt.grid()
    plt.xlabel("n")
    plt.ylabel("En")
    plt.errorbar(list(range(19)),E2,yerr=Err2,color='sienna',ecolor='sienna',label='Nconf=100')
    plt.legend()
    plt.figure('Npar=400')
    plt.grid()
    plt.xlabel("n")
    plt.ylabel("En")
    plt.errorbar(list(range(19)),E3,yerr=Err3,color='orange',ecolor='orange',label='Nconf=400')
    plt.legend()
    plt.figure('Npar=1000')
    plt.grid()
    plt.xlabel("n")
    plt.ylabel("En")
    plt.errorbar(list(range(19)),E4,yerr=Err4,color='indigo',ecolor='indigo',label='Nconf=1000')
    plt.legend()
    plt.show()
    #print(ConG[0],ConG[1])

def Extra():
    xarray=np.zeros(Npar)
    for i in range(10*Ncor):
        xarray=Metropolis(xarray)
    ConG=[]
    for i in range(10000):
        for j in range(Ncor):
            xarray=Metropolis(xarray)
        ConG.append(xarray.copy())
    JK,bar,err=Jacknife(ConG,10000)
    E=np.zeros(Npar-1)
    for i in range(Npar-1):
        E[i]=np.log(bar[i]/bar[i+1])
    err=JackE(JK,10000)
    '''
    Em=np.zeros(Npar-2)
    for i in range(1,Npar-1):
        a=G[i]/G[i+1]
        b=G[i-1]/G[i+1]
        Em[i-1]=np.log((1+b+np.sqrt((1+b)**2-4*a**2))/(2*a))
    plt.figure("Fine")
    plt.grid()
    plt.plot(list(range(1,Npar-1)),Em,color='red',label='Fine Delta E_n')
    plt.legend()
    '''
    plt.figure('Npar=10000')
    plt.grid()
    plt.errorbar(list(range(Npar-1)),E,err,color='blue',label="Delta E_n")
    plt.legend()
    plt.show()
    
'''
def expS(x,xlist):
    S=(x-xlist[0])**2+(x-xlist[-1])**2+0.25*(x**2)+0.25*(xlist[-1])**2
    for i in range(len(xlist)-1):
        S+=(xlist[i+1]-xlist[i])**2+0.25*(xlist[i]**2)
    return np.exp(-S)

def expSlat(x,xlist):#xlist被积变量
    S=0.25*(np.log((1+xlist[0])/(1-xlist[0]))-2*x)**2\
        +0.25*(x**2)+0.0625*(np.log((1+xlist[-1])/(1-xlist[-1]))**2)\
            +0.25*(2*x-np.log((1+xlist[-1])/(1-xlist[-1])))**2
    for i in range(len(xlist)-1):
        S+=0.25*(np.log((1+xlist[i+1])/(1-xlist[i+1]))-np.log((1+xlist[i])/(1-xlist[i])))**2+0.0625*(np.log((1+xlist[i])/(1-xlist[i]))**2)
    p=np.exp(-S)
    for i in xlist:
        p=p/(1-i**2)
    return p

def test():
    Vol=[[-10,10] for i in range(No-1)]
    integ=vegas.Integrator(Vol)
    nstrat=(No-1)*[5]
    def PathInt(x):
        integ(lambda y:expS(x,y),nitn=10,nstrat=nstrat)
        result=integ(lambda y:expS(x,y),nitn=10,nstrat=nstrat)
        #print(result.summary())
        return Ao*result
    print(PathInt(1))
'''

def TriK():
    A=pd.read_excel("Possibility-a.xlsx")
    B=pd.read_excel("Possibility-b.xlsx")
    Anum=0
    Bnum=0
    for key in A.columns:
        if A[key][0]==0:
            del A[key]
    for key in B.columns:
        if B[key][0]==0:
            del B[key]
    print("Number of triangle kinds in 1a:",len(A.keys()))
    print("Number of triangle kinds in 1a:",len(B.keys()))
    A.sort_values(by=0,axis=1,ascending=False,inplace=True)
    B.sort_values(by=0,axis=1,ascending=False,inplace=True)
    f1=open("texA.txt",'w+')
    f2=open("texB.txt",'w+')
    for i in range(7):
        f1.write("{}&{:.4}&{:.2}&{}&{:.4}&{:.2}&{}&{:.4}&{:.2}\\\\\n".format(\
            A.columns[i],100*A[A.columns[i]][0],A[A.columns[i]][1]/A[A.columns[3*i]][0],\
            A.columns[7+i],100*A[A.columns[7+i]][0],A[A.columns[7+i]][1]/A[A.columns[7+i]][0],\
            A.columns[14+i],100*A[A.columns[14+i]][0],A[A.columns[14+i]][1]/A[A.columns[14+i]][0]))
        f2.write("{}&{:.4}&{:.2}&{}&{:.4}&{:.2}&{}&{:.4}&{:.2}\\\\\n".format(\
            B.columns[i],100*B[B.columns[i]][0],B[B.columns[i]][1]/B[B.columns[i]][0],\
            B.columns[7+i],100*B[B.columns[7+i]][0],B[B.columns[7+i]][1]/B[B.columns[7+i]][0],\
            B.columns[14+i],100*B[B.columns[14+i]][0],B[B.columns[14+i]][1]/B[B.columns[14+i]][0]))
    f1.close()
    f2.close()
    f3=open("texC.txt",'w+')
    C=pd.read_excel("Ratio.xlsx",header=None)
    C=C.drop(index=[0])
    C.columns=['TriangleKind','Ratio']
    C=C.sort_values(by='Ratio',ascending=False)
    C.reset_index(drop=True,inplace=True)
    lenc=C.shape[0]
    groupc=lenc//5
    restc=lenc-5*groupc
    L=[0,groupc,2*groupc,3*groupc,4*groupc]
    for i in range(1,restc+1):
        for j in range(i,5):
            L[j]+=1
    for i in range(groupc):
        s=''
        for j in range(4):
            s+='{}&{:.4}&'.format(C['TriangleKind'][L[j]+i],C['Ratio'][L[j]+i])
        s+='{}&{:.4}\\\\\n'.format(C['TriangleKind'][L[4]+i],C['Ratio'][L[4]+i])
        f3.write(s)
    for j in range(restc):
        s+='{}&{:.4}&'.format(C['TriangleKind'][L[j+1]-1],C['Ratio'][L[j+1]-1])
    for j in range(restc,5):
        s+='~&~&'
    s+='\\\\\n'
    f3.write(s)
    f3.close()
#homework4_1_a()
#homework4_1_b()
#homework4_1_c()
#homework4_2_a()
#homework4_2_bc()
#TriK()
Extra()