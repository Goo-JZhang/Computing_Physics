import LinearEqs
import matplotlib.pyplot as plt
import math
import numpy as np
#常数
eps=1e-14
ieps=1e-7
#多项式
def norm(a:list)->list :
    N=len(a)
    if N==1:
        return a
    else:
        for i in range(N-1,0,-1):
            if abs(a[i])<=eps:
                a.pop(i)
            else:
                return a#从最高项开始到一次，去除0项直到找到一个非零项
        return a

def PolyAdd(a,b)->list:
    if isinstance(a,list):
        if isinstance(b,list):
            if len(a)>=len(b):
                L=a
                S=b
            else:
                L=b
                S=a
            N1=len(S)
            N2=len(L)
            newcoef=[]
            for i in range(N1):
                newcoef.append(L[i]+S[i])
            for i in range(N1,N2):
                newcoef.append(L[i])
            return newcoef
        else:
            newcoef=a[:]
            newcoef[0]+=b
            return newcoef
    else:
        if isinstance(b,list):
            newcoef=b[:]
            newcoef[0]+=a
            return newcoef
        else:
            return [a+b]

def PolySub(a,b)->list:
    if isinstance(a,list):
        if isinstance(b,list):
            if len(a)>=len(b):
                N2=len(a)
                N1=len(b)
                p=True
            else:
                N2=len(b)
                N1=len(a)
                p=False
            newcoef=[]
            for i in range(N1):
                newcoef.append(a[i]-b[i])
            if p:
                for i in range(N1,N2):
                    newcoef.append(a[i])
            else:
                for i in range(N1,N2):
                    newcoef.append(-b[i])
            return newcoef
        else:
            newcoef=a[:]
            newcoef[0]-=b
            return newcoef
    else:
        if isinstance(b,list):
            newcoef=b[:]
            newcoef[0]-=a
            return newcoef
        else:
            return [a-b]

def PolyMul(a,b)->list:
    if isinstance(a,list):
        if isinstance(b,list):
            N1=len(a)
            N2=len(b)
            newcoef=[]
            for i in range(N1+N2-1):
                p=0
                for j in range(max(0,i-N2+1),min(i+1,N1)):
                    p+=a[j]*b[i-j]
                newcoef.append(p)
            return newcoef
        else:
            return [i*b for i in a]
    else:
        if isinstance(b,list):
            return [i*a for i in b]
        else:
            return [a*b]

def PolyCal(b,x):
    if len(b)==0:
        return 0
    else:
        return b[0]+x*PolyCal(b[1:],x)

def SortPoint(x,y):
    if len(x)!=len(y):
        raise IndexError("Length not fit")
    N=len(x)
    for i in range(N):
        temp=x[i]
        for j in range(i+1,N):
            if x[j]<temp:
                x[j],x[i]=x[i],x[j]
                y[j],y[i]=y[i],y[j]
                temp=x[i]
    return x,y

def FindInterval(xlist,x,start=0):
    if x<xlist[0] or x>xlist[-1]:
        raise ValueError("Out of range")
    else:
        for i in range(start,len(xlist)-1):
            if x>=xlist[i] and x<=xlist[i+1]:
                return i


class Polynomial():
    def __init__(self,coefficient=[0]):
        self.coef=norm(coefficient)
        self.degree=len(self.coef)-1

    def __add__(self,other):
        if isinstance(other,Polynomial):
            return Polynomial(norm(PolyAdd(self.coef,other.coef)))
        else:
            return Polynomial(norm(PolyAdd(self.coef,other)))

    def __radd__(self,other):
        if isinstance(other,Polynomial):
            return Polynomial(norm(PolyAdd(self.coef,other.coef)))
        else:
            return Polynomial(norm(PolyAdd(self.coef,other)))

    def __sub__(self,other):
        if isinstance(other,Polynomial):
            return Polynomial(norm(PolySub(self.coef,other.coef)))
        else:
            return Polynomial(norm(PolySub(self.coef,other)))

    def __rsub__(self,other):
        if isinstance(other,Polynomial):
            return Polynomial(norm(PolySub(self.coef,other.coef)))
        else:
            return Polynomial(norm(PolySub(self.coef,other)))

    def __mul__(self,other):
        if isinstance(other,Polynomial):
            return Polynomial(norm(PolyMul(self.coef,other.coef)))
        else:
            return Polynomial(norm(PolyMul(self.coef,other)))

    def __rmul__(self,other):
        if isinstance(other,Polynomial):
            return Polynomial(norm(PolyMul(self.coef,other.coef)))
        else:
            return Polynomial(norm(PolyMul(self.coef,other)))

    def __repr__(self):
        return str(self.coef)

    def str(self,N=15):
        s=str(round(self.coef[0],N))
        for i in range(1,len(self.coef)):
            if self.coef[i]>eps:
                s+=" +"+str(round(self.coef[i],N))+'x**{}'.format(i)
            elif self.coef[i]<-eps:
                s+=" -"+str(round(-self.coef[i],N))+'x**{}'.format(i)
            else:
                pass
        return s

    def cal(self,x):
        return PolyCal(self.coef,x)

    def dif(self):
        p=[]
        if len(self.coef)>1:
            for i in range(1,len(self.coef)):
                p.append(i*self.coef[i])
            return Polynomial(norm(p))
        else:
            return Polynomial()
    
    def difcal(self,x):
        return self.dif().cal(x)



#拉格朗日内插
def LagrangeInterpolation(x:list,y:list)->Polynomial:
    if len(x)!=len(y):
        raise IndexError("Length not fit")
    N=len(x)
    L=Polynomial()
    for i in range(N):
        tempL=Polynomial([1])
        for j in range(i):
            tempL=tempL*Polynomial([-x[j],1])*(1/(x[i]-x[j]))
        for j in range(i+1,N):
            tempL=tempL*Polynomial([-x[j],1])*(1/(x[i]-x[j]))
        L+=y[i]*tempL
    return L

def NewtonAdd(x:list,y:list,a:list=[],ns:list=[]):
    if len(x)!=len(y)+1:
        raise IndexError("Length not fit")
    elif len(x)==1:
        N=Polynomial()
        for i in range(len(a)):
            N+=a[i]*ns[i]
        return a,ns,N
    else:
        ns.append(Polynomial([-x[0],1])*ns[-1])
        aend=y[0]
        for i in range(len(a)):
            aend-=a[i]*ns[i].cal(x[1])
        a.append(aend/ns[-1].cal(x[1]))
        return NewtonAdd(x[1:],y[1:],a,ns)

def NewtonInterpolation(x:list,y:list):
    if len(x)==0:
        raise IndexError("No point")
    else:
        return NewtonAdd(x,y[1:],[y[0]],[Polynomial([1])])

def AddNeville(x:list,y:list,xlist,Nmap):#Nmap[i][j]=T_{i+j,i}
    if len(x)!=len(y):
        raise IndexError("Length not fit")
    elif len(x)==0:
        return xlist,Nmap
    else:
        Nmap[0].append(Polynomial([y[0]]))
        xlist.append(x[0])
        L=len(xlist)-1
        for i in range(1,L):
            Nmap[i].append((Polynomial([-xlist[L-i],1])*Nmap[i-1][-1]-Polynomial([-xlist[L],1])*Nmap[i-1][-2])*(1/(xlist[L]-xlist[L-i])))
        if L!=0:
            Nmap.append([(Polynomial([-xlist[0],1])*Nmap[L-1][-1]-Polynomial([-xlist[L],1])*Nmap[L-1][-2])*(1/(xlist[L]-xlist[0]))])
        return AddNeville(x[1:],y[1:],xlist,Nmap)

def Neville(x:list,y:list):
    if len(x)==0:
        raise IndexError("No point")
    else:
        return AddNeville(x,y,[],[[]])

def Spline3Interpolation(x,y,arg=1,dy=[0,0]):#arg=1 means boundary condition with S''=0, arg=2 means Periodic boundary condition, arg=3 means boundary condition with S'=y' 
    x,y=SortPoint(x,y)
    N=len(x)
    if N<2:
        raise IndexError("Number of points not enough")
    H=[[0 for i in range(N)] for j in range(N)]
    d=[0 for i in range(N)]
    for i in range(1,N-1):
        H[i][i-1]+=(x[i]-x[i-1])/(x[i+1]-x[i-1])
        H[i][i]+=2
        H[i][i+1]+=(x[i+1]-x[i])/(x[i+1]-x[i-1])
        d[i]+=6*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))/(x[i+1]-x[i-1])
    if arg==2:#S''(x[0]+0)=S''(x[-1]-0), S'(x[0]+0)=S'(x[-1]-0), y[0]=y[-1]
        #if y[0]!=y[-1]:
        #    raise ValueError("Not periodic")
        #else:
            H[0][0]+=1
            H[0][-1]+=-1
            d[0]+=0#M[0]=M[-1]
            H[-1][1]+=(x[1]-x[0])/((x[1]-x[0])+(x[-1]-x[-2]))
            H[-1][-2]+=(x[-1]-x[-2])/((x[1]-x[0])+(x[-1]-x[-2]))
            H[-1][-1]+=2
            d[-1]+=6*((y[1]-y[0])/(x[1]-x[0])-(y[-1]-y[-2])/(x[-1]-x[-2]))/((x[1]-x[0])+(x[-1]-x[-2]))
    elif arg==3:#S'(x[0]+0)=dy[0], S'(x[-1]-0)=dy[1]
        H[0][0]+=2
        H[0][1]+=1
        d[0]+=-6*(dy[0]-(y[1]-y[0])/(x[1]-x[0]))/(x[1]-x[0])
        H[-1][-1]+=2
        H[-1][-2]+=1
        d[-1]+=6*(dy[1]-(y[-1]-y[-2])/(x[-1]-x[-2]))/(x[-1]-x[-2])
    else:#S''(x[0]+0)=S''(x[-1]-0)=0
        H[0][0]+=1
        d[0]+=0#HM=d and we need M[0]=0
        H[-1][-1]+=1
        d[-1]+=0#HM=d and we need M[-1]=0
    M=LinearEqs.GaussianElimination(H,d)
    #print(M)
    B=[y[i]-M[i]*(x[i+1]-x[i])*(x[i+1]-x[i])/6 for i in range(N-1)]
    A=[(y[i+1]-y[i])/(x[i+1]-x[i])-(x[i+1]-x[i])*(M[i+1]-M[i])/6 for i in range(N-1)]
    result=[]
    for i in range(N-1):
        P1=Polynomial([x[i+1],-1])*Polynomial([x[i+1],-1])*Polynomial([x[i+1],-1])*(M[i]/(6*(x[i+1]-x[i])))
        P2=Polynomial([-x[i],1])*Polynomial([-x[i],1])*Polynomial([-x[i],1])*(M[i+1]/(6*(x[i+1]-x[i])))
        P3=Polynomial([-x[i],1])*A[i]
        result.append(P1+P2+P3+B[i])
    return result

def Spline3Cal(xlist,spinelist,inputlist):
    if isinstance(inputlist,list):
        j=0
        result=[]
        for i in range(len(inputlist)):
            j=FindInterval(xlist,inputlist[i],j)
            result.append(spinelist[j].cal(inputlist[i]))
        return result
    else:
        j=FindInterval(xlist,inputlist)
        return spinelist[j].cal(inputlist)

def ChebyshevPolynomial(n):
    if n==0:
        return Polynomial([1])
    elif n==1:
        return Polynomial([0,1])
    else:
        return Polynomial([0,2])*ChebyshevPolynomial(n-1)-ChebyshevPolynomial(n-2)

def ChebyshevList(n):
    if n==0:
        return [Polynomial([1])]
    elif n==1:
        return [Polynomial([1]),Polynomial([0,1])]
    else:
        result=[Polynomial([1]),Polynomial([0,1])]
        for i in range(2,n+1):
            result.append(Polynomial([0,2])*result[i-1]-result[i-2])
        return result

def ChebyshevInterpolation(n,y):
    if n==0:
        raise IndexError("No point")
    CL=ChebyshevList(n-1)
    CP=Polynomial()
    temp=0
    for j in range(n):
        temp+=y(math.cos(math.pi*(0.5+j)/n))
    CP+=temp*CL[0]*(1/n)
    for i in range(1,n):
        temp=0
        for j in range(n):
            temp+=y(math.cos(math.pi*(0.5+j)/n))*math.cos(math.pi*i*(j+0.5)/n)
        CP+=temp*CL[i]*(2/n)
    return CP

def TrapezoidalIntegral(floor,ceil,f,N):
    x=np.linspace(floor,ceil,N+1)
    return (np.sum(f(x))-0.5*(f(floor)+f(ceil)))*(ceil-floor)/N

def TrapInt(floor,ceil,f):
    N=4
    if floor==ceil:
        return 0
    tempint=TrapezoidalIntegral(floor,ceil,f,2)
    while abs((ceil-floor)/N)>eps:
        nextint=TrapezoidalIntegral(floor,ceil,f,N)
        if abs(nextint-tempint)<=ieps:
            return nextint
        else:
            N=N*2
            tempint=nextint
    raise ValueError("Do not converge")

def NevilleInt(x,y,xlist,vlist):
    if len(x)!=len(y):
        raise IndexError("Length not fit")
    elif len(x)==0:
        return xlist,vlist
    else:
        vlist[0].append(y[0])
        xlist.append(x[0])
        L=len(xlist)-1
        for i in range(1,L):
            vlist[i].append((-xlist[L-i]*vlist[i-1][-1]+xlist[L]*vlist[i-1][-2])/(xlist[L]-xlist[L-i]))
        if L!=0:
            vlist.append([(-xlist[0]*vlist[L-1][-1]+xlist[L]*vlist[L-1][-2])/(xlist[L]-xlist[0])])
        return NevilleInt(x[1:],y[1:],xlist,vlist)

def ExtraPolationIntegral(floor,ceil,f,N):
    if floor==ceil:
        return 0
    h=[((ceil-floor)/2**i)**2 for i in range(N)]
    y=[TrapezoidalIntegral(floor,ceil,f,2**i) for i in range(0,N)]
    return NevilleInt(h,y,[],[[]])[1][-1][0]

def EPInt(floor,ceil,f):
    if floor==ceil:
        return 0
    tempint=[TrapezoidalIntegral(floor,ceil,f,1)]
    xlist,Nev=NevilleInt([((ceil-floor))**2],tempint,[],[[]])
    N=2
    h=((ceil-floor)/N)**2
    while h>eps:
        tempint=[TrapezoidalIntegral(floor,ceil,f,N)]
        xlist,Nev=NevilleInt([h],tempint,xlist,Nev)
        if abs(Nev[-1][0]-Nev[-2][0])<ieps:
            return Nev[-1][0]
        else:
            N=N*2
            h=((ceil-floor)/N)**2
    raise ValueError("Do not converge")

def GaussianInt(x,w,f):
    if len(x)!=len(w):
        raise IndexError("Length not fit")
    return sum([f(x[i])*w[i] for i in range(len(x))])

x=[1,2,3]
y=[6,1,6]

#print(TrapInt(-10,10,lambda y:np.exp(-y*y)*np.cos(y)))
#print(ExtraPolationIntegral(-10,10,lambda y:np.exp(-y*y)*np.cos(y),8))
#print(EPInt(-10,10,lambda y:np.exp(-y*y)*np.cos(y)))
#xlist,Nev=Neville([400.0, 100.0, 44.44444444444445, 25.0, 16.0, 11.111111111111112, 8.16326530612245, 6.25],[-6.242819674960758e-43, 10.0, -0.0001956193376184662, 5.000000000039395, -0.0609759614630906, 3.333235523664524, 0.1052342418452217, 2.4922671449879723])
#print(Nev[-1][0].coef[0])

'''
print(LagrangeInterpolation(x,y))
print(NewtonInterpolation(x,y))
print(Neville(x,y))

P=Spine3Interpolation(x,y,2)
print(P)
print(P[0].dif().cal(1),P[1].dif().cal(3))
plt.plot(Spine3Cal(x,P,list(np.linspace(1,3))))
plt.show()
'''










    