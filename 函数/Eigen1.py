import math
from typing import ForwardRef
import numpy as np

eps=1e-14
lambdaG=0.61803399
lambdaGb=1-lambdaG

def is0(a):
    if abs(a)<eps:
        return 0
    else:
        return a

def MatAdd(A,B):
    R=[]
    for i in range(len(A)):
        c=[]
        for j in range(len(A[0])):
            c.append(is0(A[i][j]+B[i][j]))
        R.append(c)
    return R

def MatSub(A,B):
    R=[]
    for i in range(len(A)):
        c=[]
        for j in range(len(A[0])):
            c.append(is0(A[i][j]-B[i][j]))
        R.append(c)
    return R
def MatMul(A,B):
    if isinstance(A,list) or isinstance(A,np.ndarray):
        if isinstance(B,list) or isinstance(B,np.ndarray):
            R=[]
            for i in range(len(A)):
                c=[]
                for j in range(len(B[0])):
                    c.append(is0(sum([A[i][k]*B[k][j] for k in range(len(B))])))
                R.append(c)
            return R
        else:
            return [[is0(B*A[i][j]) for j in range(len(A[0]))] for i in range(len(A))]
    else:
        return [[is0(A*B[i][j]) for j in range(len(B[0]))] for i in range(len(B))]

class Matrix():
    def __init__(self,Ac):#
        if isinstance(Ac,list) or isinstance(Ac,np.ndarray):
            L=len(Ac)
        else:
            raise ValueError("Invalid input")
        if isinstance(Ac[0],list) or isinstance(Ac[0],np.ndarray):
            N=len(Ac[0])
            for i in Ac:
                if len(i)!=N:
                    raise IndexError("Invalid input")
        self.column=N#列的数目
        self.row=L#行的数目
        self.shape=(N,L)
        self._v=Ac.copy()
    def __getitem__(self,index):
        if isinstance(index,tuple):
            return self._v[index[0]][index[1]]
        else:
            raise IndexError("Indalid input")
    def __setitem__(self,index,value):
        if isinstance(index,tuple):
            self._v[index[0]][index[1]]=value
        else:
            raise IndexError("Indalid input")
    def __add__(self, other):
        if isinstance(other,Matrix):
            if self.shape==other.shape:
                return Matrix(MatAdd(self._v,other._v))
            else:
                raise ValueError("Dimension not match")
        else:
            raise ValueError("Not a matrix added")
    def __radd__(self, other):
        if isinstance(other,Matrix):
            if self.shape==other.shape:
                return Matrix(MatAdd(self._v,other._v))
            else:
                raise ValueError("Dimension not match")
        else:
            raise ValueError("Not a matrix added")
    def __sub__(self,other):
        if isinstance(other,Matrix):
            if self.shape==other.shape:
                return Matrix(MatSub(self._v,other._v))
            else:
                raise ValueError("Dimension not match")
        else:
            raise ValueError("Not a matrix added")
    def __rsub__(self,other):
        if isinstance(other,Matrix):
            if self.shape==other.shape:
                return Matrix(MatSub(other._v,self._v))
            else:
                raise ValueError("Dimension not match")
        else:
            raise ValueError("Not a matrix added")
    def __mul__(self,other):
        if isinstance(other,Matrix):
            if self.column==other.row:
                return Matrix(MatMul(self._v,other._v))
            else:
                raise ValueError("Dimension not match")
        elif isinstance(other,int) or isinstance(other,float):
            return Matrix(MatMul(self._v,other))
        else:
            raise ValueError("Invalid input")
    def __rmul__(self,other):
        if isinstance(other,Matrix):
            if self.column==other.row:
                return Matrix(MatMul(other._v,self._v))
            else:
                raise ValueError("Dimension not match")
        elif isinstance(other,int) or isinstance(other,float):
            return Matrix(MatMul(self._v,other))
        else:
            raise ValueError("Invalid input")
    def getcolumn(self,column):
        return Matrix([[self._v[i][column]] for i in range(self.row)])
    def getrow(self,row):
        return Matrix([self._v[row]])
    def __eq__(self,other):
        if isinstance(other,Matrix):
            return False
        else:
            for i in range(self.row):
                for j in range(self.column):
                    if self._v[i][j]!=other[i,j]:
                        return False
        return True
    def copy(self):
        return Matrix([[i for i in j] for j in self._v])
    def trans(self):
        return Matrix([[self._v[i][j] for i in range(self.row)] for j in range(self.column)])
    def __repr__(self):
        s=''
        for i in self._v:
            for j in i:
                s+=str(j)+'\t'
            s+='\n'
        return s
    def __str__(self):
        s=''
        for i in self._v:
            for j in i:
                s+=str(j)+'\t'
            s+='\n'
        return s

def Ones(n):
    o=[]
    for i in range(n):
        c=[]
        for j in range(i):
            c.append(0)
        c.append(1)
        for j in range(i+1,n):
            c.append(0)
        o.append(c)
    return Matrix(o)

def UniVec(dim,pos=0):
    u=[[0] for i in range(dim)]
    u[pos]=[1]
    return Matrix(u)

def FrobeniusMode(m:Matrix):
    c=0
    for i in m._v:
        for j in i:
            c+=j*j.conjugate()
    return pow(c,0.5)

def GivensMat(i,j,x,y,n):
    G=Ones(n)
    c,s=GivensValue(x,y)
    G[i,i]=c
    G[i,j]=s
    G[j,i]=-s
    G[j,j]=c
    return G

def GivensValue(x,y):
    if y!=0:
        if abs(x/y)<=1:
            t=x/abs(y)
            deno=pow(1+t**2,0.5)
            c=t/deno
            s=1/deno if y>0 else -1/deno
        else:
            t=y/abs(x)
            deno=pow(1+t**2,0.5)
            c=1/deno if x>0 else -1/deno
            s=t/deno
    elif x!=0:
        if abs(y/x)>=1:
            t=x/abs(y)
            deno=pow(1+t**2,0.5)
            c=t/deno
            s=1/deno if y>0 else -1/deno
        else:
            t=y/abs(x)
            deno=pow(1+t**2,0.5)
            c=1/deno if x>0 else -1/deno
            s=t/deno
    else:
        raise ValueError("(x,y) could not be (0,0)")
    return c,s#A 2*2 Givens matrix=[[c,s],[-s,c]]

def SparseMul(A,B):#稀疏乘法
    R=[]
    for i in range(A.row):
        c=[]
        for j in range(B.column):
            s=0
            for k in range(A.column):
                if A[i,k]==0 or B[k,j]==0:
                    pass
                else:
                    s+=A[i,k]*B[k,j]
            c.append(is0(s))
        R.append(c)
    return Matrix(R)

def HouseholderQRD(Ac):
    if isinstance(Ac,Matrix):
        A=Ac.copy()
    else:
        A=Matrix(Ac.copy())
    if A.row!=A.column:
        raise IndexError("Not a square matrix")
    N=A.column
    if N==1:
        return Matrix([[1]]),A
    else:
        v=A.getcolumn(0)
        v[0,0]=v[0,0]-FrobeniusMode(v)
        fbm=FrobeniusMode(v)
        if is0(fbm)==0:
            Q=Ones(N)
            R=A
        else:
            v=v*(1/fbm)
            Q=Ones(N)-(2*v)*v.trans()
            R=A-(2*v)*v.trans()*A
        for i in range(1,N-1):#不需要最后一列进行变换
            v=R.getcolumn(i)
            for j in range(i):
                v[j,0]=0
            v[i,0]=v[i,0]-FrobeniusMode(v)
            fbm=FrobeniusMode(v)
            if is0(fbm)==0:
                continue
            else:
                v=v*(1/fbm)
                Q=Q-(Q*v)*(2*v.trans())
                R=R-(2*v)*(v.trans()*R)
        return Q,R

def GivensQRD(Ac):
    if isinstance(Ac,Matrix):
        A=Ac.copy()
    else:
        A=Matrix(Ac.copy())
    if A.row!=A.column:
        raise IndexError("Not a square matrix")
    N=A.column
    if N==1:
        return Matrix([[1]]),A
    else:
        #按列进行Givens变换
        G=Ones(N)
        for i in range(N-1):#从第一列开始
            for j in range(N-1,i,-1):#从最后一行开始
                if is0(A[j,i])==0:
                    continue
                c,s=GivensValue(A[j-1,i],A[j,i])
                #G=SparseMul(G,GivensMat(j-1,j,A[j-1,i],-A[j,i],N))
                for k in range(N):
                    temp=G[k,j-1]
                    G[k,j-1]=is0(c*G[k,j-1]+s*G[k,j])
                    G[k,j]=is0(-s*temp+c*G[k,j])
                A[j-1,i]=pow(A[j-1,i]**2+A[j,i]**2,0.5)
                A[j,i]=0
                for k in range(i+1,N):
                    temp=A[j-1,k]
                    A[j-1,k]=is0(c*A[j-1,k]+s*A[j,k])
                    A[j,k]=is0(-s*temp+c*A[j,k])
        return G,A

def Hessenberg(Ac):
    if isinstance(Ac,Matrix):
        A=Ac.copy()
    else:
        A=Matrix(Ac.copy())
    if A.row!=A.column:
        raise IndexError("Not a square matrix")
    N=A.column
    if N==1:
        return Matrix([[1]]),A
    elif N==2:
        return Matrix([[1,0],[0,1]]),A
    else:
        v=A.getcolumn(0)
        v[0,0]=0
        v[1,0]=v[1,0]-FrobeniusMode(v)
        fbm=FrobeniusMode(v)
        if is0(fbm)==0:
            Q=Ones(N)
            H=A
        else:
            v=v*(1/fbm)
            Q=Ones(N)-(2*v)*v.trans()
            H=A-(2*v)*(v.trans()*A)-(A*v)*(v.trans()*2)+(v*4)*(((v.trans()*A)*v)*v.trans())
        #H=Q*A*Q
        for i in range(1,N-2):
            v=H.getcolumn(i)
            for j in range(i+1):
                v[j,0]=0
            v[i+1,0]=is0(v[i+1,0]-FrobeniusMode(v))
            fbm=FrobeniusMode(v)
            if is0(fbm)==0:
                continue
            else:
                v=v*(1/fbm)
                Q=Q-(Q*v)*(2*v.trans())
                H=H-(2*v)*(v.trans()*H)-(H*v)*(v.trans()*2)+(v*4)*(((v.trans()*H)*v)*v.trans())
        return Q,H

def GivensHessenberg(H):
    N=H.column
    T=H.copy()
    G=Ones(N)
    angl=[]#记录每次变换的givens矩阵
    for i in range(N-1):
        if is0(T[i+1,i])==0:
            angl.append((1,0))
            continue
        c,s=GivensValue(T[i,i],T[i+1,i])
        angl.append((c,s))
        #G=SparseMul(G,GivensMat(i,i+1,T[i,i],-T[i+1,i],N))
        '''
        for k in range(N):
            temp=G[k,i]
            G[k,i]=is0(c*G[k,i]+s*G[k,i+1])
            G[k,i+1]=is0(-s*temp+c*G[k,i+1])
        '''
        T[i,i]=is0(pow(T[i,i]**2+T[i+1,i]**2,0.5))
        T[i+1,i]=0
        for j in range(i+1,N):#左乘
            temp=T[i,j]
            T[i,j]=is0(c*T[i,j]+s*T[i+1,j])
            T[i+1,j]=is0(-s*temp+c*T[i+1,j])
    for i in range(N-1):#右乘
        c,s=angl[i]
        for j in range(i+2):
            temp=T[j,i]
            T[j,i]=is0(c*T[j,i]+s*T[j,i+1])
            T[j,i+1]=is0(-s*temp+c*T[j,i+1])
    #G=RotMat([angl])
    #print(G)
    #print(G*T*G.trans())
    return angl,T#G.T*H=R,T=RG=G.T*H*G, G.T=Rot(n-1,n)...Rot(1,2)*Rot(0,1); Rot=[[c,s],[-s,c]]

def RotMat(Angl):
    N=len(Angl[0])+1
    G=Ones(N)
    for angl in Angl:
        for i in range(N-1):
            c,s=angl[i]
            for k in range(N):#G=Rot.T(0,1)*Rot.T(1,2)...Rot.T(n-1,n) Rot.T=[[c,-s],[s,c]]
                temp=G[k,i]
                G[k,i]=is0(c*G[k,i]+s*G[k,i+1])
                G[k,i+1]=is0(-s*temp+c*G[k,i+1])
    return G
 
def HessenbergQR(Ac,times=10):
    Q,H=Hessenberg(Ac)
    angl,T=GivensHessenberg(H)
    Angl=[angl]
    for i in range(times-1):
        angl,T=GivensHessenberg(T)
        Angl.append(angl)
    return Q*RotMat(Angl),T#T=G.T*Ac*G=G.T*Q*H*Q*G

def SolveMaxEigen(Ac):
    if isinstance(Ac,Matrix):
        A=Ac.copy()
    else:
        A=Matrix(Ac.copy())
    if A.column!=A.row:
        raise ValueError("Not a square matrix")
    N=A.column
    if N==1:
        return A[0,0], UniVec(1,0)
    x=UniVec(N,0)
    q=A*x
    q=q*(1/FrobeniusMode(q))
    v=q.trans()*A*q
    while is0(FrobeniusMode(x-q)):
        x=q
        q=A*q
        q=q*(1/FrobeniusMode(q))
        v=(q.trans()*A*q)[0,0]
    #vl=[v]
    #ql=[q]
    maxv=v
    maxq=q
    for i in range(1,N):
        x=UniVec(N,0)
        q=A*x
        q=q*(1/FrobeniusMode(q))
        v=q.trans()*A*q
        while is0(FrobeniusMode(x-q)):
            x=q
            q=A*q
            q=q*(1/FrobeniusMode(q))
            v=(q.trans()*A*q)[0,0]
        if v>maxv:
            maxv=v
            maxq=q
    return maxv,maxq

def Minimize1D(f,x0,x1,x3):
    assert (x0<x1 and x1<x3),'Invalid section'
    assert (f(x0)>f(x1) and f(x1)<f(x3)),'Function value of middle is bigger than the edge'
    #print(x3-x0,f(x1))
    if (x3-x0)**2<eps*(abs(x0)+abs(x3))**2:
        return f(x1)
    if x3+x0>2*x1:
        x2=lambdaG*x1+lambdaGb*x3
        if f(x2)>f(x1):
            return Minimize1D(f,x0,x1,x2)
        else:
            return Minimize1D(f,x1,x2,x3)
    else:
        x2=lambdaG*x1+lambdaGb*x0
        if f(x2)>f(x1):
            return Minimize1D(f,x2,x1,x3)
        else:
            return Minimize1D(f,x0,x2,x1)

alpha=1.0
beta=2.0
gamma=0.5
sigma=0.5

def IterMinimize(f,args:list):
    n=len(args)-1
    ixm=0#最小值坐标index
    ixmu=0#次大值坐标index
    ixM=0#最大值坐标index
    fm=f(*args[0])#最小值
    fmu=fm#次大值
    fM=fm#最大值
    xc=np.zeros(n)#中心点坐标
    fbar=0
    #delta=sum((f(xi)-fbar)**2)/n=sum(f(xi)**2-fbar**2)/n
    for i in range(n+1):
        fv=f(*args[i])
        if fv<fm:
            fm=fv
            ixm=i
        elif fv>fM:
            fM=fv
            ixM=i
        elif fv<fM and fv>fmu:
            fmu=fv
            ixmu=i
        fbar+=fv
        xc+=args[i]
    fbar=fbar/(n+1)
    delta=sum([(f(*arg)-fbar)**2 for arg in args])/n#避免大数相消
    if delta<eps:
        lowlim=args[0]
        uplim=args[0]
        for arg in args:
            for i in range(1,n):
                if arg[i]<lowlim[i]:
                    lowlim[i]=arg[i]
                if arg[i]>uplim[i]:
                    uplim[i]=arg[i]
        er=uplim-lowlim
        return args[ixm],er,f(*args[ixm])#小于临界则返回最小值
    else:#反射
        fmu=f(*args[ixmu])
        fM=f(*args[ixM])
        xc=(xc-args[ixM])/n
        xr=(1+alpha)*xc-alpha*args[ixM]
        fr=f(*xr)
        #扩展
        if fr<fmu and fr>fm:#去掉xM,换上xr
            args[ixM]=xr
        elif fr<fm:
            xe=beta*xr+(1-beta)*xc
            if f(*xe)<fm:
                args[ixm]=xe
            else:
                args[ixM]=xr
        else:#收缩
            if fr<fM:
                xcon=gamma*xr+(1-gamma)*xc
            else:
                xcon=gamma*args[ixM]+(1-gamma)*xc
            fcon=f(*xcon)
            if fcon<fM and fcon<fr:
                args[ixM]=xcon
            else:#收缩
                for i in range(ixm):
                    args[i]=sigma*args[i]+(1-sigma)*args[ixm]
                for i in range(ixm+1,n+1):
                    args[i]=sigma*args[i]+(1-sigma)*args[ixm]
        return IterMinimize(f,args)
            


def SimplexMinimize(f,args):#args是list
    assert len(args)-1==len(args[0])
    n=len(args)-1
    args=list(args)
    for i in range(n+1):
        args[i]=np.array(args[i])
    #规范输入格式为args=[np.ndarray,...]
    return IterMinimize(f,args)


"""f=lambda x,y:np.exp(x**2+x*y)+y**2
ags=[[1,2],[3,4],[-1,-2]]
print(SimplexMinimize(f,ags))"""
    




#t=[lambda x: x**2+2*x+1,-5,1,3]
#print(Minimize1D(*t))


#print(SolveMaxEigen([[1,2],[2,1]]))

C=Matrix([[1,2,3],[4,5,6],[7,8,9]])
''''
print(C)
print(GivensQRD(C)[0])
print(GivensQRD(C)[1])

print(HouseholderQRD(C)[0])
print(HouseholderQRD(C)[1])

print(GivensQRD(C)[0]*GivensQRD(C)[1])
print(HouseholderQRD(C)[0]*HouseholderQRD(C)[1])

print(Hessenberg(C)[0])
print(Hessenberg(C)[1])
print(Hessenberg(C)[0]*Hessenberg(C)[1]*Hessenberg(C)[0])
'''


P,R=HessenbergQR(Matrix([[1,2,3],[4,5,6],[7,8,9]]),30)
print(P)
print(R)
print(P*R*P.trans())


#print(GivensHessenberg(Matrix([[1,2,3],[4,5,6],[0,7,8]])))