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

def UniVec(dim,pos=0):
    u=[[0] for i in range(dim)]
    u[pos]=[1]
    return np.matrix(u)

def FrobeniusMode(m:np.ndarray):
    c=0
    N=m.shape[0]#行数
    M=m.shape[1]#列数
    for i in range(N):
        for j in range(M):
            c+=m[i,j]*m[i,j].conjugate()
    return pow(c,0.5)

def GivensMat(i,j,x,y,n):
    G=np.matrix(np.identity(n))
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

def HouseholderQRD(Ac):
    if isinstance(Ac,np.matrix):
        A=Ac.copy()
    else:
        A=np.matrix(Ac)
    if A.shape[0]!=A.shape[1]:
        raise IndexError("Not a square matrix")
    N=A.shape[1]#列
    if N==1:
        return np.matrix([[1]]),A
    else:
        v=A[:,0].copy()
        v[0,0]=v[0,0]-FrobeniusMode(v)
        fbm=FrobeniusMode(v)
        if is0(fbm)==0:
            Q=np.matrix(np.identity(N))
            R=A
        else:
            v=v*(1/fbm)
            Q=np.matrix(np.identity(N))
            R=A
            for j in range(N):
                    #s:Q*v,p:v.T*R
                    s=0
                    p=0
                    for k in range(N):
                        s+=Q[j,k]*v[k,0]#计算Q*v的第j个元素
                        p+=v[k,0]*R[k,j]#计算v.t*R的第j个元素
                    for k in range(N):
                        Q[j,k]=is0(Q[j,k]-2*s*v[k,0])#Q=Q-2*Q*v*v.t-->Q[j,k]=Q[j,k]-2*(Q*v)[j,0]*v.t[0,k]
                        R[k,j]=is0(R[k,j]-2*p*v[k,0])#R=R-2*v*v.t*R-->R[k,j]=R[k,j]-2*v[k,0]*(v.t*R)[0,j]
        #print(Q,'\n',R,'\n',Q*R)
        for i in range(1,N-1):#不需要最后一列进行变换
            v=R[:,i].copy()
            for j in range(i):
                v[j,0]=0
            v[i,0]=v[i,0]-FrobeniusMode(v)
            fbm=FrobeniusMode(v)
            if is0(fbm)==0:
                continue
            else:
                v=v*(1/fbm)
                for j in range(i,N):
                    p=0#v.T*R
                    for k in range(i,N):
                        p+=v[k,0]*R[k,j]#计算v.t*R的第j个元素
                    p=2*p
                    for k in range(i,N):
                        R[k,j]=is0(R[k,j]-p*v[k,0])#R=R-2*v*v.t*R-->R[k,j]=R[k,j]-2*v[k,0]*(v.t*R)[0,j]
                for j in range(N):
                    s=0#Q*v
                    for k in range(i,N):
                        s+=Q[j,k]*v[k,0]#计算Q*v的第j个元素
                    s=s*2
                    for k in range(i,N):
                        Q[j,k]=is0(Q[j,k]-s*v[k,0])#Q=Q-2*Q*v*v.t-->Q[j,k]=Q[j,k]-2*(Q*v)[j,0]*v.t[0,k]
            #print(Q,'\n',R,'\n',Q*R)
        return Q,R

def GivensQRD(Ac):
    if isinstance(Ac,np.matrix):
        A=Ac.copy()
    else:
        A=np.matrix(Ac)
    if A.shape[0]!=A.shape[1]:
        raise IndexError("Not a square matrix")
    N=A.shape[1]
    if N==1:
        return np.matrix([[1]]),A
    else:
        #按列进行Givens变换
        G=np.matrix(np.identity(N))
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
            #print(G,'\n',A)
        return G,A



def Hessenberg(Ac):
    if isinstance(Ac,np.matrix):
        A=Ac.copy()
    else:
        A=np.matrix(Ac)
    if A.shape[0]!=A.shape[1]:
        raise IndexError("Not a square matrix")
    N=A.shape[1]
    if N==1:
        return np.matrix([[1]]),A
    elif N==2:
        return np.matrix([[1,0],[0,1]]),A
    else:
        v=A[:,0].copy()
        v[0,0]=0
        v[1,0]=v[1,0]-FrobeniusMode(v)
        fbm=FrobeniusMode(v)
        if is0(fbm)==0:
            Q=np.matrix(np.identity(N))
            H=A
        else:
            v=v*(1/fbm)
            Q=np.matrix(np.identity(N))
            for j in range(N):
                s=0#Q*v
                for k in range(1,N):
                    s+=Q[j,k]*v[k,0]#计算Q*v的第j个元素
                s=2*s
                for k in range(1,N):
                    Q[j,k]=is0(Q[j,k]-s*v[k,0])#Q=Q-2*Q*v*v.t-->Q[j,k]=Q[j,k]-2*(Q*v)[j,0]*v.t[0,k]
            H=A
            for j in range(N):#H'=(1-2*v*v.T)H, H[i,i-1]不为0
                p=0#v.T*H
                for k in range(1,N):
                    p+=v[k,0]*H[k,j]#p=(v.T*H)[j]
                p=2*p
                for k in range(1,N):
                    H[k,j]=is0(H[k,j]-p*v[k,0])#H[k,j]=(H-2*v*v.T*H)[k,j]=H[k,j]-2*v[k](v.T*H)[j]
            for j in range(N):#H=H'(1-2*v*v.T)
                p=0#H'*v
                for k in range(1,N):
                    p+=H[j,k]*v[k,0]#p=(H'*v)[j]
                p=2*p
                for k in range(1,N):
                    H[j,k]=H[j,k]-p*v[k,0]
        #H=Q*A*Q
        for i in range(1,N-2):
            v=H[:,i].copy()
            for j in range(i+1):
                v[j,0]=0
            v[i+1,0]=is0(v[i+1,0]-FrobeniusMode(v))
            fbm=FrobeniusMode(v)
            if is0(fbm)==0:
                continue
            else:
                v=v*(1/fbm)
                #Q=Q-(Q*v)*(2*v.T)
                for j in range(N):
                    s=0#Q*v
                    for k in range(i,N):
                        s+=Q[j,k]*v[k,0]#计算Q*v的第j个元素
                    s=2*s
                    for k in range(i,N):
                        Q[j,k]=is0(Q[j,k]-s*v[k,0])#Q=Q-2*Q*v*v.t-->Q[j,k]=Q[j,k]-2*(Q*v)[j,0]*v.t[0,k]
                #H=(1-2*v*v.T)*H*(1-2*v*v.T)
                #H=H-(2*v)*(v.T*H)-(H*v)*(v.T*2)+(v*4)*(((v.T*H)*v)*v.T)
                for j in range(i,N):#H'=(1-2*v*v.T)H, H[i,i-1]不为0
                    p=0#v.T*H
                    for k in range(i+1,N):
                        p+=v[k,0]*H[k,j]#p=(v.T*H)[j]
                    p=2*p
                    for k in range(i+1,N):
                        H[k,j]=is0(H[k,j]-p*v[k,0])#H[k,j]=(H-2*v*v.T*H)[k,j]=H[k,j]-2*v[k](v.T*H)[j]
                for j in range(N):#H=H'(1-2*v*v.T)
                    p=0#H'*v
                    for k in range(i+1,N):
                        p+=H[j,k]*v[k,0]#p=(H'*v)[j]
                    p=2*p
                    for k in range(i+1,N):
                        H[j,k]=H[j,k]-p*v[k,0]
        return Q,H

def GivensHessenberg(H):
    N=H.shape[0]
    T=H.copy()
    G=np.matrix(np.identity(N))
    angl=[]#记录每次变换的givens矩阵
    for i in range(N-1):
        if is0(T[i+1,i])==0:
            angl.append((1,0))
            continue
        c,s=GivensValue(T[i,i],T[i+1,i])
        angl.append((c,s))
        #G=SparseMul(G,GivensMat(i,i+1,T[i,i],-T[i+1,i],N))
        for k in range(N):
            temp=G[k,i]
            G[k,i]=is0(c*G[k,i]+s*G[k,i+1])
            G[k,i+1]=is0(-s*temp+c*G[k,i+1])
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
    return angl,T#R=TG

def RotMat(Angl):
    N=len(Angl[0])+1
    G=np.matrix(np.identity(N))
    for angl in Angl:
        for i in range(N-1):
            c,s=angl[i]
            for k in range(N):
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
    return Q*RotMat(Angl),T

def SolveMaxEigen(Ac):
    if isinstance(np.matrix):
        A=Ac.copy()
    else:
        A=np.matrix(Ac)
    if A.shape[0]!=A.shape[1]:
        raise ValueError("Not a square matrix")
    N=A.shape[0]
    if N==1:
        return A[0,0], UniVec(1,0)
    x=UniVec(N,0)
    q=A*x
    q=q*(1/FrobeniusMode(q))
    v=q.T*A*q
    while is0(FrobeniusMode(x-q)):
        x=q
        q=A*q
        q=q*(1/FrobeniusMode(q))
        v=(q.T*A*q)[0,0]
    #vl=[v]
    #ql=[q]
    maxv=v
    maxq=q
    for i in range(1,N):
        x=UniVec(N,0)
        q=A*x
        q=q*(1/FrobeniusMode(q))
        v=q.T*A*q
        while is0(FrobeniusMode(x-q)):
            x=q
            q=A*q
            q=q*(1/FrobeniusMode(q))
            v=(q.T*A*q)[0,0]
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
        lowlim=args[0].copy()
        uplim=args[0].copy()
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

def SeekRoots(f,a,b,ep=1e-14):
    A=f(a)
    B=f(b)
    C=f((a+b)/2)
    #print(A,B,C)
    #print(abs(C)<ep)
    assert A*B<=0,print(A,B)
    if abs(C)<ep:
        #print(True)
        return (a+b)/2
    elif A*C<=0:
        return SeekRoots(f,a,(a+b)/2,ep)
    elif B*C<=0:
        return SeekRoots(f,(a+b)/2,b,ep)
    else:
        raise ValueError("No Roots")

#print(SeekRoots(lambda x:x*x-2*x-5,0,6))


"""f=lambda x,y:np.exp(x**2+x*y)+y**2
ags=[[1,2],[3,4],[-1,-2]]
print(SimplexMinimize(f,ags))"""
    




#t=[lambda x: x**2+2*x+1,-5,1,3]
#print(Minimize1D(*t))


#print(SolveMaxEigen([[1,2],[2,1]]))

#C=np.matrix([[1.0,2.0,3.0],[4.0,5.0,6.0],[7.0,8.0,9.0]])
'''
print(C)
#x,y=GivensQRD(C)
x,y=HouseholderQRD(C)
print(x)
print(y)

print(x*y)
'''
'''
print(Hessenberg(C)[0])
print(Hessenberg(C)[1])
print(Hessenberg(C)[0]*Hessenberg(C)[1]*Hessenberg(C)[0])

'''
'''
Q,H=HessenbergQR(C,30)
print(Q)
print(H)
print(Q*H*Q.T)
'''
