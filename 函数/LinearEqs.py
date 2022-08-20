import math

eps=1e-14#非0判断
###高斯消元
def maxindex(a):#输入一维list，输出最大值对应的索引号
    count=0
    index=0
    maxabs=0
    for i in a:
        if abs(i)>maxabs:
            count=index
            maxabs=abs(i)
        index+=1
    if maxabs==0:
        print(a)
        raise IndexError("Not a full map")
    else:
        return count

def GaussianElimination(Ac,bc):#Ac:二维list，作为系数矩阵;bc为一维list
    A=Ac.copy()
    b=bc.copy()
    N=len(A)#A的行数
    L=len(b)#b的长度
    if N!=L:#判断长度是否匹配
        raise IndexError("Dimension does not match")
    else:#判断列数是否匹配
        for i in A:
            if len(i)!=L:
                #print("Not a square matrix")
                raise IndexError("Not a square matrix")
    
    for i in range(N):
        t=i+maxindex([A[k][i] for k in range(i,N)])
        A[t], A[i] = A[i], A[t]
        b[t], b[i] = b[i], b[t]
        for j in range(i+1,N):
            l=A[j][i]/A[i][i]
            b[j]=b[j]-l*b[i]
            for k in range(i,N):
                if A[i][k]!=0:
                    A[j][k]=A[j][k]-l*A[i][k]
    for i in range(N-1,0,-1):
        l=b[i]/A[i][i]
        for j in range(i):
            b[j]=b[j]-l*A[j][i]
    return [b[i]/A[i][i] for i in range(N)]


###Cholestky分解
def SolveUp(Uc,hc):#解上三角，Ac:二维list，作为系数矩阵;bc为一维list
    U=Uc.copy()
    h=hc.copy()
    N=len(U)
    M=len(h)
    if N!=M:#方阵判断
        raise IndexError("Dimension does not match")
    else:
        for i in U:
            if len(i)!=N:
                raise IndexError("Not a square matrix")

    for i in range(N):
        if abs(U[i][i])<eps:#判断对角元是否为0
            raise ValueError("Not a full map")
        for j in range(i):#判断是否上三角
            if U[i][j]!=0:
                raise ValueError("Not a upper triangle matrix")
    #解方程
    for i in range(N-1,0,-1):#遍历列
        l=h[i]/U[i][i]
        for j in range(i):#行相减
            h[j]=h[j]-l*U[j][i]
    return [h[i]/U[i][i] for i in range(N)]

def SolveLow(Lc,hc):#解下三角，Ac:二维list，作为系数矩阵;bc为一维list
    L=Lc.copy()
    h=hc.copy()
    N=len(L)
    M=len(h)
    if N!=M:#方阵判断
        raise IndexError("Dimension does not match")
    else:
        for i in L:
            if len(i)!=N:
                raise IndexError("Not a square matrix")
    for i in range(N):
        if abs(L[i][i])<eps:#判断对角元是否为0
            raise ValueError("Not a full map")
        for j in range(i+1,N):#判断是否下三角
            if L[i][j]!=0:
                raise ValueError("Not a lower triangle matrix")
    #解方程
    for i in range(N-1):#遍历列
        l=h[i]/L[i][i]
        for j in range(i+1,N):#行相减
            h[j]=h[j]-l*L[j][i]
    return [h[i]/L[i][i] for i in range(N)]

def Cholesky(A):#给出A的Cholesky分解的上三角矩阵，输入A为一个二维list
    N=len(A)
    for i in range(N):
        if len(A[i])!=N:#方阵判断
            raise IndexError("Not a square matrix")
        else:
            for j in range(i,N):
                if abs(A[i][j]-A[j][i].conjugate())>eps:#判断是否厄米
                    #print(i,j,A[i][j],A[j][i])
                    raise ValueError("Not an hermitian matrix")
    H=[[0 for i in range(N)] for j in range(N)]#创建空矩阵
    if A[0][0]<0:
        raise ValueError("Not a positive defined matrix")
    else:
        H[0][0]=pow(A[0][0],0.5)
    for i in range(1,N):#从2*2开始算到N*N
        temph=[A[i][j] for j in range(i)]#取长为i的第i行
        #H^dagger v=h,两边同时取共轭得到H^T v^*=h^*，其中h^*[j]又是A[i][j]
        for j in range(i):#解下三角
            l=temph[j]/H[j][j]
            for k in range(j+1,i):
                temph[k]=temph[k]-l*H[j][k]
        tempv=[temph[j]/H[j][j] for j in range(i)]#解得v^*
        tempv2=0
        for j in range(i):
            H[j][i]=tempv[j].conjugate()#扩展H
            tempv2+=abs(tempv[j])*abs(tempv[j])#计算tempv的模方
        if A[i][i]-tempv2<=0:
            raise ValueError("Not a positive defined matrix")
        else:
            H[i][i]=pow(A[i][i]-tempv2,0.5)
    return H

def CholeskySolve(Ac,bc):#利用Cholesky分解求解线性方程组，Ac二维list，bc一维list
    A=Ac.copy()
    b=bc.copy()
    N=len(A)
    H=Cholesky(A)
    Hdagger=[[H[k][m].conjugate() for k in range(N)] for m in range(N)]#得到H的共轭转置
    return SolveUp(H,SolveLow(Hdagger,b))
