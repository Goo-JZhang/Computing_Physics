import LinearEqs
import sys

def HilbertMatrixGenerator(N):#生成Hilbert矩阵
    H=[[1/(i+j+1) for i in range(N)] for j in range(N)]
    return H

def homework():
    data=open("HilbertSolution.txt","w",encoding="utf-8")#输出结果到文件
    for i in range(1,11):
        b=[1 for j in range(i)]
        H=HilbertMatrixGenerator(i)#生成矩阵Hilbert矩阵
        print("N=",i,file=data)
        print("Solution of Cholesky Decomposition is",LinearEqs.CholeskySolve(H,b),file=data)
        print("Solution of Gaussian Elimination is",LinearEqs.GaussianElimination(H,b),"\n",file=data)
        
homework()
###以下为测试代码，忽略
#print(LinearEqs.CholeskySolve(HilbertMatrixGenerator(3),[1 for j in range(3)]))
#print(LinearEqs.GaussianElimination(HilbertMatrixGenerator(3),[1 for j in range(3)]))