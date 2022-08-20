import math
import numpy as np
import time
from scipy import integrate

zero=0

def SumFormResult(N):
    q2=0.5
    q=pow(q2,0.5)
    pi=math.pi
    axis=0#x正半轴上求和值（不含原点）
    plane=0#xy平面上求和值（不含x,y轴）
    volume=0#第一象限内求和的值（不含xy,xz,yz平面）00
    fc=-np.exp(q2)/q2+integrate.quad(lambda x:(np.exp(q2*x)-1)*pow(pi/x,1.5),zero,1)[0]-2*pow(pi,1.5)#与求和无关的常数项
    #print(fc)
    for i in range(1,N):#算一个象限内，不含轴上的点
        x2=i*i#减少乘法次数
        axis+=np.exp(q2-x2)/(x2-q2)+integrate.quad(lambda x:np.exp(q2*x-pi*pi*x2/x)*pow(pi/x,1.5),zero,1)[0]
        for j in range(1,N):
            s2=x2+j*j#减少乘法次数
            plane+=np.exp(q2-s2)/(s2-q2)+integrate.quad(lambda x:np.exp(q2*x-pi*pi*s2/x)*pow(pi/x,1.5),zero,1)[0]
            for k in range(1,N):
                r2=s2+k*k
                volume+=np.exp(q2-r2)/(r2-q2)+integrate.quad(lambda x:np.exp(q2*x-pi*pi*r2/x)*pow(pi/x,1.5),zero,1)[0]
                
    return 6*axis+12*plane+8*volume+fc#返回所有值

def homework4c():
    for i in range(1,11):
        print("When N= ",i," , the result is ",SumFormResult(i))

homework4c()


