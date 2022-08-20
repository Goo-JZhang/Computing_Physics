# -*-coding:utf-8 -*-
import Interpolation as itp
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

pd.set_option("display.max_column",None)
pd.set_option("max_colwidth",100)
pd.set_option("precision",16)

def homework2_1a():
    x1=np.linspace(0,math.pi/2,5)
    y1=np.sin(x1)
    L=itp.LagrangeInterpolation(list(x1),list(y1))
    L=itp.Neville(list(x1),list(y1))[1][-1][0]
    #print(L)
    x2=np.linspace(0,math.pi/2)
    y2=L.cal(x2)
    plt.figure(1)
    plt.grid()
    plt.plot(x2,y2,label="Polynomial")
    plt.plot(x2,np.sin(x2),color='red',linestyle='--',label="sin(x)")
    plt.legend()
    plt.figure(2)
    plt.grid()
    plt.plot(x2,y2-np.sin(x2),label="Error=Polynomial-sin(x)")
    print("Max Error =",max(abs((y2)-np.sin(x2))))
    plt.legend()
    plt.show()

def homework2_1b():
    x1=np.array([0,0.1,1.4,1.5,math.pi/2])
    y1=np.sin(x1)
    L=itp.LagrangeInterpolation(list(x1),list(y1))
    #print(L)
    x2=np.linspace(0,math.pi/2)
    y2=L.cal(x2)
    plt.figure(1)
    plt.grid()
    plt.plot(x2,y2,label="Polynomial")
    plt.plot(x2,np.sin(x2),color='red',linestyle='--',label="sin(x)")
    plt.legend()
    plt.figure(2)
    plt.grid()
    plt.plot(x2,y2-np.sin(x2),label="Error=Polynomial-sin(x)")
    print("Max Error =",max(abs((y2)-np.sin(x2))))
    plt.legend()
    plt.show()

def homework2_1c():
    x1=np.linspace(0,math.pi/2,10)
    y1=np.sin(x1)
    x2=np.linspace(0,math.pi/2,14)
    y2=np.sin(x2)
    #L1=itp.LagrangeInterpolation(list(x1),list(y1))
    #L2=itp.LagrangeInterpolation(list(x2),list(y2))
    #L1=itp.NewtonInterpolation(list(x1),list(y1))[-1]
    #L2=itp.NewtonInterpolation(list(x2),list(y2))[-1]
    L1=itp.Neville(list(x1),list(y1))[1][-1][0]
    L2=itp.Neville(list(x2),list(y2))[1][-1][0]
    #print(L2)
    #print(L)
    x3=np.linspace(0,math.pi/2,50)
    y4=L2.cal(x3)
    y3=L1.cal(x3)
    plt.figure(1)
    plt.grid()
    plt.plot(x3,y3-np.sin(x3),label="n=10,Error=Polynomial-sin(x)")
    plt.legend()
    plt.figure(2)
    plt.grid()
    plt.plot(x3,y4-np.sin(x3),label="n=14,Error=Polynomial-sin(x)")
    print("n=10, Max Error =",max(abs((y3)-np.sin(x3))))
    print("n=14, Max Error =",max(abs((y4)-np.sin(x3))))
    plt.legend()
    plt.show()

def homework2_2a():
    x1=np.linspace(-1,1,21)
    x2=np.linspace(-1,1,41)
    x3=np.linspace(-1,1,200)
    f=lambda y:1/(1+25*y*y)
    P=itp.Neville(list(x1),list(f(x1)))[1][-1][0]
    y1=f(x2)
    y11=f(x3)
    y2=P.cal(x2)
    y3=P.cal(x3)
    data={'x':x2,'f(x)':y1,'P(x)':y2,"Error":abs(y2-y1)}
    print(pd.DataFrame(data))
    plt.figure(1)
    plt.grid()
    plt.plot(x3,y11,label="f(x)")
    plt.plot(x3,y3,color='red',linestyle='--',label="P(x)")
    plt.plot(x2,y2,color='green',linestyle=':',label="P(x),n=41")
    plt.legend()
    plt.figure(2)
    plt.grid()
    plt.plot(x3,y3-y11,label="Error=P(x)-f(x)")
    plt.plot(x2,y2-y1,color='red',linestyle=':',label="Error=P(x)-f(x),n=41")
    plt.legend()
    plt.show()

def homework2_2b():
    x2=np.cos(np.linspace(0,math.pi,41))
    x3=np.linspace(-1,1,200)
    f=lambda y:1/(1+25*y*y)
    P=itp.ChebyshevInterpolation(20,f)
    y1=f(x2)
    y11=f(x3)
    y2=P.cal(x2)
    y3=P.cal(x3)
    data={'x':x2,'f(x)':y1,'C(x)':y2,"Error":abs(y2-y1)}
    print(pd.DataFrame(data))
    plt.figure(1)
    plt.grid()
    plt.plot(x3,y11,label="f(x)")
    plt.plot(x3,y3,color='red',linestyle='--',label="C(x)")
    plt.plot(x2,y2,color='green',linestyle=':',label="C(x),n=41")
    plt.legend()
    plt.figure(2)
    plt.grid()
    plt.plot(x3,y3-y11,label="Error=C(x)-f(x)")
    plt.plot(x2,y2-y1,color='red',linestyle=':',label="Error=C(x)-f(x),n=41")
    plt.legend()
    plt.show()

def homework2_2c():
    x1=np.linspace(-1,1,21)
    x2=np.linspace(-1,1,41)
    x3=np.linspace(-1,1,200)
    f=lambda y:1/(1+25*y**2)
    #P=itp.Spine3Interpolation(list(x1),list(f(x1)))
    P=itp.Spline3Interpolation(list(x1),list(f(x1)),2)
    #P=itp.Spine3Interpolation(list(x1),list(f(x1)),3,[25/338,-25/338])
    y1=f(x2)
    y11=f(x3)
    y2=np.array(itp.Spline3Cal(list(x1),P,list(x2)))
    y3=np.array(itp.Spline3Cal(list(x1),P,list(x3)))
    data={'x':x2,'f(x)':y1,'S(x)':y2,"Error":abs(y2-y1)}
    print(pd.DataFrame(data))
    plt.figure(1)
    plt.grid()
    plt.plot(x3,y11,label="f(x)")
    plt.plot(x3,y3,color='red',linestyle='--',label="S(x)")
    plt.plot(x2,y2,color='green',linestyle=':',label="S(x),n=41")
    plt.legend()
    plt.figure(2)
    plt.grid()
    plt.plot(x3,y3-y11,label="Error=S(x)-f(x)")
    plt.plot(x2,y2-y1,color='red',linestyle=':',label="Error=S(x)-f(x),n=41")
    plt.legend()
    plt.show()

def homework2_3a():
    phi=np.linspace(0,2*math.pi,9)
    x=(1-np.cos(phi))*np.cos(phi)
    y=(1-np.cos(phi))*np.sin(phi)
    data={'x':x,'y':y}
    print(pd.DataFrame(data))

def homework2_3b():
    phi=np.linspace(0,2*math.pi,9)
    x=(1-np.cos(phi))*np.cos(phi)
    y=(1-np.cos(phi))*np.sin(phi)
    t=list(range(9))
    S1=itp.Spline3Interpolation(t,list(x),2)
    S2=itp.Spline3Interpolation(t,list(y),2)
    data={'S(X;t)':[i.str(6) for i in S1],'S(Y;t)':[i.str(6) for i in S2]}
    print(pd.DataFrame(data))

def homework2_3c():
    phi=np.linspace(0,2*math.pi,9)
    x=(1-np.cos(phi))*np.cos(phi)
    y=(1-np.cos(phi))*np.sin(phi)
    t=list(range(9))
    #S1=itp.Spine3Interpolation(t,list(x))
    #S2=itp.Spine3Interpolation(t,list(y))
    #S1=itp.Spine3Interpolation(t,list(x),2)
    #S2=itp.Spine3Interpolation(t,list(y),2)
    S1=itp.Spline3Interpolation(t,list(x),3)
    S2=itp.Spline3Interpolation(t,list(y),3)
    theta=np.linspace(0,2*math.pi,200)
    n=np.linspace(0,8,200)
    x1=(1-np.cos(theta))*np.cos(theta)
    y1=(1-np.cos(theta))*np.sin(theta)
    x2=np.array(itp.Spine3Cal(t,S1,list(n)))
    y2=np.array(itp.Spine3Cal(t,S2,list(n)))
    plt.figure(1)
    plt.grid()
    plt.scatter(x,y,color="green",label="Scatter")
    plt.plot(x1,y1,label="Standard")
    plt.plot(x2,y2,color="red",linestyle="--",label="Spline")
    plt.legend()
    plt.figure(2)
    plt.grid()
    plt.plot(x2-x1,y2-y1)
    '''
    plt.figure(3)
    plt.grid()
    plt.plot(n,x2)
    plt.plot(n,x1,color="red")
    plt.figure(4)
    plt.grid()
    plt.plot(n,x2-x1)
    plt.figure(5)
    plt.grid()
    plt.plot(n,y2)
    plt.plot(n,y1,color="red")
    plt.figure(4)
    plt.grid()
    plt.plot(n,y2-y1)
    '''
    plt.show()

def homework2_4a():
    f=lambda x: np.exp(-x*x)*np.cos(x)
    x=[0.1**i for i in range(8)]
    data1={'步长':x}
    data2={"N":range(1,13)}
    for scale in [1,2,5,10]:
        data1[str([-scale,scale])]=[itp.TrapezoidalIntegral(-scale,scale,f,int(scale/i)) for i in x]
        data2[str([-scale,scale])]=[itp.ExtraPolationIntegral(-scale,scale,f,i) for i in range(1,13)]
    print(pd.DataFrame(data1).round(16))
    print(pd.DataFrame(data2).round(16))

def homework2_4b():
    f=lambda y: np.cos(y)
    df = pd.read_excel (r'GaussionPoint.xlsx', sheet_name='Sheet1')
    name=list(df.columns)
    for i in range(4):
        x=df[name[2*i]].tolist()[1:5*i+6]
        w=df[name[2*i+1]].tolist()[1:5*i+6]
        print("N={}, result is {}".format((i+1)*5,itp.GaussianInt(x,w,f)))

homework2_4b()
