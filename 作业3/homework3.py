import numpy as np
import matplotlib.pyplot as plt
import math

from pandas.core.frame import DataFrame
import Eigen as eg
import time
import pandas as pd

pd.set_option("display.max_column",None)
pd.set_option("max_colwidth",100)
pd.options.display.float_format='{:.6e}'.format

#数据文档路径
P='2pt_etac_Nconf50'

def homework3_1_d():
    rd = np.random.RandomState(0)
    n=20
    N=6
    A = np.zeros((n,N,N))
    for i in range(n):
        A[i] = rd.uniform(-1,1,(N,N))
    A=[np.matrix(i) for i in A]
    '''
    print('A[0]=')
    print(A[0])
    print("HouseholderQRD result:")
    print("Q=")
    print(eg.HouseholderQRD(A[0])[0])
    print("R=")
    print(eg.HouseholderQRD(A[0])[1])
    print("GivensQRD result:")
    print("Q=")
    print(eg.GivensQRD(A[0])[0])
    print("R=")
    print(eg.GivensQRD(A[0])[1])
    '''
    Ht=[]#Householder用时
    Gt=[]#Givens用时
    #Hstart=time.time()
    #for m in range(9):
    for i in A:
        Hstart=time.time()
        for j in range(10):
            eg.HouseholderQRD(i)
        Hend=time.time()
        Ht.append((Hend-Hstart)/10)
    #print("Householder time used:",Hstart)
    #Gstart=time.time()
    for i in A:
        Gstart=time.time()
        for j in range(10):
            eg.GivensQRD(i)
        Gend=time.time()
        Gt.append((Gend-Gstart)/10)
    #print("Givens time used:",Gend-Gstart)
    data=pd.DataFrame({'Householder Time':Ht,'Givens Time':Gt})
    print(data)
    print("First A:",A[0])
    print('Householder QRD of First A:\n')
    QH,RH=eg.HouseholderQRD(A[0])
    print('QH:\n',QH)
    print('RH:\n',RH)
    print('Givens QRD of First A:\n')
    QG,RG=eg.GivensQRD(A[0])
    print('QH:\n',QG)
    print('RH:\n',RG)
    #print(QH*QH.T,QG*QG.T,QH*RH-A[0],QG*RG-A[0])

def homework3_2_b():
    N=10
    A=np.zeros((N,N))
    for i in range(N):
        A[i][(i+N-1)%N]=-1
        A[i][i]=2
        A[i][(i+1)%N]=-1
    print(eg.SolveMaxEigen(A))

def homework3_3():
    N=50
    Nt=48
    datas=[pd.read_csv(r"{}\\2pt_etac_conf_1{:0>2d}0.txt".format(P,i),names=['time','Re','Im'],sep='\t') for i in range(N)]
    #对称化
    #print(datas[0])
    for data in datas:
        data.loc[Nt]=['t=48',data.iloc[0,1],data.iloc[0,2]]
        for j in range(1,Nt//2+1):
            data.iloc[j,1:3]=0.5*(data.iloc[j,1:3]+data.iloc[Nt-j,1:3])
            data.iloc[Nt-j,1:3]=data.iloc[j,1:3]
    #print(datas[0])
    #第一问
    #求平均
    Cbar=pd.DataFrame({'time':range(49),'Re':np.zeros(49),'dRe':np.zeros(49),'Im':np.zeros(49),'dIm':np.zeros(49)})
    for data in datas:
        Cbar['Re']+=data['Re']
        Cbar['dRe']+=data['Re']*data['Re']
        Cbar['Im']+=data['Im']
        Cbar['dIm']+=data['Im']*data['Im']
    Cbar['Re']=Cbar['Re']/N
    Cbar['Im']=Cbar['Im']/N
    Cbar['dRe']=np.sqrt((Cbar['dRe']/N-Cbar['Re']*Cbar['Re'])/(N-1))
    Cbar['dIm']=np.sqrt((Cbar['dIm']/N-Cbar['Im']*Cbar['Im'])/(N-1))
    dC=pd.DataFrame({'time':range(49),'dRe/Re':Cbar['dRe']/Cbar['Re'],'dIm/Im':Cbar['dIm']/Cbar['Im'],'dMode/Mode':np.sqrt((Cbar['dRe']*Cbar['dRe']+Cbar['dIm']*Cbar['dIm'])/(Cbar['Re']*Cbar['Re']+Cbar['Im']*Cbar['Im']))})
    print(dC)
    DeltaC=lambda x: dC['dRe/Re'][x]
    #第一问作图
    plt.figure('(a)')
    plt.grid()
    plt.xlabel('t')
    plt.ylabel('dC/C*100%')
    plt.plot(dC['time'],100*dC['dMode/Mode'],label='d|C|/|C|',color="red",linestyle='--')
    plt.plot(dC['time'],100*dC['dRe/Re'],label='dRe(C)/Re(C)')
    plt.legend()
    R=np.array([Cbar['Re'][i]/Cbar['Re'][i+1] for i in range(Nt//2)])
    #print(R)
    Meff=pd.DataFrame({'time':range(Nt//2),'meff':np.log(R)})
    #print(Meff)
    #Jackknife
    JKC=pd.DataFrame({'time':range(Nt+1)})
    for i in range(N):#剔除组态
        JKC['C{:0>2d}'.format(i)]=(N*Cbar['Re']-datas[i]['Re'])/(N-1)
    JKM=pd.DataFrame({'time':range(Nt//2)})
    for i in range(N):#去掉各组态的meff表
        JKM['meff{0:2>d}'.format(i)]=[np.log(JKC['C{:0>2d}'.format(i)][j]/JKC['C{:0>2d}'.format(i)][j+1]) for j in range(Nt//2)]
    JKM['meff']=np.zeros(Nt//2)#均值
    for i in range(N):
        JKM['meff']+=JKM['meff{0:2>d}'.format(i)]
    JKM['meff']=JKM['meff']/N
    Meff['dmeff']=np.zeros(Nt//2)#计算误差
    for i in range(N):
        Meff['dmeff']+=(JKM['meff{0:2>d}'.format(i)]-JKM['meff'])**2
    Meff['dmeff']=np.sqrt((N-1)*Meff['dmeff']/N)
    print(Meff)

    #meff修正
    FineMeff=np.array([np.log((1+R[i]*R[i-1]+pow((1+R[i]*R[i-1])**2-4*R[i]**2,0.5))/(2*R[i])) for i in range(2,Nt//2-1)])
    #R[0]过小导致方程无实根
    #print(FineMeff)

    #拟合
    #确定时间片起点终点
    Cfit=lambda A,m,t:A*(np.exp(-m*t)+np.exp(m*(t-Nt)))
    Chimin=float('inf')
    for tmin in range(Nt//2):
        for tmax in range(tmin+3,Nt//2):
            Chi2=lambda A,m:np.sum(((Cbar['Re'][tmin:tmax+1]-Cfit(A,m,np.arange(tmin,tmax+1)))/(Cbar['dRe'][tmin:tmax+1]))**2)/(tmax-tmin-1)
            #初始范围
            m1=Meff['meff'][tmin]
            m2=Meff['meff'][tmax]
            A1=Cbar['Re'][0]/3
            A2=Cbar['Re'][0]*3
            ags=[np.array([A1,m1]),np.array([A2,m1]),np.array([A1,m2])]
            point,er,mvalue=eg.SimplexMinimize(Chi2,ags)
            if mvalue<Chimin:
                Chimin=mvalue
                minpoint=point
                miner=er
                Tmin=tmin
                Tmax=tmax
    print("scale: tmin={}, tmax={}".format(Tmin,Tmax))
    print("Minimize point and value of Chi2 are A={} m={} and Chi2min={}".format(minpoint[0],minpoint[1],Chimin))
    #errm=miner[1]+np.max(Meff['dmeff'][Tmin:Tmax+1])#误差暴力估计，直接把极值偏差与拟合范围内meff的最大偏差相加
    #print("Error of m is {}".format(errm))
    Chi2=lambda A,m:np.sum(((Cbar['Re'][Tmin:Tmax+1]-Cfit(A,m,np.arange(Tmin,Tmax+1)))/(Cbar['dRe'][Tmin:Tmax+1]))**2)/(Tmax-Tmin-1)
    #求m误差,自由度d=3，3sigma置信度
    chi2=14.2
    A=minpoint[0]
    A0=minpoint[0]
    m=minpoint[1]
    m0=minpoint[1]
    #1.0007和0.999的参数不能乱改
    #正向
    for i in range(50):
        #print(Chi2(A,1.0007*m)-chi2)
        m=eg.SeekRoots(lambda x:Chi2(A,x)-chi2,m,1.0007*m,ep=1e-4)
        A=eg.SeekRoots(lambda x:Chi2(x,m)-chi2,A,1.0007*A,ep=1e-4)
    errm=miner[1]+m-m0
    #反向
    m=minpoint[1]
    A=minpoint[0]
    for i in range(50):
        #print(Chi2(A,1.0007*m)-chi2)
        m=eg.SeekRoots(lambda x:Chi2(A,x)-chi2,m,0.999*m,ep=1e-4)
        A=eg.SeekRoots(lambda x:Chi2(x,m)-chi2,A,0.999*A,ep=1e-4)
    errm=max(miner[1]+m0-m,errm)
    print("Error of m is {}".format(errm))
    '''
    figb=plt.figure('(b)')
    axb1=figb.add_subplot(111)
    plt.grid()
    axb1.set_xlabel('t')
    axb1.set_ylim(0,1.65)
    axb1.set_ylabel('meff')
    axb1.scatter(Meff['time'],Meff['meff'],label='meff')
    axb1.plot([Tmin,Tmax],[point[1],point[1]],color='red',linestyle='--',label='mfit')
    axb1.legend(loc=2)
    axb2=axb1.twinx()
    axb2.set_ylabel('dmeff')
    axb2.set_ylim(0,1e-3)
    axb2.scatter(Meff['time'],Meff['dmeff'],label='dmeff',color='red',marker='x')
    axb2.scatter(Meff['time'][Tmin:Tmax+1],abs(minpoint[1]-Meff['meff'][Tmin:Tmax+1]),label='|mfit-meff|',color='green',marker='v')
    axb2.legend(loc=1)
    '''
    plt.figure('(b)')
    plt.grid()
    plt.xlabel('t')
    plt.ylabel('meff')
    plt.errorbar(Meff['time'],Meff['meff'],yerr=Meff['dmeff'],fmt='o',ms=3,ecolor='green',label='meff')
    plt.plot([0,24],[minpoint[1],minpoint[1]],linestyle='-',label='mfit')
    plt.plot([Tmin,Tmax],[minpoint[1]-errm,point[1]-errm],color='red',linestyle='--')
    plt.plot([Tmin,Tmax],[minpoint[1]+errm,point[1]+errm],color='red',linestyle='--')
    plt.axvline(x=Tmin,ymin=minpoint[1]-0.5,ymax=minpoint[1]-0.3,color='dodgerblue',label='tmin')
    plt.axvline(x=Tmax,ymin=minpoint[1]-0.5,ymax=minpoint[1]-0.3,color='aqua',label='tmax')
    xt=np.linspace(0,23,200)
    plt.plot(xt,np.log(Cfit(minpoint[0],minpoint[1],xt)/Cfit(minpoint[0],minpoint[1],xt+1)),color='brown',linestyle='-',label='mfit(t)')
    plt.legend()

    plt.figure('(c)')
    plt.grid()
    plt.xlabel('t')
    plt.ylabel('meff')
    plt.errorbar(Meff['time'][Tmin-1:Tmax+2],Meff['meff'][Tmin-1:Tmax+2],yerr=Meff['dmeff'][Tmin-1:Tmax+2],fmt='o',ecolor='green',label='meff')
    #plt.errorbar(Meff['time'][Tmin:Tmax+1],[point[1] for i in range(Tmax+1-Tmin)],yerr=Meff['dmeff'][Tmin:Tmax+1],fmt='-',ecolor='aquamarine',label='mfit')
    plt.plot([Tmin,Tmax],[minpoint[1],minpoint[1]],linestyle='-',label='mfit')
    plt.plot([Tmin,Tmax],[minpoint[1]-errm,minpoint[1]-errm],color='red',linestyle='--')
    plt.plot([Tmin,Tmax],[minpoint[1]+errm,minpoint[1]+errm],color='red',linestyle='--')
    xt=np.linspace(Tmin-1,Tmax+1)
    plt.plot(xt,np.log(Cfit(minpoint[0],minpoint[1],xt)/Cfit(minpoint[0],minpoint[1],xt+1)),color='brown',linestyle='-',label='mfit(t)')
    plt.legend()
    #print(minpoint[1]-Meff['meff'][Tmax])
    #print(tmin,tmax)
    #meff修正
    #FineMeff=np.array([0.5*np.log((np.exp(1)-R[i+1])*(np.exp(1)*R[i]-1)/((np.exp(1)-R[i])*(np.exp(1)*R[i+1]-1))) for i in range(Nt//2-1)])

    plt.figure('(d)')
    plt.grid()
    plt.xlabel('t')
    plt.ylabel('meff')
    plt.errorbar(Meff['time'],Meff['meff'],yerr=Meff['dmeff'],fmt='o',ms=3,ecolor='green',label='meff')
    plt.scatter(range(2,Nt//2-1),FineMeff,color='slateblue',marker='v',label='fine meff')
    plt.plot([0,24],[minpoint[1],minpoint[1]],linestyle='-',label='mfit')
    xt=np.linspace(0,23,200)
    plt.plot(xt,np.log(Cfit(minpoint[0],minpoint[1],xt)/Cfit(minpoint[0],minpoint[1],xt+1)),color='brown',linestyle='-',label='mfit(t)')
    plt.legend()

    plt.show()
    #print(dC)


#homework3_1_d()
#homework3_2_b()
#homework3_3()



