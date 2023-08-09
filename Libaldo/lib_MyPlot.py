import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from  lib_MyEq import *
from  lib_MyEqEq import *
from lib_MyIntegral import *
from matplotlib.axis import Axis
from lib_tools import *
import mplcursors


def vecspace(x1,x2,x3):

    tt=np.linspace(float(x1),float(x2),x3)
    return tt 
def dataplotXY(obj,x1,x2,x3=200):
    obj=obj2str(obj)
    vecx=vecspace(x1,x2,x3)
    xx,yy=dataplot(vecx,obj)
    xx1,xx2=min(xx),max(xx)
    vecx=getinterx(obj,var=x)
    if len(vecx)>0:
        ee=obj2MyEq(obj)
        xmin,xmax=min(vecx),max(vecx)
        if xmin<xx1:
            y1=float(ee(x=xmin))
            xx=[xmin]+xx
            yy=[y1]+yy
        if xmax>xx2:
            y2=float(ee(x=xmax))
            xx=xx+[xmax]
            yy=yy+[y2]
        
    if len(xx)<100:
        return dataplotXY(obj,xx1,xx2,x3=200)
    else:    
        return xx,yy,xx1,xx2
def dataplot(vecx,sfunc):
    X=[]
    Y=[]
    ee=MyEq(parse_expr(sfunc),'ee',var=x,kshow=False)
    for i in vecx:
        kres=ee(x=i)
        if not 'I' in str(kres):
            Y.append(ee(x=i))
            X.append(i)
    return X,Y
def obj2str(obj):
    if type(obj)==str:
        return obj
    elif type(obj)==MyEq:
        return str(obj.ksym)
    elif type(obj)==MyEqEq:
        return str(simplify(expand(obj.L-obj.R)))
    else:
        return str(obj)
def obj2func(obj):
    if type(obj)==str:
        return parse_expr(obj)
    elif type(obj)==MyEq:
        return  obj.ksym 
    elif type(obj)==MyEqEq:
        return  simplify(expand(obj.L-obj.R)) 
    else:
        return obj
def obj2MyEq(obj,var=x):
    obj=obj2func(obj)
    ee=MyEq(obj,'ee',var=var,kshow=False)
    return ee

def getinterx(obj,var=x):
    ee=obj2MyEq(obj,var=var)
    vx=ee.roots(kshow=False)
    if 'I' in str(vx) and ee.degree()==3:
        vec3=ee.coef_list()
        vx= list(solve3(*vec3))
    vx2=[]
    for i in vx:
        if not 'I' in str(i):
            vx2.append(i)
            
    return vx2

def getinflexpoint(obj,var=x):
    ee=obj2MyEq(obj,var=var)
    dee=ee.diff('dee',kshow=False)
    de2=dee.diff('dee',kshow=False)
    vx=de2.roots(kshow=False)
    if 'I' in str(vx) and de2.degree()==3:
        vec3=ee.coef_list()
        vx= list(solve3(*vec3))
    vx2=[]
    for i in vx:
        if not 'I' in str(i):
            vx2.append(i) 
            

def getintery(obj,var=x):
    ee=obj2MyEq(obj,var=var)
    try:
        yy=ee(x=0)
        if not 'I' in str(yy):
            return yy
        return 0
    except:
        return 0


def get_maxminXY(obj,var=x):
    vecy=[]
    vecx=[]
    ee=obj2MyEq(obj,var=var)
    dee=ee.diff('dee',kshow=False)
    vx=dee.roots(kshow=False)
    if 'I' in str(vx) and ee.degree()==3:
        vec3=ee.coef_list()
        vx= list(solve3(*vec3))
    if type=='float' :
        kres=[float(i) for i in vx]
        vx=kres
    for i in vx:
        yy=ee(x=i)
        if not 'I' in str(yy):
            vecx.append(float(i))
            vecy.append(float(yy))
            
    return vecx,vecy

def get_diferfunc(obj1,obj2):
    return obj2func(obj1)-obj2func(obj2) 
    
def get_interfunc(*args):
    vecx,vecy,vecname=[],[],[]
    qq=len(args)
    vecee=[]
    for i in range(0,qq-1):
        for j in range(i+1,qq):
            kres=get_diferfunc(args[i],args[j])
            kres=expand(kres)
            kres=simplify(kres)
             
            ee=obj2MyEq(kres)
            if ee.degree()==3:
                vec3=ee.coef_list()
                xx= list(solve3(*vec3)) 
            else:
                xx=ee.roots(kshow=False)
                
            
            for k in xx:
                if not 'I' in str(k):
                    ee2=obj2MyEq(args[i])
                    if not k in vecx:
                        yy=ee2(x=k)
                        if not 'I' in str(yy):
                            vecx.append(k)
                            vecy.append(yy)
                            vecname.append(str(i)+str(j))
    return vecx,vecy,vecname         

def setlimitg(xx1,yy1,veclimit):
    if veclimit==[]:
        x1,x2,y1,y2=xx1,xx1,yy1,yy1
    else:    
        x1=veclimit[0]
        x2=veclimit[1]
        y1=veclimit[2]
        y2=veclimit[3]
    
    if xx1<x1:
        x1=flot2int(float(xx1))
    if xx1>x2:
        x2=flot2int(float(xx1))
    if yy1<y1:
        y1=flot2int(float(yy1))
    if yy1>y2:
        y2=flot2int(float(yy1))
        
    return [x1,x2,y1,y2]
def flot2int(expr):
    sexp=str(expr)
    p=sexp.find('.')
    if p==-1:
        return int(expr)
    else:
        expr=round(float(expr),2)
        pd=expr-int(expr)
        if pd==0:
            return int(expr)
        else:
            return expr    

def MyPlot(*args):
    vecf=[]
    x1=''
    x2=''
    XX1=''
    XX2=''
    for i in args:
        if type(i)==MyEq:
            vecf.append(i)
        elif type(i)==MyIntg:
            vecf.append(MyEq(i.ksym,'ee',kshow=False))
        else:
            if x1=='':
                x1=i
                XX1=i    
            elif x2=='':
                x2=i
                XX2=i
            else:
                x3=i
             
                
            
         
    veccolor=['blue','green','cyan','magenta']
    X,Y,L,kC,name=[],[],[],[],[]
    Ox=False
    Oy=False
    veccolor=['blue','green','cyan','magenta']
    cc=0
    for obj in vecf:
        vx=getinterx(obj)
        for i in vx:
            X.append(i)
            Y.append(0)
            L.append('Axis(X)')
            kC.append('black')
            name.append(obj.name)
            Ox=True
        i= getintery(obj)
        if i!=0:
            X.append(0)
            Y.append(i)
            L.append('Axis(Y)')
            kC.append('black')
            name.append(obj.name)
            
        vx,vy=get_maxminXY(obj)
        if len(vx)>0:
            ee=obj2MyEq(obj) 
            for i2,j2 in zip(vx,vy):
                X.append(i2)
                Y.append(j2)
                ksim='Max'
                y2=float(ee(i2+0.01))
                if y2>j2:
                    ksim='Min'
                L.append(ksim)
                kC.append(veccolor[cc])
                name.append(obj.name)
        
        df=obj.diff('df',kshow=False)
        df2=df.diff('df2',kshow=False)
        vecx=df2.solve('x','all',kshow=False)
        vecy=[float(obj(kk)) for kk in vecx]
 
        for ii,jj in zip(vecx,vecy):

            X.append(ii)
            Y.append(jj)
            L.append('Inflex')
            kC.append('red')
            name.append(obj.name)
        
        cc+=1
 
    cc=0
    for obj in vecf:
        vx=getinterx(obj)
        for i in vx:
            X.append(i)
            Y.append(0)
            L.append('Axis(X)')
            kC.append('black')
            name.append(obj.name)
            Ox=True
        i= getintery(obj)
        if i!=0:
            X.append(0)
            Y.append(i)
            L.append('Axis(Y)')
            kC.append('black')
            name.append(obj.name)
            
        vx,vy=get_maxminXY(obj)
        if len(vx)>0:
            ee=obj2MyEq(obj) 
            for i2,j2 in zip(vx,vy):
                X.append(i2)
                Y.append(j2)
                ksim='Max'
                y2=float(ee(i2+0.01))
                if y2>j2:
                    ksim='Min'
                L.append(ksim)
                kC.append(veccolor[cc])
                name.append(obj.name)
        df=obj.diff('df',kshow=False)
        df2=df.diff('df2',kshow=False)
        vecx=df2.solve('x','all',kshow=False)
        vecy=[float(obj(kk)) for kk in vecx]
        for ii,jj in zip(vecx,vecy):

            X.append(ii)
            Y.append(jj)
            L.append('Inflex')
            kC.append('red')
            name.append(obj.name)
        
        cc+=1

    
    vx,vy,vn=get_interfunc(*vecf)
    for i,j,k in zip(vx,vy,vn):
        X.append(i)
        Y.append(j)
    
        L.append('Intersc')
        kC.append('red')
        name.append(k)
    datap=[]
    for i in range(len(X)):
        datap.append((name[i],L[i],X[i],Y[i]))

    xmax=max(X)
    xmin=min(X)
    ymax=max(Y)
    ymin=min(Y)
    fcx5=(xmax-xmin)/5
    fcx10=(xmax-xmin)/10
    fcy5=(ymax-ymin)/5
    fcy10=(ymax-ymin)/10
    xmax2=xmax+fcx10
    xmin2=xmin-fcx10
    ymax2=ymax+fcy10
    ymin2=ymin-fcy10
    fig, bx = plt.subplots(figsize=(7, 3))
    veccolor=['blue','green','cyan','magenta']
    for i in range(len(X)):
        x=X[i]
        y=Y[i]
        Lab=L[i]
        kc=kC[i]
        kmarker="o"
        kmarkersize=5
        if Lab=='Min':
            kmarker="v"
        if Lab=='Max':
            kmarker="^"
        if Lab=='Intersc':
            kmarker="*"
        if Lab=='Axis(Y)':
            kmarker="x"
        if Lab=='Axis(X)':
            kmarker="x"
        if Lab=='Inflex':
            kmarker="_"
            kc="black"
            kmarkersize=25            
        bx.plot(x,y, marker=kmarker, markersize=kmarkersize, markeredgecolor=kc, markerfacecolor=kc,label=Lab)  
    mplcursors.cursor(bx)
     
     
    if x1!='':
        xmin=x1
    if x2!='':
        xmax=x2        
    
    
        
    bx.plot(xmin,0, marker=".", markersize=0.1, markeredgecolor='white')
    bx.plot(xmax,0, marker=".", markersize=0.1, markeredgecolor='white')
    bx.plot(0,ymin, marker=".", markersize=0.1, markeredgecolor='white')
    bx.plot(0,ymax, marker=".", markersize=0.1, markeredgecolor='white')
    if x1=='':
        x1=xmin2
        x2=xmax2
    y1=ymin2
    y2=ymax2
    cc=0
    for ff in vecf:
     
        xx,yy,xx1,xx2=dataplotXY(ff,x1,x2,x3=100)     
        bx.plot(xx,yy,linewidth=1.5,label="$"+latex(obj2func(ff))+"$",zorder=-1.0,color=veccolor[cc]) # plot the function     
        cc+=1
    
    if XX1!='':
        x1,x2=XX1,XX2
    else:    
        x1,x2=bx.get_xlim()
    y1,y2=bx.get_ylim()
    if x1<=0 and x2>=0:
        bx.plot([0,0],[y1,y2],linestyle='dashdot',color='red',linewidth=1.8,zorder=-1.0)
        Oy=True
    if y1<=0 and y2>=0:
        bx.plot([x1,x2],[0,0],linestyle='dashdot',color='red',linewidth=1.8,zorder=-1.0)
        Ox=True
    if Ox and Oy:
        bx.plot(0, 0, marker="o", markersize=3, markeredgecolor='black', markerfacecolor='black',label='(0,0)') # plot inter X
        
    plt.grid()
    plt.show()
