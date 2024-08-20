import sys
sys.path.insert(0, 'Libaldo/')

from sympy import *
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import Image, display



from lib_Variables import *
from lib_Mathematica import *
from lib_Algorith import * 
from lib_Mathbasic import *


from  lib_MyEq import *
from  lib_MyEqEq import *

from lib_MyFunctions import *
from lib_MyDiff import *
from lib_MyIntegral import *
from lib_MyMatrix import * 
 
 
 
init_printing()

def gradienteDiff(obj,*args):
    '''
    if 
        F=MyEq(2*x*y+x-3*y*y,'F')
    input     
        V=gradienteDiff(F,'V',x,y)
    output:
        [2y+1,2x-6y]
 
    '''
    expr=obj
    if type(obj)==MyEq:
        expr=obj.ksym
        
    name=''
    vecvar=[]
    vecdiff=[]
    for data in args:
        if type(data)==str:
            name=data
        if type(data)==Symbol:
            vecvar.append(data)        
    for var in vecvar:
        vecdiff.append(diff(expr,var))
    if name!='':
        return MyVecx(name,vecdiff)
    else:
        return vecdiff 
        
dx,dy,dz,dt,dw,dv,du,da=symbols('d_x d_y d_z d_t d_w d_v d_u  d_alpha')
def diff_respect(obj,*args):
    VECV=[x,y,z,t,alpha]
    VECD=[dx,dy,dz,dt,da]
    ops=[]
    vecvv=[]
    for data in args[1::]:
        if type(data)==str:
            ops.append(data)
        else:
            vecvv.append(data)
        
    lvar=args[0]
    rvar=[data for data in vecvv]
    ldvar=VECD[VECV.index(lvar)]
    rdvar=[VECD[VECV.index(data)] for data in rvar]
    expr1=obj.e1.ksym
    '''
    if len(rvar)==1:
        dexpr1=diff(expr1,lvar)*(ldvar/rdvar[0])
    else:
    '''
    dexpr1=diff(expr1,lvar)*ldvar
    expr2=obj.e2.ksym
    '''
    if len(rvar)==1:
        dexpr2=diff(expr2,rvar[0])
    else:
    '''
    kres=0
    for vv,vd in zip(rvar,rdvar):
        kres=kres+diff(expr2,vv)*vd
    dexpr2=kres
    if 'update' in ops:
        obj.e1.ksym=dexpr1
        obj.e2.ksym=dexpr2
        obj.s()
    else:
        return MQ(dexpr1,dexpr2)        
        
def get_difvar(*args):
    VECV=[x,y,z,w,v,u,t,alpha]
    VECD=[dx,dy,dz,dw,dv,du,dt,da]
    kres=[VECD[VECV.index(data)] for data in args]
    return kres
def get_derivadaf(expr,*args):
    ops=[x,y,z,w,v,u,t,alpha]
     
    kres=0
    varv=[data for data in args]
    vard=get_difvar(*varv)
    for vv,dv in zip(varv,vard):
        kres+=diff(expr,vv)*dv
    return kres
def implicitDiff(obj,*args,alone=False):
    '''
    Q=Type: MyEqEq
    Q= f(x,y)=f(z,v)
        Q2=diff_implicit(Q,(x,y),z,v)
        Q2=diff_implicit(Q,(x),z)
        Q=MyEqEq, (v1,v2,v3).. variables in Left side to diferential, ,z,v,..  variables in R side                   
    '''
 
    lvar=args[0]
    expr1=obj.e1.ksym
     
    if expr1==0:
        dexpr1=0
    else:    
        if type(lvar)==tuple:
            if lvar==():
                dexpr1=0
            else:
                lvar=list(lvar)
                dexpr1=get_derivadaf(expr1,*lvar)
        else:
            if lvar==0:
                dexpr1=0
            else:    
                lvar=[lvar]
                dexpr1=get_derivadaf(expr1,*lvar)
     
    try:
        rvar=args[1::]
        expr2=obj.e2.ksym
        if expr2==0:
            dexpr2=0
        else:    
            dexpr2=get_derivadaf(expr2,*rvar)
         
    except:
        dexpr2=0
    if alone:
          
        dfac=get_difvar(rvar[0])
        dexpr1=cfrac(dexpr1,dfac[0])
        dexpr2=simplify(dexpr2/dfac[0])
    kres= MQ(dexpr1,dexpr2,kshow=False)       
    return kres
    
def diffS2(p1,p2):
    # input variables example V,t*gg
    # return symbols with name dV2/dt2
    sp1='d'+alphaname(p2)+'^{2}'
    sp2='d^{2}'+alphaname(p1) 
    name='\\frac{'+sp2+'}{'+sp1+'}'
    return symbols(name)


def diffS(p1,p2,op=''):
    # input variables example V,t*gg
    # return symbols with name dV/dt
    sp2='d'+alphaname(p1)
    sp1='d'+alphaname(p2)
    name='\\frac{'+sp2+'}{'+sp1+'}'
    dp1=symbols(name)
    if op==2:
        dp2=diffS2(p1,p2)
        return dp1,dp2
    else:
        return dp1     