from sympy import *
 
from IPython.display import Math  # ,display
from matplotlib.pyplot import ylabel, plot, show, xlabel, title
from libaldo_math2 import *
from libaldo_algorith import *
 
from lib_toolsEq import *
from lib_tools import * 
from mathbasic import *
from lib_MyEq import *
from lib_MyEqEq import *
from sympy.solvers.inequalities import reduce_rational_inequalities
class MyIneq3(MyEqEq):
    def __init__(self, *args, var=x, kshow=True,ktype='iq3'):

        self.var=var
         
        if type(args[0])==str:
            self.p1=parse_expr(args[0])
        else:
            self.p1=args[0]
        self.s1=args[1]
        self.px=args[2]
        self.s2=args[3]
         
        if type(args[4])==str:
            self.p2=parse_expr(args[4])
        else:
            self.p2=args[4]
        self.type=ktype 
        self.Iq1=MQ(self.p1,self.s1,self.px,var=var,kshow=False)
        self.Iq2=MQ(self.px,self.s2,self.p2,var=var,kshow=False)
        ineqsynbo=['>=','<=','>','<']
        self.s1n=nice_iq_symbol(self.s1)
        self.s2n=nice_iq_symbol(self.s2)
        if kshow:
            display(Math(latex(self.p1)+' '+self.s1n+' '+latex(self.px)+' '+self.s2n+latex(self.p2)))

    def getdata(self):
        p1=self.Iq1.L
        s1n=self.Iq1.symb
        px=self.Iq1.R
        s2n=self.Iq2.symb
        p2=self.Iq2.R
        return p1,s1n,px,s2n,p2
    def s(self):
        p1,s1n,px,s2n,p2=self.getdata()
        display(Math(latex(p1)+' '+s1n+' '+latex(px)+' '+s2n+latex(p2)))
    def getIqs(self):
        return self.Iq1,self.Iq2
    @property
    def L(self):
        Iq1,Iq2=self.getIqs()
        return Iq1.L
    @property
    def R(self):
        Iq1,Iq2=self.getIqs()
        return Iq2.R    
    def getdata(self):
        Iq1,Iq2=self.getIqs()
        p1=Iq1.L
        s1n=Iq1.symb
        px=Iq1.R
        s2n=Iq2.symb
        p2=Iq2.R
        return p1,s1n,px,s2n,p2
    def Add(self,kval):
        Iq1,Iq2=self.getIqs()
        Iq1.Add(kval,'noshow')
        Iq2.Add(kval,'noshow')
        self.s()
    def Substrac(self,kval):
        Iq1,Iq2=self.getIqs()
        Iq1.Substrac(kval,'noshow')
        Iq2.Substrac(kval,'noshow')
        self.s()
    def Mul(self,kval):
        Iq1,Iq2=self.getIqs()
        Iq1.Mul(kval,'noshow')
        Iq2.Mul(kval,'noshow')
        if signo(kval)==-1:
            s1=Iq1.psymb
            sn1=Iq1.symb
            s2=Iq2.psymb
            sn2=Iq2.symb
            p1=Iq1.L
            p2=Iq2.R
            Iq1.psymb=s2
            Iq1.symb=sn2
            Iq1.e1.ksym=p2
            Iq2.psymb=s1
            Iq2.symb=sn1
            Iq2.e2.ksym=p1
            
            
        self.s()
    def Div(self,kval):
        Iq1,Iq2=self.getIqs()
        Iq1.Div(kval,'noshow')
        Iq2.Div(kval,'noshow')
        self.s()
    def symplify(self):
        Iq1,Iq2=self.getIqs()
        Iq1.symplify(kshow=False)
        Iq2.symplify(kshow=False)
        self.s()
    def expand(self):
        Iq1,Iq2=self.getIqs()
        Iq1.expand(kshow=False)
        Iq2.expand(kshow=False)
        self.s()
    def apart(self):
        Iq1,Iq2=self.getIqs()
        Iq1.apart(kshow=False)
        Iq2.apart(kshow=False)
        self.s()     
    def factor(self):
        Iq1,Iq2=self.getIqs()
        Iq1.factor(kshow=False)
        Iq2.factor(kshow=False)
        self.s()

    def isolve(self,*args):
        Iq1,Iq2=self.getIqs()
        s1=get_str_eq(Iq1)
        s2=get_str_eq(Iq2)
        var=self.var
        pexpr="[["+s1+","+s2+"]]"
        kres=reduce_rational_inequalities(parse_expr(pexpr),var)
        if 'Iq' in args:
            p1,p2,p3,p4=getsimpleIq(kres)
            if "." in p1 or "/" in p1:
                p1=Rational(p1)
            elif p1=='oo':
                p1=oo
            elif p1=='-oo':
                p1=-oo    
            else:
                p1=Integer(p1)
            if "." in p2 or "/" in p2:
                p2=Rational(p2)
            elif p2=='oo':
                p2=oo
            elif p2=='-oo':
                p2=-oo    
            else:
                p2=Integer(p2)    
                
            Iq1,Iq2=self.getIqs()
            Iq1.psymb=p2
            Iq1.symb=nice_iq_symbol(p2)
            Iq1.e1.ksym=p1
            Iq2.psymb=p3
            Iq2.symb=nice_iq_symbol(p3)
            Iq2.e2.ksym=p4
            self.s()
        else:    
            return kres
    def inverse(self):
        Iq1,Iq2=self.getIqs()
        Iq1.inverse(kshow=False)
        Iq2.inverse(kshow=False)
        s1=Iq1.psymb
        sn1=Iq1.symb
        s2=Iq2.psymb
        sn2=Iq2.symb
        p1=Iq1.L
        p2=Iq2.R
        Iq1.psymb=s2
        Iq1.symb=sn2
        Iq1.e1.ksym=p2
        Iq2.psymb=s1
        Iq2.symb=sn1
        Iq2.e2.ksym=p1
        self.s()
        
    def str(self):
        p1,s1,px,s2,p2=self.getdata()     
        return str(p1)+s1+str(px)+s2+str(p2)
        
        
def get_str_eq(Obj):
    p1=str(Obj.L)
    sx=Obj.psymb
    p2=str(Obj.R)
    sexpr=p1+' '+sx+' '+p2
    return sexpr
    
def iinverse(expr):
    if Is_Div(expr):      
        p1,p2=fraction(expr)
        return cfrac(p2,p1)        
    elif expr==0:         
        return oo
    else:       
        return cfrac(1,expr)   
    return cfrac(1,expr)

def isolve(*args,var=x):
    mm=[]
    for i in args:
        if type(i)==MyIneq3:
            Iq1,Iq2=i.getIqs()
            mm.append(Iq1)
            mm.append(Iq2)
        else:
            mm.append(i)
     
    ms=[]
    for i in mm:
        ms.append(get_str_eq(i)) 
    cc=1
    pexpr="[["
    for i in ms:
        if cc==1:
            pexpr=pexpr+i
            cc=cc+1
        else:
            pexpr=pexpr+","+i
            cc=cc+1
            
    pexpr=pexpr+"]]"
    kres=reduce_rational_inequalities(parse_expr(pexpr),var)
    return kres
    
def getsimpleIq(Obj):
    if type(Obj)!=str:
        ss=str(Obj)
    else:
        ss=Obj
    ee1=ss.find(' ')
    p1=ss[1:ee1]
    s2=ss[ee1+1::]
    ee2=s2.find(' ')
    p2=s2[0:ee2]
    s3=s2[len(p2)+9::]
    ee3=s3.find(' ')
    p3=s3[0:ee3]
    s4=s3[len(p3)+1::]
    ee4=s4.find(')')
    p4=s4[0:ee4]
    if "." in p1 or "/" in p1:
        p1=Rational(p1)
    else:
        p1=Integer(p1)
    if "." in p4 or "/" in p4:
        p4=Rational(p4)
    else:
        p4=Integer(p4)

    return p1,p2,p3,p4


def getisymb(strs):
    veci=['>','<','=']
    msim=''
    for i in strs:
        if i in veci:
            msim=msim+i
    return msim
def getinum(strs,var=x):
    sx=str(var)
    veci=['>','<','=',sx,' ']
    msim=''
    for i in strs:
        if i not in veci:
            msim=msim+i
    if '.' in msim or '/' in msim:
        return Rational(msim)
    elif msim=='-oo':
        return -oo
    elif msim=='oo':
        return oo    
    else:
        return Integer(msim)
def Iorden(obj,var=x):
    if type(obj)!=str:
        obj=str(obj)
    pp=obj.find(str(var))
    kres=1
    if pp>0:
        kres=1
    return kres
def getIargs(obj,var=x):
    if type(obj)!=str:
        obj=str(obj)
    ps=getisymb(obj)    
    px=getinum(obj,var=var)
    po=Iorden(obj,var=var)
    return ps,px,po

def partIq(obj,var=x):
    s1,s2=obj.args
    sg1,v1,o1=getIargs(s1)
    sg2,v2,o2=getIargs(s2)
    
    if v2<v1:    
        return v2,sg2,var,sg1,v1
    else:
        return v1,sg1,var,sg2,v2

def isolu2Iq(*args,var=x,kshow=True):
    Obj=args[0]
    try:
        Obj=args[0] 
        p1,sg1,var,sg2,p2 =partIq(Obj,var=var)
        return MyIneq3(p1,sg1,var,sg2,p2,kshow=kshow)
         
    except:
        obj=args[0]
        pp=''
        mm=obj.args
        cc=1
        for i in mm:
            Iqq=isolu2Iq(i,var=x,kshow=False)
            p1,s1,px,s2,p2=Iqq.getdata()
            if cc==1:
                pp='['+latex(p1)+s1+latex(px)+s2+latex(p2)+']'
                cc=cc+1
            else:
                pp=pp+ ' U ['+latex(p1)+s1+latex(px)+s2+latex(p2)+ ']' 
             
        display(Math(pp))       
