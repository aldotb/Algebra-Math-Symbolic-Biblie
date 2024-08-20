from sympy import *
from lib_Mathematica import *
from lib_Mathbasic import *
from lib_Algorith import *

from lib_Exponencial import *
from lib_tools import *
from lib_MyEq import *
from lib_MyEqEq import *
 
from lib_MyMatrix import MyMat
from  copy import deepcopy
from numpy import *
from sympy import sin,cos,tan,exp,log,pi,cosh,sinh
t1,t2,t3,t4,t5,t6=symbols('t1 t2 t3 t4 t5 t6')
# Convertidores
########################
def Obj2GeometryObj(obj):
    
    if type(obj)==MyPoint3D:
        return obj.Obj
    elif type(obj)==MyLine3D:
        return obj.Obj
    elif type(obj)==MyPlane3D:
        return obj.Obj
    elif type(obj)==MyVector:
        return obj.Obj
    else:
        return obj
        
def Obj2Point3D(obj):
    if type(obj)==Point3D:
        return obj
    elif type(obj)==MyPoint3D:
        return Point3D(obj())
        
    else:
        return Point(obj)

        
 

        
 

def Obj2Vector(obj):
    if type(obj)==MyPoint3D:
        kres= Matrix(obj.P) 
    elif  type(obj)==Point3D:
        kres= Matrix(list(obj))
    elif type(obj)==MyVector:
        return Matrix(obj.Obj)
    else:
        kres= Matrix(obj) 
    return kres.T

 
    
###################

X,Y=symbols('X Y')
cx,cy,cz=symbols('cx,cy,cz')
def obj2list(obj):
    if type(obj)==Point3D:
        return [obj.x,obj.y,obj.z]
    elif type(obj)==MyPoint3D:
        return [obj.P.x,obj.P.y,obj.P.z]
    else:
        return list(obj)

def  list2Point(*args):
    if len(args)==1:
        P=args[0]
        if type(P)==Point or type(P)==Point3D:
            return P
        elif type(P)==list:
            x1,y1,z1=P
            return Point3D(x1,y1,z1)
         
    else:
        return Point3D(args[0],args1[1],args[2])


def obj2Point3D(*args):
    return Obj2Array(args)
    
class MyVector():
    def __init__(self, *args,var=t,kshow=True):
        if len(args)>1 and type(args[0])!=tuple:
            vec=[]
            for i in args:
                vec.append(i)
            self.Obj=Array(vec)
        elif len(args)==1:
            self.Obj=Obj2Array(args[0])           
        else:
            P1=Obj2Array(args[0])
            P2=Obj2Array(args[1])
            P=P2-P1
            self.Obj=Obj2Array(P)
        if kshow:
            display(Math(latex(self.Obj)))

    def __call__(self):
        kres=self.Obj
        return  Array(kres)
 
    def __repr__(self):
        kres = str(self.Obj)
        return kres
    def _latex(self, obj):
        return latex(self.Obj)

    def __str__(self):
        return self.__repr__()

    @property
    def modulo(self):
        L=Obj2List(self.Obj)
        kres=0
        for i in L:
            kres=kres+i*i
        return simplify(rsimplify(sqrt(kres)))
        
    @property
    def unit(self):
        return MyVector(unitary(self.Obj),kshow=False)
    
    def dot(self,Obj2):
        obj2=Obj2Matrix(Obj2)
        obj1=Obj2Matrix(self.Obj)
        kres=obj1.dot(obj2)
        return kres
        
    def cross(self,Obj2):
        obj2=Obj2Matrix(Obj2)
        obj1=Obj2Matrix(self.Obj)
        kres=obj1.cross(obj2)
        return MyVector(Obj2Array(kres),kshow=False)
        
    def cosangle(self,Obj2):
         
        V2=Obj2Array(Obj2)
        vec2=MyVector(V2,kshow=False)
        p1=self.dot(vec2) 
        p2=self.modulo*vec2.modulo
        return p1/p2
        
    def angle(self,Obj2):
        return acos(self.cosangle(Obj2))
    
    def reduce(self):
        ksym=self.Obj
        kres=simplifydirection(ksym)
        return MyVector(kres,kshow=False)
    
    @property    
    def L(self):
        kres=self.Obj
        return Obj2List(kres)
    
    @property    
    def M(self):
        kres=self.Obj
        return Obj2Matrix(kres)    
    
def proyection(vec1,vec2):
    V1=MyVector(Obj2Array(vec1),kshow=False)
    V2=MyVector(Obj2Array(vec2),kshow=False)
    Ca=V1.cosangle(V2)
    Mo=V1.modulo
    kres=Ca*Mo*V2.unitary()
    return kres        
        
class MyPoint3D(Point3D):
    def __init__(self, *args,var=t):
        if len(args)==1:
            P=obj2list(args[0])
            self.P=Point3D(Point(*P))
            self.Obj=Point3D(Obj2GeometryObj(self.P))
            self.L=self.Obj
        else:
              
            self.P=Point3D(Point(*args))

        self.var=var
        self.Obj=Obj2GeometryObj(self.P)    
    def __call__(self, *args, **kwargs):
        if len(args)==1:
            val=obj2expr(args[0])
            x1=unisymbols(self.x.subs(self.var,val))
            y1=unisymbols(self.y.subs(self.var,val))
            z1=unisymbols(self.z.subs(self.var,val))
            return Point3D(x1,y1,z1)
            
        elif len(kwargs)>0:
            P=self.P
            x1=real_subs(P.x,**kwargs)
            y1=real_subs(P.y,**kwargs)
            z1=real_subs(P.z,**kwargs)
            return Point3D(x1,y1,z1)
        else:
            return Point3D(self.Obj)
    def midpoint(self,obj):
        P1=self.P
        if type(obj)==Point3D or type(obj)==MyPoint3D:
            P2=obj.P
        else:
            P2=Point3D(obj)
        return P1.midpoint(P2)
            
    
    def __repr__(self):
        kres = str(self.P)
        return kres
    def _latex(self, obj):
        return latex(self.P)

    def __str__(self):
        return self.__repr__()
    def distance(self,Obj):
        obj1=Obj2GeometryObj(self.P)
        obj2=Obj2GeometryObj(Obj)
        return obj1.distance(obj2)
    def Projection(self,Obj):
        if type(Obj)==MyLine3D:
            PP=Obj.arbitrarypoint(kshow=False)
            P1=Array(PP(0))
            P2=Array(PP(1))
            P3=Array(self.P)
            kres=ClosestPointOnLine(P1, P2, P3)
            return MyPoint3D(kres)
        else:
            
            Vo=list(Obj.ortogonal())
           
            P=tuple(self.Obj)
            L=MyLine3D(P,Vo,kshow=False)
            PL2=Obj.Obj
            L2=L.Obj
             
            kres= L2.intersection(PL2)
            if len(kres)==1:
                kres=kres[0]
                return MyPoint3D(kres)
            else:
                kres2=[MyPoint3D(i) for i in kres]
                return kres2
    def OrthogonalLine(self,Obj):
        P1=tuple(self.Obj)
        P2=list(self.Projection(Obj))
        return MyLine3D(P1,P2)
        
def MyLine3Direc(P,direc):
    P1=Point3D(P)
    P2=Point3D(P1.x+direc[0],P1.y+direc[1],P1.z+direc[2])
    return MyLine3D(P1,P2)

'''
class MyLine3D(Line3D):
    def __init__(self, *args,var=t,kshow=True,**kwargs):
        self.var=var
        if len(args)==1 and len(kwargs)==0:
            self.Obj=Obj2GeometryObj(args[0])
            print(1)
        elif  len(args)==2 and len(kwargs)==0 and type(args[0])==tuple and type(args[1])==list:
            self.Obj=Line3D(args[0],direction_ratio=args[1])
            print<(2)
        else:
            self.Obj=Line3D(*args,**kwargs)
            print(3)
        self.L=self.Obj    
        if kshow:
            display(Math(latex(self.Obj)))
'''    
class MyLine3D():
    def __init__(self,*args,var=t,kshow=True):
        self.var=var
        P1=Obj2Point3D(args[0])
        P2=args[1]
        if type(P2)!=list:
            vec=P1-Obj2Point3D(P2)
            self.Obj=Line3D(P1,direction_ratio=vec)
        else:
            vec=P2
        self.Obj=Line3D(P1,direction_ratio=vec)
 
        if kshow:
            L=self.Obj
            v1=L.direction_ratio
            var=self.var
            P1=list(L.p1)
            display(Math(latex(P1)+'\:'+'\:'+'+'+str(var)+latex(v1)))
            
        


            
    def __call__(self, *args, **kwargs):
        if len(args)==0 and len(kwargs)==0:
            return self.Obj
        if len(args)==1:
            L=self.Obj
            P=L.p1
            P1=MyPoint3D(P)
            P11=P1(args[0])
            P=L.p2
            P2=MyPoint3D(P)
            P22=P2(args[0])
            return Line3D(P11,P22)
            
             
    def __repr__(self):
        kres = str(self.Obj)
        return kres
    def _latex(self, obj):
        return latex(self.Obj)

    def __str__(self):
        return self.__repr__()

    def arbitrarypoint(self,var='',kshow=True):
        if var=='':
            var=self.var
        L=self.Obj
        PP=L.arbitrary_point(var)
        kres= Obj2Point3D(PP)
        kres=MyPoint3D(kres)
        if kshow:
            display(kres)
        return kres
    @property    
    def p1(self):
        L=self.Obj
        return L.p1
    @property    
    def p1(self):
        L=self.Obj
        return L.p1 
    def set(self,*args,**kwargs):
        var=self.var
        L=self.Obj
        P1=L.p2
        P2=L.p2
        P11=subsArray(Obj2Array(P1),*args,var=t,**kwargs)
        
        
    def dot(self,Obj,kname=''):
        Obj1=list(self.direction)
        if type(Obj)==Line3D or type(Obj)==MyLine3D:
            Obj=list(Obj.direction)
        else:
            Obj=list(Obj)
        kres= SumProVec(Obj1,Obj)
        if kname!='':
            return MyEq(kres,kname=kname,var=self.var)
        else:
            return kres
    
    def parametric(self,var='',kshow=True):
        if var=='':
            var=self.var
        vec=self.arbitrarypoint(var=var,kshow=False)
        if kshow:
            display(Math(latex(vec)))
        return Array(vec)

        
    def eQvectorial(self):
        
        sdr=str(tuple(simplifydirection(self.direction())))
        vecp=str(tuple(self.p1))
        svar=str(self.var)
        display(Math(vecp+'\:'+ svar+sdr)) 
    
    def vecdirection(self,*args):
        L=self.Obj
        return Obj2Array(L.direction)
    
    def direction(self):
        L=self.Obj
        return L.direction

    def eQrecta(self):
        vp=Array(self.p1)
        vpara=Array(self.arbitrarypoint(kshow=False))
        var=self.var
        vt=vpara-vp
        vt2=vt.subs(var,1)
        vxyz=[x,y,z]
        vecxyz=[]
        for i,j,k in zip(vp,vt2,vxyz):
            vecxyz.append(factor((k-i)/j))
        return vecxyz
        
    @property
    def vectorDirecction(self):
        L=self.Obj
        kres=list(L.direction)
         
        return kres

    def angle(self,Obj,*args):
        v1=self.vectorDirecction
        if type(Obj)==MyLine3D or type(Obj)==Line3D:
            if type(Obj)==MyLine3D:
                v2=Obj.vectorDirecction 
            else:
                Obj2=Obj.Obj
                v2=list(Obj2.direction)
            kres=angle2vector(v1,v2,*args)
        else:
            if type(Obj)==MyPlane3D:
                v2=Obj.vectorOrtogonal
            else:
                Obj2=Obj.Obj
                v2=list(Obj2.normal_vector)
            kres=angle2vector(v1,v2,*args,type='sin')
        return kres    
            
    def distance(self,Obj):
        obj2=Obj2GeometryObj(Obj)
        obj1=self.Obj
        D=obj1.distance(obj2)
        return D
            
    
def eQ2MyPlane3D(expr):
    P,V=func2PointNormal(unisymbols(expr))
    return MyPlane3D(P,V)

    
class MyPlane3D(Plane):
    def __init__(self,p1,a=None,b=None,var=t,kshow=True,**kwargs): #Plane(p1, a=None, b=None, **kwargs)  
        self.var=var
        if a and b:
            P1=Obj2List(Point3D(p1))
            P2=Obj2List(Point3D(a))
            P3=Obj2List(Point3D(b))
            self.PL=Plane(P1,P2,P3)
            self.Obj=Plane(P1,P2,P3)

        else:
            P1=Point3D(Obj2List(p1)) 
            self.PL=Plane(p1,a)
            self.Obj=Plane(p1,a)
            
        if kshow:
            display(Math(latex(self.PL)))
    def __call__(self):
         return self.Obj       
    def __repr__(self):
        kres = self.PL
        return kres
    def _latex(self, obj):
        return latex(self.PL)

    def __str__(self):
        return self.__repr__()

    def ortogonal(self,*args):
        PL=self.PL
        kres=PL.normal_vector
        return MyVector(list(kres))
        
         
    def vecNormal(self,*args):
        PL=self.PL
        kres=list(PL.normal_vector) 
         
        return MyVector(kres,kshow=False)  
        
    def intersection(self,Obj):
        return self.intersecta(Obj)

        
    def intersecta(self,Obj):
        PL1=self.PL
        PL2=Obj.PL
        PL3=PL1.intersection(PL2)
        if len(PL3)==1:
            PL3=PL3[0]
        if type(PL3)==Line3D:
            P1=PL3.p1
            direct=PL3.direction_ratio
            if signo(direct[0])==-1:
                direc2=[-1*i for i in direct]
                direct=direc2
            return MyLine3D(P1,list(direct))
        else:
            return PL3
    def arbitrarypoint(self,*args,kshow=True):
        L=self.Obj
        var=self.var
        val=''

        for i in args:
            if type(i)==Symbol:
                var=i
            if type(i)!=Symbol:
                val=i
        P=L.arbitrary_point(var)
        if val!='':
            P=P.subs(var,val)
        return P    
        
    def eQ(self,kshow=True,**kwargs):
        PL=self.Obj
        kres=unisymbols(PL.equation())
        if len(kwargs)>0:
            kres=real_subs(kres,**kwargs)
        if kshow:
            display(Math(latex(kres)))
        return kres
    def equation(self,*args):
        obj=self.Obj
        kres=obj.equation()
        if 'simplify'in args:
            kres=numer(factor(kres))
        return kres    
    
class  MyParametric():
    def __init__(self, *args,var=t,var2=y,L1='',L2='',npoint=100,kshow=True):
         
        vecc=[]
        self.var=var
        self.name=''
        if type(args[0])==str:
            self.name=args[0]
        self.F=[unisymbols(i) for i in args if type(i)!=str]
        self.L1=L1
        self.L2=L2
        self.npoint=npoint 
        self.A=Array(self.F)
        self.Obj=Array(self.F)
        if kshow:
            display(Math(self.name+'= '+latex(self.F)))
    def __call__(self, *args, **kwargs):
        var=self.var
        kres=self.Obj
        if len(args)==1 and len(kwargs)==0:
            kres=kres.subs(var,args[0])
            return kres
        elif len(args)==0 and len(kwargs)>0:
            Lkres=list(kres)
            kres=Array([real_subs(i,**kwargs) for i in Lkres])
            return kres
        else:
            return kres
            
    def __repr__(self):
        kres = str(self.Obj)
        return kres
    def _latex(self, obj):
        return latex(self.Obj)

    def __str__(self):
        return self.__repr__()

    
    
        
    def unitary(self):
        ''' return unitary vector of f type: Array '''
        vecL=list(self.Obj)
        dres=0
        for i in vecL:
            dres=dres+i*i
        kres=simplify(dres)
        kres=sqrt(kres)
        kres=rsimplify(kres)
        vecr=[]
        for i in vecL:
            vecr.append(i/kres)
        return Array(vecr)    
    
    def function(self,*args,**kwargs):
        ''' return  list(f.vector, *args, **kwargs)  type: list '''
        F=self.Obj
        var=self.var
        kres=subsArray(F,*args,var=var,**kwargs)
        return kres
           
    def args(self,pos):
        ''' return  f.vectors[pos]  type: expr '''
        return self.Obj[pos]

    def diff(self, *args,kshow=True):
        ''' f.diff(): return f.diff(var)
            f.diff('df'): return new MyParametric Obj wht name 'df' '''
            
        kname=''
        vecvar=[]
        var=self.var
        val=''
        for i in args:
            if type(i)==str:
                kname=i
            elif i==var:
                vecvar.append(i)
            else:
                val=i
        kres=self.Obj
        if vecvar==[]:
            vecvar=[var]
        dkres=kres.diff(*vecvar)
        if val!='':
            dkres=dkres.subs(var,val)
        if kname!='':
            nP=[kname]
            for i in dkres:
                nP.append(i)
               
            nPP= MyParametric(*nP,kshow=kshow)
             
            return nPP    
        else:
            return dkres
            
   
        
    def dot(self,v2,**kwargs):
        V1= Obj2Matrix(self.Obj)
        V2= Obj2Matrix(v2)
        kres=V1.dot(V2)
        if len(kwargs)>0:
            kres=real_subs(kres,**kwargs)
        return kres
         

    def limit(self,*args):
        var=self.var
        lm=0
        dir=''
        for i in args:
            if type(i)==Symbol:
                var=i
            elif i=='+' or i=='-':
                dir=i
            else:
                lm=i
        if dir=='':
           vecres=[limit(i,var,lm) for i in self.A]
        else:
           vecres=[limit(i,var,lm,dir=dir) for i in self.A]

        return Array(vecres)    
    
    def __mul__(self, other):
        """ Returns the vector addition of self and other """
         
        AA=Matrix(self.Obj)
        AA=AA.T
        done=True
        if type(other) ==  MyParametric:
 
            BB=Matrix(other.Obj)
            BB=BB.T
             
        elif type(other)== ImmutableDenseNDimArray:
            BB=Matrix(other)
            BB=BB.T
             
        elif type(other)== list:
            BB=Matrix(other)
            BB=BB.T
             
        else:
            done=False
        
        
        if done:
            kres=list(AA.cross(BB).T)
            return Array(kres)
        else:
            return other*self.A

        
    def __rmul__(self, other):
        """ Returns the vector addition of self and other """
        """ Returns the vector addition of self and other """
        AA=Matrix(self.A)
        AA=AA.T
        done=True
        if type(other) == vFunction:
            BB=Matrix(other.A)
            BB=BB.T
            
        elif type(other)== ImmutableDenseNDimArray:
            BB=Matrix(other)
            BB=BB.T
             
        elif type(other)== list:
            BB=Matrix(other)
            BB=BB.T
             
        else:
            done=False
        if done:
            kres=list(BB.cross(AA).T)
            return Array(kres)
        else:
            return other*self.A  
    
    @property
    def list(self):
        return list(self.A)

    def parametric(self,*args):
        vec=self.list
        var=self.var
        if len(args)==0:
            args=[x,y,z]
        P = [factor(simplesolve(vec[i]-args[i],var, kshow=False)) for i in range(len(vec))]
        return P
        
    def squareSum(self,kname='SC2',kshow=True):
        kres=0
        for i in self.F:
            kres=kres+i*i
        return MyEq(kres,kname=kname,var=self.var,kshow=kshow)    
        
    def simplify(self):
        kres=[]
        for i in self.F:
            kres.append(simplify(i))
        self.F=kres
        self.A=Array(self.F)
        display(Math(latex(self.F)))
    
    def factor(self):
        kres=[]
        for i in self.F:
            kres.append(factor(i))
        self.F=kres
        self.A=Array(self.F)
        display(Math(latex(self.F)))
    def expand(self):
        kres=[]
        for i in self.F:
            kres.append(expand(i))
        self.F=kres
        self.A=Array(self.F)
        display(Math(latex(self.F)))
    def tsimplify(self):
        kres=[]
        for i in self.F:
            kres.append(tsimplify(i))
        self.F=kres
        self.A=Array(self.F)
        display(Math(latex(self.F)))

    def lengthCurve(self,*args):
        t=symbols('t',real=True,positive=True)
        kname=''
        x1=''
        x2=''
        for i in args:
            if type(i)==str:
                kname=i
            else:
                if x1=='':
                    x1=i
                else:
                    x2=i
                    
        
        kres=self.diff('dC',kshow=False)
        kres2=kres.squareSum(kshow=False)
        kres2.ksym=sqrt(kres2.ksym)
        return kres2.integral(kname,x1,x2,self.var)
        
    def Pow(self,*args,kshow=True):
        mm=Obj2List(self.Obj)
        mm2=[i**args[0] for i in mm]
        mm3=Obj2Array(mm2)
        self.Obj=mm3
        if kshow:
            self.s()
    def s(self):
        kname=self.name
        display(Math(kname+': \:'+'\:'+latex(self.Obj))) 

    def module(self,*args):
        kk=0 
        for i in self.F:
            kk=kk+(i*i)
        kres=tsimplify(simplify(kk))
        kres=rsimplify(sqrt(kres))
         
        if len(args)==1:
            return kres.subs(self.var,args[0])
        else: 
            return kres
            
    
       
    def subs(self,v1,v2,kshow=True):
        A=self.A
        sv1=str(v1)
        vec=A.subs(sv1,v2)
        self.F=list(vec)
        self.A=vec
        if kshow:
            display(Math(self.name+'= '+latex(self.F)))

    def getPoint(self,*args):
        if len(args)==0:
            P= Point3D(self.F)
        else:
            vec=self(args[0])
            vec=list(vec)
            P= Point3D(vec)
        display(Math(latex(P)))
        return P

    def vtangente(self,*args):
        vec=self.Obj
        var=self.var
        dvec=vec.diff(var,*args)
        dmodulo=modulo(dvec)
        kres=dvec/dmodulo
        if len(args)==1:
            kres=simplify(kres.subs(var,args[0]))
        
        return kres
    def vNormal(self,*args):
        vec=self.vtangente()
        var=self.var
        vec1=vec.diff(var,*args)
        vec1=simplify(vec1)
        dmodulo=simplify(modulo(vec1))
        kres=vec1/dmodulo
        if len(args)==1:
            kres=simplify(kres.subs(var,args[0]))
        return kres

    def vBinormal(self,*args):
        vecT=Obj2Matrix(self.vtangente(*args))
        vecN=Obj2Matrix(self.vNormal(*args))
        kres=crossProduct(vecT,vecN)
        return kres
        
    

        

    def tangente(self,*args):
        vec=self.Obj
        var=self.var
        dvec=vec.diff(var)
        if len(args)==1:
            return dvec.subs(self.var,args[0])
        else:
            return dvec 

    def curvature(self,*args):
        var=self.var
        if len(args)==1:
            p1=self.diff(var,args[0])
            p2=self.diff(var,var,args[0])
        else:
            p1=self.diff(var)
            p2=self.diff(var,var)            
        p3=crossProduct(p1,p2)
        mod1=modulo(p3)
        mod2=modulo(p1)
        kres1=simplify(rsimplify(mod1))
        kres2=simplify(rsimplify(mod1))
        kres2=kres2**3
        return kres1/kres2

        
    def veccurvature(self,*args):
        kname=''
        valor=''
        for i in args:
            if type(i)==str:
                kname=i
            else:
                valor=i
                
        if kname!='':
            if valor!='':
                d2F=self.diff(self.var,self.var,valor,kname)
            else:
                d2F=self.diff(self.var,self.var,kname)
        else:
            if valor!='':
                d2F=self.diff(self.var,self.var,valor)
            else:
                d2F=self.diff(self.var,self.var)
                
        return d2F
        
          
    def torsion(self,*args):
        var=self.var
        if len(args)>0:
            p1=self.diff(var,args[0])
            p2=self.diff(var,var,args[0])
            p3=self.diff(var,var,var,args[0])
        else:
            p1=self.diff(var)
            p2=self.diff(var,var)
            p3=self.diff(var,var,var)
        p4=crossProduct(p1,p2)
        m1=Obj2Matrix(p4)
        m2=Obj2Matrix(p3)
        m3=modulo(p4)
        m4=m1.dot(m2)
        kres=simplify(rsimplify(m4/m3**2))
        return kres 

    def acosangle(self,Obj,*args,**kwargs):
        kres= angle2vector(self.function,Obj,type='angle')
        if len(args)==1:
            var=self.var
            kres= kres.subs(var,args[0])
        if len(kwargs)>0:
            kres=real_subs(kres,**kwargs)
        return kres    
            
            
    def arbitrarypoint(self,*args):
        pp=Obj2List(self.Obj)
        P=MyPoint3D(pp)
        if len(args)==1:
            var=self.var
            P=P.subs(var,args[0])
        return P

    def integralf(self,*args,x1='',x2=''):
        cc=0
        C1,C2,C3=symbols('C1 C2 C3')
        vecc=[C1,C2,C3]
        kres=Obj2List(self.Obj)
        vec=[]
        var=self.var
        for i in kres:
            newI=creteIntegral(i,var,x1=x1,x2=x2).doit()
            newI=newI+vecc[cc]
            cc=cc+1
            vec.append(newI)
        if len(args)==1 and type(args[0])==str:
            arg2=[args[0]]
            for i in vec:
                arg2.append(i)
                
            return MyParametric(*arg2)
        else:
            self.Obj=Obj2Array(arg2)
        if kshow:
            if self.name!='':
                display(Math(self.name+'='+'\;'+latex(self.Obj)))
            else:    
                display(Math(latex(self.Obj)))
        else:
            return self.Obj
    def set(self, kshow=True,**kwargs):
        vec=Obj2List(self.Obj)
        vec2=[]
        for i in vec:
            vec2.append(real_subs(i,**kwargs))
        self.Obj=Obj2Array(vec2)
        if kshow:
            if self.name!='':
                display(Math(self.name+'='+'\;'+latex(self.Obj)))
            else:    
                display(Math(latex(self.Obj)))
        else:
            return self.Obj
    def doit(self,kshow=True):
        cc=0
        C1,C2,C3,C4=symbols('C1 C2 C3 C4')
        vecC=[C1,C2,C3,C4]
        L=Obj2List(self.Obj)
        vec=[]
        for i in L:
            try:
                vec.append(i.doit() )
            except:
                vec.append(i)
        self.Obj=Array(vec)
        if kshow:
            if self.name!='':
                display(Math(self.name+'='+'\;'+latex(self.Obj)))
            else:    
                display(Math(latex(self.Obj)))
        else:
            return self.Obj
    @property        
    def List(self):
        return Obj2List(self.Obj)         
    
    def L(self,*args):
        kres=Obj2Array(self.Obj)
        if len(args)==1:
            var=self.var
            kres=kres.subs(var,args[0])
        return Obj2List(kres)    
    @property
    def Array(self):
        return Obj2Array(self.Obj)
    
    def A(self,*args):
        kres=Obj2Array(self.Obj)
        if len(args)==1:
            var=self.var
            kres=kres.subs(var,args[0])
        return kres
    
    
    @property
    def Matrix(self):
        vec=Obj2Matrix(self.Obj)
        return Matrix(vec)
        
    def M(self,*args):
        kres=Obj2Array(self.Obj)
        if len(args)==1:
            var=self.var
            kres=kres.subs(var,args[0])
        kres=Obj2Matrix(kres)
        return kres.T
        
    @property
    def T(self):
        M=self.Matrix
        return M.T
        
    def cross(self,obj2,*args):
 
        vec4=crossProduct(self.Obj,obj2)
        if len(args)==1:
            val=args[0]
            var=self.var
            vec5=vec4.subs(var,val)
            return simplify(rsimplify(vec5))
        else:
            return simplify(rsimplify(vec4))
    def Ucross(self,obj2,*args):
        vec1=self.cross(obj2,*args)
        kmod=modulo(vec1)
        return vec1/kmod
        
class MyPlane2():
    def __init__(self,*args):
         
        if len(args)==3:
            xx,yy,zz=args
            vecp=[]
            for i in args:
                
                self.PL=Plane(xx,yy,zz)
                self.P=xx
                self.vecNorm=NormalVec(obj2list(yy),obj2list(zz))
        if len(args)==2:
            P=obj2Point3D(args[0])

            self.PL=Plane(P,obj2list(args[1]))
            self.P=P
            self.vecNorm=obj2list(args[1])
             
        else:
            P,VecN=func2PointNormal(args[0])
            self.PL=Plane(P,VecN)
            self.P=P
            self.vecNorm=VecN

        
    def equation(self,kname=''):
        Pe=unisymbols(self.PL.equation())
        if kname!='':
            ee=MyEq(Pe,kname=kname)
            return ee
        else:    
            return Pe

     
    def evaluatePoint(self,*args,**kwargs):
        
        eQQ=unisymbols(self.equation())
        eQQ=real_subs(eQQ,**kwargs)
        ee=MyEq(eQQ,'ee',kshow=False)
        
        if len(args)==1:
            P=args[0]
            if type(P)==Point3D:
                x1=P.x
                y1=P.y
                z1=P.z
            else:
                x1,y1,z1=args[0]
            kres=ee(x=x1,y=y1,z=z1)
            return kres
        return ee.ksym
        
        

              
def obj2MyGeoObj(vec):
    if type(vec)==list:
        kres=[]
        for i in vec:
            kres.appenc(obj2Myobj(i))
        if len(kres)==1:
            return kres[0]
    elif type(vec)==Point2D:
        return MyPoint(vec)
    elif type(vec)==Line2D:
        return MyLine(vec)
    elif type(vec)==Circle:
        return MyCircle(vec)
    else:
        return vec
def dataargs(*args):
    kk=[]
    for i in args:
        if type(i)!=str:
            kk.append(i)
    return kk
def dataMyGeo(*args):
    args2=[]
    for i in args:
        if type(i)==MyPoint:
            args2.append(i.P)
        else:
            args2.append(i)
    return args2        
def Pointxy(*args):
    if len(args)==1:
        return args[0].x,args[0].y
    elif len(args)==4:
        return args[0],args[1],args[2],args[3]
    elif len(args)==2:
        if type(args[0])==Point2D:
            return args[0].x,args[0].y,args[1].x,args[1].y
        else:
            return args[0],args[1]
    else:
        try:
            return args[0].x,args[0].y,args[1],args[2]
        except:
            return args[0],args[1],args[2].x,args[2].y

def eQrecta(*args,slope=''):
    vecargs=['eQ','Eq','flat','noshow','Line']
    kname=''
    for i in args:
        if type(i)==str and i not in vecargs:
            kname=i

    args2=dataargs(*args)
    if type(args2[0])==Line2D:
        L=args2[0]
    elif len(args2)==1 and slope=='' and (str(x) in str(args2[0]) or str(y) in str(args2[0])):
        cx,cy,cc=getcoefrecta(args2[0])
        L=Line(cx*x+cy*y+cc)    
    else:
        if slope=='':
            x1,y1,x2,y2=Pointxy(*args2)
            L=Line((x1,y1),(x2,y2))
        else:
            x1,y1=Pointxy(*args2)
            L=Line((x1,y1),slope=slope)
    eqr=L.equation()
    if 'flat' in args:
        eqr=reducecero(factor(eqr))
    if 'Line' in args:
            return L    
    if kname!='':
        if kname=='y':
            Q = MQ(0,eqr,kshow=False)
            Q.alone(y,kshow=False)
             
            return MyEq(Q.R,kname='y',var=x)
         
        else:
            ee = MyEq(eqr,kname,var=x)
            return ee 
    else:
        return eqr
        
def midlepoint(*args):
    x1,y1,x2,y2=Pointxy(*args)
    return Point((x2+x1)/2,(y2+y1)/2)      

def vec_idempoint(*args):
    x1,y1,x2,y2=Pointxy(*args)
    return [x2-x1,y2-y1] 
    

def angle(obj1,obj2):
    return obj1.angle_between(obj2)

def lenght(*args):
    x1,y1,x2,y2=Pointxy(*args)
    return sqrt((x2-x1)**2+(y2-y1)**2)
    
def distanceP(*args):
    x1,y1,x2,y2=Pointxy(*args)
    return sqrt((x2-x1)**2+(y2-y1)**2)    


def distance(*args):
    vecC,vecL,vecP,vecxy,vecs=obj2gobj(*args)
    if len(vecL)==1:
        LL=vecL[0]
        if len(vecP)==1:
            P=vecP[0]
        if len(vecxy)==2:
            P=Point(vecxy[0],vecxy[1])
        kres= P.distance(LL)   
    else:
        if len(vecxy)==4:
            x1,y1,x2,y2=vecxy
        if len(vecxy)==2 and len(vecP)==1:
            x1,y1,x2,y2=vecxy[0],vecxy[0],vecP[0].x,vecP[0].y
        if len(vecP)==2:
            x1,y1,x2,y2=vecP[0].x,vecP[0].y,vecP[1].x,vecP[1].y
         
        kres= sqrt((x2-x1)**2+(y2-y1)**2)
    if 'simplify' in args:    
        return rsimplify(factor(expand(kres)))
    else:
        return kres
    
def distancesquare(*args):
    kres=distance(*args)
    kres=kres**2
    kres=rsimplify(kres)
    return kres
def distance2(*args):
    kres=distance(*args)
    kres=kres**2
    kres=rsimplify(kres)
    return kres    
def getcoefrecta(expr):
    expryc=expr.subs(x,0)
    exprc=expryc.subs(y,0)
    exprx=expr-expryc
    cx=exprx.subs(x,1)
    expry=expryc-exprc
    cy=expry.subs(y,1)
    return cx,cy,exprc
   

def obj2Line(obj):
    if type(obj)==Line2D:
        return obj
    elif type(obj)==MyLine:
        return obj.L
    else:
        if type(obj)==MyEq:
            expr=obj.ksym
        else:
            expr=obj
        cx,cy,cc=getcoefrecta(expr)
        lexpr=Line(cx*x+cy*y+cc)
        return lexpr
        
            
     
        

def Lintersec(L1,L2):
    kres=L1.intersection(L2)
    if type(kres)==list and len(kres)==1:
        return kres[0]
    else:
        return kres
def obj2sympy(obj):
    if type(obj)==MyLine:
        return obj.L
    if type(obj)==MyCircle:
        return obj.O
    return obj

def Lcircle1(*args):
    kname=''
    for i in args:
        if type(i)==Point2D:
            P=i
        if type(i)==MyPoint:
            P=i.P
        if type(i)==MyCircle:
            O=i
        if type(i)==str:
            kname=i
    px,py=P.x,P.y
    cx,cy=O.center()
    Lm=MyLine(cx,cy,px,py)
    mm=Lm.slope()
    mm2=-1/mm
    return MyLine(P,kname,slope=mm2)
 

def vecanswer(kres):
    try:
        if type(kres)==list:
            if len(kres)==1:
                return kres[0]
            else:
                return kres
    except:
        return kres
    
def obj2gobj(*args):
    vecxy=[]
    vecp=[]
    vecs=[]
    vecL=[]
    vecC=[]
    vecF=[]
    for i in args:
        if type(i)==str:
            vecs.append(i)
        elif type(i)==Point2D:
            vecp.append(i)
        elif type(i)==MyPoint:
            vecp.append(i.P)    
        elif type(i)==Line2D:
            vecL.append(i)
        elif type(i)==MyLine:
            vecL.append(i.L)    
        elif type(i)==Circle:
            vecC.append(i)
        elif type(i)==MyCircle:
            vecL.append(i.O) 
        elif type(i)==MyCurve:
            vecF.append(i.F)
        else:
            vecxy.append(i)
    return vecC,vecL,vecp,vecxy,vecs,vecF 
def objgetP(*args):
    vecC,vecL,vecp,vecxy,vecs,vecF=obj2gobj(*args)
    return vecp
def objgetL(*args):
    vecC,vecL,vecp,vecxy,vecs,vecF=obj2gobj(*args)
    return vecL
def objgetC(*args):
    vecC,vecL,vecp,vecxy,vecs,vecF=obj2gobj(*args)
    return vecC
def objgetxy(*args):
    vecC,vecL,vecp,vecxy,vecs,vecF=obj2gobj(*args)
    return vecxy
def objgets(*args):
    vecC,vecL,vecp,vecxy,vecs,vecF=obj2gobj(*args)
    return vecs    
def objgetF(*args):
    vecC,vecL,vecp,vecxy,vecs,vecF=obj2gobj(*args)
    return vecF   
def geoObj(expr):
    if type(expr)==MyPoint:
        return expr.P
    if type(expr)==MyPoint3D:
        return expr.P    
    elif type(expr)==MyLine:
        return expr.L
    elif type(expr)==MyCircle:
        return expr.O
    elif type(expr)==MyPlane:
        return expr.PL
    elif type(expr)==MyLine3D:
        return expr.L3    
    else:
        return expr
def point_over(P,gobj):
        P=obj2Point(P)
         
        eq=unisymbols(gobj.eq)
        x1=P.x
        y1=P.y
        eq=eq.subs(x,x1)
        eq=eq.subs(y,y1)
        eq=simplify(eq)
        if eq==0:
            return True
        else:
            return False       
def get_gobjfunction(obj):
    try:
        return unisymbols(obj.equation())
    except:
        return obj   
def gobjFunc(obj):
    try:
        return unisymbols(obj.equation())
    except:
        return obj 


def eQ2SphereFormat(expr,kshow=True):
    cc,rr=completesquarexyz(expr,'center',kshow=False)
    ss=MySphere(center=cc,radio=rr)
    kres=ss.eQ()
    return kres
    
    
def completesquarexyz(QQQ,*args,var=t,kshow=True):
    '''
    input MyEq=A*x*x+B*y*y+C*z*z+D*x+E*y+F*z*G
    return Cx,Cy,Cx,R2
    '''

    QQ=deepcopy(QQQ)
    if type(QQ)!=MyEq:
        QQ=MyEq(QQ,'QQ',kshow=False)
    R=QQ(x=0,y=0,z=0)
    QQ.Substrac(R,kshow=False)
    ex2=x*x*vecsubs(QQ.ksym,[x**2,y**2,z**2,x,y,z],[1,0,0,0,0,0])
    ey2=y*y*vecsubs(QQ.ksym,[x**2,y**2,z**2,x,y,z],[0,1,0,0,0,0])
    ez2=z*z*vecsubs(QQ.ksym,[x**2,y**2,z**2,x,y,z],[0,0,1,0,0,0])
    ex=x*vecsubs(QQ.ksym,[x**2,y**2,z**2,x,y,z],[0,0,0,1,0,0])
    ey=y*vecsubs(QQ.ksym,[x**2,y**2,z**2,x,y,z],[0,0,0,0,1,0])
    ez=z*vecsubs(QQ.ksym,[x**2,y**2,z**2,x,y,z],[0,0,0,0,0,1])
    X,ax=completesquare(ex2+ex+a,var=x)
    Y,ay=completesquare(ey2+ey+a,var=y)
    Z,az=completesquare(ez2+ez+a,var=z)
    P=factor(simplify(expand(ax+ay+az-R)))
     
    if  'center' in  args:
        if kshow:
            display(Math('CenterPoint=('+latex(X)+','+latex(Y)+','+latex(Z)+')'))
            display(Math('R2='+latex(P)))
        return Point3D(X,Y,Z),P
    elif 'parametric' in args:
        vecee=[]
        ee=MQ((x-ax)/X,cos(var)**2)
        vecee.append(ee)
        ee=MQ((y-ay)/Y,sin(var)**2)
        vecee.append(ee)
        return (vecee)
            
    else:
        if kshow:
            display(Math('Cx='+latex(X)))
            display(Math('Cy='+latex(Y)))
            display(Math('Cz='+latex(Z)))
            display(Math('R2='+latex(P)))
        return (X,Y,Z,P)
def distancePoinPlane(px,py,pz,PP):
    kres1=PP(x=px,y=py,z=pz)
    cx1,cy1,cz1,R1=eQPlanoCoef(PP,kshow=False)
    r=cx1*cx1+cy1*cy1+cz1*cz1
    return kres1*kres1/r   
def coef_PlaneEq(QQ,kshow=True):
    expr=obj2expr(QQ)
    expr2=expr.subs(x,0)
    expr2=expr2.subs(y,0)
    expr2=expr2.subs(z,0)
    
    R=expr2
    expr3=expr-R
     
    ex=vecsubs(expr3,['x','y','z'],[1,0,0])
    ey=vecsubs(expr3,['x','y','z'],[0,1,0])
    ez=vecsubs(expr3,['x','y','z'],[0,0,1])
    if kshow:
        display(Math('X='+latex(ex)))
        display(Math('Y='+latex(ey)))
        display(Math('Z='+latex(ez)))
        display(Math('R='+latex(R)))
    return (ex,ey,ez,R)

def obj2Point3D(*args):

    if len(args)==1:
        obj=args[0]
        if type(obj)==Point3D:
            return obj
        elif type(obj)==tuple:
            p1,p2,p3=list(obj)
            return Point3D(p1,p2,p3)
        elif  type(obj)==list:
            p1,p2,p3=obj
            return Point3D(p1,p2,p3)
        else:
            return obj
    else:
        p1,p2,p3=args[0],args[1],args[2]
        return Point3D(p1,p2,p3)
        
        
        

def func2PointNormal(expr):
    if type(expr)==MyEq:
        expr=expr.ksym
    z1,X,Y,Z=symbols('z1,X,Y,Z')
    P2= Plane(Point3D(1, 1, z1),(X, Y, Z))
    e2=P2.equation()
    kk=simplify(expand(e2-expr))
    R=vecsubs(kk,[x,y,z],[0,0,0])
    kk2=kk-R
    p1=factor(vecsubs(kk2,[x,y,z],[1,0,0]))
    p2=factor(vecsubs(kk2,[z,y,x],[1,0,0]))
    p3=factor(vecsubs(kk2,[z,y,x],[0,1,0]))
    X,Z,Y,z1=simplesolve(p1,p2,p3,R,X,Z,Y,z1,'noshow',kshow=False)
    P=Point3D(1,1,z1)
    vN=(X,Y,Z)
    return P,vN

def function2Plane(fu):
    PP,Vn=func2PointNormal(unisymbols(fu))
    return MyPlane3D(PP,Vn)

    
def Point_from_eQplane(expr):
    P,vN=get_Plane_Poin_Normal(expr)
    return P
def NormalVec_from_eQplane(expr):
    P,vN=get_Plane_Poin_Normal(expr)
    return vN
    
    
class MySphere():
    def __init__(self, *args,center='',radio=''):
        if center!='':
            if type(center)==list or type(center)==tuple:
                self.C=Point3D(*center)
            else:
                self.C=center
            if radio!='':
                self.R=radio
            else:
                self.R=r
            P=self.C   
            self.eQq=(x-P.x)**2+(y-P.y)**2+(z-P.z)**2-(self.R)**2    
        elif len(args)==1:
            
            expr=args[0]
            if type(expr)==MyEq:
                expr=expr.ksym
            if type(expr)==MyEqEq:
                expr=expr.L-expr.R
            self.eQq=expr
            C1,R1= completesquarexyz(expand(expr),'center',kshow=False)
            self.C=C1
            self.R=sqrt(R1)
            x1,y1,z1,rr= completesquarexyz(expand(expr),kshow=False)
            self.eQq=(x-x1)**2+(y-y1)**2+(z-z1)**2-rr
    def set(self,**kwargs):
        Se=self.eQq
        Se=real_subs(Se,**kwargs)
        self.eQq=Se
        C1=self.C
        C1=real_subs(C1,**kwargs)
        self.C=C1
        R1=self.R
        R1=real_subs(R1,**kwargs)
        self.R=R1
    
    
    @property
    def getradio(self):
        return self.R
    @property    
    def getcenter(self):
        return self.C
        
    def equation(self,*args):
        kres=unisymbols(self.eQq)
        kname=''
        op=['expand','simplify','factor']
        for i in args:
            if 'expand' in args:
                kres=expand(kres)
            if 'simplify' in args:
                kres=expand(kres)
            if 'factor' in args:
                kres=factor(kres)    
            if type(i)==str and i not in op:
                kname=i
        if kname!='':
            ee=MyEq(unisymbols(kres),kname=kname)
            return ee
        else:    
            return kres 

     
    def evaluatePoint(self,*args):
        
        if len(args)>1:
            vec=[]
            for i in args:
                vec.append(i)
            x1,y1,z1=vec
        else:    
            x1,y1,z1=Obj2List(*args)
        kres =self.eQq
        kres=kres.subs(x,x1)
        kres=kres.subs(y,y1)
        kres=kres.subs(z,z1)
        return kres
    @property     
    def Fs(self):
        return self.equation()
    def getEquation(self,kname='eQs'):
        expr=reducecero(unisymbols(self.equation()))
        ee=MyEq(expr,kname=kname)
        return ee
    def evaleq(self,**kwargs):
        ee=self.eQq
        
        kres=real_subs(ee,**kwargs)
        return kres
    def eQ(self,kshow=True):
        kres=self.eQq
        if  kshow:
            display(Math(latex(kres)))
        return kres    
        
def ProdVectorial(vec1,vec2):
    return list( NormalVec(vec1,vec2))
def pvector(vec1,vec2):
    return list( NormalVec(vec1,vec2))
def vproduc(vec1,vec2):
    return list( NormalVec(vec1,vec2))

def reducefactores(vec):
    mm=[i for i in vec if i!=0]
    mcd=gcd(mm)
    return [i/mcd for i in vec]
def vecNormal(vec1,vec2):
    return NormalVec(vec1,vec2)
def NormalVec(vec1,vec2):
    x1,y1,z1=vec1
    x2,y2,z2=vec2
    p1=y1*z2 - y2*z1
    p2=-x1*z2 + x2*z1
    p3=x1*y2 - x2*y1
    kres=[p1,p2,p3]
    kres=simplifydirection(kres)
     
    return kres

def function2Plane(expr,var=t):
    xx,yy,zz=0,0,0
    ee=MyEq(expr,'ee',kshow=False)
    vi=ee(x=0,y=0,z=0,kshow=False)
    ee.Substrac(vi,kshow=False)
    x1=ee(x=1,y=0,z=0)
    y1=ee(x=0,y=1,z=0)
    z1=ee(x=0,y=0,z=1)
    ee.Add(vi,kshow=False)
    zz=ee.solve(z,kshow=False,x=0,y=0)
     
    return MyPlane3D(Point3D(0,0,zz),[x1,y1,z1],var=var)
    
def Func2Plane(*args):
    if len(args)==1:
        return func2Plane(args[0])
    else:
        return Plane(*args)
    
def eQ2Plane(expr):
    '''
    return Plane form equation plane
    '''
    eP=MyEq(expr,'eP',kshow=False)
    xx=[0,1,0]
    yy=[0,0,1]
    zz=[]
    for i,j in zip(xx,yy):
        zz.append(eP.solve(z,kshow=False,x=i,y=j))
    vecP=[]
    for i,j,k in zip(xx,yy,zz):
        vecP.append(Point3D(i,j,k))
    P1,P2,P3=vecP
    return MyPlane3D(P1,P2,P3)


def coef0(obj):
    '''
    return coefficient zero in f(x,y)
    '''
    
    expr=obj2expr(obj)
    w=list(expr.free_symbols)
    kres=0
    for i in expr.args:
        done=True
        for j in w:
            if str(j) in str(i):
                done=False
        if done:
            kres=kres+i        
    return kres 

def squarefillxy(obj,*args):
    expr=obj2expr(obj)
    if 'factor' in args:
        expr=squarefill(expr,x,'factor')
        expr=squarefill(expr,y,'factor')
    return expr


def gintersec(obj1,obj2):
    try:
        kres=obj1.intersection(obj2)
        if type(kres)==list and len(kres)==1:
            return kres[0]
        return kres
    except:
        pass

def gobj2list(gobj):
    if type(gobj)==Array:
        return list(gobj)
    elif  type(gobj)==vFunction:
        return gobj.F
    else:
        return gobj
        
def gobj2array(gobj):
    if type(gobj)==list:
        return Array(gobj)
    elif type(gobj)==vFunction:
        return Array(gobj.F)
    else:
        return gobj

def modulo(expr,**kwargs):
    vec=Obj2List(expr)
    kres=0
    for i in vec:
        kres=kres+i*i
    kres=sqrt(kres)
    if len(kwargs)>0:
        kres=real_subs(kres,**kwargs)
    return simplify(kres)
    
def unitary(expr):
    vec=Obj2List(expr)
    kres=0
    for i in vec:
        kres=kres+i*i
    kres= rsimplify(sqrt(kres)) 
    return Obj2Array([i/kres for i in expr])

def LineXYZ(*args,var=t,kshow=True):
    vecv=[x,y,z]
    vecr=[]
    kres=[]
    for i in args:
        vecr.append(i-var)
    for i,j in zip(vecr,vecv):
         kres.append(simplesolve(i,j,'value',kshow=False ))
     
    R1=[i.subs(var,0) for i in kres]
    P1=Point3D(R1)
    R2=[i.subs(var,1) for i in kres]
    P2=Point3D(R2)
    kres=MyLine3D(P1,P2,kshow=False)
    if kshow:
        display(Math(latex(kres)))
    return kres
def simetric2Line3D(*args,var=t):
    vec=[]
    for i in args:
        vec.append(i-var)
    vecv=[x,y,z]
    for i in vecv:
        vec.append(i)
    R1=Array(simplesolve(*vec,kshow=False))
    R1=Array(R1)
    P1=R1.subs(str(var),0)
    P11=Point3D(unisymbols(P1[0]),unisymbols(P1[1]),unisymbols(P1[2]))
    R2=R1-P1
    P2=R2.subs(str(var),1)
    P22=Point3D(unisymbols(P2[0]),unisymbols(P2[1]),unisymbols(P2[2]))
    kres=Line3D(Point3D(P1),Point3D(P2))
    display(Math(latex(kres)))
    return MyLine3D(kres)
    
def SumProVec(v1,v2):
    if type(v1)==Array:
        v1=list(v1)
    if type(v2)==Array:
        v2=list(v2)
    kres=0    
    for i,j in zip(v1,v2):
        kres=kres+i*j
    return unisymbols(kres)

def Points2Array(P1,P2):
    return Array(P1-P2)


def Ldirec(Line3):
    return Array(Line3.direction)

def PointDiff(obj1,obj2):
    obj1=Array(obj1)
    obj2=Array(obj2)
    V1,V2=[],[]
    V1=Array(unisymbols(obj1))
    V2=Array(unisymbols(obj2))
    V3=V1-V2
    return list(V3)
    
def ClosestPointOnLine(a, b, p):
    ap = p-a
    ab = b-a
    result = a + dot(ap,ab)/dot(ab,ab) * ab
    return result
from sympy import sqrt
def LenghtVector(V):
    V1=list(V)
    kres=0
    for i in V1:
        kres=kres+i*i
    return sqrt(kres) 
def angle2vector(v1,v2,*args,type='cos',kname=''):
    SL=SumProVec(v1,v2)
    Lv1=LenghtVector(v1)
    Lv2=LenghtVector(v2)
    Lv=Lv1*Lv2
    if 'arg' in args:
        kres=SL/Lv
    else:
        if type=='angle':
            return SL/Lv
            
        elif type=='sin':
            kres=asin(SL/Lv)
        else:
            kres=acos(SL/Lv)
        
    if kname!='':
        return MyEq(kres,kname=kname)
    else:    
        return kres 

        
def simplifydirection(vecd):
    from numpy import lcm
    dd=[denom(i) for i in vecd]
    kres=dd[0]
    for i in range(1,len(dd)):
        kres=lcm(kres,dd[i])
        
    rr=[simplify(i*kres) for i in vecd]
    sexpr=str(rr[0])
    if sexpr[0]=='-':
        rr2=[-1*i for i in rr]
        rr=rr2

    kres=rr[0]
    for i in range(1,len(rr)):
        kres=np.gcd(kres,rr[i])
    rr2=[simplify(i/kres) for i in rr] 
    
    
    
    return rr2 

def Points2Vector(P1,P2):
    kres=P2-P1
    kres=list(kres)
    display(Math(latex(kres)))
    return list(kres)

def Projection(self,Obj):
    if type(self)==MyVector or type(self)==Array:
        V1=MyVector(Obj2Array(self),kshow=False)
        V2=MyVector(Obj2Array(Obj),kshow=False)
        Ca=V1.cosangle(V2)
        Mo=V1.modulo
        kres=Ca*Mo*V2.unitary()
        return kres
    else:
        if type(Obj)==MyLine3D:
            PP=Obj.arbitrarypoint(kshow=False)
            P1=Array(PP(0))
            P2=Array(PP(1))
            P3=Array(self.P)
            kres=ClosestPointOnLine(P1, P2, P3)
            return MyPoint3D(kres)
        else:
            
            Vo=list(Obj.ortogonal())
           
            P=tuple(self.Obj)
            L=MyLine3D(P,Vo,kshow=False)
            PL2=Obj.Obj
            L2=L.Obj
             
            kres= L2.intersection(PL2)
            if len(kres)==1:
                kres=kres[0]
                return MyPoint3D(kres)
            else:
                kres2=[MyPoint3D(i) for i in kres]
                return kres2

def subsArray(A,*args,var=t,**kwargs):
    args2=[i for i in args if type(i)!=str]
    if len(args2)==1:
        kres=A.subs(var,args2[0])
        return kres
    elif len(args2)==0 and len(kwargs)>0:
        Lkres=list(A)
        kres=Array([real_subs(i,**kwargs) for i in Lkres])
        return kres
    else:
        return A

def modulovector(expr):
    L=list(expr)
    k=[i*i for i in L]
    Ss=sum(k)
    return sqrt(Ss)

def cos2vectors(V1,V2):
    U1=SumProVec(V1,V2)
    U2=modulovector(V1)
    U3=modulovector(V2)
    U4=U2*U3
    return U1/U4

def modifyArray(Obj,var,*args,**kwargs):
    obj=Obj2List(Obj)
    if len(args)==1:
        expr=args[0]
        expr=obj2expr(expr)
        obj=[i.subs(str(var),expr) for i in obj]
    if len(kwargs)>0:
        obj=[real_subs(i,**kwargs) for i in obj]
    return obj 

def creadiffvar(var):
    if var==alpha:
        return dalpha
        
    svar=str(var)
    sdvar='d'+svar
    return symbols(sdvar)
def creteIntegral(expr,var,x1='',x2=''):
    expr=obj2expr(expr)
    dvar=creadiffvar(var)
    
    try:
        expr=expr.subs(dvar,1)
    except:
        pass
    if x1=='': 
        kres=Integral(expr,var)
    else:
        kres=Integral(expr, (var,x1,x2))
    return kres


def crossProduct(p1,p2):
    if type(p1)!=Array:
        p1=Obj2Array(p1)
    if type(p2)!=Array:
        p2=Obj2Array(p2)
    return Obj2Array(cross(p1,p2))
    
def comparePoints(p1,p2):
    return comparevec(p1,p2)

    
def comparevec(v1,v2):
    V1=Obj2Array(v1)
    V2=Obj2Array(v2)
    V3=V1-V2
    V4=list(V3)
    vecvar=[]
    for i in V4:
        kres=list(i.free_symbols)
        for j in kres:
            if j not in vecvar:
                vecvar.append(j)
    L1=len(vecvar)
    L2=len(V4)
    if L1==L2:
        vec=V4+vecvar
        return simplesolve(*vec)
    else:
        for i in range(0,L1):
            vec2=V4[i:i+L1]
            try:
                vec=vec2+vecvar
                kres=simplesolve(*vec,'value')
                return kres
                break
            except:
                pass
       
def MyGeoQ(Obj,kname):
    if type(Obj)==MyPlane3D:
        return MyEq(Obj.eQ(kshow=False),kname)
    elif type(Obj)==MySphere:
        return MyEq(Obj.eQ(kshow=False),kname)
    else:
        return Obj

class MySurface(MyEq):
    def __init__(self, *args,kshow=True):
        kname=''
        for i in args:
            if type(i)==str:
                kname=i
            if type(i)!=str:
                self.Obj=i
        self.name=kname        
        self.ee=MyEq(self.Obj,kname)
    def __call__(self,**kwargs):
        ksym=self.Obj
        if len(kwargs)>0:
            ksym=real_subs(ksym,**kwargs)
        return ksym
    def __repr__(self):
        kres = str(self.Obj)
        return kres
    def _latex(self, obj):
        return latex(self.Obj)

    def __str__(self):
        return self.__repr__()        
        
    def vecdiff(self,**kwargs):
        ksym=self.Obj
        kres=Obj2Matrix((ksym.diff(x,**kwargs),ksym.diff(y,**kwargs),ksym.diff(z,**kwargs)))
        
    def vecNormal(self,**kwargs):
        return Func2diffVector(self.Obj,**kwargs)
    


        
def Func2diffVector(eQ,**kwargs):
    kres=MyVector((eQ.diff(x,**kwargs),eQ.diff(y,**kwargs),eQ.diff(z,**kwargs)),kshow=False)
    return kres


def coef2Line2D(a,b,c,xx1=0,xx2=1):
    '''
        L=ax+by+c
        coef2Line2D(a,b,c)
        return L whit (a,b,c) coef line
        '''
    y=cfrac(-x*a,b)-cfrac(c,a)
    x1=xx1
    y1=y.subs(x,x1)
    P1=Point(x1,y1)
    x2=xx2
    y2=y.subs(x,x2)
    P2=Point(x2,y2)
    L=Line2D(P1,P2)
    return L


######### 2024 
def Eq2Plane(expr,*args):
    e1=MyEq(expr,'e1',kshow=False)

    if not 'x' in str(expr):
        y1=e1.simplesolve(y,z=0)
        z1=e1.simplesolve(z,y=0)
        P1=Point(0,y1,0)
        P2=Point(0,0,z1)
        P3=Point(1,y1,0)
    
          
    elif not 'y' in str(expr):
        x1=e1.simplesolve(x,z=0)
        z1=e1.simplesolve(z,x=0)
        P1=Point3D(x1,0,0)
        P2=Point3D(0,0,z1)
        P3=Point3D(x1,1,0) 
 
         
    elif not 'z' in str(expr):
        x1=e1.simplesolve(x,y=0)
        y1=e1.simplesolve(y,x=0)
        P1=Point3D(x1,0,0)
        P2=Point3D(0,y1,0)
        P3=Point3D(0,y1,1)   
          
    else:
        x1=e1.simplesolve(x,y=0,z=0)
        y1=e1.simplesolve(y,x=0,z=0)
        z1=e1.simplesolve(z,x=0,y=0)
        P1=Point3D(x1,0,0)
        P2=Point3D(0,y1,0)
        P3=Point3D(0,0,z1)
    if 'dots' in args:
        return P1,P2,P3
    else:
        return Plane(P1,P2,P3)