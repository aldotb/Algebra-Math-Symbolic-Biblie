from sympy import *
from lib_Mathematica import *
from lib_Mathbasic import *
from lib_Algorith import *

from lib_MyEq import *
from lib_MyEqEq import *
 
import copy
def obj2mathexpr(*args):
    vecres=[]
    for i in args:
        if type(i)==MyEq:
            vecres.append(i.ksym)
        elif type(i)==MyEqEq:
            vecres.append(i.e1.ksym-i.e2.ksym) 
        else:
            vecres.append(i)
    if len(vecres)==1:
        return vecres[0]
    else:
        return vecres
        
        
def getdataQ(*args):
    ops=['float','value','eQ','equation']
    vecexp=[]
    vecvar=[]
    vecnam=[]
    vecslist=[]
    args2=[]
    for data in args:
        if type(data)==list:
            for data2 in data:
                args2.append(data2)
        else:
            args2.append(data)
    args=args2        
    for data in args:
        if Is_String(data):
            if not data in ops:
                vecnam.append(data)
        elif Is_Symbol(data):
            vecvar.append(data)
            vecnam.append(str(data))
        elif type(data)==list:
            for  item in data:
                vecexp.append(item)        
        else:
            vecp=list(data.free_symbols)
            if vecp!=[]:
                if not vecp in vecslist:
                    vecslist.append(vecp)
                    vecexp.append(obj2mathexpr(data))
    return vecexp,vecvar,vecnam
    
def predataQ(*args):
    vecexp,vecvar,vecnam=getdataQ(*args)
    if vecvar==[]:
        vecvar2=[]
        for data in vecexp:
            vres=list(data.free_symbols)
            for qvar in vres:
                if not qvar in vecvar2:
                    vecvar2.append(unisymbols(qvar))
        vecvar=vecvar2
    for data in vecvar:
        vecnam.append(str(data))
    return vecexp,vecvar,vecnam

def ndegree(vecexp,vecvar):
    etype=1
    if len(vecvar)>1:
        try:
            var=vecvar[0]
            lvar=vecvar[1::]
            qq=len(vecvar)
            for nv in lvar:
                nvec=[]
                for data in vecexp:
                    nvec.append(data.subs(nv,var))
                vecexp=nvec
        except:
            pass
    for dataq in vecexp:
        for datav in vecvar:
            kres=degree(dataq,datav)
            if kres>etype:
                etype=kres
    return etype 

    

def simplesolve2(*args,kshow=True,**kwargs):
    vecvar=[]
    vecexp=[]
    vecnam=[]
    qq2=0
    ops=['float','value','eq','equation']
    for data in args:
        if not data in ops:
            if Is_Symbol(data):
                vecvar.append(data)
            elif type(data)==list:
                for  item in data:
                    vecexp.append(item)        
            else:
                vecexp.append(obj2mathexpr(data))
    for data in vecvar:
         vecnam.append(str(data))
    if len(kwargs)>0:
        vecexp2=[]
        for data in vecexp:
            vecexp2.append(real_subs(data,**kwargs))
        vecexp=vecexp2    
    if vecvar==[]:
        for data in vecexp:
            vecv=list(data.free_symbols)
            for dvar in vecv:
                if not dvar in vecvar:
                    vecvar.append(dvar)
        vecname=[str(data) for data in vecvar]            

    kres=solve(vecexp,vecvar)
     
    if type(kres)==dict: 
        sval,vvar=list(kres.keys()),list(kres.values())
        sname=[str(data) for data in sval]
        
    elif type(kres)==list:    
        if type(kres[0])==tuple:
            if len(kres)>1:
                qq=len(kres)
                qq2=len(list(kres[0]))
            else:
                vvar=kres[0]
                sval=vecvar[0]
                sname=str(sval)
                qq2=0
        else:
            if len(kres)==1:
                vvar=kres[0]
                sval=vecvar[0]
                sname=str(sval)
                qq2=0
        
    if qq2>0:
        nvecval=[]
        nvecname=[]
        for i in range(0,qq2):
            for cc in range(qq):
                lvalor=list(kres[i])
                nvecval.append(lvalor[i])
                nvecname.append(vecname[i]+str(i+1))
        sname=nvecname
        vvar=nvecval
 
    if 'float' in args:
        vecflo=[]
        for data in vvar:
            try:
                vecflo.append(float(data))
            except:
                vecflo.append(data)
        vvar=vecflo        
    veceq=[]
    for kname,kval in zip(sname,vvar):
        e1=MyEq(kval,kname)
        veceq.append(e1)
    if  'eq' in args or 'equation' in args:
        return veceq
    else:
        if type(vvar)==list and len(vvar)==1:
            return  vvar[0]
        else:
            return vvar     


def kunpack(**kwargs):
    sval=[]
    vval=[]
    for key, value in kwargs.items():
        sval.append(key)
        vval.append(value)
    return sval,vval    
        
    
def completesquare(expr,var=x):
    '''
    expr=x*x+2*x
    if x*x+2*x+a=(x+b)**2
    return a,b
    ''' 
    a, b=symbols('a b')
     
    vec1=coef_list(var)
    vec2=[]
    for i in vec1:
        if 'a' in str(i) or 'b' in str(i):
            vec2.append(i)
    a,b=simplesolve(*vec2,a,b)
    return -b,a
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
def squarecompletexy(*args,var1=x,var2=y):
    var1=symbols(str(var1))
    var2=symbols(str(var2))
    exprx=args[0]
    var=var1
    Mx=exprx.subs(var2,0)
    Ry=exprx-Mx
    Mx2=squarecomplete(Mx,x)
    R=coef0(Mx2)
    P1=Mx2-R
    My=Ry+R
    var=y
    My2=squarecomplete(My,var)
    R=coef0(My2)
    P2=My2-R
    if 'parts' in args:
        return P1,P2,R
    else:
        return P1+P2+R

def squarefill(*args):
    '''
    squarefill(x*x+2*x+1,x)  return (x+1)**2
    squarefill(2*x*x+4*x+1,x,'factor') return 2*(x+1)**2-2    factor
    squarefill(2*x*x+4*x+1,x,'factor',parts) return (x+1)**2,-2,2     
 
    '''
    var=x
    vecop=['factor','parts']
    fdone=False
    cf=1
    nexpr=args[0]
    for i in args:
        if type(i)==MyEq:
            expr=i.ksym
        elif type(i)==Symbol:
            var=i
        elif i=='factor':
            fdone=True 
        elif type(i)!=str:
            expr=i
        else:
           pass  
    sresto=0
    w=list(expr.free_symbols)
    if len(w)>1:
        expr8=expr
        for i in w:
            if i!=var:
                expr8=expr8.subs(i,0)
        nexpr=expr8        
        sresto=expr-nexpr
        
    x1,x2,x3=symbols('x1,x2,x3')
    vecx=[x1,x2,x3]
    
    if fdone:
         
        mm=nexpr.args
        mm2=[i.subs(var,j) for i,j in zip(mm,vecx)]
        expr2=factor(sum(mm2))
        expr3=reducecero(expr2)
        cf=abs(simplify(expr2/expr3))
        nexpr=sum(mm)/cf
    clist=coef_list(nexpr,var)
    a=sqrt(clist[0])
    b=(clist[1])/(2*a)
    c=b*b
    d=clist[2]-c
    
    if 'parts' in args:
        return (a*var+c),d*cf+sresto,cf
    else:    
        if fdone:
            return  cf*(a*var+c)**2+d*cf +sresto
        else:
            return (a*var+c)**2+d +sresto


def square2polar(obj,var=t,kshow=True):
    expr=obj2expr(obj)
    mm=expr.args
    mx=sum([i for i in mm if 'x' in str(i)])
    my=sum([i for i in mm if 'y' in str(i)])
    mz=[mx,my]
    Nm=[numer(i) for i in mz]
    Dn=[denom(i) for i in mz]
    R=[rsimplify(sqrt(i)) for i in Nm]
    ss=[Nm[0].subs(x,0),Nm[1].subs(y,0)]
    kres=[unisymbols(ss[0]+sqrt(Dn[0])*cos(var)),unisymbols(ss[1]+sqrt(Dn[1])*sin(var))]
    if kshow:
        display(Math(latex(kres)))
    return kres    
   
def completesquarexyz(QQQ,*args,var=t,kshow=True):
    '''
    input MyEq=A*x*x+B*y*y+C*z*z+D*x+E*y+F*z*G
    return Cx,Cy,Cx,R2
    '''

    QQ=copy.deepcopy(QQQ)
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
        
def obj2float(expr):
    if type(expr)==Add:
        kres=0
        for data in expr.args:
            kres+=expr2float(data)
        return kres
    elif type(expr)==Mul:
        kres=1
        for data in expr.args:
            kres=kres*expr2float(data)
        return kres    
    else:
        try:
            kres=float(expr)
        except:
            kres=expr
        return kres

# matematic expr to sympy srags map symbols to analize equation structure

# input: 3*x+3 outout : 'A(M(I,S),I)'
# input: x**3+x+3*log(y)  outout : 'A(P(S,I),S,M(I,L(S)))'
def math2map(expr):
    mm=srepr(expr)
    done=True
    cc=0
    while done:
        try:
                
            p1,p2=find_simetric_end_par(mm,cc)
            if not'(' in mm[p1::]:
                done=False
            elif p2-p1<6 or onlynumber(mm[p1:p2]):
                mm=mm.replace(mm[p1:p2+1],'*')
                cc=cc+1
            else:
                cc=p1+1
        
        except:
            done=False
    vecs=['Symbol*','Integer*','Mul','Add','Pow','Rational*',' ','log']
    vece=['S','I','M','A','P','R','','L']
    for v1,v2 in zip(vecs,vece):
        qq=mm.count(v1)
        for i in range(qq):
            mm=mm.replace(v1,v2)
    return mm

def onlynumber(sexpr):
    vec='0123456789, ()'
    done=True
    for data in sexpr:
        if not data in  vec:
            return False
    return True

def find_simetric_end_par(expr,p1):
    sexpr=str(expr)
    p2=find_first_par_from(str(expr),p1)
    done=True
    qbal=1
    p3=p2+1
    while done:
        if sexpr[p3]=='(':
            qbal=qbal+1
        if sexpr[p3]==')':
            qbal=qbal-1
        if qbal==0:
            return p2,p3
        else:
            p3=p3+1
def find_first_par_from(sexpr,p1):
    sres='('
    done=True
    while done:
        svar=sexpr[p1]
        if svar=='(':
            return p1
        else:
            p1=p1+1  


def solve_coefcompare(obj,var):
    p1=obj.e1.ksym
    p2=obj.e2.ksym
    L1=coef_list(p1,var)
    L2=coef_list(p2,var)
    vece=[]
    for data1,data2 in zip(L1,L2):
        vece.append(data1-data2)  
    return simplesolve(*vece)


def simplesolve(*args):
    
    vecexpr,vecvar=eQandSymbols(*args)
     
    vecname,vecval,qq,numvar=eqpresolve(vecexpr,vecvar)
    if 'factor' in args:
        vecval2=[]
        for data in vecval:
            try:
                vecval2.append(factor(data))
            except:
                vecval2.append(data)
        vecval=vecval2
    if 'simplify' in args:
        vecval2=[]
        for data in vecval:
            try:
                vecval2.append(simplify(data))
            except:
                vecval2.append(data)
        vecval=vecval2        
        vecval=[factor(data) for data in vecval]
    if 'eQ' in args:
        return resolve2Eq(vecname,vecval,qq,numvar)
    else:
        showresolve2Eq(vecname,vecval,qq,numvar)
        return vecval
    
def eQandSymbols(*args):
    vece1=[]
    vece2=[]
    vece3=[]
    vece4=[]
    vecsymbols=[]
    for data in args:
        if Is_Symbol(data):
            vecsymbols.append(data)
            
    for data in args:
        if type(data)!=str:
            if type(data)!=Symbol:
                vece1.append(data)
    for data in vece1:
        if type(data)==list:
            for data2 in data:
                vece2.append(data2)
        else:
            vece2.append(data)
    vecgrup=[]        
    for data in vece2:
        if type(data)==MyEqEq:
            vece3.append(data.e1.ksym-data.e2.ksym)
        elif type(data)==MyEq:
            vece3.append(data.ksym)
        else:
            vece3.append(data)
    for data in vece3:
        try:
            veri=list(data.free_symbols)
            if len(veri)>0:
                if not data in vece4:
                    vecgrup.append(veri)
                    vece4.append(data)
        except:
            pass
    if vecsymbols==[]:
        vecvar=[]
        for data in vece4:
            try:
                vecv=list(data.free_symbols)
                for sdata in vecv:
                    if not sdata in vecvar:  
                        vecvar.append(sdata)
            except:
                pass
        vecvar.sort(key=str)
        vecsymbols=vecvar
    return vece4,vecsymbols 
    
def eqpresolve(vecexpr,vecvar):
    kres=solve(vecexpr,vecvar,dict=True)
    numvar=len(vecvar)
    if type(kres)==dict:
        svar,ssol=list(kres.keys()),list(kres.values())
        snam=[str(data) for data in svar]
        knames=snam
        kvalue=ssol
        qq=1         
    elif type(kres)==list:
        ssvar,sssol,ssnam=[],[],[]
        qq=len(kres)
        cc=1
        for data in kres:
            svar,ssol=list(data.keys()),list(data.values())
            snam=[str(data) for data in svar]
            ssvar=ssvar+svar
            sssol=sssol+ssol
            ssnam=ssnam+snam
            cc+=1
        knames=ssnam
        kvalue=sssol
 
    else:
        knames=[str(data) for data in vecvar]
        kvalue=kres
        qq=1
    return knames,kvalue,qq,numvar    
    
def resolve2Eq(vecname,vecval,qq,numvar):
    veceq=[]
    if qq==1:
        for kname,kval in zip(vecname,vecval):
            veceq.append(MyEq(kval,kname))
        return veceq
    else:
        cont=0
        for i in range(numvar):
            for ss in range(1,qq+1):
                veceq.append(MyEq(vecval[cont],vecname[cont]+str(i+1)))
                cont+=1
        return veceq

def showresolve2Eq(vecname,vecval,qq,numvar):
    veceq=[]
    if qq==1:
        for kname,kval in zip(vecname,vecval):
            veceq.append(MyEq(kval,kname))
         
    else:
        cont=0
        for i in range(numvar):
            for ss in range(1,qq+1):
                veceq.append(MyEq(vecval[cont],vecname[cont]+str(i+1)))
                cont+=1
                  