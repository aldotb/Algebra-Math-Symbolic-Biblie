from sympy import *
 
from IPython.display import Math  # ,display
from matplotlib.pyplot import ylabel, plot, show, xlabel, title
from libaldo_math2 import *
from libaldo_algorith import *
 

 
from libaldo_show import *
from lib_tools import *
from lib_MyEq import *  
from lib_MyEqEq import *
import copy 
x,t,k,y=symbols('x t k y')

inesymbol=['<','<=','=','>=','>']
inesymbol2=['<','≤','=','≥','>']
C1,C2,C3,C4,C5=symbols('C1 C2 C3 C4 C5')

 
 
class MyCelule():
    def __init__(self, *args,kshow=True):
    op=['MQ','EQ']
    kname=''

    for i in args:
        if type(i)==str and i not in op:
            kname=i
    self.p1=args[0]
    self.p2=''    
    if kname!='':         
        self.p2=args[1]
        self.expr=MQ(p1,p2)
        self.type='MQ'
    else:
        self.expr=MyEq(p1,kname)
        self.type='EQ'
    
     