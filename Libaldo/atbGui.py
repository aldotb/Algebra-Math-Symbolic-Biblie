import ipywidgets as widgets
from ipywidgets import HBox, VBox
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display


import sys
sys.path.append('imghelp')
from sympy import *
#from advancemathlib import * #this cell will be initializing your variables
from helplib import helplib # load this for help , some time run two time the same cell to HQ view
from MyDiff  import * # Ode Lib
from sympy.parsing.mathematica import mathematica
from ttools import *


def newimg(mfile):
    path="Libaldo/imghelp/"
    sfile=path+mfile
    file = open(sfile, "rb")
    image = file.read()
    imggui=widgets.Image(value=image,format='png',width=800,height=400)
    return imggui    

###
from ipywidgets import Layout, Button, Box, FloatText, Textarea, Dropdown, Label, IntSlider

form_item_layout = Layout(
    display='flex',
    flex_flow='row',
    justify_content='space-between')
    
    
    
input1 = widgets.Text( description='Math expr:')
#label1 = widgets.Label(value='Pow' ,height='200px')
input2= widgets.Text( description='Pow')  
bb   = widgets.Button( description='done' )
out1  = widgets.Label(value='out')
mfoot1=widgets.HBox([input1,input2,bb],layout=form_item_layout)
mfoot2=widgets.HBox([out1],layout=form_item_layout)
vb=widgets.VBox([mfoot1,mfoot2])


def on_click_compute_stress(b):

    sbase=input1.value
    base=eval(sbase)
    sbase2=str(base)
    expo=input2.value
    sexpr='('+sbase2+')**('+expo+')'
     
    sexpr2="$"+latex(expand(parse_expr(sexpr)))+"$"
    out1.value= sexpr2
bb.on_click(on_click_compute_stress)
    
 

###
 
F00=['MyEq.png','MyEqset.png','MyEqEq.png']
T00=['MyEq','MyEq set','MyEqEq or MQ']
S00=[newimg(file) for file in F00]
TAB00 =[ VBox(children=[i]) for i in S00]
tab0 = widgets.Tab(children=TAB00) 
for i in range(len(S00)):
    tab0.set_title(i, T00[i])

F11=['Add.png','cfrac.png',('Pow.png',vb),'rpow.png']
T11=['Basic', 'Fraction','Pow','Root']
S11=[]
for i in F11:

    if type(i)==str:
        S11.append(newimg(i))
    
    else:
        vvbox=widgets.VBox([newimg(i[0]),i[1]])
        S11.append(vvbox)

        
        
    #S11.append(
    #S11=[newimg(file) for file in F11]
TAB11 =[ VBox(children=[i]) for i in S11]
tab1 = widgets.Tab(children=TAB11) 
for i in range(len(S11)):
    tab1.set_title(i, T11[i])


F22=['expand.png','factor.png','rsimplify.png','factorize.png','texpand.png','tsimplify.png']
T22=['Expand', 'Factor','Rsimplify','Factorize','Trig Expand','Trig Simplif' ]
S22=[newimg(file) for file in F22]
TAB22 =[ VBox(children=[i]) for i in S22]
tab2 = widgets.Tab(children=TAB22) 
for i in range(len(S22)):
    tab2.set_title(i, T22[i])


F33=['getbase.png','expandbase.png','insidepow.png' ]
T33=['Get Base/Expo','Exd,Fact (Base,Exponent)', 'Inside Root/Pow' ]
S33=[newimg(file) for file in F33]
TAB33 =[ VBox(children=[i]) for i in S33]
tab3 = widgets.Tab(children=TAB33) 
for i in range(len(S33)):
    tab3.set_title(i, T33[i])
    
F44=['unfog.png','transforma.png','squarecomplete.png','primefactor.png' ]
T44=['descomp func','trans Poly','SqrtComplete','Prime facts' ]
S44=[newimg(file) for file in F44]
TAB44 =[ VBox(children=[i]) for i in S44]
tab4 = widgets.Tab(children=TAB44) 
for i in range(len(S44)):
    tab4.set_title(i, T44[i]) 
    
F55=['getexpo.png','simplifyexpo.png','pow2powpow.png' ]
T55=['Basic Expon','classic','transf exponent' ]
S55=[newimg(file) for file in F55]
TAB55 =[ VBox(children=[i]) for i in S55]
tab5 = widgets.Tab(children=TAB55) 
for i in range(len(S55)):
    tab5.set_title(i, T55[i])     

    
##  MATRIX
F66=['MyMat1.png','MyMat2.png','MyMat3.png','MyMat4.png','MyMat5.png','MyMat6.png']
T66=['MyMat ','Basic Mat ','Inv,Det,Trans ','Simplesolve ','Matrix Tools ','Symb Matrix']
S66=[newimg(file) for file in F66]
TAB66 =[ VBox(children=[i]) for i in S66]
tab6 = widgets.Tab(children=TAB66) 
for i in range(len(S66)):
    tab6.set_title(i, T66[i])
    
F77=['factorinte.png']
T77=['Integer Factor']
S77=[newimg(file) for file in F77]
TAB77 =[ VBox(children=[i]) for i in S77]
tab7 = widgets.Tab(children=TAB77) 
for i in range(len(S77)):
    tab7.set_title(i, T77[i])     

   
F88=['complex.png','myplot.png']
T88=['Complex Num','Plot Funct']
S88=[newimg(file) for file in F88]
TAB88 =[ VBox(children=[i]) for i in S88]
tab8 = widgets.Tab(children=TAB88) 
for i in range(len(S88)):
    tab8.set_title(i, T88[i])     
    
helpgui = widgets.Tab(children=[tab0,tab1, tab2, tab3, tab4, tab5, tab6,tab7,tab8])
helpgui.set_title(0,'Class MyEq')
helpgui.set_title(1,'Basic')
helpgui.set_title(2,'Transformat')
helpgui.set_title(3,'Trans Monomie')
helpgui.set_title(4,'Polynomies')
helpgui.set_title(5,'Exponents')
helpgui.set_title(6,'Matrix') 
helpgui.set_title(7,'Diferential')
helpgui.set_title(8,'More..')  


##  MATRIX
MF00=['MyMat1.png','MyMat2.png','MyMat3.png','MyMat4.png','MyMat5.png','MyMat6.png']
MT00=['MyMat','Basic Mat ','Inv,Det,Trans','Simplesolve ','Matrix Tools ','Symb Matrix']
MS00=[newimg(file) for file in MF00]
MTAB00 =[ VBox(children=[i]) for i in MS00]
MatTab = widgets.Tab(children=MTAB00) 
for i in range(len(MS00)):
    MatTab.set_title(i, MT00[i])

def helpMyMat():
    display(MatTab)
    

 


    
def helplib(kimg=''):
    if kimg=='':
        return get_filesh()
    else:
        try:
            kfile=kimg+'.png'
            camino='Libaldo/imghelp/'+kfile
            baseA = np.asarray(Image.open(camino))
            plt.imshow(baseA)
            plt.axis('off')
            plt.gcf().set_facecolor("black")
            plt.rcParams['figure.dpi'] = 800
            plt.show()
        except:
            print( "don't find function, please type helplib() for more options")
    


def get_filesh():
    res = []
    for path in os.listdir('Libaldo/imghelp/'):
        res.append(path)
    mm=''
    for i in res:
        kk=i.replace('.png','(), ')
        mm=mm+kk
    return mm[1::] 
