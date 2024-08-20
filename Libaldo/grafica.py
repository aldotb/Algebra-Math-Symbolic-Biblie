# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 09:48:01 2024

@author: aldotb
"""
 
import numpy as np
from PIL import Image
import IPython.display
 

 
from IPython.display import clear_output  #  clear_output(wait=True)
 

# this cell will be initializing your variables
#from lib_MyFunctions import *
 

XX = np.asarray(Image.open('xfondo.png'))
YY = np.asarray(Image.open('yfondo.png'))
X0 = np.asarray(Image.open('D:/Libaldo/x0.png'))
X1 = np.asarray(Image.open('D:/Libaldo/x1.png'))
X2 = np.asarray(Image.open('D:/Libaldo/x2.png'))
X3 = np.asarray(Image.open('D:/Libaldo/x3.png'))
X4 = np.asarray(Image.open('D:/Libaldo/x4.png'))
X5 = np.asarray(Image.open('D:/Libaldo/x5.png'))
X6 = np.asarray(Image.open('D:/Libaldo/x6.png'))
X7 = np.asarray(Image.open('D:/Libaldo/x7.png'))  

Y0 = np.asarray(Image.open('D:/Libaldo/y0.png'))
Y1 = np.asarray(Image.open('D:/Libaldo/y1.png'))
Y2 = np.asarray(Image.open('D:/Libaldo/y2.png'))
Y3 = np.asarray(Image.open('D:/Libaldo/y3.png'))

Gpos=[X0,X1,X2,X3,X4,X5,X6,X7]
Ypos=[Y0,Y1,Y2,Y3]

def vecset2img(xlist):
    ndata=[int(str(data)[1]) for data in xlist]
    xf=XX
    for data in ndata:
        xf=xf+Gpos[data]
    return xf 

def vecset2img0(xlist):
    ndata=[int(str(data)[1]) for data in xlist]
    xf=YY
    for data in ndata:
        xf=xf+Ypos[data]
    return xf    
    
def render3(obj):
    xf=vecset2img(obj)
    img=Image.fromarray(xf)
    img.resize(size=(300,300))
    IPython.display.display(img)
    
def render(obj):
    xf=vecset2img0(obj)
    img=Image.fromarray(xf)
    img.resize(size=(300,300))
    IPython.display.display(img) 
    
    
    