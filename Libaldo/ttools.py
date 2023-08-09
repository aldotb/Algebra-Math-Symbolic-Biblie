# -*- coding: utf-8 -*-
"""
Created on Wed May  3 10:47:45 2023

@author: aldot
"""


import sys, os, requests, latex 
from sympy import *
from sympy.printing.mathml import mathml
from PIL import Image
from IPython.display import Image, display
from IPython.lib.latextools import latex_to_png
import cv2  
import numpy as np
from sympy import preview
from io import BytesIO
from PIL import Image

import sys
from lib_simpleintegral import *
 

init_printing(use_unicode=True) # allow LaTeX printing



import numpy as np
import matplotlib.pyplot as plt

def implot(img):
    if type(img)==str:
        img = cv2.imread(img)
    plt.imshow(img)
    plt.show()
    
def Pil2cv2img(img):
    return np.asarray(img)[:,:,::-1].copy()

def Cv2toPilimg(img):
    return img.asarray(im.convert('L'))

def imgnormal(img):
    if type(img)==str:
        return cv2.imread(img)
    else:
        return Pil2cv2img(img)
    
    
def imshow(img):
    if type(img)==str:
        img = cv2.imread(img)
    plt.figure(figsize=(10,1))
    plt.axis('off')
    plt.imshow(img)

    plt.show()     

def formula2file( formula, file, negate=False ):
    if type(formula)!=str:
        formula=latex(formula,mode="inline")
    tfile = file
    if negate:
        tfile = 'tmp.png'
    r = requests.get( 'http://latex.codecogs.com/png.latex?\dpi{300} \huge %s' % formula )
    f = open( tfile, 'wb' )
    f.write( r.content )
    f.close()
    if negate:
        os.system( 'convert tmp.png -channel RGB -negate -colorspace rgb %s' %file )
        
        
        
def testformula():
    
    formula2file((x*x+3*x*y+y*y)*dx-x*x*dy,'formula1.png')
    formula2file(-u*x**2*dx + du*x + dx*(u**2*x**2 + 3*u*x**2 + x**2),'formula2.png')
    formula2file(Eq(-u*x**2*dx + du*x + dx*(u**2*x**2 + 3*u*x**2 + x**2),6*z*(dx*(u**2*x**2 + 3*u*x**2 + x**2))),'formula3.png')
    formula2file(Integral(du/(x**2*(4*u + 1)), x) + Integral(-dx*(u**2 + 2*u + 1)/(x*(4*u + 1)), u),'formula4.png')   
    formula2file((x*x+3*x*y+y*y)*dx-x*x*dy+z*(-u*x**2*dx + du*x + dx*(u**2*x**2 + 3*u*x**2 + x**2)),'formula5.png')
    formula2file(z**((y + z)**((y + z)**z)),'formula6.png')




def whiteimg(yy,xx):
    img = np.zeros([yy,xx,3],dtype=np.uint8)
    img.fill(255) # or img[:] = 255
    return img



def appliedbode(img):
    if type(img)==str:
        img=cv2.imread(img)
    yy,xx,zz=img.shape
    yy2=int(yy/8)
    nxx=2*yy2+xx
    nyy=2*yy2+yy
    imgempty=whiteimg(nyy,nxx)
    x2=nxx-yy2
    x1=yy2
    y2=yy+yy2
    y1=yy2
    imgempty[y1:y2,x1:x2,:]=img
    return imgempty

    
    
    


def tranimg2file(img,kname):
    img=imgnormal(img)
    img1=appliedbode(img)
    img2=completeimg(img1)
    fimg=cv2.imwrite(kname,img2)
    
    


def getshape(img):
    if type(img)==str:
        img=cv2.imread(img)
    return get_shape_img2d(img)    

def get_shape_img2d(img):
    img=imgnormal(img)
    yy,xx,zz=img.shape
    return xx,yy
    

def imsave(img,kname=''):
     
    img=imgnormal(img)
    cv2.imwrite(kname,img)
     
        
    
def getminimos():
    img = cv2.imread('wfor.png')
    xx,yy=get_shape_img2d(img)
    xx2=min(xx,1000)
    yy2=min(yy,50)
    return (xx2,yy2) 

def showMathExpr(expr,dpi=200):
    img=MySympy2Img(expr,dpi=dpi)
    imshow(img)

def MySympy2Img(expr,dpi=200):
    obj = BytesIO()
    preview(expr, output='png', viewer='BytesIO',outputbuffer=obj,euler=False, dvioptions=['-D',str(dpi)])
    new = Image.open(obj)
    return new

def resobresize(img):
    img=imgnormal(img)
    xx,yy=get_shape_img2d(img)
    XX=580
    YY=140
    fx=XX/xx
    fy=YY/yy
    reducec = False
    if fy<1 or fx<1:
        reducec = True
        ff=fx
        if fy<fx:
            ff = fy
        nxx=int(xx*ff)
        nyy=int(yy*ff)
        img2=resizeimg(img,W=nxx,H=nyy)
        return (img2,nxx,nyy)
    else:
        return (img,xx,yy)


def resizeimg(*args,W='',H=''):
    '''
    input:
        img,factor
        img, W value 
         
        img,H,vale
        img,Wvale,Hvale

    return 
        img resized

    '''
    img=args[0]
    img=imgnormal(img)
    #ratio=1
    factor=0
    x,y=get_shape_img2d(img)
    dim=(x,y)
     
    if len(args)==2 and  W=='' and H=='':
        
        
        x,y=get_shape_img2d(img)
        factor=args[1]
        if factor<3:
            factor=int(factor*100)
        ny=int(y*factor/100)
        nx=int(x*factor/100)
        dim=(nx,ny)
        img2=cv2.resize(img, dim, interpolation = cv2.INTER_LINEAR)  
         
        return img2
    if len(args)==1 :
        
        nx,ny= get_shape_img2d(img)
        if W!='':
            nx=W
        if H!='':
            ny=H

        dim=(nx,ny)
        img2=cv2.resize(img, dim, interpolation = cv2.INTER_LINEAR)  
         
        return img2


# def generaf():
#     vecname=['gra'+str(i)+'.png' for i in range(6)]
#     A=MyEq(x*x+y*y,'A',kshow=False)
#     r=[]
#     r.append(A())    
#     r.append(A()**A())    
#     r.append(A()*12345678912345698)    
#     r.append(A()**A()**A())    
#     r.append((A()**A()**A())/((A()**A()**A())+1))    
#     r.append(r[-1]**r[-1])    
#     return vecname  