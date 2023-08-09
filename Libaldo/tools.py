# -*- coding: utf-8 -*-
"""
Created on Wed May  3 10:47:45 2023

@author: aldot
"""


import sys, os, requests, latex 
from sympy import *
from sympy.printing.mathml import mathml
init_printing(use_unicode=True) # allow LaTeX printing


def formula2file( formula, file, negate=False ):
    if type(formula)!=str:
        formula=latex(formula)
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


from PIL import Image
from IPython.display import Image, display
from IPython.lib.latextools import latex_to_png
import cv2  
import numpy as np

def whiteimg(yy,xx):
    img = np.zeros([yy,xx,3],dtype=np.uint8)
    img.fill(255) # or img[:] = 255
    return img

def showimg(img):
    ''' 
    input image name  or cv2.img
    outut display img
    '''
    if type(img)==str:
        img=cv2.imread(img)
    img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    display(Image.fromarray(img))
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

#formula1A.png
def completeimg(img):
    if type(img)==str:
        img=cv2.imread(img)
    yy,xx,zz=img.shape 
    if xx<3138:
        imgempty=whiteimg(yy,3137)
        x2=xx
        x1=0
        y2=yy
        y1=0
        imgempty[y1:y2,x1:x2,:]=img
        return imgempty
        
    else:
        return img

def tranimg2file(file,kname):
    img1=appliedbode(img)
    img2=completeimg(img1)
    fimg=cv2.imwrite(kname,img2)
    
    
            
# minifullimg=cv2.imread('formula1.png')
# yy,xx,zz=img.shape
# xx=3091
# yy=187
# yy2=int(yy/8)
# nxx=2*yy2+xx
# nyy=2*yy2+yy
 
# x2=nxx-yy2
# x1=yy2
# y2=yy+yy2
# y1=yy2

# imgempty=whiteimg(nyy,nxx)
# minifullimg=cv2.imread('formula5.png')
# imgempty[y1:y2,x1:x2,:]=minifullimg

# cv2.imwrite('formula5A.png', imgempty) 



  