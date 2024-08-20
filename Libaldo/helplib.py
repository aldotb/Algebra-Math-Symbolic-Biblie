from PIL import Image               # to load images
from IPython.display import display # to display images
import numpy as np
from matplotlib import pyplot as plt


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
    
import os

def get_filesh():
    res = []
    for path in os.listdir('Libaldo/imghelp/'):
        res.append(path)
    mm=''
    for i in res:
        kk=i.replace('.png','(), ')
        mm=mm+kk
    return mm[1::]  