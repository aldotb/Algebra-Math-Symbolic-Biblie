import copy
import itertools as it
import random

from IPython.display import clear_output, Math
from lib_Mathematica import cfrac, sE
from sympy import *
from lib_MyFunctions import *

# from IPython.display import clear_output
# clear_output(wait=True)
z1, z2, z3, z4, z5, z6, z7, z8, z9, z10 = symbols('z1 z2 z3 z4 z5 z6 z7 z8 z9 z10', positive=True)
z11, z12, z13, z14, z15, z16, z0 = symbols('z11 z12 z13 z14 z15 z16 z0', positive=True)
varmain = [z0, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16]

def obj2list(obj):
    if type(obj)==list:
        return obj 
    else:
        return obj()


def smul(*args):
    vec = [str(data) for data in args]
    kres = vec[0]
    for data in vec[1::]:
        kres = kres + '*' + data
    return parse_expr(kres, evaluate=False)


def helpstr2var():
    sE(' example: ')
    sE(" your variables name , sexpr=['A','B','C'] ")
    sE(" run this script '-->  var(','.join(sexpr)) ")
    sE("-------------------------")


def obj2set(*args):
    if len(args) == 1:
        if type(args[0]) == set:
            return args[0]
        if type(args[0]) == list:
            vec = []
            for data in args[0]:
                vec.append(data)
            return set(vec)
        else:
            return {args[0]}
    else:
        vec = []
        for data in args:
            vec.append(data)
        return set(vec)


def traducecarta(carta):
    global spalo
    valor = carta[0]
    palo = carta[1]
    if palo == 'C':
        spalo = u'♥'
    elif palo == 'T':
        spalo = u'♣'
    elif palo == 'D':
        spalo = u'♦'
    else:
        var = spalo == u'♠'
    return valor + spalo


def realc(carta):
    suits = u'♥♠♦♣'
    if not carta[1] in suits:
        return traducecarta(carta)
    else:
        return carta


def getnumber(sexpr):
    slet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    if type(sexpr) == list:
        vec = []
        for ldata in sexpr:
            vec.append(getnumber(ldata))
        return vec

    for data in slet:
        if data in sexpr:
            sexpr = sexpr.replace(data, '')
    return int(sexpr)


def unProb(event, space):
    return 1 - Prob(event, space)


def Prob(event, space):
    favorable = set.intersection  # Outcomes that are in the event and in the sample space
    cases = len
    event = obj2set(event)
    space = obj2set(space)
    numer = cases(favorable(event, space))
    denom = cases(space)
    return simplify(cfrac(numer, denom))


def fus(vec):
    vec2 = ''
    for data in vec:
        vec2 = vec2 + data
    return vec2


def vec_producto(*args):
    '''
    vec_product('AB','CD') = [('C', 'A'), ('C', 'B'), ('D', 'A'), ('D', 'B')]
    vec_product('AB','CD','flat')=['CA', 'CB', 'DA', 'DB']
    vec_product('AB','CD','count')=4
    '''
    global ranks
    ops = ['count', 'flat']
    vecin = []
    for i in args:
        if not i in ops:
            vecin.append(i)
    suits = list(vecin[0])
    p1 = ''
    p2 = ''
    if len(vecin) == 2:
        if type(vecin[1]) == str:
            ranks = list(vecin[1])
        elif type(vecin[1]) == list:
            ranks = vecin[1]
        else:
            p1 = args[1]
            ranks = [str(x) for x in range(p1)]

    if len(vecin) == 3:
        p1 = vecin[1]
        p2 = vecin[2]
        ranks = [str(x) for x in range(p1, p2 + 1)]
    kres = list(it.product(ranks, suits))

    if 'count' in args:
        return len(kres)
    else:
        if 'flat' in args:
            return fuslist(kres)
        return kres


def vec_combina(*args):
    '''
    vec_combina('ABC',2) = [('A', 'B'), ('A', 'C'), ('B', 'C')]
    vec_combina('ABC',2,'repeat')=
            [('A', 'A'), ('A', 'B'), ('A', 'C'), ('B', 'B'), ('B', 'C'), ('C', 'C')]
    vec_combina('ABC',2,'flat')=['AB', 'AC', 'BC']=
            ['AA', 'AB', 'AC', 'BB', 'BC', 'CC']
    vec_combina('ABC',2,'count')=3
    vec_combina('ABC',2,'count','repeat')=6
    '''
    ops = ['count', 'flat', 'repeat']
    suits = list(args[0])
    qq = args[1]
    if 'repeat' in args:
        kres = it.combinations_with_replacement(suits, qq)
    else:
        kres = it.combinations(suits, qq)

    if 'count' in args:
        return len(list(kres))
    else:
        kres = list(kres)
    if 'flat' in args:
        kres = fuslist(kres)
        return kres
    else:
        return kres


def vec_permute(*args):
    '''
    vec_permuta('ABC',2) = 
            [('A', 'B'), ('A', 'C'), ('B', 'A'), ('B', 'C'), ('C', 'A'), ('C', 'B')]         
    vec_combina('ABC',2,'flat')=['AB', 'AC', 'BA', 'BC', 'CA', 'CB']
    vec_combina('ABC',2,'count')=6

    '''

    ops = ['count', 'flat', 'repeat', 'norepeat']
    suits = list(args[0])
    qq = args[1]

    kres = it.permutations(suits, r=qq)
    if 'norepeat' in args:
        kres2 = []
        done = True
        while done:
            item = list(it.islice(kres, 1))
            if item == []:
                done = False
                break
            else:
                item = item[0]
                if not item in kres2:
                    kres2.append(item)

        kres = kres2

    kres = list(kres)

    if 'count' in args:
        return len(list(kres))
    else:
        kres = list(kres)
    if 'flat' in args:
        kres = fuslist(kres)
        return kres
    else:
        return kres


def vec_repart(*args):
    sexpr = args[0]
    qq = 1
    ops = ['list', 'repart', 'flat', 'norepeat', 'repeat']
    for i in args:
        if not i in ops:
            if type(i) == str:
                sexpr = i
            if type(i) == int:
                qq = i
    if qq == 1:
        kres = [x for x in sexpr]
        return kres

    vec1 = sexpr
    vec3 = []
    vec2 = []
    for data in vec1:
        vec2.append(data)
    for L in Range(qq - 1):
        vec3 = []
        for data1 in vec2:
            for data2 in vec1:
                kdata = data1 + data2
                if not 'norepeat' in args:
                    if not kdata in vec3:
                        vec3.append(kdata)
                else:
                    vec3.append(kdata)
        vec2 = vec3
    if 'list' in args:
        vec = []
        for data in vec3:
            vec2 = []
            for data2 in data:
                vec2.append(data2)
            vec.append(vec2)
        vec3 = vec
    return vec3

def crossvec(vec1,vec2):
    kres=[]
    for data1 in vec1:
        for data2 in vec2:
            if type(data1)==list:
                if type(data2)==list:
                    kres.append(data1+data2)
                else: 
                    kres.append(data1+[data2])
            else:
                if type(data2)==list:
                    kres.append([data1]+data2)
                else: 
                    kres.append([data1,data2])
    return kres

def vec_mult(*args,sexpr='',cond=''):
    ops=[]
 
    args2=[]
     
    for data in args:
        if type(data)==str:
            ops.append(data)
        else:
            args2.append(data)
     
    args=args2        
         
    if len(args)==1:
        vec=[args[0],args[0]] 
    elif  len(args) and  type(args[1])==int:
        vec=args[1]*[args[0]]
    else:
        vec=[]
        for data in args:
            vec.append(data)        
    vec1=vec[0]
    for i in range(1,len(vec)):
        vec1=crossvec(vec1,vec[i])
    kres=vec1
    if sexpr!='':
        sexpr2='['+sexpr+' for data in kres]'
        kres=eval(sexpr2)
    if cond!='':
        sexpr2='[data for data in kres if '+ cond+' ]'
        kres=eval(sexpr2)        
    if 'tuple' in ops:
        kres=[tuple(data) for data in kres]
    if 'set' in ops:
        kres= set(kres)
    return kres
    
cevn = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']


def extraecondision(sexpr):
    vecL = []
    veca = []
    p = 0
    for data in sexpr:
        if not data in cevn:
            vecL.append(p)
            veca.append(sexpr[p])
        p += 1
    vecP = []
    vecq = []
    for pos in vecL:
        if len(vecP) == 0:
            vecP.append([0, pos])
        else:
            vecP.append([vecP[-1][1] + 1, pos])
    for data in vecP:
        x1, x2 = data
        vecq.append(int(sexpr[x1:x2]))
    return veca, vecq


def check_cond(sexpr, vecL, vecN):
    done = True
    for datL, datN in zip(vecL, vecN):
        if sexpr.count(datL) != datN:
            return False
    return done


def fuslist(vec):
    vec2 = []
    for data in vec:
        vec2.append(fus(data))
    return vec2


def fcreatelist(vec, scond):
    list1 = vec
    list2 = fuslist(list1)
    vecL, vecN = extraecondision(scond)
    vec3 = []
    for data1, data2 in zip(list2, list1):
        if check_cond(data1, vecL, vecN):
            vec3.append(data2)
    return vec3


def str2num(sexpr):
    valfa = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ><='
    sexpr2 = ''
    for data in sexpr:
        if not data in valfa:
            sexpr2 += data
    if sexpr2 == '':
        return ''
    else:
        return int(sexpr2)


def getalpha(sexpr):
    vnum = '0123456789'
    sexpr2 = ''
    for data in sexpr:
        if not data in vnum:
            sexpr2 += data
    return sexpr2


def xfilter(mvec, sexpr):
    if sexpr == 'odd':
        return (xfilter(mvec, 'str2num(data)%2!=0'))
    elif sexpr == 'even':
        return (xfilter(mvec, 'str2num(data)%2==0'))
    else:

        if type(mvec) == fullevents:
            mvec = mvec.main
        for i in range(len(sexpr)):
            if sexpr[i] != '*':
                vec = []
                for data in mvec:
                    if sexpr[i] == data[i]:
                        vec.append(data)
                mvec = vec
        return mvec


def check_only(mexpr, sexpr):
    cc = 0
    qq = len(sexpr)
    for i in range(len(mexpr)):
        if sexpr == mexpr[i:i + qq]:
            cc += 1
            if cc > 1:
                return False
    if cc == 0:
        return False
    return True


# noinspection PyTypeChecker
def xfilter2(self: object, cond: object) -> object:
    if cond == 'odd':
        return [data for data in self.main if str2num(data) % 2 != 0]
    elif cond == 'even':
        return [data for data in self.main if str2num(data) % 2 == 0]
    else:
        if self.mode == 'short':
            if not '*' in cond:
                newvec = [data for data in self.main if check_only(data, cond)]
                return newvec
            else:
                return xfilter3(self.main, cond)
        else:
            return xfilter4(self.main, cond)


def xfilter3(mvec, cond):
    vecN = '1234567890'
    adone = True
    for data in vecN:
        if data in cond:
            adone = False
            break
    if adone and not '*' in cond:
        vec3 = []
        for data in mvec:
            if cond in data:
                vec3.append(data)
        return vec3

    elif '*' in cond:
        veck = mvec
        for i in range(len(cond)):
            if cond[i] != '*':
                vec0 = [data for data in veck if data[i] == cond[i]]
                veck = vec0

        return veck
    else:
        list1 = mvec
        vecL, vecN = extraecondision(cond)
        vec3 = []
        for data1 in list1:
            if check_cond(data1, vecL, vecN):
                vec3.append(data1)
        return vec3


def xfilter4(mvec2, cond):
    mvec = fuslist(mvec2)
    vec3 = xfilter3(mvec, cond)
    vec4 = [tuple(data) for data in vec3]
    return vec4


ops1 = ('=', '>', '>=', '<', '<=')
ops2 = ('&', '|')
ops3 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
ops4 = '0123456789'


def listdrop(list1, list2):
    vec = []
    for data in list1:
        if not data in list2:
            vec.append(data)

    return vec


#################################


# noinspection PyPropertyDefinition,PyTypeChecker
class fullevents:
    def __init__(self, *args):
        self.cond=''
        args2=[]
        for sdata in args:
            if type(sdata)==str and 'data' in sdata:
                self.cond=sdata
            else:
                args2.append(sdata)
        args=args2        
        ops = ['deck', 'ball', 'coin', 'combine', 'producto', 'permute', 'repeat', 'repart', 'norepeat', 'zero', 'flat',
               'product', 'drop', 'noshow']
        opsb = ['combine', 'producto', 'permute', 'repeat', 'repart', 'product']
        self.main = []
        self.hmain = []
        self.type = ''
        self.cant = []
        self.qty = 0
        self.fam = ''
        self.listfam = []
        self.value = []
        self.flist = []
        self.ownp = []
        self.mode = 'short'
        self.grado = 1
        fdone100 = False
        fdone10 = False
        fdone1000 = False

        vec = []
        ops1 = ('=', '>', '>=', '<', '<=')
        ops2 = ('&', '|')
        ops3 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        ops4 = '0123456789'
        if len(args) == 1 or (len(args) == 2 and 'noshow' in args):
            fam = ''
            vec = args[0]

            for data in vec:
                if str(data) in ops3 and not str(data) in fam:
                    fam = fam + data
            self.fam = fam
            self.main = vec
            self.mode = 'short'
            self.hmain = copy.deepcopy(vec)
        else:

            for data in opsb:

                if data in args:
                    self.mode = 'long'

            if self.mode == 'short':
                self.fam = args[0]
                if type(args[1]) == list:
                    self.cant = args[1]
                else:
                    self.cant = [args[1]] * len(self.fam)

                for i in self.cant:
                    if i > 999:
                        fdone1000 = True
                        self.grado = 4
                    elif i > 99:
                        fdone100 = True
                        self.grado = 3
                    elif i > 9:
                        fdone10 = True
                        self.grado = 2
                    else:
                        pass

                for dataf, dataq in zip(self.fam, self.cant):
                    if 'zero' in args:
                        qn = 0
                    else:
                        qn = 1
                    for i in range(qn, dataq + 1):
                        if i < 10:
                            if fdone1000:
                                self.main.append(dataf + '000' + str(i))
                            elif fdone100:
                                self.main.append(dataf + '00' + str(i))
                            elif fdone10:
                                self.main.append(dataf + '0' + str(i))
                            else:
                                self.main.append(dataf + str(i))
                        elif i < 100:
                            if fdone1000:
                                self.main.append(dataf + '00' + str(i))
                            elif fdone100:
                                self.main.append(dataf + '0' + str(i))
                            else:
                                self.main.append(dataf + str(i))
                        elif i < 1000:
                            if fdone1000:
                                self.main.append(dataf + '0' + str(i))
                            else:
                                self.main.append(dataf + str(i))
                        else:
                            self.main.append(dataf + str(i))
                kres = self.main
                self.hmain = copy.deepcopy(kres)
                if 'flat' in args:
                    self.mode = 'short'
                self.qty = sum(self.cant)
            else:
                if 'combine' in args:
                    fdone = True
                    if 'repeat' in args:
                        if 'flat' in args:
                            self.mode = 'short'
                            vec = vec_combina(args[0], args[1], 'repeat', 'flat')
                        else:
                            vec = vec_combina(args[0], args[1], 'repeat')
                    else:
                        if 'flat' in args:
                            self.mode = 'short'
                            vec = vec_combina(args[0], args[1], 'flat')
                        else:
                            vec = vec_combina(args[0], args[1])
                    self.qty = qcombine(len(args[0]), args[1])
                if 'permute' in args:
                    fdone = True
                    if 'flat' in args:
                        self.mode = 'short'
                        vec = vec_permute(args[0], args[1], 'flat')
                    else:
                        vec = vec_permute(args[0], args[1])
                    self.qty = qpermute(len(args[0]), args[1])
                if 'producto' in args or 'product' in args:
                    fdone = True
                    data1 = ''
                    data2 = ''
                    for data in args:
                        if not data in ops:
                            if data1 == '':
                                data1 = data
                                data2 = data
                            else:
                                data2 = data
                    if type(data2) == int:
                        vec = vecproduct2(data1, data2)
                        self.qty = len(data1) ** data2
                    else:
                        vec = []
                        for v1 in data1:
                            for v2 in data2:
                                vec.append([v1, v2])
                        self.qty = len(data1) * len(data2)
                    if 'flat' in args:
                        self.mode = 'short'
                        try:
                            vec = fuslist(vec)
                        except:
                            pass
                if 'repart' in args:
                    if 'flat' in args:
                        self.mode = 'short'
                    fdone = True
                    vec = vec_repart(*args)
                    self.qty = len(args[0]) ** args[1]
                if 'norepeat' in args:
                    vec2 = []
                    for data in vec:
                        if not data in vec2:
                            vec2.append(data)
                    vec = vec2

                self.main = vec
                self.hmain = vec
                self.fam = args[0]
                
        if self.cond!='':
            L=self.select(self.cond)
            self.main=L 
            self.hmain=L
        
        if not  'noshow' in args:
            display(Math(latex(self.main)))

    def __call__(self, *args):
        slist = self.main
        
        if len(args) == 0:
            return slist
        elif len(args)==1 and type(args[0])==int:
            return slist[args[0]]
        else:
            sexpr = args[0]
            return self.select(*args)


    def s(self):
        display(Math(latex(self.main)))

    @property
    def clone(self):
        return copy.deepcopy(self)

    @property
    def reset(self) -> object:
        self.main = copy.deepcopy(self.hmain)

    @property
    def sizelist(self):
        return len(self.main)

    @property
    def size(self):
        return len(self.main)
    @property    
    def set(self):
        L=self.main 
        if type(L[0])==list:
            L2=[]
            for data in L:
                L2.append(tuple(data))
            return set(L2)
        else:
            return set(L)
                    
    @property
    def sizeitem(self):
        item = self.main[0]
        return len(item)

    def applyfunc(self, func):
        L1 = self.main
        self.flist = [func(data) for data in self.main]

    def dselect(self, sexpr, *args):
        vec = self.select(sexpr, *args)
        kres = listdrop(self.main, vec)
        return kres

    def unionselect(self, *args):
        global kres
        ops = ['norepeat', 'size']
        vec = []
        for data in args:
            if not data in ops:
                vec = vec + self.select(data)
            kres = vec
        if 'norepeat' in args:
            kres = set(vec)
            kres = list(kres)
        if 'size' in args:
            return len(kres)
        return kres

    def pickanddrop(self):
        pick = random.choice(self.main)
        self.drop(pick)
        return pick

    def iselect(self, *args):
        vec = []
        for data in args:
            if type(data) == list:
                vec.append(data)
            else:
                vec.append(self.select(data))

        return get_intersec(*vec)

    def uselect(self, *args):
        vec = []
        for data in args:
            if type(data) == list:
                vec.append(data)
            else:
                vec.append(self.select(data))

        return get_union(*vec)

    def and_select(self, *args):
        vec = []
        for data in args:
            if type(data) == list:
                vec.append(data)
            else:
                vec.append(self.select(data))

        return get_intersec(*vec)

    def or_select(self, *args):
        vec = []
        for data in args:
            if type(data) == list:
                vec.append(data)
            else:
                vec.append(self.select(data))

        return get_union(*vec)

    def add_select(self, *args):
        vec = []
        for data in args:
            if type(data) == list:
                vec = vec + data
            else:
                vec = vec + self.select(data)
        return vec

    def selectcomplement(self, sexpr, *args):
        vec = self.select(sexpr, *args)
        kres = [data for data in self.main if not data in vec]
        return kres

    def secselect(self, *args):
        for cdata in args:
            self.main = self.select(cdata)

        kres = self.main
        self.reset()()()
        return kres
        
    def union(self,obj):
        L2=copy.deepcopy(self.main)
        L=copy.deepcopy(obj2list(obj))
        for sdata in L2:
            if not sdata in L:
                L.append(sdata)
        return L
       
    def intersect(self,obj):
        L3=[]
        L2=copy.deepcopy(self.main)
        L=copy.deepcopy(obj2list(obj))
        for sdata in L2:
            if  sdata in L  and not sdata in L3:
                L3.append(sdata)
        return L3       
       
        
        
    '''
    def select(self, sexpr, *args):
        if type(sexpr) == list:
            return sexpr
        if 'data' in sexpr:
            p1 = "[data for data in self.main if "
            p2 = sexpr
            p3 = "]"
            pp = p1 + p2 + p3
            vecr = eval(pp)
            return vecr
        else:

            ops = ['norepeat', 'size','set']
            ops1 = ('=', '>', '>=', '<', '<=')
            ops2 = ('&', '|')
            ops3 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
            ops4 = '0123456789'

            doneM = True
            for data1 in ops1:
                if data1 in sexpr:
                    doneM = False
            for data1 in ops2:
                if data1 in sexpr:
                    doneM = False

            if not doneM:
                kres = self.filter(sexpr)
            else:
                kres = self.xfilter(sexpr)
            if 'size' in args:
                return len(kres)
            else:
                if 'norepeat'in args or 'set' in args:
                    kres=list(set(kres))
                return (kres)
    '''
    def unselect(self,*args):
        kres=self.select(*args)
        vmain=self.main
        kres=[data for data in vmain if not data in kres]
        return kres
        
     
    def select(self,*args):
        vecs=[]
        vops=[]
        ops=['set','size']
        for data in args:    
            if type(data)==str and not data in ops:
                vops.append(data)
        for datas in vops:
            if vecs==[]:
                vecs=self.select2(datas)
            else:
                vecs=vecs+self.select2(datas)
        kres=vecs
        if 'set' in args:
            kres=list(set(kres))
        if 'size' in args:
            return len(kres)
        return kres    
                 
            
        
    def select2(self,sexpr,*args):
        if type(sexpr) == list:
            return sexpr
        if 'data' in sexpr:
            p1 = "[data for data in self.main if "
            p2 = sexpr
            p3 = "]"
            pp = p1 + p2 + p3
            vecr = eval(pp)
            return vecr
        else:

            ops = ['norepeat', 'size','set']
            ops1 = ('=', '>', '>=', '<', '<=')
            ops2 = ('&', '|')
            ops3 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
            ops4 = '0123456789'

            doneM = True
            for data1 in ops1:
                if data1 in sexpr:
                    doneM = False
            for data1 in ops2:
                if data1 in sexpr:
                    doneM = False

            if not doneM:
                kres = self.filter(sexpr)
            else:
                kres = self.xfilter(sexpr)
 
            return (kres)   
    
    def xfilter(self, sexpr):
        kres = xfilter2(self, sexpr)
        return kres

    def unfilter(self, sexpr):
        kres = xfilter2(self, sexpr)
        kres2 = [data for data in self.main if not data in kres]
        return kres2

    def andfilter(self, sexpr1, sexpr2):
        kres1 = xfilter(self.main, sexpr1)
        kres2 = xfilter(kres1, sexpr2)
        return kres2

    def removelist(self, vec2):
        vec1 = self.main
        for data in vec2:
            vec1.remove(data)
        self.main = vec1

    def drop(self, sexpr):
        vecmain = self.main
        if  sexpr in self.main:
            vecmain.remove(sexpr)
            self.main = vecmain
        elif type(sexpr) == list:
            self.droplist(sexpr)

        else:
            vec = self.select(sexpr)
            self.droplist(vec)
    
    def onedrop(self, sexpr):
        vecmain = self.main
        if  sexpr in self.main:
            vecmain.remove(sexpr)
            self.main = vecmain
        elif type(sexpr) == list:
            item=random.choice(sexpr)
            vecmain.remove(item)
            self.main = vecmain

        else:
            vec = self.select(sexpr)
            item=random.choice(vec)
            vecmain.remove(item)
            self.main = vecmain
            
    def otheranydrop(self, sexpr):
        vecmain = self.main
        if type(sexpr) == str and sexpr in self.main:
            self.otheranydrop([sexpr])
        elif type(sexpr) == list:
            newl = get_complemento(sexpr, vecmain)
            self.anydrop(newl)
        else:
            vec = self.select(sexpr)
            newl = get_complemento(vec, vecmain)
            self.anydrop(newl)

    def anydrop(self, sexpr):
        vecmain = self.main
        if type(sexpr) == str and sexpr in self.main:
            vecmain.remove(sexpr)
            self.main = vecmain
        elif type(sexpr) == list:
            item = random.choice(sexpr)
            self.drop(item)

        else:
            vec = self.select(sexpr)
            item = random.choice(vec)
            self.drop(item)

    def droplist(self, vec):
        vecm = self.main
        kres = listdrop(vecm, vec)
        self.main = kres

    # noinspection PyTypeChecker
    def removexfilter(self, sexpr):

        vec1 = self.main
        vec2 = self.xfilter(sexpr)
        for data in vec2:
            vec1.remove(data)
        self.main = vec1

    def filter(self, sexpr):
        vec = self.main
        fam = self.fam
        ops1 = ('=', '>', '>=', '<', '<=')
        ops2 = ('&', '|')
        ops3 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        ddone = False
        ndone = False
        mexpr = []
        smexpr = []
        cexpr = []
        nexpr = []
        opec = ''
        if '==' in sexpr:
            sres = '[data for data in ' + str(vec) + ' if ' + sexpr + ']'
            return eval(sres)
        elif '&' in sexpr:
            p1 = sexpr.find('&')
            opec = '&'
            mexpr = [sexpr[0:p1], sexpr[p1 + 1:len(sexpr)]]
        elif '|' in sexpr:
            p1 = sexpr.find('|')
            opec = '|'
            mexpr = [sexpr[0:p1], sexpr[p1 + 1:len(sexpr)]]
        else:
            mexpr.append(sexpr)

        for data2 in mexpr:
            data2 = data2.replace(' ', '')
            if data2 in fam:
                smexpr.append("[data for data in " + str(vec) + " if " + "'" + data2 + "'" + ' in data ]')

            else:
                if len(str(str2num(str(data2)))) == len(data2):
                    smexpr.append(str('[ data for data in ' + str(vec) + ' if str2num(str(data))== ' + data2 + ']'))
                else:
                    smexpr.append(str('[ data for data in ' + str(vec) + ' if str2num(data) ' + data2 + ']'))

        if opec != '':
            slist1, slist2 = smexpr
            sset1 = str(set(eval(slist1)))
            sset2 = str(set(eval(slist2)))
            sset3 = sset1 + opec + sset2
            set4 = eval(sset3)
            return list(set4)
        else:
            return eval(smexpr[0])

    def unprob(self, sexpr, *args):
        return simplify(1 - self.prob(sexpr, *args))

    def complement(self, *args):
        vecu = obj2list(args[0])
        vec = [data for data in vecu if not data in self.main]
 
        return vec

    
    def prob(self,sexpr,*args):
     
                
         
        global xlist
        sops=[]
        for data in args:
            if not data in sops:
                sops.append(data)
            
        if type(sexpr) == list:
            xlist = sexpr
            q1 = len(sexpr)
            q2 = (len(self.main))
            kres = cfrac(q1, q2)
        elif sexpr in self.main:
            xlist = sexpr
            xmain = self.main
            q2 = (len(self.main))
            q1 = xmain.count(sexpr)
            kres = cfrac(q1, q2)
        else:
            if type(sexpr) == str and not sexpr in self.main:
                try:
                    xlist = self.filter(sexpr)
                except:
                    xlist = self.select(sexpr)

            q2 = len(self.main)
            q1 = len(xlist)
            kres = cfrac(q1, q2)

            self.ownp.append(kres)
        if 'drop' in args:
            if type(sexpr) == str and sexpr in self.main:
                self.main.remove(sexpr)
            else:
                item = random.choice(xlist)
                self.main.remove(item)
        if 'reset' in args:
            vecr = self.hmain
            self.main = vecr
        if 'invert' in args:
            kres = 1 / kres
        if 'float' in args:
            kres = float(kres)

        return kres
    '''
    def prob(self, sexpr, *args):
        global xlist
        if type(sexpr) == list:
            xlist = sexpr
            q1 = len(sexpr)
            q2 = (len(self.main))
            kres = cfrac(q1, q2)
        elif sexpr in self.main:
            xlist = sexpr
            xmain = self.main
            q2 = (len(self.main))
            q1 = xmain.count(sexpr)
            kres = cfrac(q1, q2)
        else:
            if type(sexpr) == str and not sexpr in self.main:
                try:
                    xlist = self.filter(sexpr)
                except:
                    xlist = self.select(sexpr)

            q2 = len(self.main)
            q1 = len(xlist)
            kres = cfrac(q1, q2)

            self.ownp.append(kres)
        if 'drop' in args:
            if type(sexpr) == str and sexpr in self.main:
                self.main.remove(sexpr)
            else:
                item = random.choice(xlist)
                self.main.remove(item)
        if 'reset' in args:
            vecr = self.hmain
            self.main = vecr
        if 'invert' in args:
            kres = 1 / kres
        if 'float' in args:
            kres = float(kres)

        return kres
    '''
    def secprob(self, *args,loop=1):
        PP = 1
        vecd = []
        vop = []
        opp = ['drop', 'reset', 'float']
        varg1 = []
        varg2 = []
        for data in args:
            if data in opp:
                vop.append(data)
            else:
                vecd.append(data)    
        for L in range(loop):
            for data in vecd:
                if 'drop' in vop:
                    PP = PP * self.prob(data, 'drop')
                else:
                    PP = PP * self.prob(data)
        if 'float' in vop:
            PP = float(PP)

        if 'reset' in vop:
            self.reset 

        return PP

    def mul(self, sexpr):
        vec = []
        for data2 in sexpr:
            for data1 in self.main:
                vec.append(data1 + data2)
        self.main = vec

    def Add(self, Obj):
        self.fam = get_newfam(self.fam, Obj)
        if type(Obj) == 'list':
            self.main = self.main + Obj
        else:
            self.main.append(Obj)

    def qpermute(self, mm=''):
        vec = [self.qty]
        if mm == '':
            vec.append(mm)
        return qpermute(*vec)

    def qcombine(self, mm):
        vec = [self.qty, mm]
        return qcombine(*vec)


def removelist(vec1, vec2):
    if type(vec1) == crealist:
        vec1 = vec1.main
    for data in vec2:
        vec1.remove(data)
    return vec1


class PBayes:
    def __init__(self, vec1, vec2):
        self.S = vec1
        self.V = vec2
        self.S1 = []
        self.V1 = []

    def complemento(self, obj):
        k = []
        for data in obj:
            k.append(1 - data)
        return k

    def __call__(self, *args):
        if len(args) == 1 and len(args[0]) == 1 and args[0] in self.S:
            sval = args[0]
            return get_corres(sval, self.S, self.V)
        if len(args) == 1 and len(args[0]) == 1 and not args[0] in self.S:
            data1 = self.V
            data2 = eval('self.' + args[0])
            kres = 0
            for d1, d2 in zip(data1, data2):
                kres = kres + d1 * d2
            return kres

        if len(args) == 2 and args[1] == "'" and args[0] in self.S:
            sval = args[0]
            kres = get_corres(sval, self.S, self.V)
            return 1 - kres

        if len(args) == 2 and args[1] == "'" and args[0] not in self.S:
            kres = 0
            for data1, data2 in zip(self.V, eval('self.' + args[0])):
                kres = kres + data1 * data2
            return simplify(1 - kres)

        if len(args) == 1 and len(args[0]) == 3:
            sexpr = args[0]

            if sexpr[1] == '/':
                if sexpr[2] in self.S and sexpr[0] not in self.S:
                    for S1, V1, V2 in zip(self.S, self.V, eval('self.' + sexpr[0])):
                        if S1 == sexpr[2]:
                            return V2
                else:
                    p1 = self(sexpr[0])
                    p2 = self(sexpr[2] + sexpr[1] + sexpr[0])
                    p3 = self(sexpr[2])
                    return simplify(p1 * p2 / p3)


def respeccomp(vec):
    return [1 - data for data in vec]


def get_corres(val, vec1, vec2):
    for data1, data2 in zip(vec1, vec2):
        if data1 == val:
            return data2


def vecstr2sexpr(vec, ops='+'):
    sexpr = ''
    for fdata in vec:
        sexpr = sexpr + ops + fdata
    sexpr = sexpr[1::]
    return sexpr


def sexpr2expr(sexpr):
    return parse_expr(sexpr)


def vecstr2mathexpr(vec, ops='+'):
    sexpr = vecstr2sexpr(vec, ops=ops)
    return sexpr2expr(sexpr)


def Combine(*args):
    p1 = args[0]
    p2 = 0
    if len(args) == 2:
        p2 = args[1]
    return factorial(p1) / (factorial(p2) * factorial(p1 - p2))


def rCombine(*args):
    p1 = args[0]
    p2 = 0
    if len(args) == 2:
        p2 = args[1]
    return factorial(p1 - p2 - 1) / factorial(p2)


def rPermute(*args):
    p1 = args[0]
    kres = 1
    for i in range(1, len(args)):
        kres = kres * factorial(args[i])
    return factorial(p1) / kres


class Mzets:
    def __init__(self, *args):
        self.fam = args[0]
        self.name = args[1]
        self.qty = len(self.fam)

        self.varmain = varmain
        vecs = []
        for i in range(self.qty, 0, -1):
            L0 = fullevents(self.fam, i, 'combine', 'flat')
            vecs = vecs + L0.main
        vecvar = self.varmain[0:len(vecs)]
        self.U = self.varmain[0:len(vecs) + 1]
        VV = []
        for sdata in self.fam:
            vec2 = []
            for v1, v2 in zip(vecs, vecvar):
                if sdata in v1:
                    vec2.append(v2)
            VV.append(vec2)
        self.vlist = VV

    def genersets(self):
        vec = []
        for data, name in zip(self.vlist, self.fam):
            vec.append(Pzet(name, data, self.U))
        vec.append(set(self.U))
        return vec

x0=symbols('x0')
x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15=symbols('x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15',positive=True)
vecxp=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15]
def ssetsets(*args):
    sdata=args[0]
    qq=len(args[0])
    svec=vec_combina(sdata,1,'flat')
    kres=[[data] for data in vecxp[0:qq] ]
    cc=qq
    if qq>1:
        for i in range(2,qq+1):
            vec2=vec_combina(sdata,i,'flat')
            for data in vec2:
                for L in data:
                    pos=sdata.index(L)
                    valor=kres[pos]
                    valor=valor+[vecxp[cc]]
                    kres[pos]=valor
                cc+=1
     
    ures=[x0]+vecxp[0:cc]
    return kres,ures
def createsets(fam):
    vecf,uset=ssetsets(fam)
    qq=len(fam)
    kres=[]
    for i in range(qq):
        sname=fam[i]
        nset=Pzet(sname,vecf[i],uset)
        kres.append(nset)
    Uset= Pzet('U',uset,uset)
    kres.append(Uset)
    return kres    
class Pzet:

    def __init__(self, name, lvar, uvar,x0=x0,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,x7=x7,x8=x8,x9=x9,x10=x10):
        self.name = name
        self.lvar = lvar
        self.uvar = uvar
        self.vvar = copy.deepcopy(lvar)
        self.svar = [str(data) for data in uvar]
        self.wvar = copy.deepcopy(uvar)
        self.prob=0
        self.x0=x0
        self.x1=x1
        self.x2=x2
        self.x3=x3
        self.x4=x4
        self.x5=x5
        self.x6=x6
        self.x7=x7
        self.x8=x8
        self.x9=x9
        self.x10=x10
        
    def update_val(self,vecvar,vecval):
        for sdata,vdata in zip(vecvar,vecval):
            if sdata in self.lvar:
                kpos=self.lvar.index(sdata)
                self.vvar[kpos]=vdata
            if sdata in self.uvar:
                kpos=self.uvar.index(sdata)
                self.wvar[kpos]=vdata 
                
    def set(self,*args,**kwargs):
        if len(args)==2:
            vecvar,vecval=args[0],args[1]
            self.update_val(vecvar,vecval)
        
        if len(kwargs)!=0:
            for key, value in kwargs.items():
                if key in self.lvar:
                    kpos=self.lvar.index(key)               
                    self.vvar[kpos]=value
                if key in self.uvar:
                    kpos=self.uvar.index(key)
                    self.wvar[kpos]=value
       
            
        
    def pvalue(self):
        return sum(self.vvar)
        
    def setpro(self,valor):
        self.prob=valor
    def __call__(self, *args):
        if len(args) == 0:
            return self.lvar

    def __repr__(self):
        return str(self.lvar)

    def _latex(self, obj):
        return latex(self.lvar)

    def __str__(self):
        return self.__repr__()
    
    def Prob(self,sexpr):
        if type(sexpr)==list:
            kres=sum(sexpr)
        else:
            kres= sexpr
        return self.anyval(kres)    
        
    def valuevector(self,vec):
        for sdata,vdata in zip(self.uvar,self.wvar):
            if sdata in vec:
                kpos=self.uvar.index(sdata)
                vec[kpos]=vdata
        return vec
    
    def anyval(self,sexpr):
        for sdata,vdata in zip(self.uvar,self.wvar):
            sexpr=sexpr.subs(sdata,vdata)
        return sexpr    
    def comp(self,*args):
        kres= get_complemento(self.lvar,self.uvar)
        if 'prob' in args:
            return self.anyval(sum(kres))
        else:
            return kres

    def union(self,*args):
        obj=args[0]
        if type(obj)==Pzet:
            obj=obj.lvar
        kres= get_union(self.lvar,obj)
        if 'prob' in args:
            kres= self.anyval(sum(kres))
        if 'float' in args:
            kres=float(kres)
        return kres
    
    def cunion(self,*args):
        obj=args[0]
        if type(obj)==Pzet:
            obj=obj.lvar
        kres= get_union(self.C,obj)
        if 'prob' in args:
            kres= self.anyval(sum(kres))
        if 'float' in args:
            kres=float(kres)
        return kres

    def substrac(self,*args):
        obj=args[0]
        if type(obj)==Pzet:
            obj=obj.lvar
        kres= get_difer(self.lvar,obj)
        if 'prob' in args:
            kres=  self.anyval(sum(kres))
        if 'float' in args:
            kres=float(kres)
        return kres

    def csubstrac(self,*args):
        obj=args[0]
        if type(obj)==Pzet:
            obj=obj.lvar
        kres= get_difer(self.C,obj)
        if 'prob' in args:
            kres=   self.anyval(sum(kres))
        if 'float' in args:
            kres=float(kres)
        return kres
    
    def subsubs(self,expr):
        for data1,data2 in zip(self.uvar,self.vvar):
                expr=expr.subs(data1,data2)
        return expr
        
    def intersec(self,*args):
        obj=args[0]
        if type(obj)==Pzet:
            obj=obj.lvar
        kres= get_intersec(self.lvar,obj)
        if 'prob' in args:
            kres=   self.anyval(sum(kres))
        if 'float' in args:
            kres=float(kres)
        return kres
    def cintersec(self,*args):
        obj=args[0]
        if type(obj)==Pzet:
            obj=obj.lvar
        kres= get_intersec(self.C,obj)
        if 'prob' in args:
            kres=    self.anyval(sum(kres))
        if 'float' in args:
            kres=float(kres)
        return kres    
        
    @property
    def C(self):
        return get_complemento(self.lvar,self.uvar)
    
    def talque(self,*args):
        obj=args[0]
        if type(obj)==Pzet:
            obj=obj.lvar
        Pb=cfrac(sum(obj),sum(self.uvar))
        AB=get_intersec(self.lvar,obj)
        Pab=cfrac(sum(AB),sum(self.uvar))
        return cfrac(Pab,Pb)
        
        
            
 
    @property
    def U(self):
        return self.uvar
        
    @property
    def size(self):
        return len(self.lvar)
    @property    
    def u(self):
        return sum(self.uvar)
    @property
    def P(self):
        return sum(self.lvar)
    @property    
    def p(self ):
        kres= self.anyval(sum(self.lvar))
        return kres
        
    @property   
    def c(self):
        vecc=[data for data in self.uvar if not data in self.lvar]
        return self.anyval(sum(vecc))    

def obj2set(obj):
    if type(obj) == Pzet:
        return set(obj.lvar)
    elif type(obj) == list:
        return set(obj)
    else:
        return obj


def getdatap(*args):
    nn = ''
    mm = ''
    ops = []
    for data in args:
        if type(data) == str:
            ops.append(data)
        else:
            if nn == '':
                nn = data
            else:
                mm = data
    return nn, mm, ops


def qpermute(*args):
    nn, mm, ops = getdatap(*args)
    kres = factorial(nn)
    if ops == '':
        if mm == '':
            return kres
        else:
            kres2 = factorial(nn - mm)
            kres = kres / kres2
    else:
        if type(mm) == int:
            kres = nn ** mm
        else:
            kres2 = 1
            for data in mm:
                kres2 *= factorial(data)
            kres = kres / kres2
    if 'inv' in args:
        kres = cfrac(1, kres)
    if 'float' in args:
        kres = float(kres)
    return kres


def qcombine(*args):
    mm, nn, ops = getdatap(*args)
    if not 'repeat' in ops:
        kres = cfrac(factorial(mm), factorial(nn) * factorial(mm - nn))
    else:
        kres = cfrac(factorial(mm + nn - 1), factorial(nn) * factorial(mm - 1))
    if 'inv' in args:
        kres = cfrac(1, kres)
    if 'float' in ops:
        kres = float(kres)
    return kres


def vecproduct2(L, qq):
    L1 = [[data] for data in L]
    for i in range(qq - 1):
        vec = []
        for data in L1:
            for item in L:
                vec.append(data + [item])
        L1 = vec
    return L1


U = symbols('U')


def get_universe(*args):
    vec = set(args[0])
    for data in args[1::]:
        if type(data) == list:
            vec = vec | set(data)
        else:
            vec = list(vec)
            vec.append(data)
            vec = set(vec)
    return list(vec)
    

def get_difer(*args):
    vec1 = listlist(args[0])
    vec2 = listlist(args[1])
    vec=[]
    for data in vec1:
        if not data in vec2:
            vec.append(data)
    return vec    

def listlist(conjunto):
    kres=list(conjunto)
    if type(kres[0])==tuple or type(kres[0])==set:
        kres=[list(data) for data in kres]
    return kres    
    
def get_complemento(*args):
    vec1 = listlist(args[0])
    vec2 = listlist(args[1])
    kres = [data for data in vec2 if not data in vec1]
    return kres


def get_union(*args):
    vec1 = obj2list(args[0])
    for vec2 in args[1::]:
        for data in obj2list(vec2):
            if not data in vec1:
                vec1.append(data)
    return vec1 


def get_intersec(*args):
    vec1 = obj2list(args[0])
    for vec2 in args[1::]:
        vec3=[]
        for data in vec1:
            if data in vec2 and not data in vec3:
                vec3.append(data)
        vec1= vec3    
    return  vec1 


def cprob(CC, u=U):
    return cfrac(len(CC), len(u))


def iprob(*args, u=U):
    vec = set(args[0])
    for data in args[1::]:
        vec = vec & set(data)
    vec = list(vec)
    return cprob(vec, u=u)


def uprob(*args, u=U):
    vec = set(args[0])
    for data in args[1::]:
        vec = vec | set(data)
    vec = list(vec)
    return cprob(vec, u=u)


def dueprob(CC1, CC2, u=U):
    p1 = cprob(CC1, u=u)
    p2 = iprob(CC1, CC2, u=u)
    return p1 - p2


def sumdigit(sexpr):
    kres = 0
    for data in sexpr:
        kres += int(data)
    return kres


def alphanum2alpha(sexpr):
    vec = '0123456789'
    for data in vec:
        sexpr = sexpr.replace(data, '')
    return sexpr


def getfam(sexpr):
    if alphanum2alpha(sexpr) == '':
        return sexpr
    else:
        return alphanum2alpha(sexpr)


def get_newfam(fam, obj):
    if type(obj) == list:
        for data in obj:
            fam = get_newfam(fam, data)
        return fam
    else:
        nfam = getfam(obj)
        if not nfam in fam:
            return fam + nfam
        else:
            return fam


def get_sargs(*args):
    vec = []
    vecL = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for data in args:
        if type(data) == str:
            if not data[0] in vecL:
                vec.append(data)
    return vec


def recive_data(*args):
    for data in args:
        if type(data) != str:
            return True
    return False


def subscontruye(ss, qq, tt, vecp):
    newv = []
    for data in vecp:
        for i in range(0, min([tt, qq]) + 1):
            nsa = i * ss
            wsa = nsa + data
            if len(wsa) <= tt:
                newv.append(wsa)
            else:
                if not data in newv:
                    newv.append(data)
    return newv


class playiter:
    def __init__(self, *args):

        self.thererepit = False
        ops = get_sargs(*args)
        self.fam = args[0]
        self.qty = len(self.fam) * [1]
        for i in args:
            if type(i) == list:
                self.qty = i
        self.qtotal = sum(self.qty)
        if len(self.fam) < self.qtotal:
            self.thererepit = True

    def npossibles(self):
        T = self.qtotal
        T2 = 1
        for data in self.qty:
            T2 = T2 * factorial(data)
        kres = cfrac(factorial(T), T2)
        if 'float' in self.ops:
            kres = float(kres)
        return kres

    def npermute(self, *args):
        if not recive_data(*args):
            if not self.thererepit:
                kres = qpermute(self.qtotal)
        return True

    def gpermute(self, *args):
        if not recive_data(*args):
            if not self.thererepit:
                qpermute(self.qtotal)
        return True

    def gconvine(self, qt, *args):
        ['']

        fam = self.fam
        qvec = self.qty
        vecp = ['']

        for ss, qq in zip(fam, qvec):
            vecp = subscontruye(ss, qq, qt, vecp)
        kres4 = [data for data in vecp if len(data) == qt]
        if 'qty' in args:
            return len(kres4)
        return kres4
        
        
def real_subs(QQ,**kwargs):
    ''' 
    QQ= symbols function
    ** kwargs c7=6,g4=z..etc..
    RETURN real substitucion when variable have underscore name like 'c_7' 'g_4'    
    '''
    
    vvar=list(symbolslist(QQ))
    mvar=[]
    for i in vvar:
        nname=str(i)
        nname=nname.replace('_','')
        nvar=symbols(nname)
        mvar.append(nvar)
         
    kres=QQ 
    sres=str(kres)
    sres=sres.replace('_','')
    kres=parse_expr(sres)
    mkey=[]
    vvalue=[]
    for key, value in kwargs.items():
        mkey.append(key)
        if type(value)==MyEq:
            vvalue.append(value.ksym)
        else:
            vvalue.append(value)
        
        kres=kres.subs(parse_expr(key),value)
    for i,j in zip(mvar,vvar):
        kres=kres.subs(i,j)
    return (kres) 
    
def symbolslist(vec,var=''):
    '''
        symbolslist(x+a+b+4)=[x,a,b]
        symbolslist(x+a+b+4,x)=[a,b]
        symbolslist([x+a+b+4,a,z])=[a,b,x,z]
        symbolslist([x+a+b+4,a,z],x)=[a,b,xz]

    '''
    if type(vec)==list:
        qq=len(vec)
        kres=0
        cc=10
        for i in vec:
            kres=expand(kres+i/cc)
            cc=cc*10
    else:
        kres=vec
        
    slist=kres.free_symbols
    vecs=list(slist)
    if var!='':
        kres=[]
        for i in vecs:
            if i!=var:
                kres.append(i)
    else:
        kres=slist
    return kres 

def solvesets(*args):
    '''
    solsets=solvesets(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8)
    [vecnames,vecvalues]=solsets
    vecname=x0,x1,x2,x3....
    vecvalues=0.1,0.2,0.2......
    '''
    qs=[*args]+['pack']
    return simplesolve(*qs)

def setsets(*args):
    '''
    solsets solve real values for each elements in Pzets class.
    solsets=solvesets(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8)
    [vecnames,vecvalues]=solsets
    vecname=x0,x1,x2,x3....
    vecvalues=0.1,0.2,0.2......
    setsets(A,B,C,solsets)
    
    '''

    vecs=[]
    for data in args:
        if type(data)==Pzet:
            vecs.append(data)
        else:
            vecvar,vecval=data[0],data[1]
    for ssets in vecs:
        ssets.set(vecvar,vecval)    
        
       

def finfposp(sval,svec,vecv):
    kres=svec.index(sval)
    return vecv[kres]

Pa,Pb=symbols('P(A) P(B)')  
Ua=symbols(r'P(\bar{A})')
Ub=symbols(r'P(\bar{B})') 
Pab=symbols('P({A}\mid{B})')
Pba=symbols('P({B}\mid{A})')
Uab=symbols(r'P(\bar{A}\mid{B})')
Uba=symbols(r'P(\bar{B}\mid{A})')
abU=symbols(r'P({A}\mid\bar{B})')
baU=symbols(r'P({B}\mid\bar{A})')
UabU=symbols(r'P(\bar{A}\mid\bar{B})')
UbaU=symbols(r'P(\bar{B}\mid\bar{A})')

class symbolBayes:
    def __init__(self, fam):
        self.F = [data for data in fam]
        self.q=len(fam)
        self.A=Pa
        self.B=Pb
        self.cA = Ua
        self.cB = Ub
        self.AB= Pab
        self.BA= Pba
        self.cAB= Uab
        self.cBA= Uba
        self.ABc= abU
        self.BAc= baU
        self.cABc=UabU
        self.cBAc=UbaU

         
        V1=['BA','BAc','cBA','cBAc']
        V2=['B','B','cB','cB']
        V3=['cB','cB','B','B']
        V4=['BAc','Ba','cBAc','cBA']
        
 
        
        svec=['A','B','cA','cB','AB','BA','cAB','cBA','ABc','BAc','cABc','cBAc']
        vecv=[self.A,self.B, self.cA ,self.cB,self.AB,self.BA,self.cAB,self.cBA,self.ABc,self.BAc,self.cABc,self.cBAc]

    def expand_bayes(self,obj):
        svec=['A','B','cA','cB','AB','BA','cAB','cBA','ABc','BAc','cABc','cBAc']
        vecv=[self.A,self.B, self.cA ,self.cB,self.AB,self.BA,self.cAB,self.cBA,self.ABc,self.BAc,self.cABc,self.cBAc]
        V1=['A','BA','BAc','cA']
        V2=['B','B','cB','cB']
        V3=['cB','cB','B','B']
        V4=['BAc','Ba','cBAc','cBA']
        if obj=='AB':
            sdata=V1
        elif obj=='cAB':
            sdata==V2
        elif obj=='ABc':
            sdata=V3
        elif obj=='cABc':
            sdata=V4
        else:
            pass
        partes=[finfposp(data,svec,vecv) for data in sdata]
        x1,x2,x3,x4=partes

        return cfrac((x1*x2),(x1*x2+x3*x4))

    def intersec(self,obj):
        svec=['A','B','cA','cB','AB','BA','cAB','cBA','ABc','BAc','cABc','cBAc']
        vecv=[self.A,self.B, self.cA ,self.cB,self.AB,self.BA,self.cAB,self.cBA,self.ABc,self.BAc,self.cABc,self.cBAc]
        V1=['AB','BA','cAB','cBA','ABc','BAc','cABc','cBAc']
        V2=['BA','AB','BAc','ABc','cBA','cAB','cBAc','cABc']
        V3=['A','B','cA','cB','A','B','cA','cB']
         
        x1=finfposp(obj,V1,V2)
        X1=finfposp(x1,svec,vecv)
        x2=finfposp(obj,V1,V3)
        X2=finfposp(x2,svec,vecv)

        return X1*X2 

     
def ps_overline(sexpr):
    return '\overline{'+sexpr+'}'
def ps_union(str1,str2,pair=True):
    if pair:
        return '({'+str1+'}\cup{'+str2+'})'
    else:
        return '{'+str1+'}\cup{'+str2+'}'
def ps_inter(str1,str2,pair=True):
    if pair:
        return '({'+str1+'}\cap{'+str2+'})'
    else:
        return '{'+str1+'}\cap{'+str2+'}'
def ps_or(str1,str2):
    return '({'+str1+'}\mid{'+str2+'})'

def ptraduce(sexpr):
    sexpr=str(sexpr)
    if '\\cap' in sexpr:
        type='I'
        sexpr=sexpr.replace('cap','')
    elif '\\cup' in sexpr:
        type='U'
        sexpr=sexpr.replace('cup','')
    else:
        type='O'
    fam='ABC'
    vecv=[]
    for  data in sexpr:
        if data in fam:
            vecv.append(data)
            
    if sexpr.count('overline')==2:
        vecv2=[]
        for data in vecv:
            vecv2.append(data,lower())
    elif sexpr.count('overline')==1:
        if sexpr[0:len('({\\overline')]=='({\\overline':
            vecv[0]=vecv[0].lower()
        else:
            vecv[1]=vecv[1].lower()
    if len(vecv)==1:
        return 'P',vecv[0]
    else:    
        return type,vecv[0]+vecv[1]

U3name=[]
U3value=[]
for data1 in ['A','a']:
    for data2 in ['B','b']: 
        for data3 in ['C','c']:
            s1=data1
            if data1.islower():
                s1=ps_overline(data1.upper())
            s2=data2
            if data2.islower():
                s2=ps_overline(data2.upper())
            s3=data3
            if data3.islower():
                s3=ps_overline(data3.upper())
            U3value.append(symbols(ps_union(ps_union(s1,s2,pair=False),s3)))
            U3name.append('U'+data1+data2+data3)

I3name=[]
I3value=[]
for data1 in ['A','a']:
    for data2 in ['B','b']: 
        for data3 in ['C','c']:
            s1=data1
            if data1.islower():
                s1=ps_overline(data1.upper())
            s2=data2
            if data2.islower():
                s2=ps_overline(data2.upper())
            s3=data3
            if data3.islower():
                s3=ps_overline(data3.upper())
            I3value.append(symbols(ps_inter(ps_inter(s1,s2,pair=False),s3)))
            I3name.append('I'+data1+data2+data3)

UI5name=[]
UI5value=[]
for data1 in ['A','a']:
    for data2 in ['B','b']: 
        for data3 in ['C','c']:
            s1=data1
            if data1.islower():
                s1=ps_overline(data1.upper())
            s2=data2
            if data2.islower():
                s2=ps_overline(data2.upper())
            s3=data3
            if data3.islower():
                s3=ps_overline(data3.upper())
                 
            UI5value.append(symbols(ps_union(ps_inter(s1,s2,pair=False),s3)))
            UI5name.append(data1+'U'+data2+'I'+data3)
            UI5value.append(symbols(ps_inter(ps_union(s1,s2,pair=False),s3)))
            UI5name.append(data1+'I'+data2+'U'+data3)
vec_U=[]
vec_O=[]
vec_I=[]
vec_P=[]
svec_U=[]
svec_O=[]
svec_I=[]
svec_P=[]

svec1=['A','B','C','a','b','c']
svec2=[]
for ss in svec1:
    SS=copy.deepcopy(ss)
    if ss.islower():
        ss=ps_overline(ss.upper())
        
    vec_P.append(symbols("("+ss+")",commutative=False))
    svec_P.append('self.'+SS)

for data1 in svec1:
    for data2 in svec1:
        if data1.upper()!=data2.upper():
            svec2.append(data1+data2)
for data in svec2:
    s1=data[0]
    S1=data[0]
    if s1.islower():
        s1=ps_overline(s1.upper())
    s2=data[1]
    S2=data[1]
    if s2.islower():
        s2=ps_overline(s2.upper())
    su=ps_union(s1,s2)
    vec_U.append(symbols(su,commutative=False))
    svec_U.append('self.U'+S1+S2)
    si=ps_inter(s1,s2)
    vec_I.append(symbols(si,commutative=False))
    svec_I.append('self.I'+S1+S2)
    so=ps_or(s1,s2)
    vec_O.append(symbols(so,commutative=False))
    svec_O.append('self.'+S1+S2)
def get_complementname(obj):
    if obj.islower():
        return obj.upper()
    else:
        return obj.lower()
class probabilityEq:
    def __init__(self,sname='ABC',exclude=True):
        self.fam = sname
        self.exclude=exclude
        self.q = len(sname)
        self.A, self.B, self.C, self.a, self.b, self.c = vec_P
        self.AB, self.AC, self.Ab, self.Ac, self.BA, self.BC, self.Ba, self.Bc, self.CA, self.CB, self.Ca, self.Cb, self.aB, self.aC, self.ab, self.ac, self.bA, self.bC, self.ba, self.bc, self.cA, self.cB, self.ca, self.cb = vec_O
        self.UAB, self.UAC, self.UAb, self.UAc, self.UBA, self.UBC, self.UBa, self.UBc, self.UCA, self.UCB, self.UCa, self.UCb, self.UaB, self.UaC, self.Uab, self.Uac, self.UbA, self.UbC, self.Uba, self.Ubc, self.UcA, self.UcB, self.Uca, self.Ucb = vec_U
        self.IAB, self.IAC, self.IAb, self.IAc, self.IBA, self.IBC, self.IBa, self.IBc, self.ICA, self.ICB, self.ICa, self.ICb, self.IaB, self.IaC, self.Iab, self.Iac, self.IbA, self.IbC, self.Iba, self.Ibc, self.IcA, self.IcB, self.Ica, self.Icb = vec_I
        self.pname = ['A', 'B', 'C', 'a', 'b', 'c']
        self.oprob = [self.A, self.B, self.C, self.a, self.b, self.c]
        self.oname = ['AB', 'AC', 'Ab', 'Ac', 'BA', 'BC', 'Ba', 'Bc', 'CA', 'CB', 'Ca', 'Cb', 'aB', 'aC', 'ab', 'ac',
                      'bA', 'bC', 'ba', 'bc', 'cA', 'cB', 'ca', 'cb']
        self.oif = [self.AB, self.AC, self.Ab, self.Ac, self.BA, self.BC, self.Ba, self.Bc, self.CA, self.CB, self.Ca,
                    self.Cb, self.aB, self.aC, self.ab, self.ac, self.bA, self.bC, self.ba, self.bc, self.cA, self.cB,
                    self.ca, self.cb]
        self.ounion = [self.UAB, self.UAC, self.UAb, self.UAc, self.UBA, self.UBC, self.UBa, self.UBc, self.UCA,
                       self.UCB, self.UCa, self.UCb, self.UaB, self.UaC, self.Uab, self.Uac, self.UbA, self.UbC,
                       self.Uba, self.Ubc, self.UcA, self.UcB, self.Uca, self.Ucb]
        self.ointer = [self.IAB, self.IAC, self.IAb, self.IAc, self.IBA, self.IBC, self.IBa, self.IBc, self.ICA,
                       self.ICB, self.ICa, self.ICb, self.IaB, self.IaC, self.Iab, self.Iac, self.IbA, self.IbC,
                       self.Iba, self.Ibc, self.IcA, self.IcB, self.Ica, self.Icb]
        
        self.UABC,  self.UABc,  self.UAbC,  self.UAbc,  self.UaBC,  self.UaBc,  self.UabC,  self.Uabc=U3value
        self.IABC,  self.IABc,  self.IAbC,  self.IAbc,  self.IaBC,  self.IaBc,  self.IabC,  self.Iabc=I3value

        self.U3=[self.UABC,  self.UABc,  self.UAbC,  self.UAbc,  self.UaBC,  self.UaBc,  self.UabC,  self.Uabc]
        self.I3=[self.IABC,  self.IABc,  self.IAbC,  self.IAbc,  self.IaBC,  self.IaBc,  self.IabC,  self.Iabc]

        self.AUBIC, self.AIBUC, self.AUBIc, self.AIBUc, self.AUbIC, self.AIbUC, self.AUbIc, self.AIbUc, self.aUBIC, self.aIBUC, self.aUBIc, self.aIBUc, self.aUbIC, self.aIbUC, self.aUbIc, self.aIbUc= UI5value
        self.UI5=[self.AUBIC, self.AIBUC, self.AUBIc, self.AIBUc, self.AUbIC, self.AIbUC, self.AUbIc, self.AIbUc, self.aUBIC, self.aIBUC, self.aUBIc, self.aIBUc, self.aUbIC, self.aIbUC, self.aUbIc, self.aIbUc]
        
        self.mainname = ['A', 'B', 'C', 'a', 'b', 'c'] + self.oname + ['U' + data for data in self.oname] + ['I' + data for data in self.oname]+U3name+I3name+UI5name
        self.mainvalue = [self.A, self.B, self.C, self.a, self.b, self.c] + self.oif + self.ounion + self.ointer +self.U3+self.I3+UI5value
        self.objEQ=''
    
    def update(self):
        self.A, self.B, self.C, self.a, self.b, self.c = self.mainvalue[0:6]
        self.AB, self.AC, self.Ab, self.Ac, self.BA, self.BC, self.Ba, self.Bc, self.CA, self.CB, self.Ca, self.Cb, self.aB, self.aC, self.ab, self.ac, self.bA, self.bC, self.ba, self.bc, self.cA, self.cB, self.ca, self.cb = self.mainvalue[6:30]
        self.UAB, self.UAC, self.UAb, self.UAc, self.UBA, self.UBC, self.UBa, self.UBc, self.UCA, self.UCB, self.UCa, self.UCb, self.UaB, self.UaC, self.Uab, self.Uac, self.UbA, self.UbC, self.Uba, self.Ubc, self.UcA, self.UcB, self.Uca, self.Ucb = self.mainvalue[30:54]
        self.IAB, self.IAC, self.IAb, self.IAc, self.IBA, self.IBC, self.IBa, self.IBc, self.ICA, self.ICB, self.ICa, self.ICb, self.IaB, self.IaC, self.Iab, self.Iac, self.IbA, self.IbC, self.Iba, self.Ibc, self.IcA, self.IcB, self.Ica, self.Icb = self.mainvalue[54:78]
        self.UABC,  self.UABc,  self.UAbC,  self.UAbc,  self.UaBC,  self.UaBc,  self.UabC,  self.Uabc=self.mainvalue[78:86]
        self.IABC,  self.IABc,  self.IAbC,  self.IAbc,  self.IaBC,  self.IaBc,  self.IabC,  self.Iabc=self.mainvalue[86:94]
        self.AUBIC, self.AIBUC, self.AUBIc, self.AIBUc, self.AUbIC, self.AIbUC, self.AUbIc, self.AIbUc, self.aUBIC, self.aIBUC, self.aUBIc, self.aIBUc, self.aUbIC, self.aIbUC, self.aUbIc, self.aIbUc=self.mainvalue[94::]

    def __call__(self, *args):
        if len(args)==0 and self.objEQ!='':
            return self.objEQ
            
        elif len(args) == 1:
            return self.getvalue(self.getname(args[0]))
             
    def setEQ(self,sexpr):
        self.objEQ=sexpr 
        display(Math('Eq = '+latex(sexpr)))
        
    def getname(self,obj):
        if type(obj)==list:
            return [self.getname(data) for data in obj]
        elif  type(obj)==Add:
            return [self.getname(data) for data in obj.args]
        else:
            
            if obj==1:
                return '1'
            elif obj==-1:
                return '0'
            else:    
                obj=signo(obj)*obj
                try: 
                    kpos=self.mainvalue.index(obj)
                    return self.mainname[kpos]
                except:
                    return str(obj)
    def gettype(self,obj):
        if type(obj)==list:
            return [self.gettype(data) for data in obj]
        elif type(obj)==Add:
            return [self.gettype(data) for data in obj.args]
        else:    
            pname=self.getname(obj)
            if pname =='1': 
                return '1'
            elif pname =='0':
                return '0'
            elif len(pname)==1:
                return 'P'
            elif len(pname)==2:
                return 'O'
            else:
                if pname[0]=='U':
                    return 'U'
                else:
                    return 'I'        
    def getvalue(self,sexpr):
        if type(sexpr)==Symbol:
            sexpr=self.getname(sexpr) 
        try:    
            kpos=self.mainname.index(sexpr)
        except:
            sexpr=sexpr[1:len(sexpr)-1]
            kpos=self.mainname.index(sexpr)
        return self.mainvalue[kpos]
        
    def setvalue(self,sexpr,value):
        if type(sexpr)==list:
            for data1,data2 in zip(sexpr,value):
                self.setvalue(data1,data2)
        else:        
            if type(sexpr)==Symbol:
                sexpr=self.getname(sexpr)
                 
                
            kpos=self.mainname.index(sexpr)
            self.mainvalue[kpos]=value
            self.update()
    
    def expandunion(self,sexpr):
        p1=self.getvalue(sexpr[0])
        p2=self.getvalue(sexpr[1])
        p3=self.getvalue('I'+sexpr)
        return p1+p2-p3

    def get_uvalue(self,*args):
        sname='U'
        for data in args:
            sname=sname+data
        return self.getvalue(sname)
    def get_Ivalue(self,*args):
        sname='I'
        for data in args:
            sname=sname+data
        return self.getvalue(sname) 

    def expand(self,obj=''):
        if obj=='':
            kres=self.expand(self.objEQ)
            self.objEQ=kres
            return self.objEQ
        else:    
            kkres=0
            if type(obj)==Add:
                kkres=0
                for data in obj.args:
                    kkres+=self.expand(data)
                return kkres
            elif type(obj)==Mul:
                kkres=1
                for data in obj.args:
                    kkres=kkres*self.expand(data)
                return kkres
            else:
                kres=0
                sname=self.getname(obj)
                ktype=['I','U']
                if sname[0] in ktype:
                    s1=sname[0]
                    sname=sname[1::]
                else:
                    s1='X'
                    
                if s1=='U':
                    if self.exclude:
                        for data in sname:
                            kres+=self.getvalue(data)
                        return kres    
                    else:
                        for data in sname:
                            kres=kres+self.getvalue(data)
                        if len(sname)==2:
                            kres=kres-self.getvalue('I'+sname)
                        else:    
                            for d1 in range(len(sname)-1):
                                for d2 in range(d1+1,len(sname)):
                                    nname='I'+sname[d1]+sname[d2]
                                    kres=kres-self.getvalue(nname)
                 
                            kres+=self.getvalue('I'+sname)
                        return kres
                elif s1=='I':
                    if self.exclude:
                        kres=1
                        for data in sname:
                            kres=kres*self.getvalue(data)
                        return kres
                    else:
                         
                        v1=self.getvalue(sname[0])
                        v2=self.getvalue(sname[1])
                        v4=self.getvalue('I'+sname[0]+sname[1])
                        if len(sname)==2:
                            return v1+v2-v4
                        else:
                            v3=self.getvalue(sname[2])
                            v5=self.getvalue('I'+sname[0]+sname[2])                        
                            v6=self.getvalue('I'+sname[1]+sname[2])            
                            v7=self.getvalue('I'+sname)
                            return v1+v2+v3-v4-v5-v6+v7
                elif s1=='X':
                    if self.exclude:
                        if len(sname)==2:
                            v1=self.getvalue(sname[0])
                            v2=self.getvalue(sname[1])
                            v3=self.getvalue('I'+sname[0]+sname[1])
                            return cfrac(v3,v2)
                            
                            
                    
                    else:
                        return obj
                    
                                             
                else:
                    return obj
            
            return obj
       
    def cexpand(self,obj=''):
        if obj=='':
            kres=self.cexpand(self.objEQ)
            self.objEQ=kres
            return kres
        else:    
            if type(obj)==Add:
                kres=0
                for data in obj.args:
                    kres=kres+self.cexpand(data)
                return kres
            elif type(obj)==Mul:
                kres=1
                for data in obj.args:
                    kres=kres*self.cexpand(data)
                return kres
            else:
                try:
                    sname=self.getname(obj)
                    if len(sname)==1:
                        if sname.islower():
                            return (1-self.getvalue(sname.upper()))
                        else:
                            return obj
                    else:
                        return obj
                except:
                    return obj
            return  obj
    def subs(self,sexpr1,sexpr2):
        kres=self.objEQ
        kres=kres.subs(sexpr1,sexpr2)
        self.objEQ=kres
        return kres
        
    def simplify(self,kres=''):
        if kres=='':
            kres=self.objEQ
            kres=simplify(kres) 
            self.objEQ=kres
            return kres
        else:
            return simplify(kres)
        
        
    def expandinter(self,sexpr,*args):
        p1=self.getvalue(sexpr[0])
        p2=self.getvalue(sexpr[1])
        if 'independient' in args:
            return p1*p2
        else:    
            p3=self.getvalue(sexpr[1]+sexpr[0])
        return p1*p3   

    def expandif(self,sexpr):
        p1 = self.getvalue(sexpr[1])
        p2 = self.getvalue('I'+sexpr[1]+sexpr[0])
        return p1*p2

         
    def eQ(self,*args):
        if self.objEQ=='':
            display(Math('clss.eQ are empty, please set class.eQ'))
            display(Math('Exampe :   P.eQ=3*P.A+15*P.IAB'))
        elif len(args)==0:
            e1=MyEq(self.objEQ,'eQ')
        else:
            vecvar=[]
            vecval=[]
            for data in args:
                if type(data)==str:
                    vecvar.append(data) 
                else:
                    vecval.append(data) 
            kres=self.eQ
 
            for data1,data2 in zip(vecvar,vecval):
                kres=kres.subs(self(data1),data2)
            MQ(self.objEQ,kres)    
            return kres        


 
import numpy as np
from PIL import Image
import IPython.display
 

 
from IPython.display import clear_output  #  clear_output(wait=True)
 

# this cell will be initializing your variables
#from lib_MyFunctions import *
 
import sys
sys.path.append('D:/Libaldo/img/')
XX = np.asarray(Image.open('D:/Libaldo/img/xfondo.png'))
YY = np.asarray(Image.open('D:/Libaldo/img/yfondo.png'))

X0 = np.asarray(Image.open('D:/Libaldo/img/x0.png'))
X1 = np.asarray(Image.open('D:/Libaldo/img/x1.png'))
X2 = np.asarray(Image.open('D:/Libaldo/img/x2.png'))
X3 = np.asarray(Image.open('D:/Libaldo/img/x3.png'))
X4 = np.asarray(Image.open('D:/Libaldo/img/x4.png'))
X5 = np.asarray(Image.open('D:/Libaldo/img/x5.png'))
X6 = np.asarray(Image.open('D:/Libaldo/img/x6.png'))
X7 = np.asarray(Image.open('D:/Libaldo/img/x7.png'))  

Y0 = np.asarray(Image.open('D:/Libaldo/img/y0.png'))
Y1 = np.asarray(Image.open('D:/Libaldo/img/y1.png'))
Y2 = np.asarray(Image.open('D:/Libaldo/img/y2.png'))
Y3 = np.asarray(Image.open('D:/Libaldo/img/y3.png'))

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