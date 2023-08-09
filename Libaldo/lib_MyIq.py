


class MyInequality (MyEqEq):
    def __init__(self,xl1,s1,xx,s2,xl2,kshow=True):
        self.type='Iq'
        self.ksymL=xl1
        self.ksy=xx
        self.ksymR=xl2
        self.simb1=s1    
        self.simb2=s2
        
        display(Math(latex(xl1)+s1+latex(xx+)+s2+latex(xl2))
