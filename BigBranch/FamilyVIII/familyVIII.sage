#The family SO(4n)/(Sp(1)Sp(n)), n>=2.

def family(n):
    Gstr="D"+str(2*n)
    Hstr="A1xC"+str(n)
    o1="1X[1,1"+(n-1)*",0"+"]"
    Sym2LowerBoundsBig(Gstr,Hstr,o1)
