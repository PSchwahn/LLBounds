#The family SO(2n^2+n)/Sp(n), n>=2.

def family(n):
    Hstr="C"+str(n)
    if n%2==0:
        Gstr="D"+str(n^2+n/2)
    else:
        Gstr="B"+str(n^2+n/2-1/2)
    o1="1X[2"+(n-1)*",0"+"]"
    Sym2LowerBoundsBig(Gstr,Hstr,o1)
