#The family SO((n-1)(2n+1))/Sp(n), n>=3.

def family(n):
    Hstr="C"+str(n)
    if n%2==0:
        Gstr="B"+str(n^2-n/2-1)
    else:
        Gstr="D"+str(n^2-n/2-1/2)
    o1="1X[0,1"+(n-2)*",0"+"]"
    Sym2LowerBoundsBig(Gstr,Hstr,o1)
