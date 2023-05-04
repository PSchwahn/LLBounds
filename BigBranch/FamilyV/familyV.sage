#The family SO(n^2-1)/SU(n), n>=3.

def family(n):
    if n%2==0:
        Gstr="B"+str(n^2/2-1)
    else:
        Gstr="D"+str((n^2-1)/2)
    Hstr="A"+str(n-1)
    o1="1X[1"+(n-3)*",0"+",1]"
    Sym2LowerBoundsBig(Gstr,Hstr,o1)
