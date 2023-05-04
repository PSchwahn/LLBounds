#The family SO((n-1)(n+2)/2)/SO(n), n>=5.

def family(n):
    if n%2==0:
        k=n/2
        Hstr="D"+str(k)
        if k%2==0:
            Gstr="B"+str(k^2+k/2-1)
        else:
            Gstr="D"+str(k^2+k/2-1/2)
    else:
        k=(n-1)/2
        Hstr="B"+str(k)
        if k%2==0:
            Gstr="D"+str((2*k+3)*k/2)
        else:
            Gstr="B"+str(((2*k+3)*k-1)/2)
    o1="1X[2"+(k-1)*",0"+"]"
    Sym2LowerBoundsBig(Gstr,Hstr,o1)
