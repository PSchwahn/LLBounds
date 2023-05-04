#The family SO(n(n-1)/2)/SO(n), n>=7.

def family(n):
    if n%2==0:
        k=n/2
        Hstr="D"+str(k)
        if k%2==0:
            Gstr="D"+str(k^2-k/2)
        else:
            Gstr="B"+str((2*k^2-k-1)/2)
    else:
        k=(n-1)/2
        Hstr="B"+str(k)
        if k%2==0:
            Gstr="D"+str((2*k+1)*k/2)
        else:
            Gstr="B"+str(((2*k+1)*k-1)/2)
    o1="1X[0,1"+(k-2)*",0"+"]"
    Sym2LowerBoundsBig(Gstr,Hstr,o1)
