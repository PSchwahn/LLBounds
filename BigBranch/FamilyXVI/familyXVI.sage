#The family SO(n^2)/2SO(n), n>=3.

def family(n):
    if n==3:
        Hstr="A1xA1"
        altHstr="B1xB1"
        Gstr="B4"
        o1="1X[2,2]"
    elif n==4:
        Hstr="A1xA1xA1xA1"
        altHstr="D2xD2"
        Gstr="D8"
        o1="1X[1,1,1,1]"
    elif n%2==0:
        k=n/2
        Hstr="D"+str(k)+"xD"+str(k)
        Gstr="D"+str(n^2/2)
        o1="1X[1"+(k-1)*",0"+",1"+(k-1)*",0"+"]"
    else:
        k=(n-1)/2
        Hstr="B"+str(k)+"xB"+str(k)
        Gstr="B"+str((n^2-1)/2)
        o1="1X[1"+(k-1)*",0"+",1"+(k-1)*",0"+"]"
    Sym2LowerBoundsBig(Gstr,Hstr,o1)
