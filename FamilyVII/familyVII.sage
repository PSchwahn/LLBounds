#The family SO(2n^2+n)/Sp(n), n>=2.

def family(n):
    Hstr="C"+str(n)
    if n%2==0:
        Gstr="D"+str(n^2+n/2)
    else:
        Gstr="B"+str(n^2+n/2-1/2)
    b=branching_rule(Gstr,Hstr+"(2"+(n-1)*",0"+")","plethysm")
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    print("G="+Gstr+", H="+Hstr)
    Sym2LowerBounds(G,H,b)

#for n in range(2,21): family(n)
