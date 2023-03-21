#The family SO((n-1)(2n+1))/Sp(n), n>=3.

def family(n):
    Hstr="C"+str(n)
    if n%2==0:
        Gstr="B"+str(n^2-n/2-1)
    else:
        Gstr="D"+str(n^2-n/2-1/2)
    b=branching_rule(Gstr,Hstr+"(0,1"+(n-2)*",0"+")","plethysm")
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    print("G="+Gstr+", H="+Hstr)
    Sym2LowerBounds(G,H,b)

#for n in range(3,21): family(n)
