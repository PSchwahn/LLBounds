#The family SO(n^2-1)/SU(n), n>=3.

def family(n):
    if n%2==0:
        Gstr="B"+str(n^2/2-1)
        Hstr="A"+str(n-1)
        G=WeylCharacterRing(Gstr,style="coroots")
        H=WeylCharacterRing(Hstr,style="coroots")
        b=branching_rule(Gstr,Hstr+"(1"+(n-3)*",0"+",1)","plethysm")
    else:
        Gstr="D"+str((n^2-1)/2)
        Hstr="A"+str(n-1)
        G=WeylCharacterRing(Gstr,style="coroots")
        H=WeylCharacterRing(Hstr,style="coroots")
        b=branching_rule(Gstr,Hstr+"(1"+(n-3)*",0"+",1)","plethysm")
    print("G="+Gstr+", H="+Hstr)
    Sym2LowerBounds(G,H,b)

#for n in range(3,21): family(n)
