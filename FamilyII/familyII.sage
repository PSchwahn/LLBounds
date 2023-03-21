#The family SU(n(n+1)/2)/SU(n), n>=3.

def family(n):
    Gstr="A"+str(n*(n+1)/2-1)
    Hstr="A"+str(n-1)
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    b=branching_rule(Gstr,Hstr+"(2"+(n-2)*",0"+")","plethysm")
    print("G="+Gstr+", H="+Hstr)
    Sym2LowerBounds(G,H,b)
    
#for n in range(3,21): family(n)
