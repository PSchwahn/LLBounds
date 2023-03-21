#The family SU(n(n-1)/2)/SU(n), n>=5.

def family(n):
    Gstr="A"+str(n*(n-1)/2-1)
    Hstr="A"+str(n-1)
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    b=branching_rule(Gstr,Hstr+"(0,1"+(n-3)*",0"+")","plethysm")
    print("G="+Gstr+", H="+Hstr)
    Sym2LowerBounds(G,H,b)
    
#for n in range(5,21): family(n)
