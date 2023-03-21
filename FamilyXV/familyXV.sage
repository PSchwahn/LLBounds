#The family SO(4n^2)/2Sp(n), n>=2.

def family(n):
    Hstr="C"+str(n)+"xC"+str(n)
    Gstr="D"+str(2*n^2)
    b=branching_rule(Gstr,Hstr,"tensor")
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    Sym2LowerBounds(G,H,b)

#for n in range(2,21): family(n)
