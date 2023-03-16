def family(n):
    Gstr="A"+str(n*(n+1)/2-1)
    Hstr="A"+str(n-1)
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    b=branching_rule(Gstr,Hstr+"(2"+(n-2)*",0"+")","plethysm")
    Sym2LowerBounds(G,H,b)
    
for n in range(9,21): family(n)
