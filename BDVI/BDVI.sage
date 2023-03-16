#so(4n)/(sp(1)+sp(n)), n>=2
def family(n):
    Gstr="D"+str(2*n)
    Hstr="A1xC"+str(n)
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    b=branching_rule(Gstr,Hstr+"(1,1"+(n-1)*",0"+")","plethysm")
    print("G="+Gstr)
    print("H="+Hstr)
    Sym2LowerBounds(G,H,b)

for n in range(17,31): family(n)
