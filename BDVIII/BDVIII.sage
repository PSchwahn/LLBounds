#so(n²)/2so(n), n>=3
def family(n):
    if n==3:
        Hstr="A1xA1"
        altHstr="B1xB1"
        Gstr="B4"
        b1=branching_rule("B4","B1xB1","tensor")
        b2=branching_rule("B1","A1","isomorphic")
        b=b1*branching_rule(altHstr,Hstr,[b2,b2])
    elif n==4:
        Hstr="A1xA1xA1xA1"
        altHstr="D2xD2"
        Gstr="D8"
        b1=branching_rule("D8","D2xD2","tensor")
        b2=branching_rule("D2","A1xA1","isomorphic")
        b=b1*branching_rule(altHstr,Hstr,[b2,b2])
    elif n%2==0:
        Hstr="D"+str(n/2)+"xD"+str(n/2)
        Gstr="D"+str(n^2/2)
        b=branching_rule(Gstr,Hstr,"tensor")
    else:
        Hstr="B"+str((n-1)/2)+"xB"+str((n-1)/2)
        Gstr="B"+str((n^2-1)/2)
        b=branching_rule(Gstr,Hstr,"tensor")
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    print("G="+Gstr)
    print("H="+Hstr)
    Sym2LowerBounds(G,H,b)

for n in range(9,21): family(n)
