#The family Sp(n)/(Sp(1)SO(n)), n>=3.

def family(n):
    Gstr="C"+str(n)
    if n==3:
        Hstr="A1xA1"
        b=branching_rule("C3","A1xA1(1,2)","plethysm")
    elif n==4:
        Hstr="A1xA1xA1"
        b=branching_rule("C4","A1xA1xA1(1,1,1)","plethysm")
    elif n%2==0:
        k=n/2
        Hstr="A1xD"+str(k)
        b=branching_rule(Gstr,Hstr+"(1,1"+(k-1)*",0"+")","plethysm")
    else:
        k=(n-1)/2
        Hstr="A1xB"+str(k)
        b=branching_rule(Gstr,Hstr+"(1,1"+(k-1)*",0"+")","plethysm")
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    Sym2LowerBounds(G,H,b)

#for n in range(3,41): family(n)
