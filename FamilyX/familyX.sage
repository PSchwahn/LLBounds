#The family SO((n-1)(n+2)/2)/SO(n), n>=5.

def family(n):
    if n%2==0:
        k=n/2
        Hstr="D"+str(k)
        if k%2==0:
            Gstr="B"+str(k^2+k/2-1)
        else:
            Gstr="D"+str(k^2+k/2-1/2)
        b=branching_rule(Gstr,Hstr+"(2"+(k-1)*",0"+")","plethysm")
    else:
        k=(n-1)/2
        Hstr="B"+str(k)
        if k%2==0:
            Gstr="D"+str((2*k+3)*k/2)
        else:
            Gstr="B"+str(((2*k+3)*k-1)/2)
        b=branching_rule(Gstr,Hstr+"(2"+(k-1)*",0"+")","plethysm")
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    Sym2LowerBounds(G,H,b)

#for n in range(5,21): family(n)
