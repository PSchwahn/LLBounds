#The family SO(n(n-1)/2)/SO(n), n>=7.

def family(n):
    if n%2==0:
        k=n/2
        Hstr="D"+str(k)
        if k%2==0:
            Gstr="D"+str(k^2-k/2)
        else:
            Gstr="B"+str((2*k^2-k-1)/2)
        b=branching_rule(Gstr,Hstr+"(0,1"+(k-2)*",0"+")","plethysm")
    else:
        k=(n-1)/2
        Hstr="B"+str(k)
        if k%2==0:
            Gstr="D"+str((2*k+1)*k/2)
        else:
            Gstr="B"+str(((2*k+1)*k-1)/2)
        b=branching_rule(Gstr,Hstr+"(0,1"+(k-2)*",0"+")","plethysm")
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    Sym2LowerBounds(G,H,b)

#for n in range(7,21): family(n)
