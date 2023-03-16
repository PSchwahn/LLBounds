#SO(3n+2)/(U(n+1)SO(n)), n>=3.

def family(n):
    print("n="+str(n))
    if n%2==1:
        rankG=int((3*n+1)/2)
        G=WeylCharacterRing("B"+str(rankG),style="coroots")
        if n==3:
            Hss=WeylCharacterRing("A"+str(n)+"xA1",style="coroots")
        else:
            Hss=WeylCharacterRing("A"+str(n)+"xB"+str((n-1)/2),style="coroots")
    else:
        rankG=int(3*n/2+1)
        G=WeylCharacterRing("D"+str(rankG),style="coroots")
        if n==4:
            Hss=WeylCharacterRing("A"+str(n)+"xA1xA1",style="coroots")
        else:
            Hss=WeylCharacterRing("A"+str(n)+"xD"+str(n/2),style="coroots")
    hbasis=[list(x) for x in matrix.identity(rankG)]
    hbasis.pop(n)
    LiEG=LiEgroupfromWCR(G)
    LiEhbasis=lie(hbasis)
    rm=[list(x) for x in lie.res_mat(LiEhbasis,LiEG).sage()]
    Sym2LowerBoundsWithTorus(G,Hss,rm)

for n in range(17,41): family(n)
