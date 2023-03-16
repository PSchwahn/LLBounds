#Sp(3n-1)/(U(2n-1)Sp(n)), n>=1.

def family(n):
    print("n="+str(n))
    G=WeylCharacterRing("C"+str(3*n-1),style="coroots")
    if n==1:
        Hss=WeylCharacterRing("A1",style="coroots")
    else:
        Hss=WeylCharacterRing("A"+str(2*n-2)+"xC"+str(n),style="coroots")
    hbasis=[list(x) for x in matrix.identity(3*n-1)]
    hbasis.pop(2*n-2)
    LiEG=LiEgroupfromWCR(G)
    LiEhbasis=lie(hbasis)
    rm=[list(x) for x in lie.res_mat(LiEhbasis,LiEG).sage()]
    Sym2LowerBoundsWithTorus(G,Hss,rm)

for n in range(10,21): family(n)
