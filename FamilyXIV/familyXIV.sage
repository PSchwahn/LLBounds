#The family Sp(3n-1)/(U(2n-1)Sp(n)), n>=1.

def family(n):
    Gstr="C"+str(3*n-1)
    G=WeylCharacterRing(Gstr,style="coroots")
    if n==1:
        Hssstr="A1"
    else:
        Hssstr="A"+str(2*n-2)+"xC"+str(n)
    Hss=WeylCharacterRing(Hssstr,style="coroots")
    hbasis=[list(x) for x in matrix.identity(3*n-1)]
    hbasis.pop(2*n-2)
    LiEG=LiEgroupfromWCR(G)
    LiEhbasis=lie(hbasis)
    rm=[list(x) for x in lie.res_mat(LiEhbasis,LiEG).sage()]
    print("G="+Gstr+", Hss="+Hssstr)
    Sym2LowerBoundsWithTorus(G,Hss,rm)

#for n in range(1,21): family(n)
