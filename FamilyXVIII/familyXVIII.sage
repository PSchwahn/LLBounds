#The family SO(3n+2)/(U(n+1)SO(n)), n>=3.

def family(n):
    print("n="+str(n))
    if n%2==1:
        rankG=int((3*n+1)/2)
        Gstr="B"+str(rankG)
        if n==3:
            Hssstr="A"+str(n)+"xA1"
        else:
            Hssstr="A"+str(n)+"xB"+str((n-1)/2)
    else:
        rankG=int(3*n/2+1)
        Gstr="D"+str(rankG)
        if n==4:
            Hssstr="A"+str(n)+"xA1xA1"
        else:
            Hssstr="A"+str(n)+"xD"+str(n/2)
    G=WeylCharacterRing(Gstr,style="coroots")
    Hss=WeylCharacterRing(Hssstr,style="coroots")
    hbasis=[list(x) for x in matrix.identity(rankG)]
    hbasis.pop(n)
    LiEG=LiEgroupfromWCR(G)
    LiEhbasis=lie(hbasis)
    rm=[list(x) for x in lie.res_mat(LiEhbasis,LiEG).sage()]
    print("G="+Gstr+", Hss="+Hssstr)
    Sym2LowerBoundsWithTorus(G,Hss,rm)

#for n in range(3,41): family(n)
