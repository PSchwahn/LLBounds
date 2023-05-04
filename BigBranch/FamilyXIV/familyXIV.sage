#The family Sp(3n-1)/(U(2n-1)Sp(n)), n>=1.

def family(n):
    Gstr="C"+str(3*n-1)
    if n==1:
        Hssstr="A1"
    else:
        Hssstr="A"+str(2*n-2)+"xC"+str(n)
    hbasis=[list(x) for x in matrix.identity(3*n-1)]
    hbasis.pop(2*n-2)
    LiEG=lie(Gstr)
    LiEH=lie(Hssstr.replace("x","")+"T1")
    LiEhbasis=lie(hbasis)
    print("G: "+Gstr+", H: "+Hssstr+"T1")
    rm=[list(x) for x in lie.res_mat(LiEhbasis,LiEG).sage()]
    LiErm=listtoliebyhalf(rm)
    LiErep=listtoliebyhalf([1]+(3*n-2)*[0])
    o1=str(lie.branch(LiErep,LiEH,LiErm,LiEG)).replace(" ","").replace("\n","")
    print("o1: "+o1)
    Sym2LowerBoundsBigWithTorus(Gstr,Hssstr,o1,rm)
