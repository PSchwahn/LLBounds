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
    hbasis=[list(x) for x in matrix.identity(rankG)]
    hbasis.pop(n)
    LiEG=lie(Gstr)
    LiEH=lie(Hssstr.replace("x","")+"T1")
    LiEhbasis=lie(hbasis)
    print("G: "+Gstr+", H: "+Hssstr+"T1")
    rm=[list(x) for x in lie.res_mat(LiEhbasis,LiEG).sage()]
    LiErm=listtoliebyhalf(rm)
    LiErep=listtoliebyhalf([1]+(rankG-1)*[0])
    o1=str(lie.branch(LiErep,LiEH,LiErm,LiEG)).replace(" ","").replace("\n","")
    print("o1: "+o1)
    Sym2LowerBoundsBigWithTorus(Gstr,Hssstr,o1,rm)
