#SU(kn)/S(kU(n)), k>=3, n>=2.

def family(k,n):
    G=WeylCharacterRing("A"+str(k*n-1),style="coroots")
    Hss=WeylCharacterRing("A"+str(n-1)+(k-1)*("xA"+str(n-1)),style="coroots")
    hbasis=[list(x) for x in matrix.identity(k*n-1)]
    for i in range(1,k): hbasis.pop((k-i)*n-1)
    LiEG=LiEgroupfromWCR(G)
    LiEhbasis=lie(hbasis)
    rm=[list(x) for x in lie.res_mat(LiEhbasis,LiEG).sage()]
    Sym2LowerBoundsWithTorus(G,Hss,rm)

r=20
table=[[[k,n,k*n] for k in range(3,r+1)] for n in range(2,r+1)]
iterate=sorted(flatten(table,max_level=1),key=lambda x: x[2])
truncated=[x for x in iterate if x[2]>11]
for x in truncated: family(x[0],x[1])
