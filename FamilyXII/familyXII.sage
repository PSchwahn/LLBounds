#The family SU(l+pq)/S(U(l)U(q)U(q)), 2<=p<=q, lpq=p^2+q^2+1.

def family(p,q,l):
    Gstr="A"+str(p*q+l-1)
    G1str="A"+str(l-1)+"xA"+str(p*q-1)
    G2str="A"+str(p*q-1)
    G3str="A"+str(p-1)+"xA"+str(q-1)
    Hssstr="A"+str(l-1)+"xA"+str(p-1)+"xA"+str(q-1)
    G=WeylCharacterRing(Gstr,style="coroots")
    Hss=WeylCharacterRing(Hssstr,style="coroots")
    #rm1 branches from SU(pq) to SU(p)SU(q)
    G2=WeylCharacterRing(G2str,style="coroots")
    G3=WeylCharacterRing(G3str,style="coroots")
    b=branching_rule(G2str,G3str,"tensor")
    rm1=Matrix(resmatrix(G2,G3,b))
    #rm2 branches from U(l)SU(pq) to U(l)SU(p)SU(q)
    rm2=block_matrix([[matrix.identity(int(l-1)),0,0],[0,rm1,0],[0,0,matrix.identity(1)]])
    #rm3 branches from SU(l+pq) to S(U(l)U(pq))
    basis=[list(x) for x in matrix.identity(int(p*q+l-1))]
    basis.pop(int(l-1))
    LiEbasis=lie(basis)
    LiEbasis=lie("id("+str(p*q+l-1)+")-"+str(l))
    LiEG=LiEgroupfromWCR(G)
    LiErm3=lie.res_mat(LiEbasis,LiEG)
    rmat=LiErm3.sage()*rm2
    rm=[list(x) for x in rmat]
    print("G="+Gstr+", Hss="+Hssstr)
    Sym2LowerBoundsWithTorus(G,H,rm)

#r=40
#table=[[[p,q,(p^2+q^2+1)/(p*q),p*q+(p^2+q^2+1)/(p*q)] for q in range(p,r+1)] for p in range(2,r+1)]
#truncated=[x for x in flatten(table,max_level=1) if x[2].is_integer()]
#iterate=sorted(truncated,key=lambda x: x[3])
#for x in iterate: family(x[0],x[1],x[2])

#the first three spaces:
#family(2,5,3)
#family(5,13,3)
#family(13,34,3)
