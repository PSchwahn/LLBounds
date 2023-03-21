#The family Sp(kn)/kSp(n), n>=1, k>=3.

def block0(n):
    return zero_matrix(n)
def block1(n):
    return identity_matrix(n)
def block2(n):
    row=[0 for i in range(0,n-1)]+[1]
    return Matrix([row for i in range(0,n)])
def rm(k,n):
    bmlist=[]
    for i in range(0,k):
        bmlist+=[[block0(n) for j in range(0,k-i-1)]+[block1(n)]+[block2(n) for j in range(k-i,k)]]
    return [list(x) for x in block_matrix(bmlist)]
def family(k,n):
    Gstr="C"+str(k*n)
    if n==1:
        Hstr="A1"+(k-1)*"xA1"
        altHstr="C1"+(k-1)*"xC1"
        b1=branching_rule("C1","A1","isomorphic")
        b=branching_rule(Gstr,altHstr,"orthogonal_sum")*branching_rule(altHstr,Hstr,[b1 for i in range(0,k)])
    else:
        Hstr="C"+str(n)+(k-1)*("xC"+str(n))
        b=rm(k,n)
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    print("k="+str(k)+", n="+str(n))
    print("G="+Gstr+", H="+Hstr)
    Sym2LowerBounds(G,H,b)

#r=20
#table=[[[k,n,k*n] for k in range(3,r+1)] for n in range(1,r+1)]
#iterate=sorted(flatten(table,max_level=1),key=lambda x: x[2])
#for x in iterate: family(x[0],x[1])