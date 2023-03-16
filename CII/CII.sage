#sp(nk)/nsp(k), k>=1, n>=3
def block0(k):
    return zero_matrix(k)
def block1(k):
    return identity_matrix(k)
def block2(k):
    row=[0 for i in range(0,k-1)]+[1]
    return Matrix([row for i in range(0,k)])
def rmCII(k,n):
    bmlist=[]
    for i in range(0,n):
        bmlist+=[[block0(k) for j in range(0,n-i-1)]+[block1(k)]+[block2(k) for j in range(n-i,n)]]
    return [list(x) for x in block_matrix(bmlist)]
def family(k,n):
    Gstr="C"+str(n*k)
    if k==1:
        Hstr="A1"+(n-1)*"xA1"
        altHstr="C1"+(n-1)*"xC1"
        b1=branching_rule("C1","A1","isomorphic")
        b=branching_rule(Gstr,altHstr,"orthogonal_sum")*branching_rule(altHstr,Hstr,[b1 for i in range(0,n)])
    else:
        Hstr="C"+str(k)+(n-1)*("xC"+str(k))
        b=rmCII(k,n)
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    print("k="+str(k)+", n="+str(n))
    print("G="+Gstr)
    print("H="+Hstr)
    Sym2LowerBounds(G,H,b)


r=20
table=[[[k,n,k*n] for k in range(1,r+1)] for n in range(3,r+1)]
iterate=sorted(flatten(table,max_level=1),key=lambda x: x[2])
truncated=[x for x in iterate if x[2]>=11]
for x in truncated: family(x[0],x[1])
