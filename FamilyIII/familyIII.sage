#The family SU(pq)/(SU(p)SU(q)), 2<=p<=q, (p,q)!=(2,2).

def family(p,q):
    Gstr="A"+str(p*q-1)
    Hstr="A"+str(p-1)+"xA"+str(q-1)
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    b=branching_rule(Gstr,Hstr,"tensor")
    print("G="+Gstr+", H="+Hstr)
    Sym2LowerBounds(G,H,b)
        
#r=20
#table=[[[p,q,p*q-1] for q in range(p,r+1)] for p in range(2,r+1)]
#iterate=sorted(flatten(table,max_level=1),key=lambda x: x[2])
#for x in iterate: family(x[0],x[1])
