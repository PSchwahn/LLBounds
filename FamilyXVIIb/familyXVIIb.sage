#The family SO(kn)/kSO(n), k,n>=3

def family(k,n):
    if n==3:
        Hstr="A1"+(k-1)*"xA1"
        altHstr="B1"+(k-1)*"xB1"
        if k%2==0:
            Gstr="D"+str(3*k/2)
        else:
            Gstr="B"+str((3*k-1)/2)
        b1=branching_rule(Gstr,altHstr,"orthogonal_sum")
        b2=branching_rule("B1","A1","isomorphic")
        b3=branching_rule(altHstr,Hstr,[b2 for i in range(0,k)])
        b=b1*b3
    elif n==4:
        Hstr="A1"+(2*k-1)*"xA1"
        Gstr="D"+str(2*k)
        altHstr="D2"+(k-1)*"xD2"
        b1=branching_rule(Gstr,altHstr,"orthogonal_sum")
        b2=branching_rule("D2","A1xA1","isomorphic")
        b3=branching_rule(altHstr,Hstr,[b2 for i in range(0,k)])
        b=b1*b3
    elif n%2==0:
        l=n/2
        Hstr="D"+str(l)+(k-1)*("xD"+str(l))
        Gstr="D"+str(k*l)
        b=branching_rule(Gstr,Hstr,"orthogonal_sum")
    else:
        l=(n-1)/2
        Hstr="B"+str(l)+(k-1)*("xB"+str(l))
        if k%2==0:
            Gstr="D"+str(n*k/2)
        else:
            Gstr="B"+str((n*k-1)/2)
        b=branching_rule(Gstr,Hstr,"orthogonal_sum")
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    print("k="+str(k)+", n="+str(n))
    print("G="+Gstr+", H="+Hstr)
    Sym2LowerBounds(G,H,b)

#r=20
#table=[[[k,n,k*n] for k in range(3,r+1)] for n in range(3,r+1)]
#iterate=sorted(flatten(table,max_level=1),key=lambda x: x[2])
#for x in iterate: family(x[0],x[1])