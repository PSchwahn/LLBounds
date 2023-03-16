#so(nk)/nso(k), k,n>=3
def family(k,n):
    if k==3:
        Hstr="A1"+(n-1)*"xA1"
        altHstr="B1"+(n-1)*"xB1"
        if n%2==0:
            Gstr="D"+str(3*n/2)
        else:
            Gstr="B"+str((3*n-1)/2)
        b1=branching_rule(Gstr,altHstr,"orthogonal_sum")
        b2=branching_rule("B1","A1","isomorphic")
        b3=branching_rule(altHstr,Hstr,[b2 for i in range(0,n)])
        b=b1*b3
    elif k==4:
        Hstr="A1"+(2*n-1)*"xA1"
        Gstr="D"+str(2*n)
        altHstr="D2"+(n-1)*"xD2"
        b1=branching_rule(Gstr,altHstr,"orthogonal_sum")
        b2=branching_rule("D2","A1xA1","isomorphic")
        b3=branching_rule(altHstr,Hstr,[b2 for i in range(0,n)])
        b=b1*b3
    elif k%2==0:
        l=k/2
        Hstr="D"+str(l)+(n-1)*("xD"+str(l))
        Gstr="D"+str(n*l)
        b=branching_rule(Gstr,Hstr,"orthogonal_sum")
    else:
        l=(k-1)/2
        Hstr="B"+str(l)+(n-1)*("xB"+str(l))
        if n%2==0:
            Gstr="D"+str(k*n/2)
        else:
            Gstr="B"+str((k*n-1)/2)
        b=branching_rule(Gstr,Hstr,"orthogonal_sum")
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    print("k="+str(k)+", n="+str(n))
    print("G="+Gstr)
    print("H="+Hstr)
    Sym2LowerBounds(G,H,b)


r=20
table=[[[k,n,k*n] for k in range(4,r+1)] for n in range(3,r+1)]
iterate=sorted(flatten(table,max_level=1),key=lambda x: x[2])
truncated=[x for x in iterate if x[2]>40]
#family(3,10)
family(4,8)
#family(3,11)
#family(3,12)
family(4,9)
#family(3,13)
#family(3,14)
family(4,10)
for x in truncated: family(x[0],x[1])
