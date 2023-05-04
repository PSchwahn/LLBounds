#The family SO(4n^2)/2Sp(n), n>=2.

def family(n):
    Hstr="C"+str(n)+"xC"+str(n)
    Gstr="D"+str(2*n^2)
    b=branching_rule(Gstr,Hstr,"tensor")
    o1="1X[1"+(n-1)*",0"+",1"+(n-1)*",0"+"]"
    Sym2LowerBoundsBig(Gstr,Hstr,o1)
