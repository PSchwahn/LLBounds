#SO(p)/H where p is the isotropy rep. of a symm. space K/H.
#Here: up to dim p = 26.

todolist=[]

#S³xS³
todolist+=[[["A1",[2,0]],["A1",[2,0]]]]

#S⁴xS⁴
todolist+=[[["A1xA1",[1,1]],["A1xA1",[1,1]]]]

#3S³
todolist+=[3*[["A1",[2,0]]]]

#SU(3)/SO(3)xSU(3)/SO(3)
todolist+=[2*[["A1",[4,0]]]]

#S⁵xS⁵
todolist+=[[["B2",[1,0]],["B2",[1,0]]]]

#S³xSU(3)
todolist+=[[["A1",[2,0]],["A2",[1,1]]]]

#4S³
todolist+=[4*[["A1",[2,0]]]]

#3S⁴
todolist+=[3*[["A1xA1",[1,1]]]]

#S³xSO(5)
todolist+=[[["A1",[2,0]],["B2",[0,2]]]]

#S³xS³xSU(3)
todolist+=[[["A1",[2,0]],["A1",[2,0]],["A2",[1,1]]]]

#(SU(3)/SO(3))³
todolist+=[3*[["A1",[4,0]]]]

#5S³
todolist+=[5*[["A1",[2,0]]]]

#3S⁵
todolist+=[3*[["B2",[1,0]]]]

#SU(3)xSU(3)
todolist+=[[["A2",[1,1]],["A2",[1,1]]]]

#S³xS³xSO(5)
todolist+=[[["A1",[2,0]],["A1",[2,0]],["B2",[0,2]]]]

#4S⁴
todolist+=[4*[["A1xA1",[1,1]]]]

#S³xG2
todolist+=[[["A1",[2,0]],["G2",[0,1]]]]

#3S³xSU(3)
todolist+=[3*[["A1",[2,0]]]+[["A2",[1,1]]]]

#(SO(6)/2SO(3))² or (SU(4)/SO(4))² or SO(6)/2SO(3)xSU(4)/SO(4)
todolist+=[[["A1xA1",[2,2]],["A1xA1",[2,2]]]]

#S³xSU(4)
todolist+=[[["A1",[2,0]],["A3",[1,0,1]]]]

#SU(3)xSO(5)
todolist+=[[["A2",[1,1]],["B2",[0,2]]]]

#6S³
todolist+=[6*[["A1",[2,0]]]]

#S⁴xSU(6)/Sp(3)
todolist+=[[["A1xA1",[1,1]],["C3",[0,1,0]]]]

#3S³xSO(5)
todolist+=[3*[["A1",[2,0]]]+[["B2",[0,2]]]]

#S³xSU(3)xSU(3)
todolist+=[[["A1",[2,0]],["A2",[1,1]],["A2",[1,1]]]]

#(SU(3)/SO(3))⁴
todolist+=[4*[["A1",[4,0]]]]

#S³xS³xG2
todolist+=[[["A1",[2,0]],["A1",[2,0]],["G2",[0,1]]]]

#SO(5)xSO(5)
todolist+=[[["B2",[0,2]],["B2",[0,2]]]]

#4S³xSU(3)
todolist+=[4*[["A1",[2,0]]]+[["A2",[1,1]]]]

#5S⁴
todolist+=[5*[["A1xA1",[1,1]]]]

#4S⁵
todolist+=[4*[["B2",[1,0]]]]

#S³xS³xSU(4)
todolist+=[[["A1",[2,0]],["A1",[2,0]],["A3",[1,0,1]]]]

#S³xSU(3)xSO(5)
todolist+=[[["A1",[2,0]],["A2",[1,1]],["B2",[0,2]]]]

#7S³
todolist+=[7*[["A1",[2,0]]]]

#SU(3)xG2
todolist+=[[["A2",[1,1]],["G2",[0,1]]]]

#4S³xSO(5)
todolist+=[4*[["A1",[2,0]]]+[["B2",[0,2]]]]

#S³xS³xSU(3)xSU(3)
todolist+=[2*[["A1",[2,0]]]+2*[["A2",[1,1]]]]

#S⁴xS⁴xSU(6)/Sp(3)
todolist+=[[["A1xA1",[1,1]],["A1xA1",[1,1]],["C3",[0,1,0]]]]

#SU(3)xSU(4)
todolist+=[[["A2",[1,1]],["A3",[1,0,1]]]]

#3S³xG2
todolist+=[3*[["A1",[2,0]]]+[["G2",[0,1]]]]

#S³xSO(5)xSO(5)
todolist+=[[["A1",[2,0]],["B2",[0,2]],["B2",[0,2]]]]

#5S³xSU(3)
todolist+=[5*[["A1",[2,0]]]+[["A2",[1,1]]]]

#S³xSO(7)
todolist+=[[["A1",[2,0]],["B3",[0,1,0]]]]

#S³xSp(3)
todolist+=[[["A1",[2,0]],["C3",[2,0,0]]]]

#3S³xSU(4)
todolist+=[3*[["A1",[2,0]]]+[["A3",[1,0,1]]]]

#SO(5)xG2
todolist+=[[["B2",[0,2]],["G2",[0,1]]]]

#S³xS³xSU(3)xSO(5)
todolist+=[[["A1",[2,0]],["A1",[2,0]],["A2",[1,1]],["B2",[0,2]]]]

#3SU(3)
todolist+=[3*[["A2",[1,1]]]]

#8S³
todolist+=[8*[["A1",[2,0]]]]

#6S⁴
todolist+=[6*[["A1xA1",[1,1]]]]

#(SU(3)/SO(3))⁵
todolist+=[5*[["A1",[4,0]]]]

#SO(5)xSU(4)
todolist+=[[["B2",[0,2]],["A3",[1,0,1]]]]

#S³xSU(3)xG2
todolist+=[[["A1",[2,0]],["A2",[1,1]],["G2",[0,1]]]]

#5S³xSO(5)
todolist+=[5*[["A1",[2,0]]]+[["B2",[0,2]]]]

#3S³xSU(3)xSU(3)
todolist+=[3*[["A1",[2,0]]]+2*[["A2",[1,1]]]]

#5S⁵
todolist+=[5*[["B2",[1,0]]]]

#S³xSU(3)xSU(4)
todolist+=[[["A1",[2,0]],["A2",[1,1]],["A3",[1,0,1]]]]

#4S³xG2
todolist+=[4*[["A1",[2,0]]]+[["G2",[0,1]]]]

#SU(3)xSU(3)xSO(5)
todolist+=[[["A2",[1,1]],["A2",[1,1]],["B2",[0,2]]]]

#S³xS³xSO(5)xSO(5)
todolist+=[2*[["A1",[2,0]]]+2*[["B2",[0,2]]]]

#6S³xSU(3)
todolist+=[6*[["A1",[2,0]]]+[["A2",[1,1]]]]

#3S⁴xSU(6)/Sp(3)
todolist+=[3*[["A1xA1",[1,1]]]+[["C3",[0,1,0]]]]


def Orthstring(n):
    if n%2==0:
        return "D"+str(n/2)
    else:
        return "B"+str((n-1)/2)
#input: list of isotropy reps, e.g. inputlist=[["A1",[2,0]],["A2",[1,1]]]
def XIXbranch(inputlist):
    l=len(inputlist)
    ctlist=[x[0] for x in inputlist]
    groupslist=[WeylCharacterRing(x,style="coroots") for x in ctlist]
    isotropieslist=[groupslist[i](*inputlist[i][1]) for i in range(0,l)]
    SOmlist=[Orthstring(x.degree()) for x in isotropieslist]
    sumdim=sum([x.degree() for x in isotropieslist])
    multibranch=[]
    for i in range(0,l):
        multibranch+=[branching_rule(SOmlist[i],ctlist[i]+str(tuple(inputlist[i][1])),"plethysm")]
    Gstr=Orthstring(sumdim)
    Hstr="x".join([x for x in ctlist])
    G=WeylCharacterRing(Gstr,style="coroots")
    H=WeylCharacterRing(Hstr,style="coroots")
    b=branching_rule(Gstr,"x".join(SOmlist),"orthogonal_sum")*branching_rule("x".join(SOmlist),Hstr,multibranch)
    print("Defined by the symmetric space with isotropy rep. "+str(isotropieslist))
    return [G,H,b]

#for x in todolist: Sym2LowerBounds(*XIXbranch(x))
