#kleinvieh4: exceptional non-isotropy-irreducible spaces with semisimple isotropy
#left out: E6/D4T2 (unstable by Lauret) and full flags of E6,E7,E8 (G-stable)
todolist=[]

#D4G2
#unstable by Lauret
D4=WeylCharacterRing("D4",style="coroots")
G2=WeylCharacterRing("G2",style="coroots")
bD4G2=branching_rule("D4","B3","symmetric")*branching_rule("B3","G2(1,0)","plethysm")
todolist+=[[D4,G2,bD4G2]]

#D13A1C5D3
#unstable by Lauret
D13=WeylCharacterRing("D13",style="coroots")
A1C5D3=WeylCharacterRing("A1xC5xD3",style="coroots")
b1=branching_rule("D13","D10xD3","orthogonal_sum")
b2=branching_rule("D10","C1xC5","tensor")
b3=branching_rule("C1","A1","isomorphic")
b4=branching_rule("C1xC5","A1xC5",[b3,"identity"])
b6=branching_rule("D3","A3","isomorphic")
b7=branching_rule("A3","D3","isomorphic")
#this seems useless, but it necessary because of a bug in branching_rule
bD13A1C5D3=b1*branching_rule("D10xD3","A1xC5xD3",[b2*b4,b6*b7])
todolist+=[[D13,A1C5D3,bD13A1C5D3]]

#F4D4
#unstable by Lauret
F4=WeylCharacterRing("F4",style="coroots")
D4=WeylCharacterRing("D4",style="coroots")
bF4D4=branching_rule("F4","B4","extended")*branching_rule("B4","D4","extended")
todolist+=[[F4,D4,bF4D4]]

#E6A1A1A1
#G-semistable by Lauret, IED unknown
E6=WeylCharacterRing("E6",style="coroots")
A1A1A1=WeylCharacterRing("A1xA1xA1",style="coroots")
b1=branching_rule("E6","A2xA2xA2","extended")
b2=branching_rule("A2","A1(2,0)","plethysm")
bE6A1A1A1=b1*branching_rule("A2xA2xA2","A1xA1xA1",[b2,b2,b2])
todolist+=[[E6,A1A1A1,bE6A1A1A1]]

#E6A1D3
#G-semistable by Lauret, has IED
A1D3=WeylCharacterRing("A1xD3",style="coroots")
b1=branching_rule("E6","A1xA5","extended")
b2=branching_rule("A5","D3(1,0,0)","plethysm")
bE6A1D3=b1*branching_rule("A1xA5","A1xD3",["identity",b2])
todolist+=[[E6,A1D3,bE6A1D3]]

#E77A1
E7=WeylCharacterRing("E7",style="coroots")
SevenA1=WeylCharacterRing("A1xA1xA1xA1xA1xA1xA1", style="coroots")
b1=branching_rule("E7","A1xD6","extended")
b2=branching_rule("D6","D2xD2xD2","orthogonal_sum")
b3=branching_rule("D2","A1xA1","isomorphic")
b4=branching_rule("D2xD2xD2","A1xA1xA1xA1xA1xA1",[b3,b3,b3])
bE77A1=b1*branching_rule("A1xD6","A1xA1xA1xA1xA1xA1xA1",["identity",b2*b4])
todolist+=[[E7,SevenA1,bE77A1]]

#E7D4
bE7D4=branching_rule("E7","A7","extended")*branching_rule("A7","D4(1,0,0,0)","plethysm")
todolist+=[[E7,D4,bE7D4]]

#E73A1D4
#unstable by Lauret
ThreeA1D4=WeylCharacterRing("A1xA1xA1xD4",style="coroots")
b1=branching_rule("E7","A1xD6","extended")
b2=branching_rule("D6","D2xD4","orthogonal_sum")
b3=branching_rule("D2","A1xA1","isomorphic")
b4=branching_rule("A1xD6","A1xD2xD4",["identity",b2])
b5=branching_rule("A1xD2xD4","A1xA1xA1xD4",["identity",b3,"identity"])
bE73A1D4=b1*b4*b5
todolist+=[[E7,ThreeA1D4,bE73A1D4]]

#E88A1
E8=WeylCharacterRing("E8",style="coroots")
EightA1=WeylCharacterRing("A1xA1xA1xA1xA1xA1xA1xA1",style="coroots")
b1=branching_rule("E8","D8","extended")
b2=branching_rule("D8","D2xD2xD2xD2","orthogonal_sum")
b3=branching_rule("D2","A1xA1","isomorphic")
b4=branching_rule("D2xD2xD2xD2","A1xA1xA1xA1xA1xA1xA1xA1",[b3,b3,b3,b3])
bE88A1=b1*b2*b4
todolist+=[[E8,EightA1,bE88A1]]

#E84A1
FourA1=WeylCharacterRing("A1xA1xA1xA1",style="coroots")
b1=branching_rule("E8","E6xA2","extended")
b2=branching_rule("E6","A2xA2xA2","extended")
b3=branching_rule("A2","A1(2,0)","plethysm")
b4=branching_rule("E6xA2","A2xA2xA2xA2",[b2,"identity"])
bE84A2=b1*b4
b5=branching_rule("A2xA2xA2xA2","A1xA1xA1xA1",[b3,b3,b3,b3])
bE84A1=bE84A2*b5
todolist+=[[E8,FourA1,bE84A1]]

#E84A2
FourA2=WeylCharacterRing("A2xA2xA2xA2",style="coroots")
todolist+=[[E8,FourA2,bE84A2]]

#E8A2A2
A2A2=WeylCharacterRing("A2xA2",style="coroots")
b1=branching_rule("E8","A8","extended")
b2=branching_rule("A8","A2xA2","tensor")
bE8A2A2=b1*b2
todolist+=[[E8,A2A2,bE8A2A2]]

#E8A4A4
A4A4=WeylCharacterRing("A4xA4",style="coroots")
bE8A4A4=branching_rule("E8","A4xA4","extended")
todolist+=[[E8,A4A4,bE8A4A4]]

#E8D4D4
D4D4=WeylCharacterRing("D4xD4",style="coroots")
bE8D4D4=branching_rule("E8","D8","extended")*branching_rule("D8","D4xD4","orthogonal_sum")
todolist+=[[E8,D4D4,bE8D4D4]]

#E8B2
B2=WeylCharacterRing("B2",style="coroots")
bE8B2=branching_rule("E8","B2","miscellaneous")
todolist+=[[E8,B2,bE8B2]]

#E8C2C2
C2C2=WeylCharacterRing("C2xC2",style="coroots")
bE8C2C2=branching_rule("E8","D8","extended")*branching_rule("D8","C2xC2(1,0,1,0)","plethysm")
todolist+=[[E8,C2C2,bE8C2C2]]

#E8B4 are in separate folders

for item in todolist:
    print("G="+str(item[0])+", H="+str(item[1]))
    Sym2LowerBounds(*item)
