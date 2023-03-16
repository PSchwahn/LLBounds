#now with correct A1 Casimir
#kleinvieh3: groups with A1 factors
todolist=[]

#C2A1
C2=WeylCharacterRing("C2",style="coroots")
A1=WeylCharacterRing("A1",style="coroots")
bC2A1=branching_rule("C2","A1",rule="symmetric_power")
todolist+=[[C2,A1,bC2A1]]

#G2A1
G2=WeylCharacterRing("G2",style="coroots")
bG2A1=branching_rule("G2","A1",rule="i")
todolist+=[[G2,A1,bG2A1]]

#F4A1G2
F4=WeylCharacterRing("F4",style="coroots")
A1G2=WeylCharacterRing("A1xG2",style="coroots")
bF4A1G2=branching_rule("F4","A1xG2","miscellaneous")
todolist+=[[F4,A1G2,bF4A1G2]]

#E7A1F4
E7=WeylCharacterRing("E7",style="coroots")
A1F4=WeylCharacterRing("A1xF4",style="coroots")
bE7A1F4=branching_rule("E7","A1xF4","miscellaneous")
todolist+=[[E7,A1F4,bE7A1F4]]

for item in todolist: Sym2LowerBounds5(*item)
