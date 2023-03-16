#no A1 factors - Casimir might not be correct
#kleinvieh2: isometry group exceptional
todolist=[]

#G2A1

#G2A2
G2=WeylCharacterRing("G2",style="coroots")
A2=WeylCharacterRing("A2",style="coroots")
bG2A2=branching_rule("G2","A2",rule="extended")
todolist+=[[G2,A2,bG2A2]]

#F4A1G2

#F4A2A2
F4=WeylCharacterRing("F4",style="coroots")
A2A2=WeylCharacterRing("A2xA2",style="coroots")
bF4A2A2=branching_rule("F4","A2xA2",rule="extended")
todolist+=[[F4,A2A2,bF4A2A2]]

#E6A2
E6=WeylCharacterRing("E6",style="coroots")
bE6A2=branching_rule("E6","A2",rule="miscellaneous")
todolist+=[[E6,A2,bE6A2]]

#E6G2
bE6G2=branching_rule("E6","G2",rule="miscellaneous")
todolist+=[[E6,G2,bE6G2]]

#E6A2G2
A2G2=WeylCharacterRing("A2xG2",style="coroots")
bE6A2G2=branching_rule("E6","A2xG2",rule="miscellaneous")
todolist+=[[E6,A2G2,bE6A2G2]]

#E63A2
A2A2A2=WeylCharacterRing("A2xA2xA2",style="coroots")
bE6A2A2A2=branching_rule("E6","A2xA2xA2",rule="extended")
todolist+=[[E6,A2A2A2,bE6A2A2A2]]

#E7A2
E7=WeylCharacterRing("E7",style="coroots")
bE7A2=branching_rule("E7","A2",rule="miscellaneous")
todolist+=[[E7,A2,bE7A2]]

#E7C3G2
G2C3=WeylCharacterRing("G2xC3",style="coroots")
bE7G2C3=branching_rule("E7","G2xC3",rule="miscellaneous")
todolist+=[[E7,G2C3,bE7G2C3]]

#E7A1F4

#E7A2A5
A5A2=WeylCharacterRing("A5xA2",style="coroots")
bE7A5A2=branching_rule("E7","A5xA2",rule="extended")
todolist+=[[E7,A5A2,bE7A5A2]]

#E8G2F4
E8=WeylCharacterRing("E8",style="coroots")
G2F4=WeylCharacterRing("G2xF4",style="coroots")
bE8G2F4=branching_rule("E8","G2xF4",rule="miscellaneous")
todolist+=[[E8,G2F4,bE8G2F4]]

#E8A8
A8=WeylCharacterRing("A8",style="coroots")
bE8A8=branching_rule("E8","A8",rule="extended")
todolist+=[[E8,A8,bE8A8]]

#E8A2E6
E6A2=WeylCharacterRing("E6xA2",style="coroots")
bE8E6A2=branching_rule("E8","E6xA2",rule="extended")
todolist+=[[E8,E6A2,bE8E6A2]]

for item in todolist: Sym2LowerBounds5(*item)
