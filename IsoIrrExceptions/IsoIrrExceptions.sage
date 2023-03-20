#Isotropy irreducible exceptions (rank G <= 39)

todolist=[]

#C2A1
C2=WeylCharacterRing("C2",style="coroots")
A1=WeylCharacterRing("A1",style="coroots")
bC2A1=branching_rule("C2","A1",rule="symmetric_power")
todolist+=[[C2,A1,bC2A1]]

#C7C3
C7=WeylCharacterRing("C7",style="coroots")
C3=WeylCharacterRing("C3",style="coroots")
bC7C3=branching_rule("C7","C3(0,0,1)",rule="plethysm")
todolist+=[[C7,C3,bC7C3]]

#D8B4
D8=WeylCharacterRing("D8",style="coroots")
B4=WeylCharacterRing("B4",style="coroots")
bD8B4=branching_rule("D8","B4(0,0,0,1)",rule="plethysm")
todolist+=[[D8,B4,bD8B4]]

#B3G2
B3=WeylCharacterRing("B3",style="coroots")
G2=WeylCharacterRing("G2",style="coroots")
bB3G2=branching_rule("B3","G2(1,0)",rule="plethysm")
todolist+=[[B3,G2,bB3G2]]

#D7G2
D7=WeylCharacterRing("D7",style="coroots")
G2=WeylCharacterRing("G2",style="coroots")
bD7G2=branching_rule("D7","G2(0,1)",rule="plethysm")
todolist+=[[D7,G2,bD7G2]]

#G2A1
G2=WeylCharacterRing("G2",style="coroots")
bG2A1=branching_rule("G2","A1",rule="i")
todolist+=[[G2,A1,bG2A1]]

#G2A2
G2=WeylCharacterRing("G2",style="coroots")
A2=WeylCharacterRing("A2",style="coroots")
bG2A2=branching_rule("G2","A2",rule="extended")
todolist+=[[G2,A2,bG2A2]]

#F4A1G2
F4=WeylCharacterRing("F4",style="coroots")
A1G2=WeylCharacterRing("A1xG2",style="coroots")
bF4A1G2=branching_rule("F4","A1xG

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
E7=WeylCharacterRing("E7",style="coroots")
A1F4=WeylCharacterRing("A1xF4",style="coroots")
bE7A1F4=branching_rule("E7","A1xF4","miscellaneous")
todolist+=[[E7,A1F4,bE7A1F4]]

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

#A15D5
A15=WeylCharacterRing("A15",style="coroots")
D5=WeylCharacterRing("D5",style="coroots")
bA15D5=branching_rule("A15","D5(0,0,0,1,0)",rule="plethysm")
todolist+=[[A15,D5,bA15D5]]

#A26E6
A26=WeylCharacterRing("A26",style="coroots")
E6=WeylCharacterRing("E6",style="coroots")
bA26E6=branching_rule("A26","E6(1,0,0,0,0,0)",rule="plethysm")
todolist+=[[A26,E6,bA26E6]]

#C10A5
C10=WeylCharacterRing("C10",style="coroots")
A5=WeylCharacterRing("A5",style="coroots")
bC10A5=branching_rule("C10","A5(0,0,1,0,0)",rule="plethysm")
todolist+=[[C10,A5,bC10A5]]

#C16D6
C16=WeylCharacterRing("C16",style="coroots")
D6=WeylCharacterRing("D6",style="coroots")
bC16D6=branching_rule("C16","D6(0,0,0,0,1,0)",rule="plethysm")
todolist+=[[C16,D6,bC16D6]]

#C28E7
C28=WeylCharacterRing("C28",style="coroots")
E7=WeylCharacterRing("E7",style="coroots")
bC28E7=branching_rule("C28","E7(0,0,0,0,0,0,1)",rule="plethysm")
todolist+=[[C28,E7,bC28E7]]

#D10A3
D10=WeylCharacterRing("D10",style="coroots")
A3=WeylCharacterRing("A3",style="coroots")
bD10A3=branching_rule("D10","A3(0,2,0)",rule="plethysm")
todolist+=[[D10,A3,bD10A3]]

#D21C4
D21=WeylCharacterRing("D21",style="coroots")
C4=WeylCharacterRing("C4",style="coroots")
bD21C4=branching_rule("D21","C4(0,0,0,1)",rule="plethysm")
todolist+=[[D21,C4,bD21C4]]

#D13F4
D13=WeylCharacterRing("D13",style="coroots")
F4=WeylCharacterRing("F4",style="coroots")
bD13F4=branching_rule("D13","F4(0,0,0,1)",rule="plethysm")
todolist+=[[D13,F4,bD13F4]]

#D26F4
D26=WeylCharacterRing("D26",style="coroots")
F4=WeylCharacterRing("F4",style="coroots")
bD26F4=branching_rule("D26","F4(1,0,0,0)",rule="plethysm")
todolist+=[[D26,F4,bD26F4]]

#D35A7
D35=WeylCharacterRing("D35",style="coroots")
A7=WeylCharacterRing("A7",style="coroots")
bD35A7=branching_rule("D35","A7(0,0,0,1,0,0,0)",rule="plethysm")
todolist+=[[D35,A7,bD35A7]]

#D39E6
D39=WeylCharacterRing("D39",style="coroots")
bD39E6=branching_rule("D39","E6(0,1,0,0,0,0)",rule="plethysm")
todolist+=[[D39,E6,bD39E6]]

for item in todolist:
    print("G="+str(item[0])+", H="+str(item[1]))
    Sym2LowerBounds(*item)
