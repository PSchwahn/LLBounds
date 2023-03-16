#no A1 factors - Casimir might not be correct
#kleinvieh: rank(G) < 35, isometry group not exceptional
todolist=[]

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

#C7C3
C7=WeylCharacterRing("C7",style="coroots")
C3=WeylCharacterRing("C3",style="coroots")
bC7C3=branching_rule("C7","C3(0,0,1)",rule="plethysm")
todolist+=[[C7,C3,bC7C3]]

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

#D8B4
D8=WeylCharacterRing("D8",style="coroots")
B4=WeylCharacterRing("B4",style="coroots")
bD8B4=branching_rule("D8","B4(0,0,0,1)",rule="plethysm")
todolist+=[[D8,B4,bD8B4]]

#D21C4
D21=WeylCharacterRing("D21",style="coroots")
C4=WeylCharacterRing("C4",style="coroots")
bD21C4=branching_rule("D21","C4(0,0,0,1)",rule="plethysm")
todolist+=[[D21,C4,bD21C4]]

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
#D35=WeylCharacterRing("D35",style="coroots")
#A7=WeylCharacterRing("A7",style="coroots")
#bD35A7=branching_rule("D35","A7(0,0,0,1,0,0,0)",rule="plethysm")
#todolist+=[[D35,A7,bD35A7]]

#D64D8

#D39E6

#B66E7

#D124E8

#C2A1

for item in todolist: Sym2LowerBounds5(*item)
