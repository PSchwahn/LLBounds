#E8B4 by B4->A8->E8
E8=WeylCharacterRing("E8",style="coroots")
B4=WeylCharacterRing("B4",style="coroots")
bE8B4byA=branching_rule("E8","A8","extended")*branching_rule("A8","B4(1,0,0,0)","plethysm")
Sym2LowerBounds(E8,B4,bE8B4byA)
