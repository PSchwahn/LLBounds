#E8B4 by B4->D8->E8
E8=WeylCharacterRing("E8",style="coroots")
B4=WeylCharacterRing("B4",style="coroots")
bE8B4bySpin=branching_rule("E8","D8","extended")*branching_rule("D8","B4(0,0,0,1)","plethysm")
Sym2LowerBounds(E8,B4,bE8B4bySpin)
