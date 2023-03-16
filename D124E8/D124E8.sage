G=WeylCharacterRing("D124",style="coroots")
H=WeylCharacterRing("E8",style="coroots")
b=branching_rule("D124","E8(0,0,0,0,0,0,0,1)",rule="plethysm")
Sym2LowerBounds(G,H,b)
