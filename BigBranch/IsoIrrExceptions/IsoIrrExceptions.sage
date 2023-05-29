#Isotropy irreducible exceptions (rank G >= 35)

todolist=[
["D35","A7","1X[0,0,0,1,0,0,0]"],
["D39","E6","1X[0,1,0,0,0,0]"],
["D64","D8","1X[0,0,0,0,0,0,1,0]"],
["B66","E7","1X[1,0,0,0,0,0,0]"],
["D124","E8","1X[0,0,0,0,0,0,0,1]"],
]

for item in todolist: Sym2LowerBoundsBig(*item)
